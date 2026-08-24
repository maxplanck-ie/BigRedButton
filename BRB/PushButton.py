import concurrent.futures
import glob
import os
import shutil
import stat
from collections import namedtuple
from pathlib import Path

import BRB.ET
import BRB.jobtrack
import BRB.misc
from BRB.jobtrack import runManagedSubprocess
from BRB.logger import log

# One entry per (pipeline, organism, libraryType) group within a project.
# `project` is the raw, unpacified Parkour project name -- the pipeline
# functions pacify it themselves (RELACS needs the unpacified name to find
# its original sample sheet). `organism` is the (org_name, org_label,
# org_yaml) 3-tuple; org_label is deliberately not duplicated as its own
# field, so the two can't drift apart.
WorkItem = namedtuple(
    "WorkItem", ["project", "group", "pipeline", "organism", "libraryType", "tuples"]
)

# Name of the per-group ownership marker written into the group's outputDir.
MARKER_NAME = "running.pid"

# Fixed for now, deliberately not .ini-configurable. Kept low because not all
# of a group's work is Slurm-offloaded: RELACS runs `demultiplex_relacs -p 10`
# locally, each snakePipes driver is a real local process, and the --sequencer
# flag means two run_brb processes can already share this host.
POOL_SIZE = 2


class GroupDispatchError(RuntimeError):
    """
    Raised by runFlowcell when one or more worker threads raised. `failures`
    is a list of (WorkItem, exception) pairs; the first worker exception is
    chained as __cause__.
    """

    def __init__(self, message, failures):
        super().__init__(message)
        self.failures = failures


def parseExternalAllowlist(configValue):
    """Split a comma-separated ini value into a clean list of entries."""
    return [x.strip() for x in configValue.split(",") if x.strip()]


def isExternallyAllowed(
    libraryType, libraryProtocol, externalLibraryTypes, externalLibraryProtocols
):
    """
    Whether a sample should still be processed for an external (non-PI)
    project, based on its Parkour library type (exact match against
    externalLibraryTypes) or library protocol (startswith match against
    externalLibraryProtocols).
    """
    if libraryType in externalLibraryTypes:
        return True
    return any(libraryProtocol.startswith(p) for p in externalLibraryProtocols)


def createPath(config, group, project, org_label, libraryType, tuples):
    """Ensures that the output path exists, creates it otherwise, and return where it is"""
    if tuples[0][3]:
        baseDir = "{}/{}/Analysis_{}".format(
            config.get("Paths", "baseData"),
            config.get("Options", "runID"),
            BRB.misc.pacifier(project),
        )
    else:
        baseDir = "{}/{}/{}/{}/Analysis_{}".format(
            config.get("Paths", "groupData"),
            BRB.misc.pacifier(group),
            BRB.misc.getLatestSeqdir(config.get("Paths", "groupData"), group),
            config.get("Options", "runID"),
            BRB.misc.pacifier(project),
        )
    os.makedirs(baseDir, mode=0o700, exist_ok=True)
    oDir = os.path.join(baseDir, f"{BRB.misc.pacifier(libraryType)}_{org_label}")
    os.makedirs(oDir, mode=0o700, exist_ok=True)
    return oDir


def linkFiles(config, group, project, odir, tuples):
    """Create symlinks in odir to fastq files in {project}. Return 1 if paired-end, 0 otherwise."""
    if tuples[0][3]:
        baseDir = "{}/{}/Project_{}".format(
            config.get("Paths", "baseData"),
            config.get("Options", "runID"),
            BRB.misc.pacifier(project),
        )
    else:
        baseDir = "{}/{}/{}/{}/Project_{}".format(
            config.get("Paths", "groupData"),
            BRB.misc.pacifier(group),
            BRB.misc.getLatestSeqdir(config.get("Paths", "groupData"), group),
            config.get("Options", "runID"),
            BRB.misc.pacifier(project),
        )
    PE = False
    for t in tuples:
        currentName = "{}/{}_R1.fastq.gz".format(
            os.path.join(baseDir, f"Sample_{t[0]}"), t[1]
        )
        newName = f"{odir}/{t[1]}_R1.fastq.gz"
        if os.path.exists(currentName) and not os.path.exists(newName):
            os.symlink(currentName, newName)
        currentName = "{}/{}_R2.fastq.gz".format(
            os.path.join(baseDir, f"Sample_{t[0]}"), t[1]
        )
        newName = f"{odir}/{t[1]}_R2.fastq.gz"
        if os.path.exists(currentName):
            if not os.path.exists(newName):
                os.symlink(currentName, newName)
            PE = True
    return PE


def removeLinkFiles(d):
    """Remove symlinks created by linkFiles()"""
    files = glob.glob(f"{d}/originalFASTQ/*_R?.fastq.gz")
    if files:
        for fname in files:
            os.unlink(fname)
    files = glob.glob(f"{d}/*_R?.fastq.gz")
    for fname in files:
        os.unlink(fname)


def relinkFiles(config, group, project, org_label, libraryType, tuples):
    """
    Generate symlinks under the snakepipes originalFASTQ folder directly from the project folder.
    At this stage the multiqc files are copied over into the bioinfocoredir, as well
    (skipped for external/ignored projects, since bioinfoCoreDir is internal-only).
    """
    # relink fqs
    outputDir = createPath(config, group, project, org_label, libraryType, tuples)
    odir = os.path.join(outputDir, "originalFASTQ")
    linkFiles(config, group, project, odir, tuples)
    if tuples[0][3]:
        log.info(
            f"Skipping bioinfoCoreDir multiqc copy for external project {project}."
        )
        return
    # Copy mqc
    mqcf = os.path.join(outputDir, "multiQC", "multiqc_report.html")
    if os.path.exists(mqcf):
        log.info(f"Multiqc report found for {group} project {project}.")
        # Keyed on libraryType + org_label as well as project: two groups of
        # the same project run concurrently under the Phase 1 thread pool and
        # would otherwise interleave writes into one destination file.
        oname = (
            f"Analysis{project}_{BRB.misc.pacifier(libraryType)}"
            f"_{org_label}_multiqc.html"
        )
        of = Path(config.get("Paths", "bioinfoCoreDir")) / oname
        log.info(f"Trying to copy mqc report to {of}.")
        shutil.copyfile(mqcf, of)
    else:
        log.info(f"no multiqc report under {mqcf}.")


def getsambaPath(lane_dir, Sequencer):
    if Sequencer == "Aviti":
        current_year = str(lane_dir)[0:4]
        year_postfix = Path("Sequence_Quality_" + current_year) / Path(
            "AVITI24_" + current_year
        )
    else:
        current_year = "20" + str(lane_dir)[0:2]
        year_postfix = Path("Sequence_Quality_" + current_year) / Path(
            "Illumina_" + current_year
        )
    return current_year, year_postfix


def copyCellRanger(config, d, ignore=False):
    """
    copy Cellranger web_summaries to sequencing facility lane subdirectory & bioinfocore qc directory.
    e.g. /seqFacDir/Sequence_Quality_yyyy/Illumina_yyyy/flowcell_xxxx_lane_1/Analysis_xxx_sample_web_summary.html

    Skipped entirely when `ignore` is true (external/non-PI project), since
    both destinations are internal-only.

          :params config: configuration parsed from .ini file
          :params d: path to subdirectory of analysis folder, .e.g.
          /data/xxx/sequencing_data/yyyy_lanes_1/Analysis_2526_zzzz/RNA-Seq
          :params ignore: whether this is an external/non-PI project
          :type config: configparser.ConfigParser
          :type d: str
          :type ignore: bool
          :return: None
          :rtype: None
    """
    if ignore:
        log.info(f"Skipping copyCellRanger for external project under {d}.")
        return

    files = glob.glob(os.path.join(d, "*/outs/", "web_summary.html"))
    sequencing_type = config.get("Options", "sequencerType")

    # /data/xxx/yyyy_lanes_1/Analysis_2526_zzzz/RNA-Seqsinglecell_mouse ->
    # yyyy_lanes_1
    lane_dir = Path(d).parents[1].stem
    _current_year, year_postfix = getsambaPath(lane_dir, sequencing_type)
    for fname in files:
        # to seqfac dir.
        nname = fname.split("/")
        nname = "_".join([nname[-5], nname[-3], nname[-1]])
        # make lane directory in seqFacDir and copy it over
        seqfac_lane_dir = (
            Path(config.get("Paths", "seqFacDir")) / year_postfix / lane_dir
        )
        os.makedirs(seqfac_lane_dir, exist_ok=True)
        # Fetch flowcell ID, in case of reseq
        short_fid = str(os.path.basename(lane_dir)).split("_")[2] + "_"
        bioinfoCoreDirPath = Path(config.get("Paths", "bioinfoCoreDir")) / Path(
            short_fid + nname
        )
        nname = seqfac_lane_dir / nname
        shutil.copyfile(fname, nname)
        # to bioinfocore dir
        shutil.copyfile(fname, bioinfoCoreDirPath)


def copyRELACS(config, d):
    """
    copy RELACS demultiplexing png files to sequencing facility lane subdirectory.
    e.g. /seqFacDir/Sequence_Quality_yyyy/Illumina_yyyy/flowcell_xxxx_lane_1/xxx_RELACS_sample_fig.png

    Runs the same for internal and external (non-PI) projects: the
    sequencing facility (seqFacDir) and bioinfo-core (bioinfoCoreDir) QC
    diagnostics for a RELACS run are shared regardless of who owns the
    project, unlike the raw/processed data itself.

          :params config: configuration parsed from .ini file
          :params d: path to subdirectory of analysis folder, .e.g.
          /data/xxx/sequencing_data/yyyy_lanes_1/Analysis_2526_zzzz/ChIP-Seq_bla/RELACS_demultiplexing
          :type config: configparser.ConfigParser
          :type d: str
          :return: None
          :rtype: None
    """
    files = glob.glob(
        os.path.join(d, "RELACS_demultiplexing", "Sample*/", "*_fig.png")
    ) + glob.glob(os.path.join(d, "multiQC", "*html"))
    sequencing_type = config.get("Options", "sequencerType")

    # /data/xxx/yyyy_lanes_1/Analysis_2526_zzzz/ChIP-Seq_mouse/RELACS_demultiplexing ->
    # Sequence_Quality_yyyy/Illumina_yyyy/yyyy_lanes_1
    lane_dir = Path(d).parents[1].stem
    _current_year, year_postfix = getsambaPath(lane_dir, sequencing_type)
    log.info(f"copyRELACS - copying over RELACS files to samba path {year_postfix}")
    # `d` is .../Analysis_<proj>/<libraryType>_<orgLabel> -- its stem is
    # already the filesystem-safe (pacified) <libraryType>_<orgLabel>
    # component. Two library-groups of the same project run concurrently
    # under the Phase 1 thread pool and would otherwise collide on the same
    # destination filename in both seqFacDir and bioinfoCoreDir (the same
    # class of bug relinkFiles's multiqc naming was already fixed for).
    groupStem = Path(d).stem
    for fname in files:
        # to seqfac dir.
        nname = fname.split("/")
        if ".html" in fname:
            # ['', 'data', PI, seqdat, fid, analysis, libtype, multiqc, mqc.html]
            nname = "_".join([nname[-4], groupStem, "RELACS_analysis.html"])
        else:
            # ['', data, PI, seqdat, fid, analysis, libtype, RELACS_demultiplexing, Sample_1, mark_fig.png]
            nname = "_".join([nname[-5], groupStem, nname[-3], nname[-1]])
        # make lane directory in seqFacDir and copy it over
        seqfac_lane_dir = (
            Path(config.get("Paths", "seqFacDir")) / year_postfix / lane_dir
        )
        os.makedirs(seqfac_lane_dir, exist_ok=True)
        # bname must be computed from the plain (relative) nname, before it's
        # turned into an absolute seqFacDir path below -- Path(a) / b drops
        # `a` entirely when `b` is already absolute, which previously made
        # this silently re-copy to the seqFacDir path a second time instead
        # of ever reaching bioinfoCoreDir.
        bname = Path(config.get("Paths", "bioinfoCoreDir")) / nname
        nname = seqfac_lane_dir / nname
        shutil.copyfile(fname, nname)
        shutil.copyfile(fname, bname)


def tidyUpABit(d):
    """
    Reduce the number of files in the analysis folder.
    """
    for _d in ["clusters_logs", ".snakemake"]:
        _ = os.path.join(d, _d)
        if os.path.exists(_):
            shutil.rmtree(_)
    (Path(d) / "config.yaml").unlink(missing_ok=True)
    (Path(d) / "multiQC" / "multiqc_data" / "multiqc.log").unlink(missing_ok=True)
    (Path(d) / "multiQC" / "multiQC.out").unlink(missing_ok=True)
    (Path(d) / "multiQC" / "multiQC.err").unlink(missing_ok=True)
    (Path(d) / "config.yaml").unlink(missing_ok=True)


def stripRights(d):
    # Strip rights.
    for r, dirs, files in os.walk(d):
        for subdir in dirs:
            os.chmod(os.path.join(r, subdir), stat.S_IRWXU)
        for f in files:
            if not os.path.islink(os.path.join(r, f)):
                os.chmod(os.path.join(r, f), stat.S_IRWXU)


def touchDone(outputDir, fname="analysis.done"):
    open(os.path.join(outputDir, fname), "w").close()


def removeDone(outputDir):
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        os.remove(os.path.join(outputDir, "analysis.done"))


def RNA(config, group, project, organism, libraryType, tuples):
    """
    Need to set --libraryType
    """
    project = BRB.misc.pacifier(project)
    _org_name, org_label, org_yaml = organism
    outputDir = createPath(config, group, project, org_label, libraryType, tuples)
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        return outputDir, 0, False
    PE = linkFiles(config, group, project, outputDir, tuples)
    CMD = "PATH={}/bin:$PATH".format(
        os.path.join(config.get("Options", "snakemakeWorkflowBaseDir"))
    )
    CMD = [
        CMD,
        "mRNAseq",
        "--DAG",
        "--trim",
        "-i",
        outputDir,
        "-o",
        outputDir,
        org_yaml,
    ]
    if tuples[0][2].startswith("Smart-Seq2") or tuples[0][2].startswith(
        "NEBNext Single Cell RNA Library Preparation"
    ):
        # SMART-seq isn't a dUTP-based method!
        CMD.extend(["--libraryType", "0"])
    elif tuples[0][2].startswith("NEBNext Low Input RNA Library"):
        # Unstranded
        CMD.extend(
            [
                "--libraryType",
                "0",
                r"--trimmerOptions '-a AGATCGGAAGAGC -A AGATCGGAAGAGC'",
            ]
        )
    log.info(f"RNA wf CMD: {CMD}")
    try:
        runManagedSubprocess(" ".join(CMD))
    except:
        return outputDir, 1, False
    removeLinkFiles(outputDir)
    relinkFiles(config, group, project, org_label, libraryType, tuples)
    tidyUpABit(outputDir)
    touchDone(outputDir)
    return outputDir, 0, False


def RELACS(config, group, project, organism, libraryType, tuples):
    """
    This is a variant of the DNA mapping pipeline that does RELACS demultiplexing in addition

    This must check for the existence of a RELACS sample sheet in the run folder.

    There better not be any duplicate RELACS sample names!
    """
    runID = config.get("Options", "runID").split("_lanes")[0]
    sequencerType = config.get("Options", "sequencerType")
    _org_name, org_label, org_yaml = organism
    ignore = tuples[0][3]
    outputDir = createPath(
        config, group, BRB.misc.pacifier(project), org_label, libraryType, tuples
    )
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        if ignore and not os.path.exists(
            os.path.join(outputDir, "external_delivery.done")
        ):
            deliverExternalRELACS(config, outputDir, BRB.misc.pacifier(project))
        return outputDir, 0, True

    project = BRB.misc.pacifier(project)

    if sequencerType == "Aviti":
        # short_runs nests an extra instrument-id directory between the
        # sequencer prefix and the run itself, e.g.
        # short_runs/AVITI/AV251009/<runID>/RELACS_Project_<project>.txt
        matches = glob.glob(
            f"/dont_touch_this/short_runs/AV*/AV*/{runID}/RELACS_Project_{project}.txt"
        )
        sampleSheet = matches[0] if matches else None
    else:
        sampleSheet = (
            f"/dont_touch_this/short_runs/{runID}/RELACS_Project_{project}.txt"
        )

    # Fallback if exact path doesn't exist
    if not sampleSheet or (
        not os.path.exists(sampleSheet)
        and not os.path.exists(os.path.join(outputDir, "RELACS_sampleSheet.txt"))
    ):
        log.critical(f"RELACS: wrong samplesheet name: {sampleSheet}")
        return None, 1, False

    if ignore:
        # External (non-PI) project: the raw pooled Sample_* directories
        # live directly under baseData/runID/Project_{project}, since there's
        # no PI group directory under groupData for an external user.
        baseDir = "{}/{}/Project_{}".format(
            config.get("Paths", "baseData"),
            config.get("Options", "runID"),
            project,
        )
    else:
        baseDir = "{}/{}/{}/{}/Project_{}".format(
            config.get("Paths", "groupData"),
            BRB.misc.pacifier(group),
            BRB.misc.getLatestSeqdir(config.get("Paths", "groupData"), group),
            config.get("Options", "runID"),
            project,
        )

    # Link in files
    if not os.path.exists(os.path.join(outputDir, "RELACS_sampleSheet.txt")):
        shutil.copyfile(sampleSheet, os.path.join(outputDir, "RELACS_sampleSheet.txt"))

    # Only re-run RELACS demultiplexing if we don't have png files generated for every sample (pre-demux)
    # Infer number of samples from relacs samplesheet.
    premuxSamples = []
    with open(os.path.join(outputDir, "RELACS_sampleSheet.txt")) as f:
        for line in f:
            _s = line.strip().split("\t")[0]
            if _s not in premuxSamples:
                premuxSamples.append(_s)
    if len(premuxSamples) != len(
        list((Path(outputDir) / "RELACS_demultiplexing").rglob("*png"))
    ):
        unlinkDirs = []
        for d in glob.glob(f"{baseDir}/Sample_*"):
            bname = os.path.basename(d)
            newName = os.path.join(outputDir, bname)
            unlinkDirs.append(newName)
            if not os.path.exists(newName):
                os.symlink(d, newName)

        # -p 10 is pretty much arbitrary
        CMD = [
            "demultiplex_relacs",
            "--umiLength",
            "4",
            "-p",
            "10",
            os.path.join(outputDir, "RELACS_sampleSheet.txt"),
            os.path.join(outputDir, "RELACS_demultiplexing"),
        ]
        log.info(f"RELACS demux wf CMD: {CMD}")
        try:
            runManagedSubprocess(" ".join(CMD), cwd=outputDir)
        except:
            return outputDir, 1, False

        # clean up
        for d in unlinkDirs:
            os.unlink(d)

    # Link in the RELACS demultiplexed files
    for fname in glob.glob(
        os.path.join(outputDir, "RELACS_demultiplexing", "*", "*.gz")
    ):
        bname = os.path.basename(fname)
        if "unknown" not in bname:
            newName = os.path.join(outputDir, bname)
            if not os.path.exists(newName):
                os.symlink(fname, newName)

    # Back to the normal DNA pipeline
    CMD = "PATH={}/bin:$PATH".format(
        os.path.join(config.get("Options", "snakemakeWorkflowBaseDir"))
    )
    CMD = [
        CMD,
        "DNAmapping",
        "--DAG",
        "--trim",
        r"--trimmerOptions '-a AGATCGGAAGAG -A AGATCGGAAGAG'",
        "--UMIDedup",
        "--mapq",
        "3",
        "-i",
        outputDir,
        "-o",
        outputDir,
        org_yaml,
    ]
    log.info(f"RELACS DNA wf CMD: {CMD}")
    try:
        runManagedSubprocess(" ".join(CMD))
    except:
        return outputDir, 1, False
    removeLinkFiles(outputDir)
    tidyUpABit(outputDir)
    copyRELACS(config, outputDir)
    # Recreate links under originalFastQ
    for fname in glob.glob(
        os.path.join(outputDir, "RELACS_demultiplexing", "*", "*.gz")
    ):
        bname = os.path.basename(fname)
        if bname.startswith("unknown"):
            continue
        if not os.path.exists(os.path.join(outputDir, "originalFASTQ")):
            os.mkdir(os.path.join(outputDir, "originalFASTQ"))
        newName = os.path.join(outputDir, "originalFASTQ", bname)
        os.symlink(fname, newName)
    stripRights(outputDir)
    touchDone(outputDir)
    if ignore:
        deliverExternalRELACS(config, outputDir, project)
    return outputDir, 0, True


def deliverExternalRELACS(config, outputDir, project):
    """
    For external (non-PI) RELACS projects: once the pipeline has succeeded,
    reclaim space by truncating the bulky, fully-reproducible intermediate
    FASTQs (RELACS_demultiplexing/, FASTQ_Cutadapt/), then package the real
    result files (BAMs, bigwigs, multiQC, and run metadata -- no symlinks,
    no intermediates) into a tar.gz and deliver it via fexsend to the
    configured recipient.

    Guarded by a separate external_delivery.done marker (distinct from
    analysis.done), so if this step is interrupted after the pipeline itself
    succeeded, the next poll retries delivery instead of the outer
    analysis.done check silently skipping it forever.

    Failure here is logged but non-fatal: the analysis itself already
    succeeded, and delivery can be retried on the next poll or by hand.
    """
    doneMarker = os.path.join(outputDir, "external_delivery.done")
    if os.path.exists(doneMarker):
        return
    try:
        for subdir in ("RELACS_demultiplexing", "FASTQ_Cutadapt"):
            pattern = os.path.join(outputDir, subdir, "**", "*.gz")
            for fname in glob.glob(pattern, recursive=True):
                open(fname, "w").close()

        includeNames = ["Bowtie2", "filtered_bam", "bamCoverage", "multiQC"]
        for fname in sorted(os.listdir(outputDir)):
            full = os.path.join(outputDir, fname)
            if os.path.islink(full) or not os.path.isfile(full):
                continue
            if fname.endswith((".config.yaml", "_organism.yaml")) or (
                fname == "RELACS_sampleSheet.txt"
            ):
                includeNames.append(fname)
        includeNames = [
            n for n in includeNames if os.path.exists(os.path.join(outputDir, n))
        ]

        archiveName = f"Analysis_{project}.tar.gz"
        fexBin = os.path.expanduser(
            config.get("external", "fexsendBin", fallback="fexsend")
        )
        fexRecipient = config.get("external", "fexRecipient", fallback="")
        comment = f"BRB external delivery: project {project}"
        CMD = "tar -czf - {} | {} -s {} -C '{}' {}".format(
            " ".join(includeNames),
            fexBin,
            archiveName,
            comment,
            fexRecipient,
        )
        log.info(f"External delivery CMD: {CMD}")
        runManagedSubprocess(CMD, cwd=outputDir)
        open(doneMarker, "w").close()
    except:
        log.error(
            f"External delivery failed for {outputDir}; the analysis itself "
            "succeeded, delivery will be retried on the next poll."
        )


def DNA(config, group, project, organism, libraryType, tuples):
    """
    Run the DNA mapping pipeline on the samples. Tweals could theoretically be made
    according to the libraryProtocol (tuple[2])

    - Make /data/{group}/{LatestSeqdir}/{runID}/Analysis_{project}/{libraryType}_{org_label} directory
    - Remove previously linked in files (if any)
    - Link requested fastq files in
    - Run appropriate pipeline
    - Remove previously linked in files
    - Clean up snakemake directory
    """
    if tuples[0][2].startswith("ChIP RELACS high-throughput"):
        return RELACS(config, group, project, organism, libraryType, tuples)

    project = BRB.misc.pacifier(project)
    _org_name, org_label, org_yaml = organism
    outputDir = createPath(config, group, project, org_label, libraryType, tuples)
    log.debug("Running snakePipes in output dir " + outputDir)
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        return outputDir, 0, False
    PE = linkFiles(config, group, project, outputDir, tuples)
    log.debug(os.listdir(outputDir))
    CMD = "PATH={}/bin:$PATH".format(
        os.path.join(config.get("Options", "snakemakeWorkflowBaseDir"))
    )
    if libraryType == "CUTandTag-seq" or libraryType == "CUTandRUN-seq":
        CMD = [
            CMD,
            "DNAmapping",
            "--DAG",
            "--trim",
            "--dedup",
            "--mapq",
            "3",
            "--cutntag",
            "-i",
            outputDir,
            "-o",
            outputDir,
            org_yaml,
        ]
    elif libraryType == "ATAC-Seq":
        CMD = [
            CMD,
            "DNAmapping",
            "--DAG",
            "--trim",
            r"--trimmerOptions '-a nexteraF=CTGTCTCTTATA -A nexteraR=CTGTCTCTTATA'",
            "--dedup",
            "--mapq 2",
            "-i",
            outputDir,
            "-o",
            outputDir,
            org_yaml,
        ]
    else:
        CMD = [
            CMD,
            "DNAmapping",
            "--DAG",
            "--trim",
            "--dedup",
            "--mapq",
            "3",
            "-i",
            outputDir,
            "-o",
            outputDir,
            org_yaml,
        ]
    log.info(f"DNA wf CMD: {CMD}")
    try:
        runManagedSubprocess(" ".join(CMD))
    except:
        return outputDir, 1, False
    removeLinkFiles(outputDir)
    relinkFiles(config, group, project, org_label, libraryType, tuples)
    tidyUpABit(outputDir)
    stripRights(outputDir)
    touchDone(outputDir)
    return outputDir, 0, False


def WGBS(config, group, project, organism, libraryType, tuples):
    """
    Run the WGBS pipeline

    TODO: set trimming according to the libraryType
    TODO: I don't think we know how to send back metrics yet
    """

    project = BRB.misc.pacifier(project)
    _org_name, org_label, org_yaml = organism
    outputDir = createPath(config, group, project, org_label, libraryType, tuples)
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        return outputDir, 0, False
    PE = linkFiles(config, group, project, outputDir, tuples)
    CMD = "PATH={}/bin:$PATH".format(
        os.path.join(config.get("Options", "snakemakeWorkflowBaseDir"))
    )
    CMD = [CMD, "WGBS", "--DAG", "--trim", "-i", outputDir, "-o", outputDir, org_yaml]
    log.info(f"WGBS wf CMD: {CMD}")
    try:
        runManagedSubprocess(" ".join(CMD))
    except:
        return outputDir, 1, False
    removeLinkFiles(outputDir)
    relinkFiles(config, group, project, org_label, libraryType, tuples)
    tidyUpABit(outputDir)
    stripRights(outputDir)
    touchDone(outputDir)
    return outputDir, 0, False


def ATAC(config, group, project, organism, libraryType, tuples):
    """
    Run the DNA mapping pipeline and then the default ATAC pipeline
    """

    project = BRB.misc.pacifier(project)
    _org_name, org_label, org_yaml = organism
    outputDir = createPath(config, group, project, org_label, libraryType, tuples)
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        return outputDir, 0, False

    if not os.path.exists(os.path.join(outputDir, "DNA.done")):
        outputDir, rv, sambaret = DNA(
            config, group, project, organism, libraryType, tuples
        )
        if rv != 0:
            return outputDir, rv, sambaret

        removeDone(outputDir)
        touchDone(outputDir, "DNA.done")
    CMD = "PATH={}/bin:$PATH".format(
        os.path.join(config.get("Options", "snakemakeWorkflowBaseDir"))
    )
    CMD = [CMD, "ATACseq", "--DAG", "-d", outputDir, org_yaml]
    log.info(f"ATAC wf CMD: {CMD}")
    try:
        runManagedSubprocess(" ".join(CMD))
    except:
        return outputDir, 1, False
    tidyUpABit(outputDir)
    stripRights(outputDir)
    touchDone(outputDir)
    return outputDir, 0, False


def scRNAseq(config, group, project, organism, libraryType, tuples):
    """
    Run one of the scRNAseq pipelines (snakePipes or 10X)

    The protocol is tuples[0][2] and we assume they're all the same...

    We currently just skip unknown protocols and don't mention that!
    """

    project = BRB.misc.pacifier(project)
    _org_name, org_label, org_yaml = organism
    outputDir = createPath(config, group, project, org_label, libraryType, tuples)
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        return outputDir, 0, True

    accepted_names = [
        "Chromium_NextGEM_SingleCell3Prime_GeneExpression_v3.1_DualIndex",
        "Chromium_GEM-X_SingleCell_3primeRNA-seq_v4",
    ]

    if tuples[0][2] in accepted_names:
        if "GRCh38" in org_yaml:
            org_yaml = "GRCh38"
        PE = linkFiles(config, group, project, outputDir, tuples)
        snakeMakePath = "{}/bin".format(
            os.path.join(config.get("Options", "snakemakeWorkflowBaseDir"))
        )
        CMD = [
            config.get("10x", "RNA"),
            outputDir,
            outputDir,
            org_yaml,
            " --snakemakePath ",
            snakeMakePath,
        ]
        log.info(f"scRNA wf CMD: {' '.join(CMD)}")
        try:
            runManagedSubprocess(" ".join(CMD))
        except:
            return outputDir, 1, False
        removeLinkFiles(outputDir)
        tidyUpABit(outputDir)
        stripRights(outputDir)
        copyCellRanger(config, outputDir, tuples[0][3])
        sambaUpdate = True
    elif tuples[0][2] == "Cel-Seq 2 for single cell RNA-Seq":
        PE = linkFiles(config, group, project, outputDir, tuples)
        CMD = "PATH={}/bin:$PATH".format(
            os.path.join(config.get("Options", "snakemakeWorkflowBaseDir"))
        )
        CMD = [
            CMD,
            "scRNAseq",
            "--DAG",
            "--myKit CellSeq384",
            "--skipVelocyto",
            "-i",
            outputDir,
            "-o",
            outputDir,
            org_yaml,
        ]
        log.info(f"scRNA wf CMD: {CMD}")
        try:
            runManagedSubprocess(" ".join(CMD))
        except:
            return outputDir, 1, False
        removeLinkFiles(outputDir)
        tidyUpABit(outputDir)
        stripRights(outputDir)
        sambaUpdate = False
    else:
        log.info(f"Unsupported protocol: {tuples[0][2]}")
        sambaUpdate = False

    touchDone(outputDir)
    return outputDir, 0, sambaUpdate


def HiC(config, group, project, organism, libraryType, tuples):
    """
    Running the HiC pipeline on the samples.

    - Make /data/{group}/{LatestSeqdir}/{runID}/Analysis_{project}/{libraryType}_{org_label} directory
    - Remove previously linked in files (if any)
    - Link requested fastq files in
    - Run appropriate pipeline
    - Remove previously linked in files
    - Clean up snakemake directory
    """

    project = BRB.misc.pacifier(project)
    _org_name, org_label, org_yaml = organism
    outputDir = createPath(config, group, project, org_label, libraryType, tuples)
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        return outputDir, 0, False
    PE = linkFiles(config, group, project, outputDir, tuples)
    CMD = "PATH={}/bin:$PATH".format(
        os.path.join(config.get("Options", "snakemakeWorkflowBaseDir"))
    )
    CMD = [
        CMD,
        "HiC",
        "--DAG",
        "--noTAD",
        "--enzyme",
        "DpnII",
        "-i",
        outputDir,
        "-o",
        outputDir,
        org_yaml,
    ]
    log.info(f"HiC wf CMD: {CMD}")
    try:
        runManagedSubprocess(" ".join(CMD))
    except:
        return outputDir, 1, False
    removeLinkFiles(outputDir)
    relinkFiles(config, group, project, org_label, libraryType, tuples)
    tidyUpABit(outputDir)
    stripRights(outputDir)
    touchDone(outputDir)
    return outputDir, 0, False


def makePairs(config, group, project, organism, libraryType, tuples):
    """
    Running makePairs pipeline.
    """
    project = BRB.misc.pacifier(project)
    _org_name, org_label, org_yaml = organism
    outputDir = createPath(config, group, project, org_label, libraryType, tuples)
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        return outputDir, 0, False
    PE = linkFiles(config, group, project, outputDir, tuples)
    CMD = "PATH={}/bin:$PATH".format(
        os.path.join(config.get("Options", "snakemakeWorkflowBaseDir"))
    )
    CMD = [CMD, "makePairs", "--DAG", "-i", outputDir, "-o", outputDir, org_yaml]
    log.info(f"makePairs wf CMD: {CMD}")
    try:
        runManagedSubprocess(" ".join(CMD))
    except:
        return outputDir, 1, False
    removeLinkFiles(outputDir)
    relinkFiles(config, group, project, org_label, libraryType, tuples)
    tidyUpABit(outputDir)
    touchDone(outputDir)
    return outputDir, 0, False


def scATAC(config, group, project, organism, libraryType, tuples):
    """
    scATAC 10x
    """

    project = BRB.misc.pacifier(project)
    _org_name, org_label, org_yaml = organism
    outputDir = createPath(config, group, project, org_label, libraryType, tuples)
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        return outputDir, 0, True
    runID = config.get("Options", "runID").split("_lanes")[0]
    if (
        tuples[0][2] == "scATAC-Seq 10xGenomics"
        or tuples[0][2] == "Next GEM Single Cell ATAC"
        or tuples[0][2] == "Chromium Next GEM Single Cell ATAC v2"
    ):
        # PE = linkFiles(config, group, project, outputDir, tuples)
        samples = " ".join(i[1] for i in tuples)
        if tuples[0][3]:
            # External (non-PI) project: no group directory under groupData.
            inDir = "{}/{}/Project_{}".format(
                config.get("Paths", "baseData"),
                config.get("Options", "runID"),
                BRB.misc.pacifier(project),
            )
        else:
            inDir = "{}/{}/{}/{}/Project_{}".format(
                config.get("Paths", "groupData"),
                BRB.misc.pacifier(group),
                BRB.misc.getLatestSeqdir(config.get("Paths", "groupData"), group),
                config.get("Options", "runID"),
                BRB.misc.pacifier(project),
            )

        snakeMakePath = "{}/bin".format(
            os.path.join(config.get("Options", "snakemakeWorkflowBaseDir"))
        )
        CMD = config.get("10x", "ATAC") + " -i " + inDir
        CMD += " -o " + outputDir
        CMD += " " + org_yaml
        CMD += " --projectID " + project + " --samples " + samples
        CMD += " --snakemakePath " + snakeMakePath

        log.info(f"scATAC wf CMD: {CMD}")
        try:
            runManagedSubprocess(CMD)
        except:
            return outputDir, 1, False
        # removeLinkFiles(outputDir)
        copyCellRanger(config, outputDir, tuples[0][3])
        stripRights(outputDir)
        tidyUpABit(outputDir)
    else:
        log.info(f"Unsupported protocol: {tuples[0][2]}, just passing the data")
        touchDone(outputDir)
        return outputDir, 0, False

    touchDone(outputDir)
    return outputDir, 0, True


def runOneGroup(config, item, registry):
    """
    Ownership-guarded dispatch of a single WorkItem, coordinated through
    `registry` (a BRB.jobtrack.JobRegistry shared by every worker for this
    flowcell) and the on-disk marker jobtrack.py maintains in the group's
    outputDir.

    The marker is deliberately not a real lock (no flock, no atomicity
    guarantee). It exists only to stop a restarted run_brb from launching a
    second snakePipes run into an outputDir that a still-alive orphan from a
    crashed run owns. `registry.aborted` guards the in-process case: once the
    flowcell's crash handler has latched teardown, no further worker should
    start a fresh pipeline run.

    Always returns a list containing exactly one message row -- a SKIPPED row
    when the group was not dispatched (owned by a live PID, an unreadable
    marker, or the flowcell already being torn down), otherwise the OK/FAILED
    row that dispatching produces. Never returns bare None and never an
    unwrapped bare row, because runFlowcell does
    `msg.extend(future.result())` on this return value.
    """
    org_name, org_label, _org_yaml = item.organism
    pipelineFn = globals()[item.pipeline]
    outputDir = createPath(
        config,
        item.group,
        item.project,
        org_label,
        item.libraryType,
        item.tuples,
    )

    def skipRow(status):
        return [
            [
                item.project,
                org_name,
                item.libraryType,
                item.pipeline,
                status,
                "not updated",
                False,
                0,
            ]
        ]

    if registry.aborted:
        log.warning(f"Flowcell teardown in progress; not dispatching {outputDir}.")
        return skipRow("SKIPPED (flowcell teardown in progress)")

    state, marker = BRB.jobtrack.markerState(outputDir)
    if state == BRB.jobtrack.MARKER_LIVE:
        log.warning(
            f"Skipping {item.project} / {item.pipeline} / {item.libraryType}: "
            f"{outputDir} is owned by live pid {marker['pid']} "
            f"(job_ids={marker.get('job_ids', [])}). Not dispatching, and not "
            "counting this as a failure."
        )
        # A distinct "SKIPPED" status (not "FAILED"): this group was never
        # actually analysed this pass, so it must not be silently absent from
        # the report, but it also must not be double-counted as a pipeline
        # failure by anything that specifically checks for "FAILED" (retry
        # logic, email.finishedEmail's FAILED count).
        return skipRow("SKIPPED (owned by live PID)")
    if state == BRB.jobtrack.MARKER_CORRUPT:
        log.warning(
            f"{outputDir} has an unreadable {BRB.jobtrack.MARKER_NAME}; "
            "skipping this group this pass. Remove the file by hand once its "
            "owner is known to be dead."
        )
        return skipRow("SKIPPED (unreadable marker)")
    if state == BRB.jobtrack.MARKER_ABANDONED:
        log.warning(
            f"Removing abandoned {BRB.jobtrack.MARKER_NAME} in {outputDir} "
            f"(pid {marker['pid']}, cancelled={marker.get('cancelled', False)}, "
            f"job_ids={marker.get('job_ids', [])})."
        )
        BRB.jobtrack.clearMarker(outputDir)

    handle = registry.register_group(outputDir)
    BRB.jobtrack.writeMarker(outputDir)
    try:
        with BRB.jobtrack.bindGroup(handle):
            reruncount = 0
            outputDir, rv, sambaUpdate = pipelineFn(
                config,
                item.group,
                item.project,
                item.organism,
                item.libraryType,
                item.tuples,
            )
            if rv != 0:
                if registry.aborted:
                    log.warning(
                        f"{outputDir} failed while the flowcell was being torn "
                        "down; not retrying."
                    )
                    return skipRow("SKIPPED (flowcell teardown in progress)")
                # Allow for one re-run
                reruncount = 1
                outputDir, rv, sambaUpdate = pipelineFn(
                    config,
                    item.group,
                    item.project,
                    item.organism,
                    item.libraryType,
                    item.tuples,
                )
            if rv == 0:
                log.debug(
                    f"BRB run for {item.pipeline} {org_label} {item.libraryType} {item.tuples} complete."
                )
                log.info(
                    f"Processed project {BRB.misc.pacifier(item.project)} with the {item.pipeline} pipeline. {item.libraryType}, {org_name}. Rerun = {reruncount}"
                )
                return [
                    BRB.ET.phoneHome(
                        config,
                        outputDir,
                        item.pipeline,
                        item.tuples,
                        org_name,
                        item.project,
                        item.libraryType,
                    )
                    + [sambaUpdate, reruncount]
                ]
            log.warning(
                f"FAILED project {BRB.misc.pacifier(item.project)} with the {item.pipeline} pipeline. {item.libraryType}, {org_name}. Rerun = {reruncount}"
            )
            return [
                [
                    item.project,
                    org_name,
                    item.libraryType,
                    item.pipeline,
                    "FAILED",
                    "not updated",
                    sambaUpdate,
                    reruncount,
                ]
            ]
    finally:
        BRB.jobtrack.clearMarker(outputDir)
        registry.unregister_group(outputDir)


def GetResults(config, project, libraries):
    """
    Project is something like '352_Grzes_PearceEd' and libraries is a dictionary with libraries as keys:
        {'18L005489': ['FAT_first_A',
                       'Other',
                       'scRNA-Seq 10xGenomics',
                       'mouse'],
         '18L005490': ['FAT_first_B',
                       'Other',
                       'scRNA-Seq 10xGenomics',
                       'mouse'],
        }
    Returns (workItems, msg): a list of WorkItem entries to be dispatched by
    runFlowcell, and the message entries GetResults produces itself (the
    telegraphHome entry for skipped libraries, and the external-skip entry).
    It's assumed that everything within a single library type can be analysed together.
    """
    ignore = False
    try:
        group = project.split("_")[-1].split("-")[0].lower()
        group = BRB.misc.pacifier(group)
        dataPath = Path(
            config.get("Paths", "groupData"),
            group,
            BRB.misc.getLatestSeqdir(config.get("Paths", "groupData"), group),
            config.get("Options", "runID"),
            "Project_" + BRB.misc.pacifier(project),
        )
        log.info(f"Processing {dataPath}")
    except:
        ignore = True
    validLibraryTypes = {
        v: i
        for i, v in enumerate(config.get("Options", "validLibraryTypes").split(","))
    }
    pipelines = config.get("Options", "pipelines").split(",")
    # Library types/protocols that should still be processed for external
    # (non-PI) projects. Split into real lists up front: `x in "a,b,c"`
    # would otherwise do substring matching against the raw config string
    # (e.g. "RNA-Seq" would false-positive against "stranded mRNA-Seq").
    externalLibraryTypes = parseExternalAllowlist(
        config.get("external", "LibraryTypes", fallback="")
    )
    externalLibraryProtocols = parseExternalAllowlist(
        config.get("external", "LibraryProtocols", fallback="")
    )
    # split by analysis type and species, since we can only process some types of this
    analysisTypes = {}
    skipList = []
    external_skipList = []
    org_dict = {}
    for library, v in libraries.items():
        (
            sampleName,
            libraryType,
            libraryProtocol,
            organism,
            _indexType,
            _requestDepth,
        ) = v
        org_name, org_label, org_yaml = organism
        # Extra checks to see where we miss out
        if libraryType in validLibraryTypes:
            log.info(f"ValidLibraryType for sample {library} = {libraryType}")
        else:
            log.info(f"Not a ValidLibraryType for sample {library} = {libraryType}")
        if not (org_label or org_yaml):
            log.info(
                f"Species label or YAML was not set for {org_name} (check Parkour DB.)"
            )
        externallyAllowed = isExternallyAllowed(
            libraryType, libraryProtocol, externalLibraryTypes, externalLibraryProtocols
        )
        if (
            libraryType in validLibraryTypes
            and (org_label or org_yaml)
            and (ignore == False or externallyAllowed)
        ):
            if org_label not in org_dict:
                org_dict[org_label] = organism
            idx = validLibraryTypes[libraryType]
            pipeline = pipelines[idx]
            if pipeline not in analysisTypes:
                analysisTypes[pipeline] = {}
            if org_label not in analysisTypes[pipeline]:
                analysisTypes[pipeline][org_label] = {}
            if libraryType not in analysisTypes[pipeline][org_label]:
                analysisTypes[pipeline][org_label][libraryType] = []
            analysisTypes[pipeline][org_label][libraryType].append(
                [library, sampleName, libraryProtocol, ignore]
            )
            log.debug(f"Considering analysis types: {analysisTypes}")
        else:
            if ignore == False:
                skipList.append([library, sampleName, libraryType])
            else:
                external_skipList.append([library, sampleName, libraryType])
    msg = []
    if len(skipList):
        for i in skipList:
            log.info(
                f"Skipping sample {i[0]}/{i[0]} ({org_name} - project {BRB.misc.pacifier(project)}).\n"
            )
        msg = msg + [
            BRB.ET.telegraphHome(
                config, group, BRB.misc.pacifier(project), skipList, org_name
            )
        ]
    log.debug(config)
    workItems = []
    for pipeline, v in analysisTypes.items():
        log.debug("Queueing pipeline " + pipeline)
        for org_label, v2 in v.items():
            log.debug("Queueing organism label " + org_label)
            organism = org_dict[org_label]
            org_name, org_label, org_yaml = organism
            log.debug(organism)
            for libraryType, tuples in v2.items():
                log.debug("Queueing libraryType " + libraryType)
                log.debug(tuples)
                # RELACS needs the unpacified project name to copy the original sample sheet to the dest dir
                # hence the pacifier is applied on the project in each pipeline separately
                workItems.append(
                    WorkItem(
                        project=project,
                        group=group,
                        pipeline=pipeline,
                        organism=organism,
                        libraryType=libraryType,
                        tuples=tuples,
                    )
                )
    # In case there is an external_skipList, there shouldn't be a skipList !
    if external_skipList:
        assert not skipList
        libTypes = ",".join({i[2] for i in external_skipList})
        msg = msg + [[project, org_name, libTypes, None, None, None, False, None]]
    return workItems, msg


def runFlowcell(config, ParkourDict, registry=None):
    """
    Flowcell-wide coordinator: build every project's work items, dispatch them
    through a bounded thread pool, and collect their message entries.

    Results come back in completion order, not submission order.
    email.finishedEmail only makes order-independent checks over msg, so this
    is safe -- don't introduce positional assumptions here.

    On a crash (a worker future raising), the executor is shut down with
    `shutdown(wait=False, cancel_futures=True)` specifically so the error
    email goes out immediately, rather than waiting for surviving workers
    whose subprocess calls (Slurm queue time) can block for hours. This
    guarantees the alert is sent promptly. It does NOT, however, mean the
    `run_brb` process itself exits promptly: `concurrent.futures.thread`
    registers an interpreter-exit hook that joins every worker thread, and
    those threads are non-daemon, so after `run.py` re-raises the exception
    from this function, the process will still hang at interpreter shutdown
    until any surviving in-flight worker threads finish -- potentially
    hours later. This is an accepted tradeoff (the alert already went out
    by then); it is not a bug to "fix" by switching to `os._exit` or daemon
    threads.
    """
    bdir = "{}/{}".format(
        config.get("Paths", "baseData"), config.get("Options", "runID")
    )
    msg = []
    workItems = []
    for k, v in ParkourDict.items():
        if not os.path.exists(f"{bdir}/Project_{BRB.misc.pacifier(k)}"):
            log.info(
                f"{bdir}/Project_{BRB.misc.pacifier(k)} doesn't exist, probably lives on another lane."
            )
            continue
        projectItems, projectMsg = GetResults(config, k, v)
        workItems.extend(projectItems)
        msg.extend(projectMsg)

    if not workItems:
        log.info("runFlowcell: no work items to dispatch.")
        return msg

    log.info(f"runFlowcell: dispatching {len(workItems)} groups, POOL_SIZE={POOL_SIZE}")
    # Shared by every worker dispatched below so a fatal crash (see the
    # exception-handling pass further down) can be latched process-wide and
    # each in-flight group's Slurm jobs/local drivers can be found and
    # cancelled. Constructed fresh per flowcell -- never a module-level
    # global, since `run_brb -s illumina` and `run_brb -s aviti` share a host
    # but must never see each other's groups. Callers (namely tests) may pass
    # an existing registry in instead, so they can assert cancelAllGroups was
    # called on that exact object.
    if registry is None:
        registry = BRB.jobtrack.JobRegistry()
    # Deliberately not a `with` block: on the crash path we must shut down
    # with wait=False so the error e-mail goes out now, rather than after the
    # slowest surviving worker's Slurm jobs finish hours later.
    executor = concurrent.futures.ThreadPoolExecutor(max_workers=POOL_SIZE)
    futures = {
        executor.submit(runOneGroup, config, item, registry): item for item in workItems
    }

    failedFuture = None
    for future in concurrent.futures.as_completed(futures):
        exc = future.exception()
        if exc is not None:
            failedFuture = future
            # Cancel every in-flight Slurm job and driver for this flowcell
            # right away, before the drain pass below collects any other
            # already-done failures -- this is what keeps that drain fast,
            # since cancelled drivers exit instead of blocking for hours of
            # Slurm queue time. Must fire exactly once per crash.
            BRB.jobtrack.cancelAllGroups(registry)
            break
        msg.extend(future.result())

    if failedFuture is None:
        executor.shutdown(wait=True)
        return msg

    executor.shutdown(wait=False, cancel_futures=True)
    failures = [(futures[failedFuture], failedFuture.exception())]
    # Don't hide unrelated failures that happened around the same time: report
    # every worker exception that is already known, without waiting for the
    # workers that are still running.
    for future, item in futures.items():
        if future is failedFuture or future.cancelled() or not future.done():
            continue
        otherExc = future.exception()
        if otherExc is not None:
            failures.append((item, otherExc))
    for item, exc in failures:
        log.critical(
            f"runFlowcell: worker for project {item.project}, pipeline "
            f"{item.pipeline}, libraryType {item.libraryType} raised {exc!r}"
        )
    summary = "; ".join(
        f"{item.project} / {item.pipeline} / {item.libraryType}"
        for item, _exc in failures
    )
    raise GroupDispatchError(
        f"{len(failures)} library-group(s) raised during dispatch: {summary}",
        failures,
    ) from failures[0][1]
