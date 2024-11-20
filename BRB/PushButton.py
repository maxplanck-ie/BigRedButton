import os
import shutil
import glob
import subprocess
import BRB.ET
import BRB.misc
from BRB.logger import log
import stat
from pathlib import Path


def createPath(config, group, project, org_label, libraryType, tuples):
    """Ensures that the output path exists, creates it otherwise, and return where it is"""
    if tuples[0][3]:
        baseDir = "{}/{}/Analysis_{}".format(config.get('Paths', 'baseData'),
                                                            config.get('Options', 'runID'),
                                                            BRB.misc.pacifier(project))
    else:
        baseDir = "{}/{}/{}/{}/Analysis_{}".format(config.get('Paths', 'groupData'),
                                                            BRB.misc.pacifier(group),
                                                            BRB.misc.getLatestSeqdir(config.get('Paths','groupData'), group),
                                                            config.get('Options', 'runID'),
                                                            BRB.misc.pacifier(project))
    os.makedirs(baseDir, mode=0o700, exist_ok=True)
    oDir = os.path.join(baseDir, "{}_{}".format(BRB.misc.pacifier(libraryType), org_label))
    os.makedirs(oDir, mode=0o700, exist_ok=True)
    return oDir


def linkFiles(config, group, project, odir, tuples):
    """Create symlinks in odir to fastq files in {project}. Return 1 if paired-end, 0 otherwise."""
    if tuples[0][3]:
        baseDir = "{}/{}/Project_{}".format(config.get('Paths', 'baseData'),
                                                            config.get('Options', 'runID'),
                                                            BRB.misc.pacifier(project))
    else:
        baseDir = "{}/{}/{}/{}/Project_{}".format(config.get('Paths', 'groupData'),
                                                           BRB.misc.pacifier(group),
                                                           BRB.misc.getLatestSeqdir(config.get('Paths','groupData'), group),
                                                           config.get('Options', 'runID'),
                                                           BRB.misc.pacifier(project))
    PE = False
    for t in tuples:
        currentName = "{}/{}_R1.fastq.gz".format(os.path.join(baseDir, "Sample_{}".format(t[0])), t[1])
        newName = "{}/{}_R1.fastq.gz".format(odir, t[1])
        if os.path.exists(currentName):
            if not os.path.exists(newName):
                os.symlink(currentName, newName)
        currentName = "{}/{}_R2.fastq.gz".format(os.path.join(baseDir, "Sample_{}".format(t[0])), t[1])
        newName = "{}/{}_R2.fastq.gz".format(odir, t[1])
        if os.path.exists(currentName):
            if not os.path.exists(newName):
                os.symlink(currentName, newName)
            PE = True
    return PE


def removeLinkFiles(d):
    """Remove symlinks created by linkFiles()"""
    files = glob.glob("{}/originalFASTQ/*_R?.fastq.gz".format(d))
    if files:
        for fname in files:
            os.unlink(fname)
    files = glob.glob("{}/*_R?.fastq.gz".format(d))
    for fname in files:
        os.unlink(fname)


def relinkFiles(config, group, project, org_label, libraryType, tuples):
    """
    Generate symlinks under the snakepipes originalFASTQ folder directly from the project folder.
    At this stage the multiqc files are copied over into the bioinfocoredir, as well.
    """
    # relink fqs
    outputDir = createPath(config, group, project, org_label, libraryType, tuples)
    odir = os.path.join(outputDir, "originalFASTQ")
    linkFiles(config, group, project, odir, tuples)
    # Copy mqc
    mqcf = os.path.join(outputDir, 'multiQC', 'multiqc_report.html')
    if os.path.exists(mqcf):
        log.info(f"Multiqc report found for {group} project {project}.")
        oname = 'Analysis' + project + '_multiqc.html'
        of = Path(config.get('Paths', 'bioinfoCoreDir')) / oname
        log.info(f"Trying to copy mqc report to {of}.")
        shutil.copyfile(mqcf, of)
    else:
        log.info(f"no multiqc report under {mqcf}.")


def copyCellRanger(config, d):
    '''
    copy Cellranger web_summaries to sequencing facility lane subdirectory & bioinfocore qc directory.
    e.g. /seqFacDir/Sequence_Quality_yyyy/Illumina_yyyy/flowcell_xxxx_lane_1/Analysis_xxx_sample_web_summary.html
   
          :params config: configuration parsed from .ini file
          :params d: path to subdirectory of analysis folder, .e.g. 
          /data/xxx/sequencing_data/yyyy_lanes_1/Analysis_2526_zzzz/RNA-Seq
          :type config: configparser.ConfigParser
          :type d: str
          :return: None
          :rtype: None
    '''

    files = glob.glob(os.path.join(d, '*/outs/', 'web_summary.html'))

    # /data/xxx/yyyy_lanes_1/Analysis_2526_zzzz/RNA-Seqsinglecell_mouse ->
    # yyyy_lanes_1
    lane_dir = Path(d).parents[1].stem
    current_year = "20" + str(lane_dir)[0:2]
    year_postfix = Path("Sequence_Quality_" + current_year) / Path("Illumina_" + current_year)
    for fname in files:
        # to seqfac dir.
        nname = fname.split('/')
        nname = "_".join([nname[-5], nname[-3],nname[-1]])
        # make lane directory in seqFacDir and copy it over
        seqfac_lane_dir = Path(config.get('Paths', 'seqFacDir')) / year_postfix / lane_dir
        os.makedirs(seqfac_lane_dir, exist_ok=True)
        # Fetch flowcell ID, in case of reseq
        short_fid = str(os.path.basename(lane_dir)).split('_')[2] + '_'
        bioinfoCoreDirPath = Path(config.get('Paths', 'bioinfoCoreDir')) / Path(short_fid + nname)
        nname = seqfac_lane_dir / nname
        shutil.copyfile(fname, nname)
        # to bioinfocore dir
        shutil.copyfile(fname, bioinfoCoreDirPath)
    


def copyRELACS(config, d):
    '''
    copy RELACS demultiplexing png files to sequencing facility lane subdirectory.
    e.g. /seqFacDir/Sequence_Quality_yyyy/Illumina_yyyy/flowcell_xxxx_lane_1/xxx_RELACS_sample_fig.png
   
          :params config: configuration parsed from .ini file
          :params d: path to subdirectory of analysis folder, .e.g. 
          /data/xxx/sequencing_data/yyyy_lanes_1/Analysis_2526_zzzz/ChIP-Seq_bla/RELACS_demultiplexing
          :type config: configparser.ConfigParser
          :type d: str
          :return: None
          :rtype: None
    '''

    files = glob.glob(os.path.join(d, "RELACS_demultiplexing", 'Sample*/', '*_fig.png')) + glob.glob(os.path.join(d, "multiQC", '*html'))

    # /data/xxx/yyyy_lanes_1/Analysis_2526_zzzz/ChIP-Seq_mouse/RELACS_demultiplexing ->
    # Sequence_Quality_yyyy/Illumina_yyyy/yyyy_lanes_1
    lane_dir = Path(d).parents[1].stem
    current_year = "20" + str(lane_dir)[0:2]
    year_postfix = Path("Sequence_Quality_" + current_year) / Path("Illumina_" + current_year)
    log.info(f"copyRELACS - copying over RELACS files to samba path {year_postfix}")
    for fname in files:
        # to seqfac dir.
        nname = fname.split('/')
        if '.html' in fname:
            # ['', 'data', PI, seqdat, fid, analysis, libtype, multiqc, mqc.html]
            nname = "_".join([nname[-4], 'RELACS_analysis.html'])
        else:
            # ['', data, PI, seqdat, fid, analysis, libtype, RELACS_demultiplexing, Sample_1, mark_fig.png]
            nname = "_".join([nname[-5], nname[-3],nname[-1]])
        # make lane directory in seqFacDir and copy it over
        seqfac_lane_dir = Path(config.get('Paths', 'seqFacDir')) / year_postfix / lane_dir
        os.makedirs(seqfac_lane_dir, exist_ok=True)
        nname = seqfac_lane_dir / nname
        bname = Path(config.get('Paths', 'bioinfoCoreDir')) / nname
        shutil.copyfile(fname, nname)
        shutil.copyfile(fname, bname)


def tidyUpABit(d):
    """
    Reduce the number of files in the analysis folder.
    """
    for _d in ['clusters_logs', '.snakemake']:
        _ = os.path.join(d, _d)
        if os.path.exists(_):
            shutil.rmtree(_)
    (Path(d) / 'config.yaml').unlink(missing_ok=True)
    (Path(d) / 'multiQC' / 'multiqc_data' / 'multiqc.log').unlink(missing_ok=True)
    (Path(d) / 'multiQC' / 'multiQC.out').unlink(missing_ok=True)
    (Path(d) / 'multiQC' / 'multiQC.err').unlink(missing_ok=True)
    (Path(d) / 'config.yaml').unlink(missing_ok=True)


def stripRights(d):
    # Strip rights.
    for r, dirs, files in os.walk(d):
        for d in dirs:
            os.chmod(os.path.join(r, d), stat.S_IRWXU)
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
    org_name, org_label, org_yaml = organism
    outputDir = createPath(config, group, project, org_label, libraryType, tuples)
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        return outputDir, 0, False
    PE = linkFiles(config, group, project, outputDir, tuples)
    CMD = "PATH={}/bin:$PATH".format(os.path.join(config.get('Options', 'snakemakeWorkflowBaseDir')))
    CMD = [CMD, 'mRNAseq', '--DAG', '--trim', '-i', outputDir, '-o', outputDir, org_yaml]
    if tuples[0][2].startswith("Smart-Seq2"):
        # SMART-seq isn't a dUTP-based method!
        CMD.extend(['--libraryType', '0'])
    elif tuples[0][2].startswith("NEBNext Low Input RNA Library"):
        # Unstranded
        CMD.extend(['--libraryType', '0'])
    log.info(f"RNA wf CMD: {CMD}")
    try:
        subprocess.check_call(' '.join(CMD), shell=True)
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
    runID = config.get('Options', 'runID').split("_lanes")[0]
    org_name, org_label, org_yaml = organism
    outputDir = createPath(config, group, BRB.misc.pacifier(project), org_label, libraryType, tuples)
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        return outputDir, 0, True

    sampleSheet = "/dont_touch_this/solexa_runs/{}/RELACS_Project_{}.txt".format(runID, BRB.misc.pacifier(project))
    if not os.path.exists(sampleSheet) and not os.path.exists(os.path.join(outputDir, "RELACS_sampleSheet.txt")):
        log.critical("RELACS: wrong samplesheet name: {}".format(sampleSheet))
        print("wrong samplesheet name!", sampleSheet)
        return None, 1, False

    project = BRB.misc.pacifier(project)
    baseDir = "{}/{}/{}/{}/Project_{}".format(config.get('Paths', 'groupData'),
                                                           BRB.misc.pacifier(group),
                                                           BRB.misc.getLatestSeqdir(config.get('Paths','groupData'), group),
                                                           config.get('Options', 'runID'),
                                                           project)

    # Link in files
    if not os.path.exists(os.path.join(outputDir, "RELACS_sampleSheet.txt")):
        shutil.copyfile(sampleSheet, os.path.join(outputDir, "RELACS_sampleSheet.txt"))
    
    # Only re-run RELACS demultiplexing if we don't have png files generated for every sample (pre-demux)
    # Infer number of samples from relacs samplesheet.
    premuxSamples = []
    with open(os.path.join(outputDir, "RELACS_sampleSheet.txt")) as f:
        for line in f:
            _s = line.strip().split('\t')[0]
            if _s not in premuxSamples:
                premuxSamples.append(_s)
    if len(premuxSamples) != len(list((Path(outputDir) / 'RELACS_demultiplexing').rglob("*png"))):
        unlinkDirs = []
        for d in glob.glob("{}/Sample_*".format(baseDir)):
            bname = os.path.basename(d)
            newName = os.path.join(outputDir, bname)
            unlinkDirs.append(newName)
            if not os.path.exists(newName):
                os.symlink(d, newName)

        # -p 10 is pretty much arbitrary
        CMD = ["demultiplex_relacs", "--umiLength", "4", "-p", "10", os.path.join(outputDir, "RELACS_sampleSheet.txt"), os.path.join(outputDir, "RELACS_demultiplexing")]
        log.info(f"RELACS demux wf CMD: {CMD}")
        try:
            subprocess.check_call(' '.join(CMD), shell=True, cwd=outputDir)
        except:
            return outputDir, 1, False

        # clean up
        for d in unlinkDirs:
            os.unlink(d)

    # Link in the RELACS demultiplexed files
    for fname in glob.glob(os.path.join(outputDir, "RELACS_demultiplexing", "*", "*.gz")):
        bname = os.path.basename(fname)
        if 'unknown' not in bname:
            newName = os.path.join(outputDir, bname)
            if not os.path.exists(newName):
                os.symlink(fname, newName)


    # Back to the normal DNA pipeline
    CMD = "PATH={}/bin:$PATH".format(os.path.join(config.get('Options', 'snakemakeWorkflowBaseDir')))
    CMD = [CMD, 'DNAmapping', '--DAG', '--trim', '--UMIDedup', '--mapq', '3', '-i', outputDir, '-o', outputDir, org_yaml]
    log.info(f"RELACS DNA wf CMD: {CMD}")
    try:
        subprocess.check_call(' '.join(CMD), shell=True)
    except:
        return outputDir, 1, False
    removeLinkFiles(outputDir)
    tidyUpABit(outputDir)
    copyRELACS(config, outputDir)
    # Recreate links under originalFastQ
    for fname in glob.glob(os.path.join(outputDir, "RELACS_demultiplexing", "*", "*.gz")):
        bname = os.path.basename(fname)
        if bname.startswith('unknown'):
            continue
        if not os.path.exists(os.path.join(outputDir, 'originalFASTQ')):
            os.mkdir(os.path.join(outputDir, 'originalFASTQ'))
        newName = os.path.join(outputDir, 'originalFASTQ', bname)
        os.symlink(fname, newName)
    stripRights(outputDir)
    touchDone(outputDir)
    return outputDir, 0, True


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
    org_name, org_label, org_yaml = organism
    outputDir = createPath(config, group, project, org_label, libraryType, tuples)
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        return outputDir, 0, False
    PE = linkFiles(config, group, project, outputDir, tuples)
    CMD = "PATH={}/bin:$PATH".format(os.path.join(config.get('Options', 'snakemakeWorkflowBaseDir')))
    if libraryType == 'CUTandTag-seq' or libraryType == 'CUTandRUN-seq':
        CMD = [CMD, 'DNAmapping', '--DAG', '--trim', '--dedup', '--mapq', '3', '--cutntag', '-i', outputDir, '-o', outputDir, org_yaml]
    elif libraryType == 'ATAC-Seq':
        CMD = [CMD, 'DNAmapping', '--DAG', '--trim', r"--trimmerOptions '-a nexteraF=CTGTCTCTTATA -A nexteraR=CTGTCTCTTATA'", '--dedup', '--mapq 2', '-i', outputDir, '-o', outputDir, org_yaml]
    else:
        CMD = [CMD, 'DNAmapping', '--DAG', '--trim', '--dedup', '--mapq', '3', '-i', outputDir, '-o', outputDir, org_yaml]
    log.info(f"DNA wf CMD: {CMD}")
    try:
        subprocess.check_call(' '.join(CMD), shell=True)
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
    org_name, org_label, org_yaml = organism
    outputDir = createPath(config, group, project, org_label, libraryType, tuples)
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        return outputDir, 0, False
    PE = linkFiles(config, group, project, outputDir, tuples)
    CMD = "PATH={}/bin:$PATH".format(os.path.join(config.get('Options', 'snakemakeWorkflowBaseDir')))
    CMD = [CMD, 'WGBS', '--DAG', '--trim', '-i', outputDir, '-o', outputDir, org_yaml]
    log.info(f"WGBS wf CMD: {CMD}")
    try:
        subprocess.check_call(' '.join(CMD), shell=True)
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
    org_name, org_label, org_yaml = organism
    outputDir = createPath(config, group, project, org_label, libraryType, tuples)
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        return outputDir, 0, False

    if not os.path.exists(os.path.join(outputDir, "DNA.done")):
        outputDir, rv, sambaret  = DNA(config, group, project, organism, libraryType, tuples)
        if rv != 0:
            return outputDir, rv, sambaret

        removeDone(outputDir)
        touchDone(outputDir, "DNA.done")
    CMD = "PATH={}/bin:$PATH".format(os.path.join(config.get('Options', 'snakemakeWorkflowBaseDir')))
    CMD = [CMD, 'ATACseq', '--DAG', '-d', outputDir, org_yaml]
    log.info(f"ATAC wf CMD: {CMD}")
    try:
        subprocess.check_call(' '.join(CMD), shell=True)
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
    org_name, org_label, org_yaml = organism
    outputDir = createPath(config, group, project, org_label, libraryType, tuples)
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        return outputDir, 0, True

    accepted_names = [
        'Chromium_NextGEM_SingleCell3Prime_GeneExpression_v3.1_DualIndex',
        'Chromium_GEM-X_SingleCell_3primeRNA-seq_v4'
    ]
    if tuples[0][2] in accepted_names:
        PE = linkFiles(config, group, project, outputDir, tuples)
        # scRNA has their own organism mapping table, just make sure no spaces are included
        # FIXME: ^^^ this
        CMD = [config.get('10x', 'RNA'), outputDir, outputDir, org_yaml]
        log.info(f"scRNA wf CMD: {' '.join(CMD)}")
        try:
            subprocess.check_call(' '.join(CMD), shell=True)
        except:
            return outputDir, 1, False
        removeLinkFiles(outputDir)
        tidyUpABit(outputDir)
        stripRights(outputDir)
        copyCellRanger(config,outputDir)
        sambaUpdate = True
    elif tuples[0][2] == "Cel-Seq 2 for single cell RNA-Seq":
        PE = linkFiles(config, group, project, outputDir, tuples)
        CMD = "PATH={}/bin:$PATH".format(os.path.join(config.get('Options', 'snakemakeWorkflowBaseDir')))
        CMD = [CMD, 'scRNAseq', '--DAG', '--myKit CellSeq384','--skipVelocyto' , '-i', outputDir, '-o', outputDir, org_yaml]
        log.info(f"scRNA wf CMD: {CMD}")
        try:
            subprocess.check_call(' '.join(CMD), shell=True)
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
    org_name, org_label, org_yaml = organism
    outputDir = createPath(config, group, project, org_label, libraryType, tuples)
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        return outputDir, 0, False
    PE = linkFiles(config, group, project, outputDir, tuples)
    CMD = "PATH={}/bin:$PATH".format(os.path.join(config.get('Options', 'snakemakeWorkflowBaseDir')))
    CMD = [CMD, 'HiC', '--DAG', '--noTAD', '--enzyme', 'DpnII', '-i', outputDir, '-o', outputDir, org_yaml]
    log.info(f"HiC wf CMD: {CMD}")
    try:
        subprocess.check_call(' '.join(CMD), shell=True)
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
    org_name, org_label, org_yaml = organism
    outputDir = createPath(config, group, project, org_label, libraryType, tuples)
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        return outputDir, 0, False
    PE = linkFiles(config, group, project, outputDir, tuples)
    CMD = "PATH={}/bin:$PATH".format(os.path.join(config.get('Options', 'snakemakeWorkflowBaseDir')))
    CMD = [CMD, 'makePairs', '--DAG', '-i', outputDir, '-o', outputDir, org_yaml]
    log.info(f"makePairs wf CMD: {CMD}")
    try:
        subprocess.check_call(' '.join(CMD), shell=True)
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
    org_name, org_label, org_yaml = organism
    outputDir = createPath(config, group, project, org_label, libraryType, tuples)
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        return outputDir, 0, True
    runID = config.get('Options', 'runID').split("_lanes")[0]
    if (
        tuples[0][2] == "scATAC-Seq 10xGenomics"
        or tuples[0][2] == "Next GEM Single Cell ATAC"
        or tuples[0][2] == "Chromium Next GEM Single Cell ATAC v2"
    ):
        # PE = linkFiles(config, group, project, outputDir, tuples)
        samples = ' '.join(i[1] for i in tuples)
        inDir = "{}/{}/{}/{}/Project_{}".format(config.get('Paths', 'groupData'),
                                               BRB.misc.pacifier(group),
                                               BRB.misc.getLatestSeqdir(config.get('Paths','groupData'), group),
                                               config.get('Options', 'runID'),
                                               BRB.misc.pacifier(project))
        CMD = config.get('10x', 'ATAC')+" -i "+inDir
        CMD += " -o "+outputDir
        CMD += " "+org_yaml
        CMD += " --projectID "+project+" --samples "+samples
        log.info(f"scATAC wf CMD: {CMD}")
        try:
            subprocess.check_call(CMD, shell=True)
        except:
            return outputDir, 1, False
        # removeLinkFiles(outputDir)
        copyCellRanger(config,outputDir)
        stripRights(outputDir)
        tidyUpABit(outputDir)
    else:
        log.info(f"Unsupported protocol: {tuples[0][2]}, just passing the data")
        touchDone(outputDir)
        return outputDir, 0, False
    
    touchDone(outputDir)
    return outputDir, 0, True


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

    This doesn't return anything. It's assumed that everything within a single library type can be analysed together.
    """
    ignore = False
    try:
        group = project.split("_")[-1].split("-")[0].lower()
        group = BRB.misc.pacifier(group)
        dataPath = Path(
            config.get('Paths', 'groupData'),
            group,
            BRB.misc.getLatestSeqdir(config.get('Paths','groupData'), group),
            config.get('Options', 'runID'),
            'Project_' + BRB.misc.pacifier(project)
        )
        log.info(f"Processing {dataPath}")
    except:
        print("external data")
        ignore = True
    validLibraryTypes = {v: i for i, v in enumerate(config.get('Options', 'validLibraryTypes').split(','))}
    pipelines = config.get('Options', 'pipelines').split(',')
    # split by analysis type and species, since we can only process some types of this
    analysisTypes = dict()
    skipList = []
    external_skipList = []
    for library, v in libraries.items():
        sampleName, libraryType, libraryProtocol, organism, indexType, requestDepth = v
        org_name, org_label, org_yaml = organism
        # Extra checks to see where we miss out
        if libraryType in validLibraryTypes:
            log.info(f"ValidLibraryType = {libraryType}")
        else:
            log.info(f"Not a ValidLibraryType = {libraryType}")
        if not (org_label or org_yaml):
            log.info(f"Species label or YAML was not set for {org_name} (check Parkour DB.)")
        if libraryType in validLibraryTypes and (org_label or org_yaml) and (ignore==False or libraryType in config.get('external','LibraryTypes')):
            idx = validLibraryTypes[libraryType]
            pipeline = pipelines[idx]
            if pipeline not in analysisTypes:
                analysisTypes[pipeline] = dict()
            # TODO: iirc labels don't have the unique constrain on our DB,
            # and using org_yaml (which can be a path) feels awkward here, even if it works.
            # I can make labels unique and add this constraint, no problem.
            if org_label not in analysisTypes[pipeline]:
                analysisTypes[pipeline][org_label] = dict()
            if libraryType not in analysisTypes[pipeline][org_label]:
                analysisTypes[pipeline][org_label][libraryType] = list()
            analysisTypes[pipeline][org_label][libraryType].append([library, sampleName, libraryProtocol, ignore])
            log.debug(f"Considering analysis types: {analysisTypes}")
        else:
            if ignore == False:
               skipList.append([library, sampleName, libraryType])
            else:
               external_skipList.append([library, sampleName, libraryType])
    msg = []
    if len(skipList):
        for i in skipList:
            log.info(f"Skipping sample {i[0]}/{i[0]} ({org_name} - project {BRB.misc.pacifier(project)}).\n")
        msg = msg + [BRB.ET.telegraphHome(config, group, BRB.misc.pacifier(project), skipList, org_name)]
    for pipeline, v in analysisTypes.items():
        for org_label, v2 in v.items():
            for libraryType, tuples in v2.items():
                reruncount = 0
                # RELACS needs the unpacified project name to copy the original sample sheet to the dest dir
                # hence the pacifier is applied on the project in each pipeline separately
                outputDir, rv, sambaUpdate = globals()[pipeline](config, group, project, organism, libraryType, tuples)
                if reruncount == 0 and rv != 0:
                    # Allow for one re-run
                    reruncount += 1
                    outputDir, rv, sambaUpdate = globals()[pipeline](config, group, project, organism, libraryType, tuples)
                if rv == 0:
                    msg = msg + [BRB.ET.phoneHome(config, outputDir, pipeline, tuples, org_name, project, libraryType) + [sambaUpdate, reruncount]]
                    log.info(f"Processed project {BRB.misc.pacifier(project)} with the {pipeline} pipeline. {libraryType}, {org_name}. Rerun = {reruncount}")
                else:
                    msg = msg + [[project, org_name, libraryType, pipeline, 'FAILED', 'not updated', sambaUpdate, reruncount]]
                    log.warning(f"FAILED project {BRB.misc.pacifier(project)} with the {pipeline} pipeline. {libraryType}, {org_name}. Rerun = {reruncount}")
    # In case there is an external_skipList, there shouldn't be a skipList !
    if external_skipList:
        assert not skipList
        libTypes = ','.join(set([i[2] for i in external_skipList]))
        msg = msg + [[project, org_name, libTypes, None, None, None, False, None]]
    return msg
