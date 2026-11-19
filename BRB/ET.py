import glob
import json
import math
import os
import re
import subprocess
import threading
from pathlib import Path

import pandas as pd
import requests

import BRB.misc
from BRB.logger import log

# Up to POOL_SIZE PushButton workers can finish at once and POST results for
# the same flowcell_id. Serialise those POSTs defensively: if Parkour's
# results view does a read-modify-write over the flowcell's library list,
# concurrent POSTs would lose updates. The lock costs no throughput -- the
# POST is fast, and it's the Slurm queue wait that the thread pool overlaps.
_parkourPostLock = threading.Lock()


def _readDupFile(path):
    """(optical_dups, total) from a splitFastq stats file, None if empty or garbled."""
    try:
        fields = open(path).read().split()
    except OSError:
        return None
    if len(fields) != 2 or not all(f.isdigit() for f in fields):
        return None
    opt, total = map(int, fields)
    return (opt, total) if opt <= total else None


def getNReads(d):
    """
    Get the number of reads and % optical dupes from a directory.

    splitFastq writes "<optical_dups> <total_reads>" to <sample>.metrics
    (<sample>.duplicate.txt before ~2026-08). Runs processed before
    dissectBCL v1.3.5 can still have the file in the un-prefixed Project_
    directory, so look there too. The first usable file wins; an empty or
    garbled one is skipped.

    No stats file at all (machines that do not mark duplicates): count R1 reads
    and report 0% optical duplicates. Stats files present but all unusable:
    the optical duplicate rate is NaN (left out of the payload), not a made-up 0.
    """
    seen = False
    for base in (d, d.replace("FASTQC_", "")):
        for pat in ("*.metrics", "*.duplicate.txt"):
            for f in sorted(glob.glob(f"{base}/{pat}")):
                seen = True
                got = _readDupFile(f)
                if got:
                    optDupes, total = got
                    opt_frac = 100.0 * optDupes / total if total else float(-1)
                    return total - optDupes, opt_frac
    log.warning(
        f"getNReads - no usable *.metrics/*.duplicate.txt in {d} "
        f"({'all empty or malformed' if seen else 'none found'}); "
        "counting R1 reads"
    )
    fastq = glob.glob("{}/*_R1.fastq.gz".format(d.replace("FASTQC_", "")))
    if not fastq:
        return float("nan"), float("nan")
    # For machines in which we don't mark duplicates
    # Just count the number of reads in R1
    c1 = subprocess.Popen(["zcat", fastq[0]], stdout=subprocess.PIPE)
    res = subprocess.check_output(["wc", "-l"], stdin=c1.stdout)
    return int(res) // 4, float("nan") if seen else 0.0


def getOffSpeciesRate(d, org_label=None) -> float:
    """
    Parses kraken report for number of reads mapping to unexpected organisms
    """
    fname = glob.glob(f"{d}/*rep")[0]
    # match Parkour Organism.label (stable, admin-curated short code -- see
    # the "label" field on library_sample_shared.Organism) to kraken db
    # group. Organism.name is a free-text display string that Parkour staff
    # rename freely (e.g. the Nov-2024 organism list overhaul); label is not.
    org_map = {
        "h38": "humangrp",
        "mm10": "mousegrp",
        "m39": "mousegrp",
        "lambda": "lambdaphage",
        "nematode": "c-elegans",
        "drosophila": "flygrp",
    }
    if org_label not in org_map:
        log.info(
            f"getOffSpeciesRate - organism label {org_label} is not in the org_map!"
        )
        return 0
    off = 0
    with open(fname) as f:
        for line in f:
            if org_map[org_label] in line:
                off = 1 - (float(line.strip().split()[0]) / 100)
    # off-species actually means fraction of non-expected organism reads !
    # off-species reads vs of-species reads ;)
    log.info(f"confident reads for {fname} = {off}")
    return off


def getBaseStatistics(config, outputDir, samples_id, org_label=None):
    """
    Return a directionary with keys lib names and values:
    (sample name, nReads, off-species rate, % optical dupes)

    Also return the mapping of sample names to library names
    """
    baseDict = {}  # output dictionary
    s2l = {}  # sample to library dictionary
    _odir, adir = os.path.split(os.path.split(outputDir)[0])
    pdir = f"FASTQC_Project_{adir[9:]}"
    for sample in samples_id:
        for d in glob.glob(
            "{}/{}/{}/Sample_{}".format(
                config.get("Paths", "baseData"),
                config.get("Options", "runID"),
                pdir,
                sample,
            )
        ):
            libName = os.path.split(d)[1][7:]
            if len(glob.glob(f"{d}/*_R1_fastqc.zip")) == 0:
                continue  # Skip failed samples
            sampleName = glob.glob(f"{d}/*_R1_fastqc.zip")[0]
            sampleName = os.path.split(sampleName)[1][:-14]
            nReads, optDupes = getNReads(d)  # opt. dup.
            offRate = getOffSpeciesRate(d, org_label)
            baseDict[libName] = [sampleName, nReads, offRate, optDupes]
            s2l[sampleName] = libName
    return baseDict, s2l


def DNA(config, outputDir, baseDict, sample2lib):
    """
    Parse an output directory to get a dictionary of libraries and their associated values.

    Add % mapped, % dupped, and insert size to baseDict. Filter it for those actually in the output
    """
    # baseDict, sample2lib = getBaseStatistics(config, outputDir)
    # If we have RELACS, the sample2lib won't match what we find here.
    # We can re-parse the sampleSheet to upload actual statistics of the RELACS demuxed samples.
    if Path(outputDir, "RELACS_sampleSheet.txt").exists():
        # RELACS is a problem for parkour (matching is in sampleID / barcode level).
        # Just return a list of dicts with the previous info
        m = []
        for k, v in baseDict.items():
            m.append(
                {
                    "barcode": k,
                    "reads_pf_sequenced": v[1],
                    "confident_reads": v[2],
                    "optical_duplicates": v[3],
                }
            )
        log.info(f"ET - DNA module detected RELACS. Returning {m}")
        return m

    # % Mapped
    for fname in glob.glob(f"{outputDir}/Bowtie2/*.Bowtie2_summary.txt"):
        sampleName = os.path.basename(fname).split(".Bowtie2_summary")[0]
        if sampleName not in sample2lib:
            log.warning(
                f"DNA: {sampleName} is not in this run's samples ({outputDir}); skipped"
            )
            continue
        with open(fname) as fh:
            lastLine = fh.read().split("\n")[-2]
        mappedPercent = lastLine.split("%")[0]
        baseDict[sample2lib[sampleName]].append(float(mappedPercent))
        # % Duplicated
        dup_info = glob.glob(
            f"{outputDir}/multiQC/multiqc_data/multiqc_samtools_flagstat.txt"
        )[0]
        dup_df = pd.read_csv(
            dup_info, sep="\t", usecols=["Sample", "total_passed", "duplicates_passed"]
        )

        # Following condition handle sample-name mismatches caused by trailing '_' or any other sign in parsed filenames
        # Unknown (NaN, left out of the payload) if the sample is not present in dup_df.
        if sampleName not in list(dup_df["Sample"].astype(str)):
            dup_rate = float("nan")
        else:
            dup_df = dup_df.loc[dup_df["Sample"].astype(str) == sampleName]
            dup_rate = (
                dup_df["duplicates_passed"].values / dup_df["total_passed"].values * 100
            )
            dup_rate = dup_rate[0]
        baseDict[sample2lib[sampleName]].append(dup_rate)
        # Median insert size
        insert_size_info = os.path.join(
            outputDir, "deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv"
        )
        insert_size_df = pd.read_csv(insert_size_info, sep="\t")
        medInsertSize = insert_size_df.loc[
            insert_size_df["Unnamed: 0"]
            == "filtered_bam/" + sampleName + ".filtered.bam"
        ]
        medInsertSize = medInsertSize["Frag. Len. Median"].values[0]
        baseDict[sample2lib[sampleName]].append(int(medInsertSize))

    log.info(f"ET - DNA module parsed {baseDict}")

    # Reformat into a matrix
    m = []
    for k, v in baseDict.items():
        m.append(
            {
                "barcode": k,
                "reads_pf_sequenced": v[1],
                "confident_reads": v[2],
                "optical_duplicates": v[3],
                "dupped_reads": v[5],
                "mapped_reads": v[4],
                "insert_size": v[6],
            }
        )
    return m


def peakQC(outputDir, sampleName):
    """
    MACS2_QC/<sample>.filtered.BAM_peaks.qc.txt -> peak_count, frip,
    peak_genome_coverage. Empty dict when the run has no peak calling.
    """
    f = os.path.join(outputDir, "MACS2_QC", f"{sampleName}.filtered.BAM_peaks.qc.txt")
    if not os.path.exists(f):
        return {}
    row = pd.read_csv(f, sep="\t").iloc[0]
    return {
        "peak_count": int(row["peak_count"]),
        "frip": float(row["FRiP"]),
        "peak_genome_coverage": float(row["peak_genome_coverage"]),
    }


FRAGMENT_COLS = {
    "fragment_q1": "Frag. Len. 1st. Qu.",
    "fragment_median": "Frag. Len. Median",
    "fragment_q3": "Frag. Len. 3rd Qu.",
    "fragment_p10": "Frag. Len. 10%",
    "fragment_p90": "Frag. Len. 90%",
}


def fragmentQC(outputDir, sampleName):
    """
    deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv -> fragment_q1,
    fragment_median, fragment_q3, fragment_p10, fragment_p90 (bp);
    multiQC samtools flagstat read1_passed -> fragments_total (pairs before
    filtering). Keys are omitted when their source is missing.
    """
    out = {}
    f = os.path.join(
        outputDir, "deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv"
    )
    if os.path.exists(f):
        df = pd.read_csv(f, sep="\t", index_col=0)
        row = f"filtered_bam/{sampleName}.filtered.bam"
        if row in df.index:
            for key, col in FRAGMENT_COLS.items():
                out[key] = int(df.loc[row, col])
    f = os.path.join(outputDir, "multiQC/multiqc_data/multiqc_samtools_flagstat.txt")
    if os.path.exists(f):
        df = pd.read_csv(f, sep="\t", usecols=["Sample", "read1_passed"])
        df = df.loc[df["Sample"].astype(str) == sampleName]
        if len(df):
            out["fragments_total"] = int(df["read1_passed"].values[0])
    return out


def chromatin(config, outputDir, baseDict, sample2lib):
    """
    ChIP/CUT&Tag/CUT&RUN/ATAC: DNA-mapping stats plus fragment-length and
    peak QC when present.
    """
    m = DNA(config, outputDir, baseDict, sample2lib)
    barcode2sample = {v: k for k, v in sample2lib.items()}
    for entry in m:
        sampleName = barcode2sample.get(entry["barcode"])
        if sampleName is not None:
            entry.update(fragmentQC(outputDir, sampleName))
            entry.update(peakQC(outputDir, sampleName))
    return m


def _multiqcRows(outputDir, fileName, cols):
    """{sample: row} from a multiQC table, first row per sample; {} when the file is missing."""
    f = os.path.join(outputDir, "multiQC", "multiqc_data", fileName)
    if not os.path.exists(f):
        return {}
    rows = {}
    for row in pd.read_csv(f, sep="\t", usecols=["Sample", *cols]).to_dict("records"):
        rows.setdefault(str(row["Sample"]), row)
    return rows


def _pct(part, total):
    """part as % of total; NaN when total is 0, which dropNaN leaves out of the payload."""
    return 100.0 * part / total if total else float("nan")


def _starMapping(fname):
    """(mapped %, uniquely mapped %, multi-mapped %) from a STAR Log.final.out."""
    uniq = multimap = 0.0
    with open(fname) as f:
        for line in f:
            if "Uniquely mapped reads %" in line:
                uniq = float(line.strip().split("\t")[-1][:-1])
            elif "% of reads mapped to multiple loci" in line:
                multimap = float(line.strip().split("\t")[-1][:-1])
    return uniq + multimap, uniq, multimap


def RNA(config, outputDir, baseDict, sample2lib):
    """
    Add STAR mapping, flagstat duplication and featureCounts assignment to each library's base stats.

    A library needs all three sources. One that lacks any of them (or has no row for its
    sample) is sent base stats only, with a warning naming it; the other libraries keep theirs.
    """
    flagstat = _multiqcRows(
        outputDir,
        "multiqc_samtools_flagstat.txt",
        ["total_passed", "duplicates_passed"],
    )
    featurecounts = _multiqcRows(
        outputDir, "multiqc_featurecounts.txt", ["Total", "Assigned"]
    )
    rna = {}
    missing = {}
    for fname in glob.glob(f"{outputDir}/STAR/*/*.Log.final.out"):
        sampleName = os.path.basename(fname)[: -len(".Log.final.out")]
        if sampleName not in sample2lib:
            log.warning(
                f"RNA: {sampleName} is not in this run's samples ({outputDir}); skipped"
            )
            continue
        lib = sample2lib[sampleName]
        dup = flagstat.get(sampleName)
        assigned = featurecounts.get(f"{sampleName}.filtered")
        if dup is None or assigned is None:
            missing[lib] = [
                name
                for name, row in (
                    ("multiQC flagstat", dup),
                    ("multiQC featureCounts", assigned),
                )
                if row is None
            ]
            continue
        mapped, uniq, multimap = _starMapping(fname)
        rna[lib] = {
            "mapped_reads": mapped,
            "uniq_mapped": uniq,
            "multi_mapped": multimap,
            "dupped_reads": _pct(dup["duplicates_passed"], dup["total_passed"]),
            "assigned_reads": _pct(assigned["Assigned"], assigned["Total"]),
        }
    m = []
    for lib, v in baseDict.items():
        entry = {
            "barcode": lib,
            "reads_pf_sequenced": v[1],
            "confident_reads": v[2],
            "optical_duplicates": v[3],
        }
        if lib in rna:
            entry.update(rna[lib])
        else:
            sources = missing.get(lib, ["STAR Log.final.out"])
            log.warning(
                f"RNA: {lib} missing {', '.join(sources)} in {outputDir}; sending base statistics only"
            )
        m.append(entry)
    log.info(f"ET - RNA module parsed {m}")
    return m


def HiC(config, outputDir, baseDict, sample2lib):
    """
    Base statistics plus the unique high quality pair rate from hicexplorer.

    multiqc_hicexplorer.txt has one row per matrix, named <sample>_bs<binsize>.cool.
    Each library takes the rate of the sample its base statistics came from
    (baseDict[lib][0]): that sample's 10kb matrix if present, else its first row.
    """
    hic = f"{outputDir}/multiQC/multiqc_data/multiqc_hicexplorer.txt"
    rate = {}  # sample -> (binsize, unique high quality pair %)
    if os.path.exists(hic):
        df = pd.read_csv(
            hic,
            sep="\t",
            usecols=["Sample", "Pairs mappable, unique and high quality_pct"],
        )
        for name, pct in zip(
            df["Sample"].astype(str), df["Pairs mappable, unique and high quality_pct"]
        ):
            mt = re.fullmatch(r"(.*)_bs(\d+[A-Za-z]*)\.cool", name)
            sample, binsize = (mt.group(1), mt.group(2)) if mt else (name, None)
            if sample not in rate or (binsize == "10kb" and rate[sample][0] != "10kb"):
                rate[sample] = (binsize, pct)
    m = []
    for lib, v in baseDict.items():
        entry = {
            "barcode": lib,
            "reads_pf_sequenced": v[1],
            "confident_reads": v[2],
            "optical_duplicates": v[3],
        }
        if v[0] in rate:
            entry["uniq_mapped"] = rate[v[0]][1]
        m.append(entry)
    return m


def pairsStats(fname):
    """{key: int} from a pairtools *.dedup.stats file (tab separated key/value lines)."""
    out = {}
    with open(fname) as f:
        for line in f:
            key, _, val = line.rstrip("\n").partition("\t")
            if val.isdigit():
                out[key] = int(val)
    return out


def makePairs(config, outputDir, baseDict, sample2lib):
    """
    Base statistics plus the pairtools dedup summary of micro-C (makePairs pipeline).

    pairs/<sample>.pairs.dedup.stats, one per sample:
    mapped_reads = total_mapped / total, dupped_reads = total_dups / total_mapped (both %).
    A library without a usable stats file gets base statistics only.
    """
    m = []
    for lib, v in baseDict.items():
        entry = {
            "barcode": lib,
            "reads_pf_sequenced": v[1],
            "confident_reads": v[2],
            "optical_duplicates": v[3],
        }
        fname = f"{outputDir}/pairs/{v[0]}.pairs.dedup.stats"
        if os.path.exists(fname):
            st = pairsStats(fname)
            if "total" in st and "total_mapped" in st and "total_dups" in st:
                entry["mapped_reads"] = _pct(st["total_mapped"], st["total"])
                entry["dupped_reads"] = _pct(st["total_dups"], st["total_mapped"])
            else:
                log.warning(f"makePairs: {fname} lacks total/total_mapped/total_dups")
        else:
            log.warning(
                f"makePairs: {lib} has no {fname}; sending base statistics only"
            )
        m.append(entry)
    return m


def dropNaN(msg):
    """json.dumps writes NaN (not valid JSON); a 0/0 rate must be left out instead."""
    return [
        {k: v for k, v in e.items() if not (isinstance(v, float) and math.isnan(v))}
        for e in msg
    ]


def _parkourFailed(ret):
    """True when Parkour did not take the results: the POST raised (ret None) or HTTP >= 400."""
    return ret is None or ret.status_code >= 400


def sendToParkour(config, msg):
    msg = dropNaN(msg)
    basePath = config.get("Paths", "baseData")
    runID = config.get("Options", "runID")
    aviti_check = glob.glob(f"{basePath}/{runID}/RunManifest.csv")
    if aviti_check:
        FCID = config.get("Options", "runID").split("_")[2]
        if "-" in FCID:
            FCID = FCID.split("-")[-1]
        d = {"flowcell_id": FCID}
        d["sequences"] = json.dumps(msg)
    else:
        FCID = config.get("Options", "runID").split("_")[3][1:]
        if "-" in FCID:
            FCID = FCID.split("-")[-1]
        d = {"flowcell_id": FCID}
        d["sequences"] = json.dumps(msg)
    log.info(f"sendToParkour: Sending {d} to Parkour")
    with _parkourPostLock:
        try:
            res = requests.post(
                config.get("Parkour", "ResultsURL"),
                auth=(config.get("Parkour", "user"), config.get("Parkour", "password")),
                data=d,
                verify=config.get("Parkour", "cert"),
                timeout=60,
            )
        except requests.RequestException as e:
            log.error(
                f"sendToParkour: POST to Parkour failed, results not updated: {e}"
            )
            return None
    log.info(f"sendToParkour return {res}")
    if res.status_code >= 400:
        log.error(
            f"sendToParkour: Parkour answered HTTP {res.status_code} {res.reason}, results not updated"
        )
    return res


def phoneHome(
    config, outputDir, pipeline, samples_tuples, organism, org_label, project, libType
):
    """
    Return metrics to Parkour, the results are in outputDir and pipeline needs to be run on them
    """
    samples_id = [row[0] for row in samples_tuples]
    baseDict, sample2lib = getBaseStatistics(config, outputDir, samples_id, org_label)
    msg = None
    parsers = {
        "DNA": chromatin,
        "ATAC": chromatin,
        "RNA": RNA,
        "HiC": HiC,
        "makePairs": makePairs,
    }
    if pipeline in parsers:
        try:
            msg = parsers[pipeline](config, outputDir, baseDict, sample2lib)
        except Exception:  # noqa: BLE001 - any parser error: still send the base stats
            log.exception(
                f"phoneHome: parsing {pipeline} results of {project} failed, "
                "sending the base statistics only"
            )
    if msg is None:
        msg = [
            {
                "barcode": k,
                "reads_pf_sequenced": v[1],
                "confident_reads": v[2],
                "optical_duplicates": v[3],
            }
            for k, v in baseDict.items()
        ]
    log.info(f"phoneHome: got msg = {msg}")
    if msg is not None:
        ret = sendToParkour(config, msg)
    else:
        ret = None

    status = "FAILED" if _parkourFailed(ret) else "success"
    return [project, organism, libType, pipeline, status, ret]


def telegraphHome(config, group, project, skipList, organism=None, org_label=None):
    """
    The skipList is a list of samples/libraries for which we don't run a pipeline, but it'd be nice to still send back sequencing metrics
    Structure of skipList:
    [library, sampleName, analysisType]
    """
    log.info(f"telegraphHome triggered for {project}")
    # make a fake output directory path
    baseDir = Path(
        config.get("Paths", "groupData"),
        BRB.misc.pacifier(group),
        BRB.misc.getLatestSeqdir(config.get("Paths", "groupData"), group),
        config.get("Options", "runID"),
        f"Analysis_{BRB.misc.pacifier(project)}",
    )
    # Mock path
    outputDir = baseDir / "DNA_mouse"
    samples_id = [row[0] for row in skipList]
    baseDict, _sample2lib = getBaseStatistics(config, outputDir, samples_id, org_label)
    # Reformat into a matrix
    m = []
    for k, v in baseDict.items():
        m.append(
            {
                "barcode": k,
                "reads_pf_sequenced": v[1],
                "confident_reads": v[2],
                "optical_duplicates": v[3],
            }
        )
    ret = sendToParkour(config, m)
    # Format the libtypes
    libTypes = ",".join({i[2] for i in skipList})
    # [project, organism, libtypes, workflow, workflow status, parkour status, sambaUpdate]
    status = "FAILED" if _parkourFailed(ret) else None
    return [project, organism, libTypes, None, status, ret, False]
