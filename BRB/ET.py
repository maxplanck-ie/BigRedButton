import requests
import subprocess
import os
import glob
import json
import BRB.misc
from BRB.logger import log
import pandas as pd
from pathlib import Path

def getNReads(d):
    """
    Get the number of reads and % optical dupes from a directory
    """
    if len(glob.glob("{}/*.duplicate.txt".format(d))):
        fname = glob.glob("{}/*.duplicate.txt".format(d))[0]
        s = open(fname).read()
        optDupes, total = s.split()
        try:
            opt_frac = 100. * float(optDupes) / float(total)
        except ZeroDivisionError:
            opt_frac = float(-1)
        return int(total) - int(optDupes), opt_frac
    else:
        # For machines in which we don't mark duplicates
        # Just count the number of reads in R1
        CMD1 = ["zcat", glob.glob("{}/*_R1.fastq.gz".format(d.replace("FASTQC_", "")))[0]]
        CMD2 = ["wc", "-l"]
        c1 = subprocess.Popen(CMD1, stdout=subprocess.PIPE)
        res = subprocess.check_output(CMD2, stdin=c1.stdout)
        return (int(res) / 4), 0.


def getOffSpeciesRate(d, organism = None) -> float:
    """
    Parses 
    """
    fname = glob.glob("{}/*rep".format(d))[0]
    if not os.path.exists(fname):
        return 0
    # match parkour org to kraken db organism/group
    org_map = {
        'Human (GRCh38)': 'humangrp',
        'Human (GRCh37 / hg19)': 'humangrp',
        'Mouse (GRCm38 / mm10)': 'mousegrp',
        'Mouse (GRCm39)': 'mousegrp',
        'Escherichia phage Lambda':'lambdaphage',
        'Caenorhabditis_elegans': 'c-elegans',
        'lamprey': 'sea-lamprey',
        'medaka': 'japanese-medaka',
        'zebrafish': 'zebrafish',
        'drosophila': 'flygrp',
    }
    if organism not in org_map:
        return 0
    with open(fname) as f:
        for line in f:
            if org_map[organism] in line:
                off = 1-(float(line.strip().split()[0])/100)
    # off-species actually means fraction of non-expected organism reads !
    # off-species reads vs of-species reads ;)
    log.info(f"confident reads for {fname} = {off}")
    return off


def getBaseStatistics(config, outputDir, samples_id, organism = None):
    """
    Return a directionary with keys lib names and values:
    (sample name, nReads, off-species rate, % optical dupes)

    Also return the mapping of sample names to library names
    """
    baseDict = dict()  # output dictionary
    s2l = dict()  # sample to library dictionary
    odir, adir = os.path.split(os.path.split(outputDir)[0])
    pdir = "FASTQC_Project_{}".format(adir[9:])
    for sample in samples_id:
        for d in glob.glob("{}/{}/{}/Sample_{}".format(config.get('Paths','baseData'),
                                                       config.get('Options', 'runID'),
                                                       pdir, sample)):
            libName = os.path.split(d)[1][7:]
            if len(glob.glob("{}/*_R1_fastqc.zip".format(d))) == 0:
                continue  # Skip failed samples
            sampleName = glob.glob("{}/*_R1_fastqc.zip".format(d))[0]
            sampleName = os.path.split(sampleName)[1][:-14]
            nReads, optDupes = getNReads(d) # opt. dup.
            offRate = getOffSpeciesRate(d,organism)
            baseDict[libName] = [sampleName, nReads, offRate, optDupes]
            s2l[sampleName] = libName
    return baseDict, s2l


def DNA(config, outputDir, baseDict, sample2lib):
    """
    Parse an output directory to get a dictionary of libraries and their associated values.

    Add % mapped, % dupped, and insert size to baseDict. Filter it for those actually in the output
    """
    # baseDict, sample2lib = getBaseStatistics(config, outputDir)

    # % Mapped
    for fname in glob.glob("{}/Bowtie2/*.Bowtie2_summary.txt".format(outputDir)):
        sampleName = os.path.basename(fname).split(".Bowtie2_summary")[0]
        lastLine = open(fname).read().split("\n")[-2]
        mappedPercent = lastLine.split("%")[0]
        baseDict[sample2lib[sampleName]].append(float(mappedPercent))
        # % Duplicated
        dup_info = glob.glob("{}/multiQC/multiqc_data/multiqc_samtools_flagstat.txt".format(outputDir))[0]
        dup_df = pd.read_csv(dup_info, sep ="\t", usecols=["Sample", "total_passed", "duplicates_passed"])
        dup_df = dup_df.loc[dup_df["Sample"] == sampleName]
        dup_rate = dup_df["duplicates_passed"].values/dup_df["total_passed"].values*100
        dup_rate = dup_rate[0]
        baseDict[sample2lib[sampleName]].append(dup_rate)
        # Median insert size
        insert_size_info = os.path.join(outputDir, "deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv")
        insert_size_df = pd.read_csv(insert_size_info, sep ="\t")
        medInsertSize = insert_size_df.loc[insert_size_df["Unnamed: 0"]=="filtered_bam/"+sampleName+".filtered.bam"]
        medInsertSize = medInsertSize["Frag. Len. Median"].values[0]
        baseDict[sample2lib[sampleName]].append(int(medInsertSize))

    # # Filter
    outputDict = {k: v for k, v in baseDict.items() if len(v) == 8}
    # Reformat into a matrix
    m = []
    for k, v in outputDict.items():
        m.append({'barcode': k,
                  'reads_pf_sequenced': v[1],
                  'confident_reads': v[2],
                  'optical_duplicates': v[3],
                  'dupped_reads': v[6],
                  'mapped_reads': v[5],
                  'insert_size': v[7]})
    return m


def RNA(config, outputDir, baseDict, sample2lib):
    """
    Parse an output directory to get a dictionary of libraries and their associated values.

    Add % mapped to baseDict. Filter it for those actually in the output
    """
    # baseDict, sample2lib = getBaseStatistics(config, outputDir)
    # % Mapped
    for fname in glob.glob("{}/STAR/*/*.Log.final.out".format(outputDir)):
        f = open(fname)
        tot = 0
        uniq = 0
        multimap = 0
        for line in f:
            if 'Uniquely mapped reads %' in line:
                uniq = float(line.strip().split("\t")[-1][:-1])
                tot += uniq
            elif '% of reads mapped to multiple loci' in line:
                multimap = float(line.strip().split("\t")[-1][:-1])
                tot += multimap
        sampleName = os.path.basename(fname).split(".")[0]
        baseDict[sample2lib[sampleName]].append(tot)
        baseDict[sample2lib[sampleName]].append(uniq)
        baseDict[sample2lib[sampleName]].append(multimap)
        #  duplication
        dup_info = glob.glob("{}/multiQC/multiqc_data/multiqc_samtools_flagstat.txt".format(outputDir))[0]
        dup_df = pd.read_csv(dup_info, sep ="\t", usecols=["Sample", "total_passed", "duplicates_passed"])
        dup_df = dup_df.loc[dup_df["Sample"] == sampleName]
        dup_rate = dup_df["duplicates_passed"].values/dup_df["total_passed"].values*100
        dup_rate = dup_rate[0]
        baseDict[sample2lib[sampleName]].append(dup_rate)
        # assigned reads
        assigned_info = glob.glob("{}/multiQC/multiqc_data/multiqc_featureCounts.txt".format(outputDir))[0]
        assigned_df = pd.read_csv(assigned_info, sep ="\t", usecols=["Sample", "Total", "Assigned"])
        assigned_df = assigned_df.loc[assigned_df["Sample"] == sampleName+".filtered"]
        assigned_rate = assigned_df["Assigned"].values/assigned_df["Total"].values*100
        assigned_rate = assigned_rate[0]
        baseDict[sample2lib[sampleName]].append(assigned_rate)



    # Filter
    outputDict = {k: v for k, v in baseDict.items() if len(v) == 10}
    # Reformat into a matrix
    m = []
    for k, v in outputDict.items():
        m.append({'barcode': k,
                  'reads_pf_sequenced': v[1],
                  'confident_reads': v[2],
                  'optical_duplicates': v[3],
                  'mapped_reads': v[5],
                  'uniq_mapped': v[6],
                  'multi_mapped': v[7],
                  'dupped_reads': v[8],
                  'assigned_reads': v[9]})
    return m


def sendToParkour(config, msg):
    FCID = config.get("Options", "runID").split("_")[3][1:]
    if '-' in FCID:
        FCID = FCID.split('-')[-1]
    d = {'flowcell_id': FCID}
    d['sequences'] = json.dumps(msg)
    log.info(f"sendToParkour: Sending {d} to Parkour")
    res = requests.post(config.get("Parkour", "ResultsURL"), auth=(config.get("Parkour", "user"), config.get("Parkour", "password")), data=d, verify=config.get("Parkour", "cert"))
    log.info(f"sendToParkour return {res}")
    return res



def phoneHome(config, outputDir, pipeline, samples_tuples, organism, project, libType):
    """
    Return metrics to Parkour, the results are in outputDir and pipeline needs to be run on them
    """
    samples_id = [row[0] for row in samples_tuples]
    baseDict, sample2lib = getBaseStatistics(config, outputDir, samples_id, organism)

    log.info("phoneHome: baseDict: {}, sample2lib: {}".format(baseDict, sample2lib))

    msg = None
    if pipeline == 'DNA':
        msg = DNA(config, outputDir, baseDict, sample2lib)
    elif pipeline == 'RNA':
        msg = RNA(config, outputDir, baseDict, sample2lib)
    else:
        m = []
        for k, v in baseDict.items():
            m.append({'barcode': k,
                      'reads_pf_sequenced': v[1],
                      'confident_reads': v[2],
                      'optical_duplicates': v[3]})
        msg = m

    if msg is not None:
        ret = sendToParkour(config, msg)
    else:
        ret = None
    
    return [project, organism, libType, pipeline, 'success', ret]


def telegraphHome(config, group, project, skipList, organism=None):
    """
    The skipList is a list of samples/libraries for which we don't run a pipeline, but it'd be nice to still send back sequencing metrics
    Structure of skipList:
    [library, sampleName, libraryType]
    """
    log.info(f"telegraphHome triggered for {project}")
    # make a fake output directory path
    baseDir = Path(
        config.get('Paths', 'groupData'),
        BRB.misc.pacifier(group),
        BRB.misc.getLatestSeqdir(config.get('Paths','groupData'), group),
        config.get('Options', 'runID'),
        BRB.misc.pacifier(project)
    )
    # Mock path
    outputDir = baseDir / "DNA_mouse"
    samples_id = [row[0] for row in skipList]
    baseDict, sample2lib = getBaseStatistics(config, outputDir, samples_id, organism)
    # Reformat into a matrix
    m = []
    for k, v in baseDict.items():
        m.append({'barcode': k,
                  'reads_pf_sequenced': v[1],
                  'confident_reads': v[2],
                  'optical_duplicates': v[3]})
    ret = sendToParkour(config, m)
    # Format the libtypes
    libTypes = ','.join(set([i[2] for i in skipList]))
    # [project, organism, libtypes, workflow, workflow status, parkour status, sambaUpdate]
    return [project, organism, libTypes, None, None, ret, False]
