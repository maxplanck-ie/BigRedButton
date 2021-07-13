import requests
import subprocess
import os
import glob
import csv
import json
import BRB.misc
import pandas as pd


def getNReads(d):
    """
    Get the number of reads and % optical dupes from a directory
    """
    if len(glob.glob("{}/*.duplicate.txt".format(d))):
        fname = glob.glob("{}/*.duplicate.txt".format(d))[0]
        s = open(fname).read()
        optDupes, total = s.split()
        return int(total) - int(optDupes), 100. * float(optDupes) / float(total)
    else:
        # For machines in which we don't mark duplicates
        # Just count the number of reads in R1
        CMD1 = ["zcat", glob.glob("{}/*_R1.fastq.gz".format(d))[0]]
        CMD2 = ["wc", "-l"]
        c1 = subprocess.Popen(CMD1, stdout=subprocess.PIPE)
        res = subprocess.check_output(CMD2, stdin=c1.stdout)
        return (int(res)/4), 0


def getOffSpeciesRate(d):
    """
    Get the percentage of off-species reads from a directory
    This is copied from the bcl2fastq pipeline
    """
    fname = glob.glob("{}/*_screen.txt".format(d))[0]
    total = 0
    species=[]
    ohol=[]
    mhol=[]
    i = 0
    maxi = 0
    for line in csv.reader(open(fname, "r"), dialect="excel-tab") :
        if(len(line) == 0) :
            break
        if(line[0].startswith("#")) :
            continue
        if(line[0].startswith("Library")) :
            continue
        if(line[0].startswith("PhiX") or line[0].startswith("Adapters") or line[0].startswith("Vectors") or line[0].startswith("rRNA")):
            continue
        species.append(line[0])
        ohol.append(float(line[5]))

        if(ohol[maxi] < ohol[i]) :
            maxi = i
        i += 1

    off = 0
    for i in range(len(ohol)) :
        if(i != maxi) :
            off += ohol[i]
    return off


def getBaseStatistics(config, outputDir):
    """
    Return a directionary with keys lib names and values:
    (sample name, nReads, off-species rate, % optical dupes)

    Also return the mapping of sample names to library names
    """
    baseDict = dict()  # output dictionary
    s2l = dict()  # sample to library dictionary
    odir, adir = os.path.split(os.path.split(outputDir)[0])
    pdir = "Project_{}".format(adir[9:])
    for d in glob.glob("{}/{}/{}/Sample_*".format(config.get('Paths','baseData'), config.get('Options', 'runID'), pdir)):
        libName = os.path.split(d)[1][7:]
        if len(glob.glob("{}/*_R1.fastq.gz".format(d))) == 0:
            continue  # Skip failed samples
        sampleName = glob.glob("{}/*_R1.fastq.gz".format(d))[0]
        sampleName = os.path.split(sampleName)[1][:-12]
        nReads, optDupes = getNReads(d)
        offRate = getOffSpeciesRate(d)
        baseDict[libName] = [sampleName, nReads, offRate, optDupes]
        s2l[sampleName] = libName
    return baseDict, s2l


def DNA(config, outputDir):
    """
    Parse an output directory to get a dictionary of libraries and their associated values.

    Add % mapped, % dupped, and insert size to baseDict. Filter it for those actually in the output
    """
    baseDict, sample2lib = getBaseStatistics(config, outputDir)

    # % Mapped
    for fname in glob.glob("{}/Bowtie2/*.Bowtie2_summary.txt".format(outputDir)):
        sampleName = os.path.basename(fname).split(".Bowtie2_summary")[0]
        lastLine = open(fname).read().split("\n")[-2]
        mappedPercent = lastLine.split("%")[0]
        try:
            baseDict[sample2lib[sampleName]].append(float(mappedPercent))
        except:
            continue
        # % Duplicated
        dup_info = glob.glob("{}/multiQC/multiqc_data/multiqc_samtools_flagstat.txt".format(outputDir))[0]
        dup_df = pd.read_csv(dup_info, sep ="\t", usecols=["Sample", "total_passed", "duplicates_passed"])
        dup_df = dup_df.loc[dup_df["Sample"] == sampleName+".markdup"]
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
    outputDict = {k: v for k, v in baseDict.items() if len(v) == 7}
    # Reformat into a matrix
    m = []
    for k, v in outputDict.items():
        m.append({'barcode': k,
                  'reads_pf_sequenced': v[1],
                  'confident_reads': v[2],
                  'optical_duplicates': v[3],
                  'dupped_reads': v[5],
                  'mapped_reads': v[4],
                  'insert_size': v[6]})
    return m


def RNA(config, outputDir):
    """
    Parse an output directory to get a dictionary of libraries and their associated values.

    Add % mapped to baseDict. Filter it for those actually in the output
    """
    baseDict, sample2lib = getBaseStatistics(config, outputDir)
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
        dup_df = dup_df.loc[dup_df["Sample"] == sampleName+".markdup"]
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
    outputDict = {k: v for k, v in baseDict.items() if len(v) == 9}
    # Reformat into a matrix
    m = []
    for k, v in outputDict.items():
        m.append({'barcode': k,
                  'reads_pf_sequenced': v[1],
                  'confident_reads': v[2],
                  'optical_duplicates': v[3],
                  'mapped_reads': v[4],
                  'uniq_mapped': v[5],
                  'multi_mapped': v[6],
                  'dupped_reads': v[7],
                  'assigned_reads': v[8]})
    return m


def sendToParkour(config, msg):
    FCID = config.get("Options", "runID").split("_")[3][1:]  # C605HACXX from 150416_SN7001180_0196_BC605HACXX
    if '-' in FCID:
        FCID = FCID.split('-')[-1]
    d = {'flowcell_id': FCID}
    d['sequences'] = json.dumps(msg)
    res = requests.post(config.get("Parkour", "ResultsURL"), auth=(config.get("Parkour", "user"), config.get("Parkour", "password")), data=d)


def phoneHome(config, outputDir, pipeline):
    """
    Return metrics to Parkour, the results are in outputDir and pipeline needs to be run on them
    """
    msg = None
    if pipeline == 'DNA':
        msg = DNA(config, outputDir)
    elif pipeline == 'RNA':
        msg = RNA(config, outputDir)
    if msg is not None:
        sendToParkour(config, msg)


def telegraphHome(config, group, project, skipList):
    """
    The skipList is a list of samples/libraries for which we don't run a pipeline, but it'd be nice to still send back sequencing metrics

    ([library, sampleName])
    """
    # make a fake output directory path
    baseDir = "{}/{}/{}/{}/Analysis_{}".format(config.get('Paths', 'groupData'),
                                                            BRB.misc.pacifier(group),
                                                            BRB.misc.getLatestSeqdir(config.get('Paths','groupData'), group),
                                                            config.get('Options', 'runID'),
                                                            BRB.misc.pacifier(project))
    outputDir = os.path.join(baseDir, "DNA_mouse") # does not matter what it is, it is just a generic name. No pipeline is run on these data
    baseDict, sample2lib = getBaseStatistics(config, outputDir)
    # Reformat into a matrix
    m = []
    for k, v in baseDict.items():
        m.append({'barcode': k,
                  'reads_pf_sequenced': v[1],
                  'confident_reads': v[2],
                  'optical_duplicates': v[3]})
    sendToParkour(config, m)
    return "Sent information on {} libraries from project {} from the {} group to Parkour\n".format(len(skipList), project, group)
