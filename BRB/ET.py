import requests
import subprocess
import os
import glob
import csv
import json


def getNReads(d):
    """
    Get the number of reads and % optical dupes from a directory
    """
    fname = glob.glob("{}/*.duplicate.txt".format(d))[0]
    s = open(fname).read()
    optDupes, nonDupes = s.split()
    return int(nonDupes), 100. * float(optDupes) / float(nonDupes)


def getOffSpeciesRate(d):
    """
    Get the percentage of off-species reads from a directory
    This is copied from the bcl2fastq pipeline
    """
    fname = glob.glob("{}/*_R1_screen.txt".format(d))[0]
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


def getBaseStatistics(outputDir):
    """
    Return a directionary with keys lib names and values:
    (sample name, nReads, off-species rate, % optical dupes)

    Also return the mapping of sample names to library names
    """
    baseDict = dict()  # output dictionary
    s2l = dict()  # sample to library dictionary
    odir, adir = os.path.split(os.path.split(outputDir)[0])
    pdir = "Project_{}".format(adir[9:])
    for d in glob.glob("{}/{}/Sample_*".format(odir, pdir)):
        libName = os.path.split(d)[1][7:]
        sampleName = glob.glob("{}/*_R1.fastq.gz".format(d))[0]
        sampleName = os.path.split(sampleName)[1][:-12]
        nReads, optDupes = getNReads(d)
        offRate = getOffSpeciesRate(d)
        baseDict[libName] = [sampleName, nReads, offRate, optDupes]
        s2l[sampleName] = libName
    return baseDict, s2l


def DNA(outputDir):
    """
    Parse an output directory to get a dictionary of libraries and their associated values.

    Add % mapped, % dupped, and insert size to baseDict. Filter it for those actually in the output
    """
    baseDict, sample2lib = getBaseStatistics(outputDir)

    # % Mapped
    for fname in glob.glob("{}/Bowtie2/*.Bowtie2_summary.txt".format(outputDir)):
        try:
            sampleName = os.path.split(fname)[1][:-20]
            lastLine = open(fname).read().split("\n")[-2]
            mappedPercent = lastLine.split("%")[0]
            baseDict[sample2lib[sampleName]].append(float(mappedPercent))
        except:
            pass

    # % Duplicated
    for fname in glob.glob("{}/Picard_qc/MarkDuplicates/*.mark_duplicates_metrics.txt".format(outputDir)):
        sampleName = os.path.split(fname)[1][:-28]
        lastLine = open(fname).read().split("\n")[-4]
        duppedPercent = lastLine.split("\t")[-2]
        baseDict[sample2lib[sampleName]].append(float(duppedPercent))

    # Median insert size
    for fname in glob.glob("{}/Picard_qc/InsertSizeMetrics/*.insert_size_metrics.txt".format(outputDir)):
        sampleName = os.path.split(fname)[1][:-24]
        lastLine = open(fname).read().split("\n")
        medInsertSize = lastLine[7].split("\t")[0]
        baseDict[sample2lib[sampleName]].append(int(medInsertSize))

    # Filter
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
        msg = DNA(outputDir)
    elif pipeline == 'RNA':
        msg = RNA(outputDir)
    if msg is not None:
        sendToParkour(config, msg)
