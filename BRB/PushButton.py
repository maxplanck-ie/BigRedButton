import os
import shutil
import glob
import subprocess
import BRB.galaxy
import BRB.ET
import BRB.misc

def createPath(config, group, project, organism, libraryType):
    """Ensures that the output path exists, creates it otherwise, and return where it is"""
    baseDir = "{}/{}/sequencing_data/{}/Analysis_{}".format(config.get('Paths', 'groupData'),
                                                            BRB.misc.pacifier(group),
                                                            config.get('Options', 'runID'),
                                                            BRB.misc.pacifier(project))
    os.makedirs(baseDir, exist_ok=True)

    oDir = os.path.join(baseDir, "{}_{}".format(BRB.misc.pacifier(libraryType), organism))
    os.makedirs(oDir, exist_ok=True)
    return oDir


def linkFiles(config, group, project, odir, tuples):
    """Create symlinks in odir to fastq files in {project}. Return 1 if paired-end, 0 otherwise."""
    baseDir = "{}/{}/sequencing_data/{}/Project_{}".format(config.get('Paths', 'groupData'),
                                                           BRB.misc.pacifier(group),
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
    files = glob.glob("{}/*_R?.fastq.gz".format(d))
    for fname in files:
        os.unlink(fname)


def organism2Org(config, organism):
    """Convert a parkour organism name to a Snakemake organism ID"""
    organisms = config.get('Options', 'validOrganisms').split(',')
    orgs = config.get('Options', 'organismNames').split(',')
    for x, y in zip(organisms, orgs):
        if organism == x:
            return y
    raise RuntimeError('An apparently valid organism doesn\'t have a matching snakemake genome ID!')


def tidyUpABit(d):
    """
    If we don't tidy up we'll have a lot of dot files to upload to Galaxy
    """
    try:
        shutil.rmtree(os.path.join(d, 'cluster_logs'))
        os.unlink(os.path.join(d, 'config.yaml'))
        shutil.rmtree(os.path.join(d, '.snakemake'))
        for f in glob.glob(os.path.join(d, '*.log')):
            os.unlink(f)
    
        for d2 in glob.glob(os.path.join(d, 'FASTQ*')):
            shutil.rmtree(d2)
    except:
        pass


def touchDone(outputDir, fname="analysis.done"):
    open(os.path.join(outputDir, fname), "w").close()


def removeDone(outputDir):
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        os.remove(os.path.join(outputDir, "analysis.done"))


def RNA(config, group, project, organism, libraryType, tuples):
    """
    Need to set --libraryType
    """

    print(tuples[0])

    outputDir = createPath(config, group, project, organism, libraryType)
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        return outputDir, 0
    PE = linkFiles(config, group, project, outputDir, tuples)
    org = organism2Org(config, organism)
    CMD = "PATH={}/bin:$PATH".format(os.path.join(config.get('Options', 'snakemakeWorkflowBaseDir')))
    CMD = [CMD, 'RNA-seq', '--DAG', '-j', config.get('Queue', 'parallelProcesses'), '-i', outputDir, '-o', outputDir, '-m', 'alignment,deepTools_qc', org]
    #if org == 'dm6':
    #    CMD.extend(['--star_options', '"--limitBAMsortRAM 60000000000"'])
    if tuples[0][2].startswith("SMART-Seq"):
        # SMART-seq isn't a dUTP-based method!
        CMD.extend(['--libraryType', '0'])
    elif tuples[0][2].startswith("NEBNext Low Input RNA Library"):
        # Unstranded
        CMD.extend(['--libraryType', '0'])
    try:
        subprocess.check_call(' '.join(CMD), shell=True)
    except:
        return outputDir, 1
    removeLinkFiles(outputDir)
    tidyUpABit(outputDir)
    touchDone(outputDir)
    return outputDir, 0


    baseDir = "{}/{}/sequencing_data/{}/Project_{}".format(config.get('Paths', 'groupData'),
                                                           BRB.misc.pacifier(group),
                                                           config.get('Options', 'runID'),
                                                           BRB.misc.pacifier(project))
def RELACS(config, group, project, organism, libraryType, tuples):
    """
    This is a variant of the DNA mapping pipeline that does RELACS demultiplexing in addition

    This must check for the existence of a RELACS sample sheet in the run folder.

    There better not be any duplicate RELACS sample names!
    """
    runID = config.get('Options', 'runID').split("_lanes")[0]

    outputDir = createPath(config, group, project, organism, libraryType)
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        return outputDir, 0

    sampleSheet = "/dont_touch_this/solexa_runs/{}/RELACS_Project_{}.txt".format(runID, BRB.misc.pacifier(project))
    if not os.path.exists(sampleSheet) and not os.path.exists(os.path.join(outputDir, "RELACS_sampleSheet.txt")):
        return None, 1

    baseDir = "{}/{}/sequencing_data/{}/Project_{}".format(config.get('Paths', 'groupData'),
                                                           BRB.misc.pacifier(group),
                                                           config.get('Options', 'runID'),
                                                           BRB.misc.pacifier(project))

    # Link in files
    if not os.path.exists(os.path.join(outputDir, "RELACS_sampleSheet.txt")):
        shutil.copyfile(sampleSheet, os.path.join(outputDir, "RELACS_sampleSheet.txt"))
    unlinkDirs = []
    for d in glob.glob("{}/Sample_*".format(baseDir)):
        bname = os.path.basename(d)
        newName = os.path.join(outputDir, bname)
        unlinkDirs.append(newName)
        if not os.path.exists(newName):
            os.symlink(d, newName)

    # -p 10 is pretty much arbitrary
    CMD = ["demultiplex_relacs.py", "--umiLength", "4", "-p", "10", os.path.join(outputDir, "RELACS_sampleSheet.txt"), os.path.join(outputDir, "RELACS_demultiplexing")]
    try:
        subprocess.check_call(' '.join(CMD), shell=True, cwd=outputDir)
    except:
        return outputDir, 1

    # clean up
    for d in unlinkDirs:
        os.unlink(d)

    # Link in the RELACS demultiplexed files
    for fname in glob.glob(os.path.join(outputDir, "RELACS_demultiplexing", "*", "*.gz")):
        bname = os.path.basename(fname)
        if bname.startswith('unknown'):
            continue
        newName = os.path.join(outputDir, bname)
        if not os.path.exists(newName):
            os.symlink(fname, newName)

    # Back to the normal DNA pipeline
    org = organism2Org(config, organism)
    CMD = "PATH={}/bin:$PATH".format(os.path.join(config.get('Options', 'snakemakeWorkflowBaseDir')))
    CMD = [CMD, 'DNA-mapping', '--DAG', '--trim', '--UMIDedup', '--mapq', '3', '-j', config.get('Queue', 'parallelProcesses'), '-i', outputDir, '-o', outputDir, org]
    try:
        subprocess.check_call(' '.join(CMD), shell=True)
    except:
        return outputDir, 1
    removeLinkFiles(outputDir)
    tidyUpABit(outputDir)
    touchDone(outputDir)
    return outputDir, 0


def DNA(config, group, project, organism, libraryType, tuples):
    """
    Run the DNA mapping pipeline on the samples. Tweals could theoretically be made
    according to the libraryProtocol (tuple[2])

    - Make /data/{group}/sequencing_data/{runID}/Analysis_{project}/{libraryType}_{organism} directory
    - Remove previously linked in files (if any)
    - Link requested fastq files in
    - Run appropriate pipeline
    - Remove previously linked in files
    - Clean up snakemake directory
    """
    print(tuples[0])
    if tuples[0][2].startswith("ChIP RELACS high-throughput"):
        return RELACS(config, group, project, organism, libraryType, tuples)

    outputDir = createPath(config, group, project, organism, libraryType)
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        return outputDir, 0
    PE = linkFiles(config, group, project, outputDir, tuples)
    org = organism2Org(config, organism)
    CMD = "PATH={}/bin:$PATH".format(os.path.join(config.get('Options', 'snakemakeWorkflowBaseDir')))
    CMD = [CMD, 'DNA-mapping', '--DAG', '--trim', '--dedup', '--mapq', '3', '-j', config.get('Queue', 'parallelProcesses'), '-i', outputDir, '-o', outputDir, org]
    try:
        subprocess.check_call(' '.join(CMD), shell=True)
    except:
        return outputDir, 1
    removeLinkFiles(outputDir)
    tidyUpABit(outputDir)
    touchDone(outputDir)
    return outputDir, 0


def WGBS(config, group, project, organism, libraryType, tuples):
    """
    Run the WGBS pipeline

    TODO: set trimming according to the libraryType
    TODO: I don't think we know how to send back metrics yet
    """
    outputDir = createPath(config, group, project, organism, libraryType)
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        return outputDir, 0
    PE = linkFiles(config, group, project, outputDir, tuples)
    org = organism2Org(config, organism)
    CMD = "PATH={}/bin:$PATH".format(os.path.join(config.get('Options', 'snakemakeWorkflowBaseDir')))
    CMD = [CMD, 'WGBS', '--DAG', '--trim', '-j', config.get('Queue', 'parallelProcesses'), '-i', outputDir, '-o', outputDir, org]
    try:
        subprocess.check_call(' '.join(CMD), shell=True)
    except:
        return outputDir, 1
    removeLinkFiles(outputDir)
    tidyUpABit(outputDir)
    touchDone(outputDir)
    return outputDir, 0


def ATAC(config, group, project, organism, libraryType, tuples):
    """
    Run the DNA mapping pipeline and then the default ATAC pipeline
    """
    outputDir = createPath(config, group, project, organism, libraryType)
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        return outputDir, 0
 
    if not os.path.exists(os.path.join(outputDir, "DNA.done")):
        outputDir, rv = DNA(config, group, project, organism, libraryType, tuples)
        if rv != 0:
            return outputDir, rv

        removeDone(outputDir)
        touchDone(outputDir, "DNA.done")
    org = organism2Org(config, organism)
    CMD = "PATH={}/bin:$PATH".format(os.path.join(config.get('Options', 'snakemakeWorkflowBaseDir')))
    CMD = [CMD, 'ATAC-seq', '--DAG', '-d', outputDir, org]
    try:
        subprocess.check_call(' '.join(CMD), shell=True)
    except:
        return outputDir, 1
    tidyUpABit(outputDir)
    touchDone(outputDir)
    return outputDir, 0


def scRNAseq(config, group, project, organism, libraryType, tuples):
    """
    Run one of the scRNAseq pipelines (snakePipes or 10X)

    The protocol is tuples[0][2] and we assume they're all the same...

    We currently just skip unknown protocols and don't mention that!
    """
    outputDir = createPath(config, group, project, organism, libraryType)
    if os.path.exists(os.path.join(outputDir, "analysis.done")):
        return outputDir, 0

    org = organism2Org(config, organism)
    if tuples[0][2] == "10xGenomics for single cell RNA-Seq":
        PE = linkFiles(config, group, project, outputDir, tuples)
        CMD = ['/data/manke/repository/scripts/snakemake_workflows/10X_snakepipe-0.1.0/10X', outputDir, outputDir, org]
        try:
            subprocess.check_call(' '.join(CMD), shell=True)
        except:
            return outputDir, 1
        removeLinkFiles(outputDir)
        tidyUpABit(outputDir)
    elif tuples[0][2] == "Cel-Seq 2 for single cell RNA-Seq":
        PE = linkFiles(config, group, project, outputDir, tuples)
        CMD = "PATH={}/bin:$PATH".format(os.path.join(config.get('Options', 'snakemakeWorkflowBaseDir')))
        CMD = [CMD, 'scRNAseq', '--DAG', '-j', config.get('Queue', 'parallelProcesses'), '--trim',  '--cellBarcodeFile /data/processing/bioinfo-core/celseq_barcodes.384.txt',  '--cellBarcodePattern NNNNNNNXXXXXXX', '-i', outputDir, '-o', outputDir, org]
        try:
            subprocess.check_call(' '.join(CMD), shell=True)
        except:
            return outputDir, 1
        removeLinkFiles(outputDir)
        tidyUpABit(outputDir)
    touchDone(outputDir)
    return outputDir, 0


def GetResults(config, project, libraries):
    """
    Project is something like '352_Grzes_PearceEd' and libraries is a dictionary with libraries as keys:
        {'18L005489': ['FAT_first_A',
                       'Other',
                       '10xGenomics for single cell RNA-Seq',
                       'mouse'],
         '18L005490': ['FAT_first_B',
                       'Other',
                       '10xGenomics for single cell RNA-Seq',
                       'mouse'],

    This doesn't return anything. It's assumed that everything within a single library type can be analysed together.
    """
    group = project.split("_")[-1].split("-")[0].lower()
    dataPath = "{}/{}/sequencing_data/{}/Project_{}".format(config.get('Paths', 'groupData'),
                                                            BRB.misc.pacifier(group),
                                                            config.get('Options', 'runID'),
                                                            BRB.misc.pacifier(project))

    validLibraryTypes = {v: i for i, v in enumerate(config.get('Options', 'validLibraryTypes').split(','))}
    pipelines = config.get('Options', 'pipelines').split(',')
    validOrganisms = config.get('Options', 'validOrganisms').split(',')

    ignore = False
    if not os.path.exists(dataPath):
        ignore = True

    # split by analysis type and organism, since we can only process some types of this
    analysisTypes = dict()
    skipList = []
    for library, v in libraries.items():
        sampleName, libraryType, libraryProtocol, organism = v
        if libraryType in validLibraryTypes and organism in validOrganisms and not ignore:
            idx = validLibraryTypes[libraryType]
            pipeline = pipelines[idx]
            if pipeline not in analysisTypes:
                analysisTypes[pipeline] = dict()
            if organism not in analysisTypes[pipeline]:
                analysisTypes[pipeline][organism] = dict()
            if libraryType not in analysisTypes[pipeline][organism]:
                analysisTypes[pipeline][organism][libraryType] = list()
            analysisTypes[pipeline][organism][libraryType].append([library, sampleName, libraryProtocol])
        else:
            print("Skipping {}/{} with type {} and organism {}".format(project, library, libraryType, organism))
            skipList.append([library, sampleName])

    msg = ""
    if len(skipList):
        msg += BRB.ET.telegraphHome(config, group, BRB.misc.pacifier(project), skipList)

    for pipeline, v in analysisTypes.items():
        for organism, v2 in v.items():
            for libraryType, tuples in v2.items():
                outputDir, rv = globals()[pipeline](config, group, BRB.misc.pacifier(project), organism, libraryType, tuples)
                if rv == 0:
                    BRB.galaxy.linkIntoGalaxy(config, group, BRB.misc.pacifier(project), outputDir)
                    BRB.ET.phoneHome(config, outputDir, pipeline)
                    msg += 'Processed project {} with the {} pipeline. The samples were of type {} from a {}.\n'.format(BRB.misc.pacifier(project), pipeline, libraryType, organism)
                else:
                    msg += "I received an error processing {}_{}_{}_{} for you.\n".format(BRB.misc.pacifier(project), pipeline, libraryType, organism)

    return msg
