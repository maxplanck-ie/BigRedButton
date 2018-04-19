import os
import shutil
import glob
import subprocess
import BRB.galaxy
import BRB.ET

def createPath(config, group, project, organism, libraryType):
    """Ensures that the output path exists, creates it otherwise, and return where it is"""
    baseDir = "{}/{}/sequencing_data/{}/Analysis_{}".format(config.get('Paths', 'groupData'),
                                                            group,
                                                            config.get('Options', 'runID'),
                                                            project)
    os.makedirs(baseDir, exist_ok=True)

    oDir = os.path.join(baseDir, "{}_{}".format(libraryType, organism))
    os.makedirs(oDir, exist_ok=True)
    return oDir


def linkFiles(config, group, project, odir, tuples):
    """Create symlinks in odir to fastq files in {project}. Return 1 if paired-end, 0 otherwise."""
    baseDir = "{}/{}/sequencing_data/{}/Project_{}".format(config.get('Paths', 'groupData'),
                                                           group,
                                                           config.get('Options', 'runID'),
                                                           project)

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
            os.unlink(d2)
    except:
        pass


def RNA(config, group, project, organism, libraryType, tuples):
    """
    Need to set --library_type and maybe --start_options
    """
    outputDir = createPath(config, group, project, organism, libraryType)
    removeLinkFiles(outputDir)
    PE = linkFiles(config, group, project, outputDir, tuples)
    org = organism2Org(config, organism)
    CMD = os.path.join(config.get('Options', 'snakemakeWorkflowBaseDir'), "RNA-seq")
    CMD = [CMD, '-i', outputDir, '-o', outputDir, org]
    rv = subprocess.check_call(CMD, shell=True)
    removeLinkFiles(outputDir)
    tidyUpABit(outputDir)
    return outputDir, rv


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
    outputDir = createPath(config, group, project, organism, libraryType)
    removeLinkFiles(outputDir)
    PE = linkFiles(config, group, project, outputDir, tuples)
    org = organism2Org(config, organism)
    CMD = os.path.join(config.get('Options', 'snakemakeWorkflowBaseDir'), "DNA-mapping")
    CMD = [CMD, '--trim', '--dedup', '--mapq', '3', '-i', outputDir, '-o', outputDir, org]
    rv = subprocess.check_call(' '.join(CMD), shell=True)
    removeLinkFiles(outputDir)
    tidyUpABit(outputDir)
    return outputDir, rv


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
                                                            group,
                                                            config.get('Options', 'runID'),
                                                            project)

    validLibraryTypes = {v: i for i, v in enumerate(config.get('Options', 'validLibraryTypes').split(','))}
    pipelines = config.get('Options', 'pipelines').split(',')
    validOrganisms = config.get('Options', 'validOrganisms').split(',')

    if not os.path.exists(dataPath):
       return

    # split by analysis type and organism, since we can only process some types of this
    analysisTypes = dict()
    for library, v in libraries.items():
        sampleName, libraryType, libraryProtocol, organism = v
        if libraryType in validLibraryTypes and organism in validOrganisms:
            idx = validLibraryTypes[libraryType]
            pipeline = pipelines[idx]
            if pipeline not in analysisTypes:
                analysisTypes[pipeline] = dict()
            if organism not in analysisTypes[pipeline]:
                analysisTypes[pipeline][organism] = dict()
            if libraryType not in analysisTypes[pipeline][organism]:
                analysisTypes[pipeline][organism][libraryType] = list()
            analysisTypes[pipeline][organism][libraryType].append([library, sampleName, libraryProtocol])

    msg = ""
    for pipeline, v in analysisTypes.items():
        for organism, v2 in v.items():
            for libraryType, tuples in v2.items():
                outputDir, rv = globals()[pipeline](config, group, project, organism, libraryType, tuples)
                if rv == 0:
                    BRB.galaxy.linkIntoGalaxy(config, group, project, outputDir)
                    BRB.ET.phoneHome(config, outputDir, pipeline)
                    msg += 'Processed project {} with the {} pipeline. The samples were of type {} from a {}.\n'.format(project, pipeline, libraryType, organism)
                else:
                    msg += "I can't process {}_{}_{}_{} for you. You should panic now.\n".format(project, pipeline, libraryType, organism)
    return msg
