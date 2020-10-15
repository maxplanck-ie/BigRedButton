import unicodedata
import os

def getLatestSeqdir(groupData, PI):
    seqDirNum = 0
    for dirs in os.listdir(os.path.join(groupData, PI)):
        if 'sequencing_data' in dirs:
            seqDirStrip = dirs.replace('sequencing_data','')
            if seqDirStrip is not '':
                if int(seqDirStrip) > seqDirNum:
                    seqDirNum = int(seqDirStrip)
    if seqDirNum == 0:
        return 'sequencing_data'
    else:
        return 'sequencing_data' + str(seqDirNum)

def pacifier(s):
    """
    Pacify a string by removing umlauts such that ö becomes o. Also remove spaces, since they break things

    This only works in python 3
    """
    s = s.replace(" ", "")
    return str(unicodedata.normalize('NFKD',s).encode('ASCII', 'ignore'), 'utf-8')
