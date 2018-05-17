import unicodedata


def pacifier(s):
    """
    Pacify a string by removing umlauts such that ö becomes o. Also remove spaces, since they break things

    This only works in python 3
    """
    s = s.replace(" ", "")
    return str(unicodedata.normalize('NFKD',s).encode('ASCII', 'ignore'), 'utf-8')
