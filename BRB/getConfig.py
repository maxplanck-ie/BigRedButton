import configparser
import os

def getConfig() :
    config = configparser.ConfigParser()
    config.read_file(open("{}/BigRedButton.ini".format(os.path.expanduser("~"))))
    if("Paths" in config.sections()) :
        return config
    return None
