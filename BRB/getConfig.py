import configparser
import os

def getConfig() :
    config = configparser.ConfigParser()
    # config.read_file(open("{}/BigRedButton.ini".format(os.path.expanduser("~"))))
    config.read_file(open("{}/test.ini".format(os.path.join("/data/processing1/leily/BigRedButton"))))
    if("Paths" in config.sections()) :
        return config
    return None
