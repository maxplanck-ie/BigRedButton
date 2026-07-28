import configparser
import os
import sys


def getConfig(configFile=None):
    if configFile is None:
        print("Error: configFile is undefined.")
        sys.exit(1)

    if not os.path.exists(configFile):
        print(f"Error: configFile {configFile} does not exists")
        sys.exit(1)

    config = configparser.ConfigParser()
    with open(configFile) as fh:
        config.read_file(fh)

    if "Paths" not in config.sections():
        print(f"Error: No Paths defined in config {configFile}")
        sys.exit(1)

    return config
