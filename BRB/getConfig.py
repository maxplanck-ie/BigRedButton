import configparser
import sys
import os

def getConfig(configFile=None) :
    if configFile is None:
        print('Error: configFile is undefined.')
        sys.exit(1)

    if not os.path.exists(configFile):
        print('Error: configFile {} does not exists'.format(configFile))
        sys.exit(1)

    config = configparser.ConfigParser()
    config.read_file(open(configFile))

    if("Paths" not in config.sections()) :
        print('Error: No Paths defined in config {}'.format(configFile))
        sys.exit(1)

    return config
