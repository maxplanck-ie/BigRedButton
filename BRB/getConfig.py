import configparser
import os
import sys

from BRB.misc import configGitInfo


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

    if config.has_section("Options") and config.has_option(
        "Options", "validLibraryTypes"
    ):
        print(
            f"Error: config key [Options] validLibraryTypes was renamed to "
            f"validAnalysisTypes (Parkour LibraryType -> AnalysisType rename). "
            f"Update {configFile}."
        )
        sys.exit(1)

    if config.has_section("external") and config.has_option("external", "LibraryTypes"):
        print(
            f"Error: config key [external] LibraryTypes was renamed to "
            f"AnalysisTypes (Parkour LibraryType -> AnalysisType rename). "
            f"Update {configFile}."
        )
        sys.exit(1)

    if "Options" not in config.sections() or not config.has_option(
        "Options", "validAnalysisTypes"
    ):
        print(f"Error: No validAnalysisTypes defined in config {configFile}")
        sys.exit(1)

    gitBin = config.get("software", "git", fallback="git")
    config.set("Options", "configCommit", configGitInfo(configFile, gitBin) or "")

    return config
