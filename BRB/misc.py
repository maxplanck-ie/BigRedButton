import os
import re
import subprocess as sp
import sys
import unicodedata
from importlib.metadata import version
from pathlib import Path


def loadUserDictionary():
    d = {}
    with open("/home/pipegrp/parkourUsers.txt") as f:
        for line in f:
            cols = line.split("\t")
            d[cols[1]] = [cols[0], cols[2]]
    return d


def getLatestSeqdir(groupData, PI):
    seqDirNum = 0
    for dirs in os.listdir(os.path.join(groupData, PI)):
        if "sequencing_data" in dirs:
            seqDirStrip = dirs.replace("sequencing_data", "")
            if seqDirStrip != "":
                seqDirNum = max(seqDirNum, int(seqDirStrip))
    if seqDirNum == 0:
        return "sequencing_data"
    else:
        return "sequencing_data" + str(seqDirNum)


_GIT_DESCRIBE_RE = re.compile(
    r"^(?P<tag>.+)-(?P<n>\d+)-g(?P<hash>[0-9a-f]+)(?P<dirty>-dirty)?$"
)


def _formatGitDescribe(describeOut):
    """
    Reformat git describe --long output from 'TAG-N-gHASH[-dirty]' to
    'TAG +N:gHASH[-dirty]', which reads less ambiguously than three
    dash-separated fields. Left as-is when there's no tag to split off
    (e.g. the repo has never been tagged, so --always fell back to a bare
    hash).
    """
    m = _GIT_DESCRIBE_RE.match(describeOut)
    if not m:
        return describeOut
    dirty = m.group("dirty") or ""
    return f"{m.group('tag')} +{m.group('n')}:g{m.group('hash')}{dirty}"


def getVersion(distName, gitBin="git"):
    """
    Live version string from the checked-out git repo (tag-count-hash,
    '-dirty' if uncommitted changes), so it reflects the branch actually
    running rather than whatever setuptools_scm baked into the installed
    metadata at install time. Falls back to the installed package metadata
    when not run from a git checkout, or when gitBin can't be run (e.g. git
    isn't on PATH in the conda env - see [software] git in the config).
    """
    try:
        out = sp.run(
            [gitBin, "describe", "--tags", "--long", "--dirty", "--always"],
            cwd=Path(__file__).resolve().parent,
            capture_output=True,
            text=True,
            check=True,
            timeout=5,
        )
        return _formatGitDescribe(out.stdout.strip())
    except (FileNotFoundError, sp.CalledProcessError, sp.TimeoutExpired):
        return version(distName)


def configGitInfo(configfile, gitBin="git"):
    """
    If configfile lives inside a git repo, refuse to run when it has
    uncommitted or untracked changes - otherwise the version/commit
    reported in emails wouldn't match the config that actually ran.
    Returns the config file's latest commit hash, or None when the config
    isn't tracked in a git repo at all (nothing to check).
    """
    configDir = Path(configfile).resolve().parent
    try:
        sp.run(
            [gitBin, "rev-parse", "--is-inside-work-tree"],
            cwd=configDir,
            capture_output=True,
            text=True,
            check=True,
            timeout=5,
        )
    except (FileNotFoundError, sp.CalledProcessError, sp.TimeoutExpired):
        return None
    status = sp.run(
        [gitBin, "status", "--porcelain", "--", str(configfile)],
        cwd=configDir,
        capture_output=True,
        text=True,
        check=True,
        timeout=5,
    )
    if status.stdout.strip():
        print(
            f"Error: configfile {configfile} has uncommitted or untracked "
            "changes in its git repo - commit it before running."
        )
        sys.exit(1)
    commit = sp.run(
        [gitBin, "log", "-1", "--format=%h", "--", str(configfile)],
        cwd=configDir,
        capture_output=True,
        text=True,
        check=True,
        timeout=5,
    )
    return commit.stdout.strip() or None


def pacifier(s):
    """
    Pacify a string by removing umlauts such that ö becomes o. Also remove spaces, since they break things

    This only works in python 3
    """
    s = s.replace(" ", "")
    s = s.replace("'", "")
    return str(unicodedata.normalize("NFKD", s).encode("ASCII", "ignore"), "utf-8")
