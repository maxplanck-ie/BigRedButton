import json
import os
import re
import subprocess as sp
import sys
import unicodedata
from importlib.metadata import version
from pathlib import Path

import requests

from BRB.logger import log


def resolveDeliverTo(config):
    """
    Query Parkour's internal_pis endpoint for {PI name: deliver_to} overrides,
    used to translate a PI's Parkour name into the periphery directory name
    IT actually uses when the two diverge (dissectBCL has the same mechanism,
    see misc.deliverDirName there).

    Unlike dissectBCL's internal/external FEX routing decision, this is only
    an enhancement on top of an already-working fallback (see GetResults() in
    PushButton.py), so a failed/empty lookup here is not fatal: it's logged
    and an empty map is returned, leaving the existing fallback in charge.
    """
    try:
        response = requests.get(
            config.get("Parkour", "InternalPIsURL"),
            params={
                "organizations": config.get(
                    "Internals", "Organizations", fallback="MPI-IE"
                )
            },
            auth=(config.get("Parkour", "user"), config.get("Parkour", "password")),
            verify=config.get("Parkour", "cert"),
        )
        response.raise_for_status()
        pis = response.json()["pis"]
    except (requests.exceptions.RequestException, ValueError, KeyError) as e:
        log.warning(
            f"resolveDeliverTo: failed to fetch deliver_to overrides from "
            f"Parkour, falling back to the surname heuristic for every PI: {e}"
        )
        return {}
    if not isinstance(pis, dict):
        # Older Parkour versions returned a bare list of names with no
        # deliver_to info at all - nothing to override with.
        return {}
    return {name.lower(): token for name, token in pis.items() if token}


def resolveGroup(config, project):
    """
    Resolve the periphery group (PI) directory name for a project string
    like "4035_Demollin_Cabezas-Wallscheid": the PI's Parkour deliver_to
    override, when config["Internals"]["deliverTo"] (populated by
    resolveDeliverTo()) has one, else a best-effort guess by truncating a
    hyphenated surname at the first hyphen (e.g. "cabezas-wallscheid" ->
    "cabezas"). The fallback is only skipped once an override exists.
    """
    PI = project.split("_")[-1].lower()
    deliverTo = json.loads(config.get("Internals", "deliverTo", fallback="{}"))
    group = deliverTo.get(PI)
    if group is None:
        group = PI.split("-")[0]
    return pacifier(group)


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


def _runGit(args, gitBin, cwd):
    """
    Run a git command, returning its CompletedProcess, or None if the call
    failed for any reason: git binary missing, the directory not being a
    repo, the command erroring, or the call being slow enough to hit the
    timeout. This metadata is best-effort - a slow or flaky filesystem
    (a hung git on shared storage, say) must never take down a pipeline
    that is otherwise fine to run.
    """
    try:
        return sp.run(
            [gitBin, *args],
            cwd=cwd,
            capture_output=True,
            text=True,
            check=True,
            timeout=5,
        )
    except (FileNotFoundError, sp.CalledProcessError, sp.TimeoutExpired):
        return None


def configGitInfo(configfile, gitBin="git"):
    """
    If configfile lives inside a git repo, refuse to run when it has
    uncommitted or untracked changes - otherwise the version/commit
    reported in emails wouldn't match the config that actually ran.
    Returns the config file's latest commit hash, or None when the config
    isn't tracked in a git repo at all (nothing to check), or when git
    couldn't be queried in time (check skipped - see _runGit).
    """
    configDir = Path(configfile).resolve().parent
    if _runGit(["rev-parse", "--is-inside-work-tree"], gitBin, configDir) is None:
        return None
    status = _runGit(
        ["status", "--porcelain", "--", str(configfile)], gitBin, configDir
    )
    if status is None:
        # Unverified, not verified-clean: say so rather than pretend, but
        # keep going - an unreachable git is not a reason to abort the run.
        print(
            f"Warning: could not run git status on {configfile} - skipping "
            "the uncommitted-changes check and omitting the commit from "
            "the report.",
            file=sys.stderr,
        )
        return None
    if status.stdout.strip():
        print(
            f"Error: configfile {configfile} has uncommitted or untracked "
            "changes in its git repo - commit it before running."
        )
        sys.exit(1)
    commit = _runGit(
        ["log", "-1", "--format=%h", "--", str(configfile)], gitBin, configDir
    )
    if commit is None:
        print(
            f"Warning: could not run git log on {configfile} - the commit "
            "will be omitted from the report.",
            file=sys.stderr,
        )
        return None
    return commit.stdout.strip() or None


def pacifier(s):
    """
    Pacify a string by removing umlauts such that ö becomes o. Also remove spaces, since they break things

    This only works in python 3
    """
    s = s.replace(" ", "")
    s = s.replace("'", "")
    return str(unicodedata.normalize("NFKD", s).encode("ASCII", "ignore"), "utf-8")
