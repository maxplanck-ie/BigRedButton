import json
import os
import unicodedata

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


def pacifier(s):
    """
    Pacify a string by removing umlauts such that ö becomes o. Also remove spaces, since they break things

    This only works in python 3
    """
    s = s.replace(" ", "")
    s = s.replace("'", "")
    return str(unicodedata.normalize("NFKD", s).encode("ASCII", "ignore"), "utf-8")
