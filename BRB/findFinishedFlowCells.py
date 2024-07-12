import sys
import glob
import requests
from rich import print
from pathlib import Path
from BRB.logger import log


def flowCellProcessed(config):
    return Path(
        config.get("Paths","baseData"),
        config.get("Options","runID"),
        'analysis.done'
    ).exists()


def markFinished(config):
    _p = Path(
        config["Paths"]["baseData"],
        config["Options"]["runID"],
        'analysis.done'
    )
    _p.touch()
    log.info(f"{_p} created, flow cell processed.")
    print(f"{_p} created, flow cell processed.")


def queryParkour(config):
    FCID = config.get("Options", "runID").split("_")[3][1:]  # C605HACXX from 150416_SN7001180_0196_BC605HACXX
    if '-' in FCID:
        FCID = FCID.split('-')[-1]
    d = {'flowcell_id': FCID}
    res = requests.get(config.get("Parkour", "QueryURL"), auth=(config.get("Parkour", "user"), config.get("Parkour", "password")), params=d, verify=config.get("Parkour", "cert"))
    if res.status_code == 200:
        return res.json()
    return dict()


def newFlowCell(config):
    dirs = glob.glob("{}/*/fastq.made".format(config.get("Paths","baseData")))
    for d in dirs :
        #Get the flow cell ID (e.g., 150416_SN7001180_0196_BC605HACXX)
        config.set('Options','runID',d.split("/")[-2])

        if config.get("Options","runID")[:4] < "1804":
            continue

        if not flowCellProcessed(config):
            print(f"Found new flow cell: [green]{config.get("Options","runID")}[/green]")
            # Query parkour to see if there's anything to be done for this
            ParkourDict = queryParkour(config)
            if len(ParkourDict) == 0:
                markFinished(config)
                config.set('Options','runID', '')
                ParkourDict = None
                continue
            return config, ParkourDict

    return config, None
