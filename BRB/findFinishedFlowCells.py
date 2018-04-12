'''
This file includes anything involved in finding new flow cells to process.

Note that this includes anything done after a flow cell has been processed,
such as marking it as having been processed and sending emails.
'''

import os.path
import sys
import glob
import requests


#Returns True on processed, False on unprocessed
def flowCellProcessed(config):
    path = "{}/{}/analysis.done".format(config.get("Paths","baseData"), config.get("Options","runID"))
    return os.path.exists(path)


def markFinished(config):
    open("{}/{}/analysis.done".format(config["Paths"]["baseData"], config["Options"]["runID"]), "w").close()
    return


def queryParkour(config):
    FCID = config.get("Options", "runID").split("_")[3][1:]  # C605HACXX from 150416_SN7001180_0196_BC605HACXX
    if '-' in FCID:
        FCID = FCID.split('-')[-1]
    d = {'flowcell_id': FCID}
    res = requests.get(config.get("Parkour", "QueryURL"), auth=(config.get("Parkour", "user"), config.get("Parkour", "password")), params=d)
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
            sys.stderr.write("Found a new flow cell: %s\n" % config.get("Options","runID"))
            # Query parkour to see if there's anything to be done for this
            ParkourDict = queryParkour(config)
            if len(ParkourDict) == 0:
                markFinished(config)
                config.set('Options','runID', '')
                ParkourDict = None
                continue
            return config, ParkourDict

    return config, None
