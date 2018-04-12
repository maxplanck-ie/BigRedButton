#!/usr/bin/env python3
import sys
import os
import datetime
import time
import syslog
import BRB.getConfig
import BRB.findFinishedFlowCells
import BRB.PushButton
import importlib
import signal
from threading import Event
import urllib3

# Disable excess warning messages if we disable SSL checks
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
gotHUP = Event()

def breakSleep(signo, _frame):
    gotHUP.set()

def sleep(config):
    gotHUP.wait(timeout=float(config['Options']['sleepTime'])*60*60)
    gotHUP.clear()

signal.signal(signal.SIGHUP, breakSleep)

while True:
    #Reimport to allow reloading a new version
    importlib.reload(BRB.getConfig)
    importlib.reload(BRB.findFinishedFlowCells)

    #Read the config file
    config = BRB.getConfig.getConfig()
    if(config is None):
        #There's no recovering from this!
        sys.exit("Error: couldn't read the config file!")

    #Get the next flow cell to process, or sleep
    config, ParkourDict = BRB.findFinishedFlowCells.newFlowCell(config)
    if(config.get('Options','runID') == '') or ParkourDict is None:
        sleep(config)
        continue

    #Process each group's data
    for k, v in ParkourDict.items():
        BRB.PushButton.GetResults(config, k, v)

    ##Upload to Galaxy
    #try :
    #    message += bcl2fastq_pipeline.galaxy.linkIntoGalaxy(config)
    #except:
    #    syslog.syslog("Got an error while uploading to Galaxy!\n")
    #    bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Got an error while uploading to Galaxy!")
    #    sleep(config)

    #transferTime = datetime.datetime.now()-startTime

    ##Update parkour, errors are non-fatal here
    #try:
    #    message += bcl2fastq_pipeline.misc.jsonParkour(config, message)
    #except:
    #    syslog.syslog("Received an error while updating Parkour!\n")
    #    bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Got an error while updating Parkour!")

    ##Email finished message
    #try :
    #    BigRedButton.email.finishedEmail(config, message, runTime, transferTime)
    #except :
    #    #Unrecoverable error
    #    syslog.syslog("Couldn't send the finished email! Quiting")
    #    sys.exit()

    #Mark the flow cell as having been processed
    BRB.findFinishedFlowCells.markFinished(config)
