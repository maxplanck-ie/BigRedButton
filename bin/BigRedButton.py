#!/usr/bin/env python3
import sys
import os
import datetime
import time
import syslog
import BRB.getConfig
import BRB.findFinishedFlowCells
import BRB.PushButton
import BRB.email
import BRB.ET
import BRB.misc
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
    importlib.reload(BRB.PushButton)
    importlib.reload(BRB.ET)
    importlib.reload(BRB.email)
    importlib.reload(BRB.misc)

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

    #Process each group's data, ignore cases where the project isn't in the lanes being processed
    bdir = "{}/{}".format(config.get('Paths', 'baseData'), config.get('Options', 'runID'))
    msg = '\n'
    for k, v in ParkourDict.items():
        if not os.path.exists("{}/Project_{}".format(bdir, BRB.misc.pacifier(k))):
            continue
        try:
            msg += BRB.PushButton.GetResults(config, k, v)
        except:
            BRB.email.errorEmail(config, sys.exc_info(), "Received an error running PushButton.GetResults() with {} and {}".format(k, v))
            sys.exit("Received an error running PushButton.GetResults() with {} and {}".format(k, v))

    #Email finished message
    try :
        BRB.email.finishedEmail(config, msg)
    except :
        #Unrecoverable error
        sys.exit()

    #Mark the flow cell as having been processed
    BRB.findFinishedFlowCells.markFinished(config)
