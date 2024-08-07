#!/usr/bin/env python
import sys
import os
import BRB.getConfig
import BRB.findFinishedFlowCells
import BRB.PushButton
import BRB.email
import BRB.ET
import BRB.misc
from BRB.logger import setLog, log
import rich_click as click
from time import sleep
from pathlib import Path
from rich import print

@click.command(
    context_settings=dict(
        help_option_names=["-h", "--help"]
    )
)
@click.option(
   "-c",
   "--configfile",
   type=click.Path(exists=True),
   required=False,
   default=os.path.expanduser('~/configs/BigRedButton.ini'),
   help='specify a custom ini file.',
   show_default=True
)
def run_brb(configfile):

    while True:
        #Read the config file
        config = BRB.getConfig.getConfig(configfile)

        #Get the next flow cell to process, or sleep
        config, ParkourDict = BRB.findFinishedFlowCells.newFlowCell(config)
        if(config.get('Options','runID') == '') or ParkourDict is None:
            sleep(60*60)
            continue

        # Open log file
        logFile = Path(
            config['Paths']['logPath'],
            config.get('Options','runID') + '.log'
        )
        print(f"Logging into: {logFile}")
        setLog(logFile)

        #Process each group's data, ignore cases where the project isn't in the lanes being processed
        bdir = "{}/{}".format(config.get('Paths', 'baseData'), config.get('Options', 'runID'))
        msg = []
        for k, v in ParkourDict.items():
            if not os.path.exists("{}/Project_{}".format(bdir, BRB.misc.pacifier(k))):
                log.info("{}/Project_{} doesn't exist, probably lives on another lane.".format(bdir, BRB.misc.pacifier(k)))
                continue
            try:
                msg = msg + BRB.PushButton.GetResults(config, k, v)
            except Exception as e:
                BRB.email.errorEmail(config, sys.exc_info(), "Received an error running PushButton.GetResults() with {} and {}".format(k, v))
                log.critical("Received an error running PushButton.GetResults() with {} and {}".format(k, v))
                print("Received an error running PushButton.GetResults() with {} and {}".format(k, v), file=sys.stderr)
                raise

        #Email finished message
        log.info('Create e-mail')
        log.info(msg)
        BRB.email.finishedEmail(config, msg)

        #Mark the flow cell as having been processed
        BRB.findFinishedFlowCells.markFinished(config)
        log.info('=== finished flowcell ===')
