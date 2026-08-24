#!/usr/bin/env python
import os
import sys
from pathlib import Path
from time import sleep

import rich_click as click
from rich import print

import BRB.email
import BRB.ET
import BRB.findFinishedFlowCells
import BRB.getConfig
import BRB.misc
import BRB.PushButton
from BRB.logger import log, setLog


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "-c",
    "--configfile",
    type=click.Path(exists=True),
    required=False,
    default=os.path.expanduser("~/configs/BigRedButton.ini"),
    help="specify a custom ini file.",
    show_default=True,
)
@click.option(
    "-s",
    "--sequencer",
    type=click.Choice(["illumina", "aviti"]),
    required=False,
    default=None,
    help="Restrict polling to a single platform's baseData/logPath.",
)
def run_brb(configfile, sequencer):
    while True:
        # Read the config file
        config = BRB.getConfig.getConfig(configfile)

        # Get the next flow cell to process, or sleep
        config, ParkourDict = BRB.findFinishedFlowCells.newFlowCell(config, sequencer)
        if (config.get("Options", "runID") == "") or ParkourDict is None:
            print("Going back to sleep for 60 minutes.")
            sleep(60 * 60)
            continue

        # Open log file
        logFile = Path(
            config["Paths"]["logPath"], config.get("Options", "runID") + ".log"
        )
        logFile.parent.mkdir(parents=True, exist_ok=True)
        print(f"Logging into: {logFile}")
        setLog(logFile)

        # Dispatch every project's library-groups through the flowcell-wide
        # thread pool. Projects that don't live on the lanes being processed
        # are skipped inside runFlowcell.
        try:
            msg = BRB.PushButton.runFlowcell(config, ParkourDict)
        except Exception as e:
            errMsg = f"Received an error running PushButton.runFlowcell(): {e}"
            BRB.email.errorEmail(config, sys.exc_info(), errMsg)
            log.critical(errMsg)
            print(errMsg, file=sys.stderr)
            raise

        # Email finished message
        log.info("Create e-mail")
        log.info(msg)
        BRB.email.finishedEmail(config, msg)

        # Mark the flow cell as having been processed
        BRB.findFinishedFlowCells.markFinished(config)
        log.info("=== finished flowcell ===")
