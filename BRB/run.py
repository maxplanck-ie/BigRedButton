#!/usr/bin/env python
import glob
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


def process_data(config, ParkourDict):
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
    
    # Email finished message
    log.info('Create e-mail')
    log.info(msg)
    BRB.email.finishedEmail(config, msg)
    
    return


def validate_fcid_with_stats(ctx, param, value):
    if ctx.params.get('stats') and not value:
        raise click.UsageError('--fcid is required when --stats standalone run is active.')
    return value


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
   help='Specify a custom ini file.',
   show_default=True
)
@click.option(
    "-s",
    "--stats",
    required=False,
    is_flag=True,
    help='Standalone run, will not run any pipelines. Requires --fcid to indicate target.'
)
@click.option('--fcid', callback=validate_fcid_with_stats, help='Flowcell ID to push stats.')
def run_brb(configfile, stats, fcid):

    while True:
        # Read the config file
        config = BRB.getConfig.getConfig(configfile)
        
        if not stats:
            # Get the next flow cell to process, or sleep
            config, ParkourDict = BRB.findFinishedFlowCells.newFlowCell(config)
            if (config.get('Options','runID') == '') or ParkourDict is None:
                sleep(60*60)
                continue

            # Open log file
            logFile = Path(
                config['Paths']['logPath'],
                config.get('Options','runID') + '.log'
            )
            print(f"Logging into: {logFile}")
            setLog(logFile)

        else:
            # Push stats on-demand
            log.info(f"Pushing stats for flowcell: {fcid}")
            d = [d for d in glob.glob("{}/*/fastq.made".format(config.get('Paths', 'baseData'))) if fcid in d]
            dual_lane = len(d) == 2
            if len(d) == 0:
                log.error(f"No fastq.made files found for {fcid}")
                return  # Exit BRB if no files found.
            elif len(d) > 2:
                log.error(f"How many lanes does {fcid} have?!")
                return  # Exit BRB this error shouldn't happen at all.
            
            config.set('Options','runID',d[0].split("/")[-2])
            ParkourDict = BRB.findFinishedFlowCells.queryParkour(config)
            
            if dual_lane:
                config1 = BRB.getConfig.getConfig(configfile)
                config1.set('Options','runID',d[1].split("/")[-2])
                ParkourDict1 = BRB.findFinishedFlowCells.queryParkour(config)

            # Open log file
            if not dual_lane:
                logFile = Path(config['Paths']['logPath'], config.get('Options','runID') + '.log')
            else:
                logFile = Path(config['Paths']['logPath'], config.get('Options','runID') + '_2' + '.log')
            print(f"Logging into: {logFile}")
            setLog(logFile)

            if dual_lane:
                log.info("Same log-file is being used for both lanes. Hopefully this is not too confusing :$")

        
        # Process each group's data, ignore cases where the project isn't in the lanes being processed
        process_data(config, ParkourDict)
        
        if stats and dual_lane:
            process_data(config1, ParkourDict1)

        
        if not stats:
            # Mark the flow cell as having been processed
            BRB.findFinishedFlowCells.markFinished(config)
            log.info('=== finished flowcell ===')


        if stats:
            return  # don't do anything else.
