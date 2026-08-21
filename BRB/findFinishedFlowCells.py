import glob
from pathlib import Path

import requests
from rich import print

from BRB.logger import log


def flowCellProcessed(config):
    return Path(
        config.get("Paths", "baseData"), config.get("Options", "runID"), "analysis.done"
    ).exists()


def markFinished(config):
    _p = Path(config["Paths"]["baseData"], config["Options"]["runID"], "analysis.done")
    _p.touch()
    log.info(f"{_p} created, flow cell processed.")
    print(f"{_p} created, flow cell processed.")


def queryParkour(config):
    basePath = config.get("Paths", "baseData")
    sequencer_type = config.get("Options", "sequencerType")
    if sequencer_type == "Aviti":
        FCID = config.get("Options", "runID").split("_")[2]
        if "-" in FCID:
            FCID = FCID.split("-")[-1]
        d = {"flowcell_id": FCID}
    else:
        FCID = config.get("Options", "runID").split("_")[3][
            1:
        ]  # C605HACXX from 150416_SN7001180_0196_BC605HACXX
        if "-" in FCID:
            FCID = FCID.split("-")[-1]
        d = {"flowcell_id": FCID}
    res = requests.get(
        config.get("Parkour", "QueryURL"),
        auth=(config.get("Parkour", "user"), config.get("Parkour", "password")),
        params=d,
        verify=config.get("Parkour", "cert"),
    )
    if res.status_code == 200:
        return res.json()
    return {}


def detect_sequencer_type(base_path: str) -> str:
    """
    Detect whether a sequencing run is Aviti or Illumina
    based on the presence of the Aviti-specific RunManifest.csv file.
    """
    aviti_check = glob.glob(f"{base_path}/*/RunManifest.csv")
    if aviti_check:
        return "Aviti"
    else:
        return "Illumina"


def newFlowCell(config, sequencer=None):
    platforms = []
    if sequencer in (None, "illumina"):
        platforms.append(
            (
                "illumina",
                config.get("Paths", "baseData_illumina"),
                config.get("Paths", "logPath_illumina"),
            )
        )
    if sequencer in (None, "aviti"):
        platforms.append(
            (
                "aviti",
                config.get("Paths", "baseData_aviti"),
                config.get("Paths", "logPath_aviti"),
            )
        )

    print("Checking for new flowcells...")
    for platform, baseData, logPath in platforms:
        # aviti's baseData is the shared parent of multiple instrument
        # directories (e.g. AV251009, AV261103), so run directories sit one
        # level deeper than under illumina's single-instrument baseData.
        pattern = (
            f"{baseData}/*/*/fastq.made"
            if platform == "aviti"
            else f"{baseData}/*/fastq.made"
        )
        dirs = glob.glob(pattern)
        found = False
        for d in dirs:
            # Get the flow cell ID (e.g., 150416_SN7001180_0196_BC605HACXX)
            run_id = Path(d).parents[0].name
            config.set("Options", "runID", run_id)

            if config.get("Options", "runID")[:4] < "1804":
                continue

            # Detect sequencer type
            base_path = str(Path(d).parents[0])
            seq_type = detect_sequencer_type(base_path)
            config.set("Options", "sequencerType", seq_type)

            if platform == "aviti":
                # Output and logs mirror the serial-ID (e.g. AV251009)
                # subdir that baseData_aviti holds this flowcell under.
                serialID = Path(d).parents[1].name
                config.set("Paths", "baseData", str(Path(baseData, serialID)))
                config.set("Paths", "logPath", str(Path(logPath, serialID)))
            else:
                config.set("Paths", "baseData", baseData)
                config.set("Paths", "logPath", logPath)

            if not flowCellProcessed(config):
                found = True
                print(
                    f"Found new flow cell: [green]{config.get('Options', 'runID')}[/green]"
                )
                # Query parkour to see if there's anything to be done for this
                ParkourDict = queryParkour(config)
                if len(ParkourDict) == 0:
                    markFinished(config)
                    config.set("Options", "runID", "")
                    ParkourDict = None
                    continue
                return config, ParkourDict
        if not found:
            print(f"  No new {platform} flowcells found.")
    return config, None
