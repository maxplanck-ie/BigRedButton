import configparser
from pathlib import Path
from unittest.mock import Mock, patch

import pytest
import urllib3

from BRB.findFinishedFlowCells import (
    detect_sequencer_type,
    flowCellProcessed,
    markFinished,
    newFlowCell,
    queryParkour,
)


@pytest.fixture(scope="session")
def ifs(tmp_path_factory):
    fp = tmp_path_factory.mktemp("flowcells")
    Path(fp, "fc_fin").mkdir()
    Path(fp, "fc_fin", "analysis.done").touch()
    Path(fp, "fc_unfin").mkdir()
    (fp / "aviti_run" / "flowcellXXX" / "RunManifest.csv").parent.mkdir(
        parents=True, exist_ok=True
    )
    (fp / "aviti_run" / "flowcellXXX" / "RunManifest.csv").touch()

    return fp


def create_conf(l=None):
    """
    sets up a config, where every list (_l) in l gets set in config as :
    config[_l[0]] = {_l[1]: _l[2]}

    Additionally some garbage values for parkour API gets set
    """
    if l is None:
        l = []
    config = configparser.ConfigParser()
    config["Parkour"] = {
        "QueryURL": "https://parkour-demo.ie-freiburg.mpg.de/nonext_api",
        "user": "jefke",
        "password": "123",
        "cert": "",
    }

    for _l in l:
        config[_l[0]] = {_l[1]: _l[2]}

    if config.get("Options", "sequencerType", fallback="Illumina") != "Aviti":
        if not config.has_option("Options", "runID"):
            config.set("Options", "runID", "150416_SN7001180_0196_BC605HACXX")
    else:
        if not config.has_option("Options", "runID"):
            config.set("Options", "runID", "20250901_AV25XXX9_250443KMND")

    return config


class TestfindFinishedFlowCells:
    def test_flowCellProcessed(self, ifs):
        config = create_conf(
            [["Paths", "baseData", ifs], ["Options", "runID", "fc_fin"]]
        )
        assert flowCellProcessed(config) == True
        config = create_conf(
            [["Paths", "baseData", ifs], ["Options", "runID", "fc_unfin"]]
        )
        assert flowCellProcessed(config) == False

    def test_markFinished(self, ifs):
        config = create_conf(
            [
                ["Paths", "baseData", ifs],
                ["Options", "runID", "fc_unfin"],
            ]
        )
        markFinished(config)
        assert Path(ifs, "fc_unfin", "analysis.done").exists()

    def test_queryParkour(self):
        urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
        config = create_conf(
            [["Paths", "baseData", ifs], ["Options", "sequencerType", "Illumina"]]
        )
        assert len(queryParkour(config)) == 0

    ## Sequencing type detection
    def test_typedetection_aviti(self, ifs):
        basedir_path = Path(ifs, "aviti_run")
        seq_type = detect_sequencer_type(str(basedir_path))
        assert seq_type == "Aviti"

    def test_typedetection_illumina(self, ifs):
        basedir_path = Path(ifs, "fc_fin")
        seq_type = detect_sequencer_type(str(basedir_path))
        assert seq_type == "Illumina"

    ## query Parkour flow cell id extraction test
    @patch("BRB.findFinishedFlowCells.requests.get")
    def test_queryParkour_illumina(self, mock_get, ifs):
        urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

        mock_resp = Mock()
        mock_resp.status_code = 200
        mock_resp.json.return_value = {}
        mock_get.return_value = mock_resp

        config = create_conf(
            [["Paths", "baseData", ifs], ["Options", "sequencerType", "Illumina"]]
        )
        result = queryParkour(config)

        FCID = "C605HACXX"
        mock_get.assert_called_once_with(
            "https://parkour-demo.ie-freiburg.mpg.de/nonext_api",
            auth=("jefke", "123"),
            params={"flowcell_id": FCID},
            verify="",
        )

        assert result == {}

    @patch("BRB.findFinishedFlowCells.requests.get")
    def test_queryParkour_aviti(self, mock_get, ifs):
        urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

        mock_resp = Mock()
        mock_resp.status_code = 200
        mock_resp.json.return_value = {"some": "data"}
        mock_get.return_value = mock_resp

        config = create_conf(
            [["Paths", "baseData", ifs], ["Options", "sequencerType", "Aviti"]]
        )
        result = queryParkour(config)

        FCID = "250443KMND"
        mock_get.assert_called_once_with(
            "https://parkour-demo.ie-freiburg.mpg.de/nonext_api",
            auth=("jefke", "123"),
            params={"flowcell_id": FCID},
            verify="",
        )
        assert result == {"some": "data"}


@pytest.fixture(scope="session")
def platform_dirs(tmp_path_factory):
    fp = tmp_path_factory.mktemp("platform_dirs")
    illu_base = fp / "illumina"
    aviti_base = fp / "aviti"
    (illu_base / "20250101_illumina_runXXX").mkdir(parents=True)
    (illu_base / "20250101_illumina_runXXX" / "fastq.made").touch()
    # aviti_base is the shared parent of multiple instruments, so run
    # directories sit one level deeper than under illumina's baseData.
    (aviti_base / "AV999999" / "20250101_AV999999_runYYY").mkdir(parents=True)
    (aviti_base / "AV999999" / "20250101_AV999999_runYYY" / "fastq.made").touch()
    return illu_base, aviti_base


def create_platform_conf(illu_base, aviti_base):
    config = configparser.ConfigParser()
    config["Parkour"] = {
        "QueryURL": "https://parkour-demo.ie-freiburg.mpg.de/nonext_api",
        "user": "jefke",
        "password": "123",
        "cert": "",
    }
    config["Paths"] = {
        "baseData_illumina": str(illu_base),
        "baseData_aviti": str(aviti_base),
        "logPath_illumina": str(illu_base / "LOG"),
        "logPath_aviti": str(aviti_base / "LOG"),
    }
    config["Options"] = {}
    return config


class TestNewFlowCellSequencerGating:
    @patch("BRB.findFinishedFlowCells.queryParkour")
    def test_illumina_only_never_touches_aviti(self, mock_query, platform_dirs):
        illu_base, aviti_base = platform_dirs
        mock_query.return_value = {"some": "data"}
        config = create_platform_conf(illu_base, aviti_base)

        config, _ = newFlowCell(config, sequencer="illumina")

        assert config.get("Options", "runID") == "20250101_illumina_runXXX"
        assert config.get("Paths", "baseData") == str(illu_base)
        assert config.get("Paths", "logPath") == str(illu_base / "LOG")
        assert not Path(
            aviti_base, "AV999999", "20250101_AV999999_runYYY", "analysis.done"
        ).exists()

    @patch("BRB.findFinishedFlowCells.queryParkour")
    def test_aviti_only_never_touches_illumina(self, mock_query, platform_dirs):
        illu_base, aviti_base = platform_dirs
        mock_query.return_value = {"some": "data"}
        config = create_platform_conf(illu_base, aviti_base)

        config, _ = newFlowCell(config, sequencer="aviti")

        assert config.get("Options", "runID") == "20250101_AV999999_runYYY"
        assert config.get("Paths", "baseData") == str(aviti_base / "AV999999")
        assert config.get("Paths", "logPath") == str(aviti_base / "LOG" / "AV999999")
        assert not Path(illu_base, "20250101_illumina_runXXX", "analysis.done").exists()

    @patch("BRB.findFinishedFlowCells.queryParkour")
    def test_aviti_finds_runs_across_multiple_instruments(self, mock_query, tmp_path):
        """Regression test: baseData_aviti is the shared parent of several
        instrument directories (e.g. AV251009, AV261103), so a run's
        fastq.made sits two levels below baseData_aviti, not one."""
        illu_base = tmp_path / "illumina"
        aviti_base = tmp_path / "aviti"
        illu_base.mkdir()
        (aviti_base / "AV251009" / "20260818_AV251009_runZZZ").mkdir(parents=True)
        (aviti_base / "AV251009" / "20260818_AV251009_runZZZ" / "fastq.made").touch()

        mock_query.return_value = {"some": "data"}
        config = create_platform_conf(illu_base, aviti_base)

        config, ParkourDict = newFlowCell(config, sequencer="aviti")

        assert config.get("Options", "runID") == "20260818_AV251009_runZZZ"
        assert ParkourDict == {"some": "data"}

    @patch("BRB.findFinishedFlowCells.queryParkour")
    def test_aviti_baseData_and_logPath_nest_under_serial_id(
        self, mock_query, tmp_path
    ):
        """Regression test: baseData_aviti/logPath_aviti are the shared
        parent of multiple instruments. flowCellProcessed/markFinished and
        run.py's logfile path all build off config["Paths"]["baseData"]/
        ["logPath"], so those must already be nested under the matching
        serial-ID (e.g. AV251009), or analysis.done and the logfile get
        written to a flat path whose parent doesn't exist."""
        illu_base = tmp_path / "illumina"
        aviti_base = tmp_path / "aviti"
        aviti_log = tmp_path / "aviti_log"
        illu_base.mkdir()
        runID = "20260818_AV251009_runZZZ"
        (aviti_base / "AV251009" / runID).mkdir(parents=True)
        (aviti_base / "AV251009" / runID / "fastq.made").touch()

        mock_query.return_value = {"some": "data"}
        config = create_platform_conf(illu_base, aviti_base)
        config.set("Paths", "logPath_aviti", str(aviti_log))

        config, _ = newFlowCell(config, sequencer="aviti")

        assert config.get("Paths", "baseData") == str(aviti_base / "AV251009")
        assert config.get("Paths", "logPath") == str(aviti_log / "AV251009")
        # flowCellProcessed must resolve to the real, nested run directory.
        assert Path(config.get("Paths", "baseData"), runID) == (
            aviti_base / "AV251009" / runID
        )
