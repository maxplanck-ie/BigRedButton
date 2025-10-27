from pathlib import Path
import pytest
import configparser
import urllib3
from BRB.findFinishedFlowCells import flowCellProcessed, markFinished, queryParkour, detect_sequencer_type
from unittest.mock import patch, Mock

@pytest.fixture(scope='session')
def ifs(tmp_path_factory):
  fp = tmp_path_factory.mktemp("flowcells")
  Path(fp,"fc_fin").mkdir()
  Path(fp, "fc_fin", "analysis.done").touch()
  Path(fp,"fc_unfin").mkdir()
  (fp / "aviti_run" / "flowcellXXX" / "RunManifest.csv").parent.mkdir(parents=True, exist_ok=True)
  (fp / "aviti_run" / "flowcellXXX" / "RunManifest.csv").touch()

  return fp


def create_conf(l = []):
    '''
    sets up a config, where every list (_l) in l gets set in config as :
    config[_l[0]] = {_l[1]: _l[2]}
    
    Additionaly some garbage values for parkour API gets set
    '''
    config = configparser.ConfigParser()
    config['Parkour'] = {
        "QueryURL": "https://parkour-demo.ie-freiburg.mpg.de/nonext_api",
        "user": "jefke",
        "password": "123",
        "cert": ""
    }

    for _l in l:
        config[_l[0]] = {_l[1]: _l[2]}

    if ( config.get("Options", "sequencerType", fallback="Illumina") != "Aviti" ):
        if not config.has_option("Options", "runID"):
            config.set("Options", "runID", "150416_SN7001180_0196_BC605HACXX")
    else:
        if not config.has_option("Options", "runID"):
            config.set("Options", "runID", "20250901_AV25XXX9_250443KMND")
    
    return config

class TestfindFinishedFlowCells:
    def test_flowCellProcessed(self, ifs):
        config = create_conf([
            ['Paths', 'baseData', ifs],
            ['Options', 'runID', 'fc_fin']
        ])
        assert flowCellProcessed(config) == True
        config = create_conf([
            ['Paths', 'baseData', ifs],
            ['Options', 'runID', 'fc_unfin']
        ])
        assert flowCellProcessed(config) == False

    def test_markFinished(self, ifs):
        config = create_conf([
            ['Paths', 'baseData', ifs],
            ['Options', 'runID', 'fc_unfin'],
        ])
        markFinished(config)
        assert Path(ifs, 'fc_unfin', 'analysis.done').exists()

    def test_queryParkour(self):
        urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
        config = create_conf([
        ['Paths', 'baseData', ifs],
        ['Options', 'sequencerType', 'Illumina']])
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

        config = create_conf([
            ['Paths', 'baseData', ifs],
            ['Options', 'sequencerType', 'Illumina']
        ])
        result = queryParkour(config)

        FCID = "C605HACXX"
        mock_get.assert_called_once_with(
            "https://parkour-demo.ie-freiburg.mpg.de/nonext_api",
            auth=("jefke", "123"),
            params={'flowcell_id': FCID},
            verify=""
        )

        assert result == {}

    @patch("BRB.findFinishedFlowCells.requests.get")
    def test_queryParkour_aviti(self, mock_get, ifs):
        urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

        mock_resp = Mock()
        mock_resp.status_code = 200
        mock_resp.json.return_value = {"some": "data"}
        mock_get.return_value = mock_resp

        config = create_conf([
            ['Paths', 'baseData', ifs],
            ['Options', 'sequencerType', 'Aviti']
        ])
        result = queryParkour(config)

        FCID = "250443KMND"
        mock_get.assert_called_once_with(
            "https://parkour-demo.ie-freiburg.mpg.de/nonext_api",
            auth=("jefke", "123"),
            params={'flowcell_id': FCID},
            verify=""
        )
        assert result == {"some": "data"}
