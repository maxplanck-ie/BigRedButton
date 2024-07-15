from pathlib import Path
import pytest
import configparser
import urllib3
from BRB.findFinishedFlowCells import flowCellProcessed, markFinished, queryParkour


@pytest.fixture(scope='session')
def ifs(tmp_path_factory):
  fp = tmp_path_factory.mktemp("flowcells")
  Path(fp,"fc_fin").mkdir()
  Path(fp, "fc_fin", "analysis.done").touch()
  Path(fp,"fc_unfin").mkdir()
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
    config['Options'] = {'runID': '150416_SN7001180_0196_BC605HACXX'}
    for _l in l:
        config[_l[0]] = {_l[1]: _l[2]}
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
        config = create_conf()
        assert len(queryParkour(config)) == 0