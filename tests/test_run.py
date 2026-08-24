import configparser

import pytest

import BRB.email
import BRB.findFinishedFlowCells
import BRB.getConfig
import BRB.PushButton
import BRB.run


class StopLoop(Exception):
    """Breaks run_brb's `while True` after one flowcell."""


def runConfig(tmp_path):
    config = configparser.ConfigParser()
    config["Paths"] = {
        "baseData": str(tmp_path / "base"),
        "logPath": str(tmp_path / "logs"),
    }
    config["Options"] = {"runID": "run1"}
    return config


def wireRunBrb(monkeypatch, tmp_path, runFlowcellImpl):
    config = runConfig(tmp_path)
    parkourDict = {
        "1_A_Foo": {"L1": ["s1", "ChIP-Seq", "proto", ("human", "hg38", "y"), "i7", 30]}
    }
    monkeypatch.setattr(BRB.getConfig, "getConfig", lambda configfile: config)
    monkeypatch.setattr(
        BRB.findFinishedFlowCells,
        "newFlowCell",
        lambda cfg, sequencer: (config, parkourDict),
    )
    monkeypatch.setattr(BRB.PushButton, "runFlowcell", runFlowcellImpl)
    return config, parkourDict


class TestRunBrbWiring:
    def test_calls_runFlowcell_once_and_emails_result(self, tmp_path, monkeypatch):
        seen = {}

        def fakeRunFlowcell(cfg, parkourDict, **kwargs):
            seen["parkourDict"] = parkourDict
            return [["1_A_Foo", "human", "ChIP-Seq", "DNA", "success", "OK", False, 0]]

        _config, parkourDict = wireRunBrb(monkeypatch, tmp_path, fakeRunFlowcell)
        monkeypatch.setattr(
            BRB.email, "finishedEmail", lambda cfg, msg: seen.setdefault("msg", msg)
        )

        def stopHere(cfg):
            raise StopLoop

        monkeypatch.setattr(BRB.findFinishedFlowCells, "markFinished", stopHere)

        with pytest.raises(StopLoop):
            BRB.run.run_brb.callback(configfile=None, sequencer=None)

        assert seen["parkourDict"] == parkourDict
        assert seen["msg"] == [
            ["1_A_Foo", "human", "ChIP-Seq", "DNA", "success", "OK", False, 0]
        ]

    def test_error_email_names_the_failing_work_items(self, tmp_path, monkeypatch):
        item = BRB.PushButton.WorkItem(
            project="1_A_Foo",
            group="foo",
            pipeline="RNA",
            organism=("human", "hg38", "y"),
            libraryType="stranded mRNA-Seq",
            tuples=[["L1", "s1", "proto", False]],
        )
        err = BRB.PushButton.GroupDispatchError(
            "1 library-group(s) raised during dispatch: 1_A_Foo / RNA / stranded mRNA-Seq",
            [(item, RuntimeError("kaboom"))],
        )

        def fakeRunFlowcell(cfg, parkourDict, **kwargs):
            raise err

        wireRunBrb(monkeypatch, tmp_path, fakeRunFlowcell)
        seen = {}
        monkeypatch.setattr(
            BRB.email,
            "errorEmail",
            lambda cfg, errTuple, msg: seen.setdefault("msg", msg),
        )

        with pytest.raises(BRB.PushButton.GroupDispatchError):
            BRB.run.run_brb.callback(configfile=None, sequencer=None)

        assert "runFlowcell" in seen["msg"]
        assert "1_A_Foo / RNA / stranded mRNA-Seq" in seen["msg"]
