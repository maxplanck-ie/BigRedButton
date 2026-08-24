import configparser
import os

import pytest

import BRB.ET
from BRB import PushButton
from BRB.PushButton import WorkItem

HUMAN = ("human", "hg38", "/yaml/hg38.yaml")
MOUSE = ("mouse", "mm10", "/yaml/mm10.yaml")


def dispatchConfig(tmp_path):
    """
    Config for the dispatch tests. groupData/<group> deliberately does NOT
    exist, so GetResults sets ignore=True and the [external] allowlist keeps
    the libraries in play -- which also makes createPath() root its output
    under baseData, with no PI group directory needed.
    """
    config = configparser.ConfigParser()
    config["Paths"] = {
        "baseData": str(tmp_path / "base"),
        "groupData": str(tmp_path / "group"),
        "seqFacDir": str(tmp_path / "seqfac"),
        "bioinfoCoreDir": str(tmp_path / "bioinfocore"),
    }
    config["Options"] = {
        "runID": "run1",
        "sequencerType": "Aviti",
        "validLibraryTypes": "ChIP-Seq,stranded mRNA-Seq",
        "pipelines": "DNA,RNA",
        "snakemakeWorkflowBaseDir": "/snakepipes",
    }
    config["external"] = {
        "LibraryTypes": "ChIP-Seq,stranded mRNA-Seq",
        "LibraryProtocols": "",
    }
    (tmp_path / "base" / "run1").mkdir(parents=True)
    (tmp_path / "group").mkdir(parents=True)
    (tmp_path / "bioinfocore").mkdir(parents=True)
    return config


class TestGetResultsReturnsWorkItems:
    def test_returns_workitems_and_does_not_dispatch(self, tmp_path, monkeypatch):
        config = dispatchConfig(tmp_path)

        def boom(*args, **kwargs):
            raise AssertionError("GetResults must not dispatch pipelines any more")

        monkeypatch.setattr(PushButton, "DNA", boom)
        monkeypatch.setattr(PushButton, "RNA", boom)
        libraries = {
            "L1": ["s1", "ChIP-Seq", "proto", HUMAN, "i7", 30],
            "L2": ["s2", "stranded mRNA-Seq", "proto", HUMAN, "i7", 30],
        }

        workItems, msg = PushButton.GetResults(config, "1_A_Foo", libraries)

        assert msg == []
        assert len(workItems) == 2
        assert {(w.pipeline, w.libraryType) for w in workItems} == {
            ("DNA", "ChIP-Seq"),
            ("RNA", "stranded mRNA-Seq"),
        }
        for w in workItems:
            assert isinstance(w, WorkItem)
            assert w.project == "1_A_Foo"
            assert w.group == "foo"
            assert w.organism == HUMAN
        dnaItem = next(w for w in workItems if w.pipeline == "DNA")
        assert dnaItem.tuples == [["L1", "s1", "proto", True]]

    def test_skiplist_message_still_produced(self, tmp_path, monkeypatch):
        config = dispatchConfig(tmp_path)
        # group dir exists => ignore=False => the invalid library type lands
        # in skipList rather than external_skipList.
        (tmp_path / "group" / "foo" / "sequencing_data").mkdir(parents=True)

        def fakeTelegraph(config, group, project, skipList, organism=None):
            return ["telegraph", organism, skipList[0][2], None, None, None, False]

        monkeypatch.setattr(BRB.ET, "telegraphHome", fakeTelegraph)
        libraries = {
            "L1": ["s1", "ChIP-Seq", "proto", HUMAN, "i7", 30],
            "L2": ["s2", "Other", "proto", HUMAN, "i7", 30],
        }

        workItems, msg = PushButton.GetResults(config, "1_A_Foo", libraries)

        assert len(workItems) == 1
        assert workItems[0].pipeline == "DNA"
        assert workItems[0].tuples == [["L1", "s1", "proto", False]]
        assert msg == [["telegraph", "human", "Other", None, None, None, False]]

    def test_external_skiplist_message_still_produced(self, tmp_path):
        config = dispatchConfig(tmp_path)
        config["external"]["LibraryTypes"] = "ChIP-Seq"
        libraries = {
            "L1": ["s1", "ChIP-Seq", "proto", HUMAN, "i7", 30],
            "L2": ["s2", "stranded mRNA-Seq", "proto", HUMAN, "i7", 30],
        }

        workItems, msg = PushButton.GetResults(config, "1_A_Foo", libraries)

        assert [w.libraryType for w in workItems] == ["ChIP-Seq"]
        assert msg == [
            ["1_A_Foo", "human", "stranded mRNA-Seq", None, None, None, False, None]
        ]


def makeWorkItem(**kwargs):
    base = {
        "project": "1_A_Foo",
        "group": "foo",
        "pipeline": "RNA",
        "organism": HUMAN,
        "libraryType": "stranded mRNA-Seq",
        "tuples": [["L1", "s1", "proto", True]],
    }
    base.update(kwargs)
    return WorkItem(**base)


def fakePhoneHome(config, outputDir, pipeline, tuples, organism, project, libType):
    return [project, organism, libType, pipeline, "success", "PARKOUR_OK"]


class TestRunOneGroup:
    def test_success_first_try(self, tmp_path, monkeypatch):
        config = dispatchConfig(tmp_path)
        calls = []

        def stub(config, group, project, organism, libraryType, tuples):
            calls.append((group, project, libraryType))
            return "/out", 0, False

        monkeypatch.setattr(PushButton, "RNA", stub)
        monkeypatch.setattr(BRB.ET, "phoneHome", fakePhoneHome)

        result = PushButton.runOneGroup(config, makeWorkItem())

        assert len(calls) == 1
        assert calls[0] == ("foo", "1_A_Foo", "stranded mRNA-Seq")
        assert result == [
            [
                "1_A_Foo",
                "human",
                "stranded mRNA-Seq",
                "RNA",
                "success",
                "PARKOUR_OK",
                False,
                0,
            ]
        ]

    def test_retries_once_then_succeeds(self, tmp_path, monkeypatch):
        config = dispatchConfig(tmp_path)
        calls = []

        def stub(config, group, project, organism, libraryType, tuples):
            calls.append(1)
            return ("/out", 1, False) if len(calls) == 1 else ("/out", 0, True)

        monkeypatch.setattr(PushButton, "RNA", stub)
        monkeypatch.setattr(BRB.ET, "phoneHome", fakePhoneHome)

        result = PushButton.runOneGroup(config, makeWorkItem())

        assert len(calls) == 2
        assert result[0][-2:] == [True, 1]
        assert result[0][4] == "success"

    def test_fails_twice_gives_FAILED_entry(self, tmp_path, monkeypatch):
        config = dispatchConfig(tmp_path)
        calls = []

        def stub(config, group, project, organism, libraryType, tuples):
            calls.append(1)
            return "/out", 1, False

        monkeypatch.setattr(PushButton, "RNA", stub)

        def mustNotBeCalled(*args, **kwargs):
            raise AssertionError("phoneHome must not run for a failed group")

        monkeypatch.setattr(BRB.ET, "phoneHome", mustNotBeCalled)

        result = PushButton.runOneGroup(config, makeWorkItem())

        assert len(calls) == 2
        assert result == [
            [
                "1_A_Foo",
                "human",
                "stranded mRNA-Seq",
                "RNA",
                "FAILED",
                "not updated",
                False,
                1,
            ]
        ]


def findDeadPid():
    """A PID that is guaranteed not to name a live process right now."""
    pid = 4000000
    while pid > 2:
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return pid
        except (PermissionError, OverflowError, ValueError):
            pass
        pid -= 1
    raise RuntimeError("could not find a dead PID")


class TestOwnershipMarker:
    def _outputDir(self, config, item):
        return PushButton.createPath(
            config,
            item.group,
            item.project,
            item.organism[1],
            item.libraryType,
            item.tuples,
        )

    def test_no_marker_dispatches_and_cleans_up(self, tmp_path, monkeypatch):
        config = dispatchConfig(tmp_path)
        item = makeWorkItem()
        calls = []

        def stub(config, group, project, organism, libraryType, tuples):
            calls.append(1)
            return "/out", 0, False

        monkeypatch.setattr(PushButton, "RNA", stub)
        monkeypatch.setattr(BRB.ET, "phoneHome", fakePhoneHome)

        result = PushButton.runOneGroup(config, item)

        assert len(calls) == 1
        assert len(result) == 1
        marker = os.path.join(self._outputDir(config, item), PushButton.MARKER_NAME)
        assert not os.path.exists(marker)

    def test_marker_written_while_dispatching(self, tmp_path, monkeypatch):
        config = dispatchConfig(tmp_path)
        item = makeWorkItem()
        seen = {}

        def stub(config, group, project, organism, libraryType, tuples):
            marker = os.path.join(self._outputDir(config, item), PushButton.MARKER_NAME)
            with open(marker) as fh:
                seen["content"] = fh.read()
            return "/out", 0, False

        monkeypatch.setattr(PushButton, "RNA", stub)
        monkeypatch.setattr(BRB.ET, "phoneHome", fakePhoneHome)

        PushButton.runOneGroup(config, item)

        assert seen["content"].split()[0] == str(os.getpid())

    def test_live_pid_marker_skips_group(self, tmp_path, monkeypatch):
        config = dispatchConfig(tmp_path)
        item = makeWorkItem()
        marker = os.path.join(self._outputDir(config, item), PushButton.MARKER_NAME)
        with open(marker, "w") as fh:
            fh.write(f"{os.getpid()} 1234.5\n")

        def mustNotBeCalled(*args, **kwargs):
            raise AssertionError("a live-PID-owned group must not be dispatched")

        monkeypatch.setattr(PushButton, "RNA", mustNotBeCalled)

        result = PushButton.runOneGroup(config, item)

        assert result == []
        assert os.path.exists(marker), "another owner's marker must not be removed"

    def test_dead_pid_marker_is_abandoned_and_dispatch_proceeds(
        self, tmp_path, monkeypatch
    ):
        config = dispatchConfig(tmp_path)
        item = makeWorkItem()
        marker = os.path.join(self._outputDir(config, item), PushButton.MARKER_NAME)
        with open(marker, "w") as fh:
            fh.write(f"{findDeadPid()} 1234.5\n")
        calls = []

        def stub(config, group, project, organism, libraryType, tuples):
            calls.append(1)
            return "/out", 0, False

        monkeypatch.setattr(PushButton, "RNA", stub)
        monkeypatch.setattr(BRB.ET, "phoneHome", fakePhoneHome)

        result = PushButton.runOneGroup(config, item)

        assert len(calls) == 1
        assert len(result) == 1
        assert not os.path.exists(marker)

    def test_marker_removed_when_pipeline_raises(self, tmp_path, monkeypatch):
        config = dispatchConfig(tmp_path)
        item = makeWorkItem()

        def stub(config, group, project, organism, libraryType, tuples):
            raise RuntimeError("boom")

        monkeypatch.setattr(PushButton, "RNA", stub)

        with pytest.raises(RuntimeError):
            PushButton.runOneGroup(config, item)

        marker = os.path.join(self._outputDir(config, item), PushButton.MARKER_NAME)
        assert not os.path.exists(marker)
