import configparser

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
