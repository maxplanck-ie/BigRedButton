import configparser
import glob as glob_mod
import os

import pytest

from BRB import PushButton

HUMAN = ("Human", "h38", "/yaml/hg38.yaml")


def make_config(overrides=None):
    config = configparser.ConfigParser()
    config["Paths"] = {
        "baseData": "/base",
        "groupData": "/data",
    }
    config["Options"] = {
        "runID": "20260101_AV999999_1234567_lanes_1_2",
        "snakemakeWorkflowBaseDir": "/snakepipes",
    }
    config["10x"] = {"RNA": "/tenx/rna_wf", "ATAC": "/tenx/atac_wf"}
    for section, kv in (overrides or {}).items():
        if section not in config:
            config[section] = {}
        config[section].update(kv)
    return config


def stubHelpers(monkeypatch, outputDir, subprocessRaises=False):
    """
    Stub every side-effecting collaborator each pipeline function calls, so
    only the function's own branching/CMD-building logic is under test.
    """
    calls = {"runManagedSubprocess": [], "touchDone": [], "copyCellRanger": []}

    monkeypatch.setattr(PushButton, "createPath", lambda *a, **k: str(outputDir))
    monkeypatch.setattr(PushButton, "linkFiles", lambda *a, **k: 1)
    monkeypatch.setattr(PushButton, "removeLinkFiles", lambda *a, **k: None)
    monkeypatch.setattr(PushButton, "relinkFiles", lambda *a, **k: None)
    monkeypatch.setattr(PushButton, "tidyUpABit", lambda *a, **k: None)
    monkeypatch.setattr(PushButton, "stripRights", lambda *a, **k: None)
    monkeypatch.setattr(PushButton, "removeDone", lambda *a, **k: None)

    def fakeTouchDone(outputDir, fname="analysis.done"):
        calls["touchDone"].append(fname)

    monkeypatch.setattr(PushButton, "touchDone", fakeTouchDone)

    def fakeCopyCellRanger(*a, **k):
        calls["copyCellRanger"].append(a)

    monkeypatch.setattr(PushButton, "copyCellRanger", fakeCopyCellRanger)

    def fakeRunManagedSubprocess(cmd, **kwargs):
        calls["runManagedSubprocess"].append(cmd)
        if subprocessRaises:
            raise RuntimeError("driver failed")

    monkeypatch.setattr(PushButton, "runManagedSubprocess", fakeRunManagedSubprocess)
    return calls


class TestDNA:
    def test_relacs_protocol_delegates_to_relacs(self, tmp_path, monkeypatch):
        seen = {}

        def fakeRELACS(*a, **k):
            seen["called"] = a
            return ("d", 0, False)

        monkeypatch.setattr(PushButton, "RELACS", fakeRELACS)
        tuples = [["lib1", "s1", "ChIP RELACS high-throughput", False]]
        result = PushButton.DNA(make_config(), "grp", "Proj", HUMAN, "ChIP-Seq", tuples)
        assert "called" in seen
        assert result == ("d", 0, False)

    def test_already_done_skips_dispatch(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        (outputDir / "analysis.done").touch()
        calls = stubHelpers(monkeypatch, outputDir)
        tuples = [["lib1", "s1", "WGS", False]]

        result = PushButton.DNA(make_config(), "grp", "Proj", HUMAN, "WGS", tuples)

        assert result == (str(outputDir), 0, False)
        assert calls["runManagedSubprocess"] == []

    def test_success_runs_pipeline_and_touches_done(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        calls = stubHelpers(monkeypatch, outputDir)
        tuples = [["lib1", "s1", "WGS", False]]

        result = PushButton.DNA(make_config(), "grp", "Proj", HUMAN, "WGS", tuples)

        assert result == (str(outputDir), 0, False)
        assert len(calls["runManagedSubprocess"]) == 1
        assert calls["touchDone"] == ["analysis.done"]

    def test_failure_returns_nonzero_and_skips_touch_done(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        calls = stubHelpers(monkeypatch, outputDir, subprocessRaises=True)
        tuples = [["lib1", "s1", "WGS", False]]

        result = PushButton.DNA(make_config(), "grp", "Proj", HUMAN, "WGS", tuples)

        assert result == (str(outputDir), 1, False)
        assert calls["touchDone"] == []

    @pytest.mark.parametrize(
        "libraryType,expectedFlag",
        [
            ("CUTandTag-seq", "--cutntag"),
            ("CUTandRUN-seq", "--cutntag"),
            ("ATAC-Seq", "nexteraF"),
            ("WGS", None),
        ],
    )
    def test_cmd_varies_by_library_type(
        self, tmp_path, monkeypatch, libraryType, expectedFlag
    ):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        calls = stubHelpers(monkeypatch, outputDir)
        tuples = [["lib1", "s1", "some protocol", False]]

        PushButton.DNA(make_config(), "grp", "Proj", HUMAN, libraryType, tuples)

        cmd = calls["runManagedSubprocess"][0]
        if expectedFlag:
            assert expectedFlag in cmd
        else:
            assert "--cutntag" not in cmd
            assert "nexteraF" not in cmd


class TestWGBS:
    def test_already_done_skips_dispatch(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        (outputDir / "analysis.done").touch()
        calls = stubHelpers(monkeypatch, outputDir)
        tuples = [["lib1", "s1", "WGBS", False]]

        result = PushButton.WGBS(make_config(), "grp", "Proj", HUMAN, "WGBS", tuples)

        assert result == (str(outputDir), 0, False)
        assert calls["runManagedSubprocess"] == []

    def test_success_runs_pipeline(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        calls = stubHelpers(monkeypatch, outputDir)
        tuples = [["lib1", "s1", "WGBS", False]]

        result = PushButton.WGBS(make_config(), "grp", "Proj", HUMAN, "WGBS", tuples)

        assert result == (str(outputDir), 0, False)
        assert "WGBS" in calls["runManagedSubprocess"][0]
        assert calls["touchDone"] == ["analysis.done"]

    def test_failure_returns_nonzero(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        stubHelpers(monkeypatch, outputDir, subprocessRaises=True)
        tuples = [["lib1", "s1", "WGBS", False]]

        result = PushButton.WGBS(make_config(), "grp", "Proj", HUMAN, "WGBS", tuples)

        assert result == (str(outputDir), 1, False)


class TestATAC:
    def test_already_done_skips_everything(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        (outputDir / "analysis.done").touch()
        calls = stubHelpers(monkeypatch, outputDir)
        tuples = [["lib1", "s1", "ATAC-Seq", False]]

        result = PushButton.ATAC(
            make_config(), "grp", "Proj", HUMAN, "ATAC-Seq", tuples
        )

        assert result == (str(outputDir), 0, False)
        assert calls["runManagedSubprocess"] == []

    def test_runs_dna_stage_then_atac_stage(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        calls = stubHelpers(monkeypatch, outputDir)
        tuples = [["lib1", "s1", "ATAC-Seq", False]]

        result = PushButton.ATAC(
            make_config(), "grp", "Proj", HUMAN, "ATAC-Seq", tuples
        )

        assert result == (str(outputDir), 0, False)
        # DNA stage's DNAmapping CMD, then ATACseq CMD
        assert len(calls["runManagedSubprocess"]) == 2
        assert "ATACseq" in calls["runManagedSubprocess"][1]
        # DNA()'s own success path touches "analysis.done" first, ATAC()
        # promotes that to "DNA.done", then the ATACseq stage's own success
        # touches "analysis.done" again.
        assert calls["touchDone"] == ["analysis.done", "DNA.done", "analysis.done"]

    def test_dna_stage_already_done_skips_straight_to_atac(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        (outputDir / "DNA.done").touch()
        calls = stubHelpers(monkeypatch, outputDir)
        tuples = [["lib1", "s1", "ATAC-Seq", False]]

        PushButton.ATAC(make_config(), "grp", "Proj", HUMAN, "ATAC-Seq", tuples)

        assert len(calls["runManagedSubprocess"]) == 1
        assert "ATACseq" in calls["runManagedSubprocess"][0]

    def test_dna_stage_failure_propagates_and_skips_atac(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        calls = stubHelpers(monkeypatch, outputDir, subprocessRaises=True)
        tuples = [["lib1", "s1", "ATAC-Seq", False]]

        result = PushButton.ATAC(
            make_config(), "grp", "Proj", HUMAN, "ATAC-Seq", tuples
        )

        assert result == (str(outputDir), 1, False)


class TestHiC:
    def test_already_done_skips_dispatch(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        (outputDir / "analysis.done").touch()
        calls = stubHelpers(monkeypatch, outputDir)
        tuples = [["lib1", "s1", "HiC", False]]

        result = PushButton.HiC(make_config(), "grp", "Proj", HUMAN, "HiC", tuples)

        assert result == (str(outputDir), 0, False)
        assert calls["runManagedSubprocess"] == []

    def test_success_runs_pipeline(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        calls = stubHelpers(monkeypatch, outputDir)
        tuples = [["lib1", "s1", "HiC", False]]

        result = PushButton.HiC(make_config(), "grp", "Proj", HUMAN, "HiC", tuples)

        assert result == (str(outputDir), 0, False)
        assert "DpnII" in calls["runManagedSubprocess"][0]

    def test_failure_returns_nonzero(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        stubHelpers(monkeypatch, outputDir, subprocessRaises=True)
        tuples = [["lib1", "s1", "HiC", False]]

        result = PushButton.HiC(make_config(), "grp", "Proj", HUMAN, "HiC", tuples)

        assert result == (str(outputDir), 1, False)


class TestMakePairs:
    def test_already_done_skips_dispatch(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        (outputDir / "analysis.done").touch()
        calls = stubHelpers(monkeypatch, outputDir)
        tuples = [["lib1", "s1", "HiC", False]]

        result = PushButton.makePairs(
            make_config(), "grp", "Proj", HUMAN, "HiC", tuples
        )

        assert result == (str(outputDir), 0, False)
        assert calls["runManagedSubprocess"] == []

    def test_success_runs_pipeline(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        calls = stubHelpers(monkeypatch, outputDir)
        tuples = [["lib1", "s1", "HiC", False]]

        result = PushButton.makePairs(
            make_config(), "grp", "Proj", HUMAN, "HiC", tuples
        )

        assert result == (str(outputDir), 0, False)
        assert "makePairs" in calls["runManagedSubprocess"][0]

    def test_failure_returns_nonzero(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        stubHelpers(monkeypatch, outputDir, subprocessRaises=True)
        tuples = [["lib1", "s1", "HiC", False]]

        result = PushButton.makePairs(
            make_config(), "grp", "Proj", HUMAN, "HiC", tuples
        )

        assert result == (str(outputDir), 1, False)


class TestScRNAseq:
    def test_already_done_skips_dispatch(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        (outputDir / "analysis.done").touch()
        calls = stubHelpers(monkeypatch, outputDir)
        tuples = [
            [
                "lib1",
                "s1",
                "Chromium_NextGEM_SingleCell3Prime_GeneExpression_v3.1_DualIndex",
                False,
            ]
        ]

        result = PushButton.scRNAseq(
            make_config(), "grp", "Proj", HUMAN, "scRNA-Seq", tuples
        )

        assert result == (str(outputDir), 0, True)
        assert calls["runManagedSubprocess"] == []

    def test_10x_protocol_runs_10x_pipeline_and_copies_cellranger(
        self, tmp_path, monkeypatch
    ):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        calls = stubHelpers(monkeypatch, outputDir)
        tuples = [
            [
                "lib1",
                "s1",
                "Chromium_NextGEM_SingleCell3Prime_GeneExpression_v3.1_DualIndex",
                False,
            ]
        ]

        result = PushButton.scRNAseq(
            make_config(), "grp", "Proj", HUMAN, "scRNA-Seq", tuples
        )

        assert result == (str(outputDir), 0, True)
        assert "/tenx/rna_wf" in calls["runManagedSubprocess"][0]
        assert len(calls["copyCellRanger"]) == 1

    def test_celseq_protocol_runs_snakepipes_and_skips_cellranger(
        self, tmp_path, monkeypatch
    ):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        calls = stubHelpers(monkeypatch, outputDir)
        tuples = [["lib1", "s1", "Cel-Seq 2 for single cell RNA-Seq", False]]

        result = PushButton.scRNAseq(
            make_config(), "grp", "Proj", HUMAN, "scRNA-Seq", tuples
        )

        assert result == (str(outputDir), 0, False)
        assert "scRNAseq" in calls["runManagedSubprocess"][0]
        assert calls["copyCellRanger"] == []

    def test_unsupported_protocol_just_passes_through(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        calls = stubHelpers(monkeypatch, outputDir)
        tuples = [["lib1", "s1", "SomeUnknownProtocol", False]]

        result = PushButton.scRNAseq(
            make_config(), "grp", "Proj", HUMAN, "scRNA-Seq", tuples
        )

        assert result == (str(outputDir), 0, False)
        assert calls["runManagedSubprocess"] == []
        assert calls["touchDone"] == ["analysis.done"]

    def test_10x_protocol_failure_returns_nonzero(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        stubHelpers(monkeypatch, outputDir, subprocessRaises=True)
        tuples = [
            [
                "lib1",
                "s1",
                "Chromium_NextGEM_SingleCell3Prime_GeneExpression_v3.1_DualIndex",
                False,
            ]
        ]

        result = PushButton.scRNAseq(
            make_config(), "grp", "Proj", HUMAN, "scRNA-Seq", tuples
        )

        assert result == (str(outputDir), 1, False)


class TestScATAC:
    def test_already_done_skips_dispatch(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        (outputDir / "analysis.done").touch()
        calls = stubHelpers(monkeypatch, outputDir)
        tuples = [["lib1", "s1", "scATAC-Seq 10xGenomics", False]]

        result = PushButton.scATAC(
            make_config(), "grp", "Proj", HUMAN, "scATAC-Seq", tuples
        )

        assert result == (str(outputDir), 0, True)
        assert calls["runManagedSubprocess"] == []

    @pytest.mark.parametrize(
        "protocol",
        [
            "scATAC-Seq 10xGenomics",
            "Next GEM Single Cell ATAC",
            "Chromium Next GEM Single Cell ATAC v2",
        ],
    )
    def test_supported_protocols_run_10x_atac_pipeline(
        self, tmp_path, monkeypatch, protocol
    ):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        calls = stubHelpers(monkeypatch, outputDir)
        tuples = [["lib1", "s1", protocol, False]]
        groupSeqDir = None
        # non-external: getLatestSeqdir hits the real filesystem via
        # BRB.misc.getLatestSeqdir(groupData, group) -- stub it out too.
        monkeypatch.setattr(
            PushButton.BRB.misc, "getLatestSeqdir", lambda *a, **k: "sequencing_data"
        )

        result = PushButton.scATAC(
            make_config(), "grp", "Proj", HUMAN, "scATAC-Seq", tuples
        )

        assert result == (str(outputDir), 0, True)
        assert "/tenx/atac_wf" in calls["runManagedSubprocess"][0]
        assert len(calls["copyCellRanger"]) == 1

    def test_external_project_roots_input_dir_under_basedata(
        self, tmp_path, monkeypatch
    ):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        calls = stubHelpers(monkeypatch, outputDir)
        tuples = [["lib1", "s1", "scATAC-Seq 10xGenomics", True]]

        result = PushButton.scATAC(
            make_config(), "grp", "Proj", HUMAN, "scATAC-Seq", tuples
        )

        assert result == (str(outputDir), 0, True)
        cmd = calls["runManagedSubprocess"][0]
        assert "-i /base/20260101_AV999999_1234567_lanes_1_2/Project_Proj" in cmd

    def test_unsupported_protocol_touches_done_without_running_pipeline(
        self, tmp_path, monkeypatch
    ):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        calls = stubHelpers(monkeypatch, outputDir)
        tuples = [["lib1", "s1", "SomeUnknownProtocol", False]]

        result = PushButton.scATAC(
            make_config(), "grp", "Proj", HUMAN, "scATAC-Seq", tuples
        )

        assert result == (str(outputDir), 0, False)
        assert calls["runManagedSubprocess"] == []
        assert calls["touchDone"] == ["analysis.done"]

    def test_failure_returns_nonzero(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        stubHelpers(monkeypatch, outputDir, subprocessRaises=True)
        tuples = [["lib1", "s1", "scATAC-Seq 10xGenomics", True]]

        result = PushButton.scATAC(
            make_config(), "grp", "Proj", HUMAN, "scATAC-Seq", tuples
        )

        assert result == (str(outputDir), 1, False)


def relacsConfig(sequencerType, runID="20260101_AV999999_9999999", **pathOverrides):
    config = make_config()
    config["Options"]["runID"] = runID
    config["Options"]["sequencerType"] = sequencerType
    for k, v in pathOverrides.items():
        config["Paths"][k] = v
    return config


def stubRelacsHelpers(monkeypatch, subprocessRaises=None):
    """
    Stub RELACS's own post-pipeline collaborators (distinct from stubHelpers,
    which targets the plain DNA/WGBS/ATAC pipelines). `subprocessRaises` is a
    predicate `cmd -> bool` deciding which runManagedSubprocess calls raise.
    """
    calls = {"runManagedSubprocess": [], "touchDone": []}

    monkeypatch.setattr(
        PushButton.BRB.misc, "getLatestSeqdir", lambda *a, **k: "seqdir1"
    )
    monkeypatch.setattr(PushButton, "removeLinkFiles", lambda *a, **k: None)
    monkeypatch.setattr(PushButton, "tidyUpABit", lambda *a, **k: None)
    monkeypatch.setattr(PushButton, "copyRELACS", lambda *a, **k: None)
    monkeypatch.setattr(PushButton, "stripRights", lambda *a, **k: None)

    def fakeTouchDone(outputDir, fname="analysis.done"):
        calls["touchDone"].append(fname)

    monkeypatch.setattr(PushButton, "touchDone", fakeTouchDone)

    def fakeRunManagedSubprocess(cmd, **kwargs):
        calls["runManagedSubprocess"].append(cmd)
        if "demultiplex_relacs" in cmd:
            # Simulate the demux driver producing its outputs.
            demuxDir = os.path.join(
                kwargs.get("cwd", ""), "RELACS_demultiplexing", "out"
            )
            os.makedirs(demuxDir, exist_ok=True)
            open(os.path.join(demuxDir, "sample1.png"), "w").close()
            open(os.path.join(demuxDir, "sample1_R1.fastq.gz"), "w").close()
        if subprocessRaises and subprocessRaises(cmd):
            raise RuntimeError("driver failed")

    monkeypatch.setattr(PushButton, "runManagedSubprocess", fakeRunManagedSubprocess)
    return calls


class TestRELACS:
    def test_missing_samplesheet_returns_error(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        monkeypatch.setattr(PushButton, "createPath", lambda *a, **k: str(outputDir))
        tuples = [["lib1", "sample1", "ChIP RELACS high-throughput", False]]

        result = PushButton.RELACS(
            relacsConfig("NovaSeq"), "grp", "Proj", HUMAN, "ChIP-Seq", tuples
        )

        assert result == (None, 1, False)

    def test_already_done_external_triggers_delivery(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        (outputDir / "analysis.done").touch()
        monkeypatch.setattr(PushButton, "createPath", lambda *a, **k: str(outputDir))
        seen = {}
        monkeypatch.setattr(
            PushButton,
            "deliverExternalRELACS",
            lambda cfg, d, project: seen.setdefault("called", (d, project)),
        )
        tuples = [["lib1", "sample1", "ChIP RELACS high-throughput", True]]

        result = PushButton.RELACS(
            relacsConfig("NovaSeq"), "grp", "Proj", HUMAN, "ChIP-Seq", tuples
        )

        assert result == (str(outputDir), 0, True)
        assert seen["called"] == (str(outputDir), "Proj")

    def test_already_done_external_with_delivery_done_skips_delivery(
        self, tmp_path, monkeypatch
    ):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        (outputDir / "analysis.done").touch()
        (outputDir / "external_delivery.done").touch()
        monkeypatch.setattr(PushButton, "createPath", lambda *a, **k: str(outputDir))
        monkeypatch.setattr(
            PushButton,
            "deliverExternalRELACS",
            lambda *a, **k: pytest.fail("delivery should not run twice"),
        )
        tuples = [["lib1", "sample1", "ChIP RELACS high-throughput", True]]

        result = PushButton.RELACS(
            relacsConfig("NovaSeq"), "grp", "Proj", HUMAN, "ChIP-Seq", tuples
        )

        assert result == (str(outputDir), 0, True)

    def test_already_done_internal_skips_delivery(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        (outputDir / "analysis.done").touch()
        monkeypatch.setattr(PushButton, "createPath", lambda *a, **k: str(outputDir))
        monkeypatch.setattr(
            PushButton,
            "deliverExternalRELACS",
            lambda *a, **k: pytest.fail("internal project must not be delivered"),
        )
        tuples = [["lib1", "sample1", "ChIP RELACS high-throughput", False]]

        result = PushButton.RELACS(
            relacsConfig("NovaSeq"), "grp", "Proj", HUMAN, "ChIP-Seq", tuples
        )

        assert result == (str(outputDir), 0, True)

    def test_aviti_success_copies_samplesheet_and_runs_demux(
        self, tmp_path, monkeypatch
    ):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        monkeypatch.setattr(PushButton, "createPath", lambda *a, **k: str(outputDir))

        runID = "20260101_AV999999_9999999"
        sampleSheetSrc = tmp_path / "RELACS_Project_Proj.txt"
        sampleSheetSrc.write_text("sample1\tACGT\n")
        realGlob = glob_mod.glob

        def fakeGlob(pattern, **kwargs):
            if pattern.startswith(f"/dont_touch_this/short_runs/AV*/AV*/{runID}"):
                return [str(sampleSheetSrc)]
            return realGlob(pattern, **kwargs)

        monkeypatch.setattr(PushButton.glob, "glob", fakeGlob)

        groupData = tmp_path / "groupdata"
        monkeypatch.setattr(
            PushButton.BRB.misc, "getLatestSeqdir", lambda *a, **k: "seqdir1"
        )
        sampleDir = (
            groupData / "grp" / "seqdir1" / runID / "Project_Proj" / "Sample_lib1"
        )
        sampleDir.mkdir(parents=True)

        calls = stubRelacsHelpers(monkeypatch)
        tuples = [["lib1", "sample1", "ChIP RELACS high-throughput", False]]

        result = PushButton.RELACS(
            relacsConfig("Aviti", runID=runID, groupData=str(groupData)),
            "grp",
            "Proj",
            HUMAN,
            "ChIP-Seq",
            tuples,
        )

        assert result == (str(outputDir), 0, True)
        assert (outputDir / "RELACS_sampleSheet.txt").read_text() == "sample1\tACGT\n"
        assert any("demultiplex_relacs" in c for c in calls["runManagedSubprocess"])
        assert any("DNAmapping" in c for c in calls["runManagedSubprocess"])
        assert calls["touchDone"] == ["analysis.done"]
        # Demuxed, non-"unknown" reads get relinked under originalFASTQ.
        assert (outputDir / "originalFASTQ" / "sample1_R1.fastq.gz").is_symlink()
        # The transient Sample_* symlink used only to feed the demux CMD is
        # cleaned back up afterwards.
        assert not (outputDir / "Sample_lib1").exists()

    def test_non_aviti_skips_demux_when_pngs_already_match_and_delivers_external(
        self, tmp_path, monkeypatch
    ):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        monkeypatch.setattr(PushButton, "createPath", lambda *a, **k: str(outputDir))
        (outputDir / "RELACS_sampleSheet.txt").write_text("sample1\tACGT\n")
        demuxDir = outputDir / "RELACS_demultiplexing" / "out"
        demuxDir.mkdir(parents=True)
        (demuxDir / "sample1.png").touch()
        (demuxDir / "sample1_R1.fastq.gz").touch()

        calls = stubRelacsHelpers(monkeypatch)
        seen = {}
        monkeypatch.setattr(
            PushButton,
            "deliverExternalRELACS",
            lambda cfg, d, project: seen.setdefault("delivered", (d, project)),
        )
        tuples = [["lib1", "sample1", "ChIP RELACS high-throughput", True]]

        result = PushButton.RELACS(
            relacsConfig("NovaSeq"), "grp", "Proj", HUMAN, "ChIP-Seq", tuples
        )

        assert result == (str(outputDir), 0, True)
        assert not any("demultiplex_relacs" in c for c in calls["runManagedSubprocess"])
        assert seen["delivered"] == (str(outputDir), "Proj")

    def test_demux_failure_returns_nonzero(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        monkeypatch.setattr(PushButton, "createPath", lambda *a, **k: str(outputDir))
        (outputDir / "RELACS_sampleSheet.txt").write_text("sample1\tACGT\n")

        stubRelacsHelpers(
            monkeypatch, subprocessRaises=lambda cmd: "demultiplex_relacs" in cmd
        )
        tuples = [["lib1", "sample1", "ChIP RELACS high-throughput", False]]

        result = PushButton.RELACS(
            relacsConfig("NovaSeq"), "grp", "Proj", HUMAN, "ChIP-Seq", tuples
        )

        assert result == (str(outputDir), 1, False)

    def test_dnamapping_failure_returns_nonzero(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        monkeypatch.setattr(PushButton, "createPath", lambda *a, **k: str(outputDir))
        (outputDir / "RELACS_sampleSheet.txt").write_text("sample1\tACGT\n")
        demuxDir = outputDir / "RELACS_demultiplexing" / "out"
        demuxDir.mkdir(parents=True)
        (demuxDir / "sample1.png").touch()
        (demuxDir / "sample1_R1.fastq.gz").touch()

        stubRelacsHelpers(monkeypatch, subprocessRaises=lambda cmd: "DNAmapping" in cmd)
        tuples = [["lib1", "sample1", "ChIP RELACS high-throughput", False]]

        result = PushButton.RELACS(
            relacsConfig("NovaSeq"), "grp", "Proj", HUMAN, "ChIP-Seq", tuples
        )

        assert result == (str(outputDir), 1, False)
