"""
Lane-split flowcells (dissectBCL --forcelanesplit): one project, one pool per
lane, demultiplexed into <flowcell>_lanes_1 and <flowcell>_lanes_2. All names
below are fake; protocol strings are the Parkour protocol names from the
parkour2 fixtures (library_sample_shared/fixtures/libraryprotocol.json).
"""

import configparser

import pytest

import BRB.ET
from BRB import PushButton

HUMAN = ("human", "hg38", "/yaml/hg38.yaml")
FLOWCELL = "261028_A00931_0938_BXXXXXXXX"
PROJECT = "9999_Foo_Bar"
GEX = "Chromium_NextGEM_SingleCell3Prime_GeneExpression_v3.1_DualIndex"
ATAC = "Chromium Next GEM Single Cell ATAC v2"

# barcode -> [name, analysisType, protocol, organism, indexType, depth]
PARKOUR = {
    PROJECT: {
        "99L001": ["gex_a", "RNA-Seq single cell", GEX, HUMAN, "i7", 30],
        "99L002": ["gex_b", "RNA-Seq single cell", GEX, HUMAN, "i7", 30],
        "99L003": ["atac_a", "ATAC-Seq single cell", ATAC, HUMAN, "i7", 30],
        "99L004": ["atac_b", "ATAC-Seq single cell", ATAC, HUMAN, "i7", 30],
    }
}
LANE1 = {"99L001": "gex_a", "99L002": "gex_b"}
LANE2 = {"99L003": "atac_a", "99L004": "atac_b"}


def laneConfig(tmp_path, lane):
    config = configparser.ConfigParser()
    config["Paths"] = {
        "baseData": str(tmp_path / "base"),
        # groupData/<group> is absent -> GetResults ignore=True, so the
        # [external] allowlist keeps the libraries and output is under baseData.
        "groupData": str(tmp_path / "group"),
        "seqFacDir": str(tmp_path / "seqfac"),
        "bioinfoCoreDir": str(tmp_path / "bioinfocore"),
    }
    config["Options"] = {
        "runID": f"{FLOWCELL}_lanes_{lane}",
        "sequencerType": "Illumina",
        "validAnalysisTypes": "RNA-Seq single cell,ATAC-Seq single cell",
        "pipelines": "scRNAseq,scATAC",
        "snakemakeWorkflowBaseDir": "/snakepipes",
    }
    config["external"] = {
        "AnalysisTypes": "RNA-Seq single cell,ATAC-Seq single cell",
        "LibraryProtocols": "",
    }
    config["10x"] = {"RNA": "/tenx/rna_wf", "ATAC": "/tenx/atac_wf"}
    (tmp_path / "group").mkdir(parents=True, exist_ok=True)
    (tmp_path / "bioinfocore").mkdir(parents=True, exist_ok=True)
    return config


def makeLane(tmp_path, lane, samples, suffixes):
    """Fake dissectBCL output: <flowcell>_lanes_<n>/Project_*/Sample_<barcode>/<name>_<sfx>.fastq.gz"""
    project = tmp_path / "base" / f"{FLOWCELL}_lanes_{lane}" / f"Project_{PROJECT}"
    for barcode, name in samples.items():
        sampleDir = project / f"Sample_{barcode}"
        sampleDir.mkdir(parents=True)
        for sfx in suffixes:
            (sampleDir / f"{name}_{sfx}.fastq.gz").touch()


@pytest.fixture
def flowcell(tmp_path, monkeypatch):
    makeLane(tmp_path, 1, LANE1, ["RNA_R1", "RNA_R2", "RNA_R1_optical_duplicates"])
    makeLane(tmp_path, 2, LANE2, ["ATAC_I1", "ATAC_I2", "ATAC_R1", "ATAC_R2"])
    monkeypatch.setattr(PushButton, "POOL_SIZE", 1)
    monkeypatch.setattr(
        BRB.ET,
        "phoneHome",
        lambda config, outputDir, pipeline, tuples, org, label, project, libType: [
            project,
            org,
            libType,
            pipeline,
            "success",
            "PARKOUR_OK",
        ],
    )
    for fn in ("copyCellRanger", "tidyUpABit", "stripRights", "touchDone"):
        monkeypatch.setattr(PushButton, fn, lambda *a, **k: None)
    return tmp_path


def recordDispatch(monkeypatch):
    """Record (pipeline, [sample names]) per dispatch, without running anything."""
    seen = []

    def make(pipeline):
        def stub(config, group, project, organism, analysisType, tuples):
            seen.append((pipeline, [t[1] for t in tuples]))
            return "/out", 0, False

        return stub

    monkeypatch.setattr(PushButton, "scRNAseq", make("scRNAseq"))
    monkeypatch.setattr(PushButton, "scATAC", make("scATAC"))
    return seen


class TestLaneFilter:
    def test_each_lane_dispatches_only_its_own_samples(self, flowcell, monkeypatch):
        seen = recordDispatch(monkeypatch)

        msg1 = PushButton.runFlowcell(laneConfig(flowcell, 1), PARKOUR)
        msg2 = PushButton.runFlowcell(laneConfig(flowcell, 2), PARKOUR)

        assert [s for p, s in seen if p == "scRNAseq"] == [["gex_a", "gex_b"]]
        assert [s for p, s in seen if p == "scATAC"] == [["atac_a", "atac_b"]]
        assert [m[3] for m in msg1] == ["scRNAseq"]
        assert [m[3] for m in msg2] == ["scATAC"]

    def test_project_without_samples_in_this_lane_is_skipped_and_logged(
        self, flowcell, monkeypatch, caplog
    ):
        # Project dir exists in lane 2, but none of lane-1's samples are in it.
        seen = recordDispatch(monkeypatch)
        lane1Only = {PROJECT: {k: v for k, v in PARKOUR[PROJECT].items() if k in LANE1}}
        caplog.set_level("INFO")

        msg = PushButton.runFlowcell(laneConfig(flowcell, 2), lane1Only)

        assert msg == []
        assert seen == []
        assert f"Project_{PROJECT}" in caplog.text
        assert "no samples" in caplog.text.lower()

    def test_scatac_gets_only_lane_samples_in_samples_arg(self, flowcell, monkeypatch):
        """Real scATAC: the --samples argument must not name the other lane's samples."""
        cmds = []
        monkeypatch.setattr(
            PushButton, "runManagedSubprocess", lambda cmd, **k: cmds.append(cmd)
        )

        PushButton.runFlowcell(laneConfig(flowcell, 2), PARKOUR)

        assert len(cmds) == 1
        samples = cmds[0].split(" --samples ")[1].split(" --snakemakePath")[0]
        assert samples.split() == ["atac_a", "atac_b"]


class TestProtocolGrouping:
    def test_mixed_protocols_in_one_analysis_type_get_separate_work_items(
        self, tmp_path
    ):
        config = laneConfig(tmp_path, 1)
        other = "Chromium_GEM-X_SingleCell_3primeRNA-seq_v4"
        libraries = {
            "99L001": ["a", "RNA-Seq single cell", GEX, HUMAN, "i7", 30],
            "99L002": ["b", "RNA-Seq single cell", other, HUMAN, "i7", 30],
        }

        items, _ = PushButton.GetResults(config, PROJECT, libraries)

        for item in items:
            assert len({t[2] for t in item.tuples}) == 1


MULTI_GEX = "Chromium_NextGEM_Multiome_GEX"
MULTI_ATAC = "Chromium_NextGEM_Multiome_ATAC"
# Real Multiome names: <base>_RNA / <base>_ATAC, FASTQs named the same way.
MULTIOME = {
    PROJECT: {
        "99L001": ["con_RNA", "RNA-Seq single cell", MULTI_GEX, HUMAN, "i7", 30],
        "99L002": ["ko_RNA", "RNA-Seq single cell", MULTI_GEX, HUMAN, "i7", 30],
        "99L003": ["con_ATAC", "ATAC-Seq single cell", MULTI_ATAC, HUMAN, "i7", 30],
        "99L004": ["ko_ATAC", "ATAC-Seq single cell", MULTI_ATAC, HUMAN, "i7", 30],
    }
}


@pytest.fixture
def multiome(tmp_path, monkeypatch):
    makeLane(tmp_path, 1, {"99L001": "con_RNA", "99L002": "ko_RNA"}, ["R1", "R2"])
    makeLane(
        tmp_path,
        2,
        {"99L003": "con_ATAC", "99L004": "ko_ATAC"},
        ["I1", "I2", "R1", "R2"],
    )
    for lane in (1, 2):
        (tmp_path / "base" / f"{FLOWCELL}_lanes_{lane}" / "fastq.made").touch()
    monkeypatch.setattr(PushButton, "POOL_SIZE", 1)
    monkeypatch.setattr(
        BRB.ET,
        "phoneHome",
        lambda config, outputDir, pipeline, tuples, org, label, project, libType: [
            project,
            org,
            libType,
            pipeline,
            "success",
            "PARKOUR_OK",
        ],
    )
    for fn in ("copyCellRanger", "tidyUpABit", "stripRights"):
        monkeypatch.setattr(PushButton, fn, lambda *a, **k: None)
    cmds = []
    monkeypatch.setattr(
        PushButton, "runManagedSubprocess", lambda cmd, **k: cmds.append(cmd)
    )
    return tmp_path, cmds


def multiomeConfig(tmp_path, lane):
    config = laneConfig(tmp_path, lane)
    config["10x"]["ARC"] = "/tenx/arc_wf"
    return config


class TestMultiome:
    def test_multiome_protocols_route_to_scMultiome_not_per_lane_pipelines(
        self, tmp_path
    ):
        items, _ = PushButton.GetResults(
            multiomeConfig(tmp_path, 2), PROJECT, MULTIOME[PROJECT]
        )

        assert {i.pipeline for i in items} == {"scMultiome"}

    def test_first_lane_with_both_halves_ready_runs_one_joint_arc(self, multiome):
        tmp_path, cmds = multiome

        msg = PushButton.runFlowcell(multiomeConfig(tmp_path, 1), MULTIOME)

        assert len(cmds) == 1
        cmd = cmds[0]
        assert cmd.startswith("/tenx/arc_wf ")
        assert (
            f"--gex-indir {tmp_path}/base/{FLOWCELL}_lanes_1/Project_{PROJECT} " in cmd
        )
        assert (
            f"--atac-indir {tmp_path}/base/{FLOWCELL}_lanes_2/Project_{PROJECT} " in cmd
        )
        assert "--pair con:con_RNA:con_ATAC" in cmd
        assert "--pair ko:ko_RNA:ko_ATAC" in cmd
        # Output lives under the lane folder that ran it.
        assert f"-o {tmp_path}/base/{FLOWCELL}_lanes_1/Analysis_{PROJECT}/" in cmd
        assert [m[3] for m in msg] == ["scMultiome"]
        assert list((tmp_path / "base" / f"{FLOWCELL}_lanes_1").rglob("analysis.done"))

    def test_pairs_done_in_a_sibling_lane_are_not_rerun(self, multiome):
        tmp_path, cmds = multiome
        PushButton.runFlowcell(multiomeConfig(tmp_path, 1), MULTIOME)
        lane1Out = next(
            (tmp_path / "base" / f"{FLOWCELL}_lanes_1").glob("Analysis_*/*")
        )
        for base in ("con", "ko"):  # what the 10x_ARC wrapper writes per pair
            (lane1Out / f"{base}.done").touch()
        cmds.clear()

        PushButton.runFlowcell(multiomeConfig(tmp_path, 2), MULTIOME)

        assert cmds == []

    def test_lane_processed_before_its_partner_is_demultiplexed_hands_over(
        self, multiome
    ):
        # BRB may pick up lanes_2 while lanes_1 is still being demultiplexed
        # (no fastq.made). run.py then marks lanes_2 analysis.done, so the
        # joint run must happen in lanes_1's later pass.
        tmp_path, cmds = multiome
        (tmp_path / "base" / f"{FLOWCELL}_lanes_1" / "fastq.made").unlink()

        PushButton.runFlowcell(multiomeConfig(tmp_path, 2), MULTIOME)

        assert cmds == []
        assert not list(
            (tmp_path / "base" / f"{FLOWCELL}_lanes_2").rglob("analysis.done")
        )

        (tmp_path / "base" / f"{FLOWCELL}_lanes_1" / "fastq.made").touch()
        PushButton.runFlowcell(multiomeConfig(tmp_path, 1), MULTIOME)

        assert len(cmds) == 1
        assert f"-o {tmp_path}/base/{FLOWCELL}_lanes_1/Analysis_{PROJECT}/" in cmds[0]

    def test_unpaired_multiome_sample_is_not_run(self, multiome):
        tmp_path, cmds = multiome
        # ko has no GEX library in Parkour: only con can be paired.
        partial = {
            PROJECT: {k: v for k, v in MULTIOME[PROJECT].items() if k != "99L002"}
        }

        PushButton.runFlowcell(multiomeConfig(tmp_path, 2), partial)

        assert len(cmds) == 1
        assert "--pair con:con_RNA:con_ATAC" in cmds[0]
        assert "ko" not in cmds[0].split("--pair ", 1)[1]
