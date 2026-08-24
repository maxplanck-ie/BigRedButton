import configparser
import os
import subprocess
from pathlib import Path
from unittest.mock import patch

import pytest

from BRB import PushButton
from BRB.PushButton import (
    RELACS,
    copyCellRanger,
    copyRELACS,
    deliverExternalRELACS,
    isExternallyAllowed,
    parseExternalAllowlist,
)


def make_config(overrides=None):
    """
    Minimal config with just enough set for the functions under test.
    `overrides` is a dict of {section: {key: value}}.
    """
    config = configparser.ConfigParser()
    config["Paths"] = {
        "baseData": "/base",
        "groupData": "/data",
        "seqFacDir": "/seqfac",
        "bioinfoCoreDir": "/bioinfocore",
    }
    config["Options"] = {
        "runID": "20260101_AV999999_1234567_lanes_1_2",
        "sequencerType": "Aviti",
        "snakemakeWorkflowBaseDir": "/snakepipes",
    }
    config["external"] = {
        "LibraryTypes": "",
        "LibraryProtocols": "",
        "fexRecipient": "pipegrp@example.org",
        "fexsendBin": "/fex/fexsend",
    }
    for section, kv in (overrides or {}).items():
        if section not in config:
            config[section] = {}
        config[section].update(kv)
    return config


class TestParseExternalAllowlist:
    def test_empty_string(self):
        assert parseExternalAllowlist("") == []

    def test_single_entry(self):
        assert parseExternalAllowlist("ChIP-Seq") == ["ChIP-Seq"]

    def test_multiple_entries_stripped(self):
        assert parseExternalAllowlist("ChIP-Seq, RNA-Seq ,WGS") == [
            "ChIP-Seq",
            "RNA-Seq",
            "WGS",
        ]

    def test_ignores_blank_entries(self):
        assert parseExternalAllowlist("ChIP-Seq,,  ,RNA-Seq") == [
            "ChIP-Seq",
            "RNA-Seq",
        ]


class TestIsExternallyAllowed:
    def test_exact_type_match(self):
        assert isExternallyAllowed("ChIP-Seq", "anything", ["ChIP-Seq"], []) is True

    def test_no_match(self):
        assert isExternallyAllowed("RNA-Seq", "anything", ["ChIP-Seq"], []) is False

    def test_protocol_startswith_match(self):
        assert (
            isExternallyAllowed(
                "ChIP-Seq",
                "ChIP RELACS high-throughput v2",
                [],
                ["ChIP RELACS high-throughput"],
            )
            is True
        )

    def test_protocol_no_match(self):
        assert (
            isExternallyAllowed(
                "ChIP-Seq", "Standard ChIP-Seq", [], ["ChIP RELACS high-throughput"]
            )
            is False
        )

    def test_no_substring_false_positive(self):
        """
        Regression test: `libraryType in "a,b,c"` (the old code) does
        substring matching against the raw string, so "RNA-Seq" would
        false-positive against "stranded mRNA-Seq". With real list
        membership, it must not.
        """
        assert (
            isExternallyAllowed("stranded mRNA-Seq", "anything", ["RNA-Seq"], [])
            is False
        )

    def test_empty_allowlists(self):
        assert isExternallyAllowed("ChIP-Seq", "ChIP RELACS", [], []) is False


class TestInternalShareSkip:
    """
    copyCellRanger reads config['Paths']['seqFacDir'] and ['bioinfoCoreDir']
    as soon as it proceeds past the ignore check. A config missing those
    keys makes that a hard failure, which we use here to prove the function
    returns immediately (without touching them) when ignore=True, and does
    NOT return early when ignore=False.

    copyRELACS has no such skip: RELACS QC diagnostics (seqFacDir +
    bioinfoCoreDir) are shared for internal and external projects alike, see
    TestCopyRELACSAlwaysRuns below.
    """

    def test_copyCellRanger_skips_when_ignored(self, tmp_path):
        config = configparser.ConfigParser()
        config["Options"] = {"sequencerType": "Aviti"}
        # No [Paths] section at all -- would raise if the function proceeded.
        copyCellRanger(config, str(tmp_path), ignore=True)

    def test_copyCellRanger_raises_when_not_ignored_and_misconfigured(self, tmp_path):
        config = configparser.ConfigParser()
        config["Options"] = {"sequencerType": "Aviti"}
        (tmp_path / "sampleA" / "outs").mkdir(parents=True)
        (tmp_path / "sampleA" / "outs" / "web_summary.html").write_text("x")
        with pytest.raises(configparser.NoSectionError):
            copyCellRanger(config, str(tmp_path), ignore=False)

    def test_relinkFiles_skips_mqc_copy_when_ignored(self, tmp_path):
        config = configparser.ConfigParser()
        config["Paths"] = {"baseData": str(tmp_path), "groupData": str(tmp_path)}
        config["Options"] = {"runID": "run1"}
        # tuples[0][3] = True (ignore); createPath/linkFiles will run fine
        # under baseData, but if the mqc-copy step ran it would raise on the
        # missing bioinfoCoreDir key.
        tuples = [["lib1", "sample1", "protocol", True]]
        PushButton.relinkFiles(config, "group", "proj", "org", "ChIP-Seq", tuples)


class TestCopyRELACSAlwaysRuns:
    """
    copyRELACS() takes no `ignore` param and ships the same seqFacDir +
    bioinfoCoreDir QC diagnostics regardless of whether the project is
    internal or external -- unlike copyCellRanger, which is still skipped
    for external projects.
    """

    def _make_output_dir(self, tmp_path):
        outputDir = (
            tmp_path
            / "20260101_AV999999_1234567_lanes_1_2"
            / "Analysis_99_Foo_Bar"
            / "ChIP-Seq_h38"
        )
        (outputDir / "RELACS_demultiplexing" / "Sample_1").mkdir(parents=True)
        (outputDir / "RELACS_demultiplexing" / "Sample_1" / "mark_fig.png").write_text(
            "x"
        )
        (outputDir / "multiQC").mkdir(parents=True)
        (outputDir / "multiQC" / "multiqc_report.html").write_text("x")
        return outputDir

    def test_copies_to_seqfac_and_bioinfocore(self, tmp_path):
        outputDir = self._make_output_dir(tmp_path)
        seqfac = tmp_path / "seqfac"
        bioinfocore = tmp_path / "bioinfocore"
        bioinfocore.mkdir()
        config = make_config(
            overrides={
                "Paths": {"seqFacDir": str(seqfac), "bioinfoCoreDir": str(bioinfocore)}
            }
        )

        copyRELACS(config, str(outputDir))

        seqfac_files = list(seqfac.rglob("*"))
        assert any(f.name.endswith("mark_fig.png") for f in seqfac_files)
        assert any(f.name.endswith("RELACS_analysis.html") for f in seqfac_files)
        bioinfocore_files = list(bioinfocore.iterdir())
        assert any(f.name.endswith("mark_fig.png") for f in bioinfocore_files)
        assert any(f.name.endswith("RELACS_analysis.html") for f in bioinfocore_files)

    def test_no_ignore_param_accepted(self):
        import inspect

        assert "ignore" not in inspect.signature(copyRELACS).parameters


class TestCopyRELACSUniqueGroupName:
    """
    Two library-groups of the same project (different libraryType or
    organism) must not copy their png/html QC files to the same destination
    filename in seqFacDir or bioinfoCoreDir -- under the Phase 1 thread pool
    they would be writing there concurrently. Mirrors
    TestRelinkFilesUniqueMultiqcName's coverage of the same class of bug in
    relinkFiles's multiqc destination naming.
    """

    def _make_output_dir(self, tmp_path, libraryType, org_label):
        outputDir = (
            tmp_path
            / "20260101_AV999999_1234567_lanes_1_2"
            / "Analysis_99_Foo_Bar"
            / f"{libraryType}_{org_label}"
        )
        (outputDir / "RELACS_demultiplexing" / "Sample_1").mkdir(parents=True)
        (outputDir / "RELACS_demultiplexing" / "Sample_1" / "mark_fig.png").write_text(
            libraryType + org_label
        )
        (outputDir / "multiQC").mkdir(parents=True)
        (outputDir / "multiQC" / "multiqc_report.html").write_text(
            libraryType + org_label
        )
        return outputDir

    def test_two_groups_get_distinct_html_and_png_filenames(self, tmp_path):
        seqfac = tmp_path / "seqfac"
        bioinfocore = tmp_path / "bioinfocore"
        bioinfocore.mkdir()
        config = make_config(
            overrides={
                "Paths": {"seqFacDir": str(seqfac), "bioinfoCoreDir": str(bioinfocore)}
            }
        )

        for libraryType, org_label in (
            ("ChIP-Seq", "hg38"),
            ("CUTandTag-seq", "mm10"),
        ):
            outputDir = self._make_output_dir(tmp_path, libraryType, org_label)
            copyRELACS(config, str(outputDir))

        bioinfocore_pngs = {
            f.name for f in bioinfocore.iterdir() if f.name.endswith(".png")
        }
        bioinfocore_htmls = {
            f.name for f in bioinfocore.iterdir() if f.name.endswith(".html")
        }
        assert len(bioinfocore_pngs) == 2, (
            f"expected two distinct png filenames, got {bioinfocore_pngs}"
        )
        assert len(bioinfocore_htmls) == 2, (
            f"expected two distinct html filenames, got {bioinfocore_htmls}"
        )

        seqfac_files = list(seqfac.rglob("*"))
        seqfac_pngs = {f.name for f in seqfac_files if f.name.endswith(".png")}
        seqfac_htmls = {f.name for f in seqfac_files if f.name.endswith(".html")}
        assert len(seqfac_pngs) == 2, (
            f"expected two distinct png filenames, got {seqfac_pngs}"
        )
        assert len(seqfac_htmls) == 2, (
            f"expected two distinct html filenames, got {seqfac_htmls}"
        )


class TestDeliverExternalRELACS:
    def _make_output_dir(self, tmp_path):
        outputDir = tmp_path / "Analysis_99_Foo_Bar" / "ChIP-Seq_h38"
        outputDir.mkdir(parents=True)
        (outputDir / "RELACS_demultiplexing" / "Sample_1").mkdir(parents=True)
        (
            outputDir / "RELACS_demultiplexing" / "Sample_1" / "a_R1.fastq.gz"
        ).write_bytes(b"not actually empty")
        (outputDir / "FASTQ_Cutadapt").mkdir()
        (outputDir / "FASTQ_Cutadapt" / "a_R1.fastq.gz").write_bytes(b"also not empty")
        (outputDir / "Bowtie2").mkdir()
        (outputDir / "Bowtie2" / "a.bam").write_bytes(b"bam")
        (outputDir / "filtered_bam").mkdir()
        (outputDir / "bamCoverage").mkdir()
        (outputDir / "multiQC").mkdir()
        (outputDir / "DNAmapping.config.yaml").write_text("x: 1")
        (outputDir / "RELACS_sampleSheet.txt").write_text("a\tb\tc")
        return outputDir

    def test_truncates_intermediate_gz_files(self, tmp_path):
        outputDir = self._make_output_dir(tmp_path)
        config = make_config()
        with patch("BRB.PushButton.runManagedSubprocess") as mock_call:
            deliverExternalRELACS(config, str(outputDir), "99_Foo_Bar")
        assert mock_call.called
        gz1 = outputDir / "RELACS_demultiplexing" / "Sample_1" / "a_R1.fastq.gz"
        gz2 = outputDir / "FASTQ_Cutadapt" / "a_R1.fastq.gz"
        assert gz1.stat().st_size == 0
        assert gz2.stat().st_size == 0

    def test_tar_command_includes_expected_dirs_and_excludes_intermediates(
        self, tmp_path
    ):
        outputDir = self._make_output_dir(tmp_path)
        config = make_config()
        with patch("BRB.PushButton.runManagedSubprocess") as mock_call:
            deliverExternalRELACS(config, str(outputDir), "99_Foo_Bar")
        cmd = mock_call.call_args[0][0]
        assert "Bowtie2" in cmd
        assert "filtered_bam" in cmd
        assert "bamCoverage" in cmd
        assert "multiQC" in cmd
        assert "DNAmapping.config.yaml" in cmd
        assert "RELACS_sampleSheet.txt" in cmd
        assert "RELACS_demultiplexing" not in cmd
        assert "FASTQ_Cutadapt" not in cmd
        assert "/fex/fexsend" in cmd
        assert "pipegrp@example.org" in cmd
        assert "Analysis_99_Foo_Bar.tar.gz" in cmd

    def test_creates_done_marker_on_success(self, tmp_path):
        outputDir = self._make_output_dir(tmp_path)
        config = make_config()
        with patch("BRB.PushButton.runManagedSubprocess"):
            deliverExternalRELACS(config, str(outputDir), "99_Foo_Bar")
        assert (outputDir / "external_delivery.done").exists()

    def test_second_call_is_a_noop(self, tmp_path):
        outputDir = self._make_output_dir(tmp_path)
        (outputDir / "external_delivery.done").touch()
        config = make_config()
        with patch("BRB.PushButton.runManagedSubprocess") as mock_call:
            deliverExternalRELACS(config, str(outputDir), "99_Foo_Bar")
        mock_call.assert_not_called()

    def test_failure_is_logged_not_raised(self, tmp_path):
        outputDir = self._make_output_dir(tmp_path)
        config = make_config()
        with patch("BRB.PushButton.runManagedSubprocess", side_effect=OSError("boom")):
            deliverExternalRELACS(config, str(outputDir), "99_Foo_Bar")
        assert not (outputDir / "external_delivery.done").exists()


class TestRELACSExternalPaths:
    """
    End-to-end (with subprocess/network mocked) coverage of RELACS()'s
    external-user code path: the sample-sheet glob and baseDir must resolve
    without a PI group directory, and delivery must be triggered.
    """

    def _base_setup(self, tmp_path, ignore):
        baseData = tmp_path / "base"
        runID = "20260101_AV999999_1234567_lanes_1_2"
        project = "99_Foo_Bar"
        outputDir = (
            baseData / runID / f"Analysis_{project}" / "ChIP-Seq_h38"
            if ignore
            else tmp_path
            / "data"
            / "foo"
            / "sequencing_data"
            / runID
            / f"Analysis_{project}"
            / "ChIP-Seq_h38"
        )
        outputDir.mkdir(parents=True)
        sampleSheetContent = "Sample_1\tACGT\tsampleA\n"

        config = make_config(
            {
                "Paths": {
                    "baseData": str(baseData),
                    "groupData": str(tmp_path / "data"),
                },
                "Options": {"runID": runID},
            }
        )
        tuples = [["lib1", "sampleA", "ChIP RELACS high-throughput", ignore]]
        return config, tuples, project, outputDir, sampleSheetContent

    @patch("BRB.PushButton.createPath")
    @patch("BRB.PushButton.copyRELACS")
    @patch("BRB.PushButton.deliverExternalRELACS")
    @patch("BRB.PushButton.runManagedSubprocess")
    @patch("BRB.PushButton.glob.glob")
    def test_external_uses_baseData_not_groupData(
        self,
        mock_glob,
        mock_check_call,
        mock_deliver,
        mock_copyRELACS,
        mock_createPath,
        tmp_path,
    ):
        config, tuples, project, outputDir, sheetContent = self._base_setup(
            tmp_path, ignore=True
        )
        mock_createPath.return_value = str(outputDir)
        (outputDir / "RELACS_sampleSheet.txt").write_text(sheetContent)
        # 1 premux sample, 1 png already present => skip the demux subprocess call
        (outputDir / "RELACS_demultiplexing").mkdir()
        (outputDir / "RELACS_demultiplexing" / "sampleA_fig.png").write_text("x")

        def glob_side_effect(pattern, recursive=False):
            if "RELACS_Project_" in pattern:
                return [
                    "/dont_touch_this/short_runs/AVITI/AV1/run/RELACS_Project_99_Foo_Bar.txt"
                ]
            if pattern.endswith(".png") or "png" in pattern:
                return [str(outputDir / "RELACS_demultiplexing" / "sampleA_fig.png")]
            if "Sample_*" in pattern:
                # This is the assertion target: baseDir must be baseData-rooted.
                assert pattern.startswith(
                    str(
                        Path(config.get("Paths", "baseData"))
                        / config.get("Options", "runID")
                        / f"Project_{project}"
                    )
                )
                return []
            if "*.gz" in pattern:
                return []
            return []

        mock_glob.side_effect = glob_side_effect
        organism = ("Human", "h38", "GRCh38")
        result = RELACS(config, "villa", project, organism, "ChIP-Seq", tuples)
        assert result[1] == 0
        mock_deliver.assert_called_once()

    @patch("BRB.PushButton.createPath")
    @patch("BRB.PushButton.deliverExternalRELACS")
    @patch("BRB.PushButton.runManagedSubprocess")
    @patch("BRB.PushButton.glob.glob")
    def test_internal_uses_groupData_and_skips_delivery(
        self, mock_glob, mock_check_call, mock_deliver, mock_createPath, tmp_path
    ):
        config, tuples, project, outputDir, sheetContent = self._base_setup(
            tmp_path, ignore=False
        )
        mock_createPath.return_value = str(outputDir)
        (outputDir / "RELACS_sampleSheet.txt").write_text(sheetContent)
        (outputDir / "RELACS_demultiplexing").mkdir()
        (outputDir / "RELACS_demultiplexing" / "sampleA_fig.png").write_text("x")

        groupSeqDir = tmp_path / "data" / "villa" / "sequencing_data"
        groupSeqDir.mkdir(parents=True)

        def glob_side_effect(pattern, recursive=False):
            if "RELACS_Project_" in pattern:
                return ["/dont_touch_this/short_runs/run/RELACS_Project_99_Foo_Bar.txt"]
            if "Sample_*" in pattern:
                assert str(tmp_path / "data" / "villa") in pattern
                return []
            # covers both the copyRELACS() png/html lookup and the
            # RELACS-demultiplexed-fastq linking globs -- nothing to do here
            return []

        mock_glob.side_effect = glob_side_effect
        organism = ("Human", "h38", "GRCh38")
        result = RELACS(config, "villa", project, organism, "ChIP-Seq", tuples)
        assert result[1] == 0
        mock_deliver.assert_not_called()

    @patch("BRB.PushButton.createPath")
    @patch("BRB.PushButton.deliverExternalRELACS")
    @patch("BRB.PushButton.glob.glob")
    def test_sample_sheet_glob_matches_nested_avX_layout(
        self, mock_glob, mock_deliver, mock_createPath, tmp_path
    ):
        """
        Regression test for the /dont_touch_this/short_runs/AV*/AV*/<runID>
        nesting fix (was AV*/<runID>, one level too shallow).
        """
        config, tuples, project, outputDir, _sheetContent = self._base_setup(
            tmp_path, ignore=True
        )
        mock_createPath.return_value = str(outputDir)
        seenPatterns = []

        def glob_side_effect(pattern, recursive=False):
            seenPatterns.append(pattern)
            return []

        mock_glob.side_effect = glob_side_effect
        organism = ("Human", "h38", "GRCh38")
        result = RELACS(config, "villa", project, organism, "ChIP-Seq", tuples)
        assert result[1] == 1  # no sample sheet found -> expected failure
        samplesheetPatterns = [p for p in seenPatterns if "RELACS_Project_" in p]
        assert len(samplesheetPatterns) == 1
        assert "AV*/AV*/" in samplesheetPatterns[0]


class TestRelinkFilesUniqueMultiqcName:
    """
    Two library-groups of the same project must not copy their multiQC
    report to the same bioinfoCoreDir filename -- under the Phase 1 thread
    pool they would be writing there concurrently.
    """

    def test_two_groups_same_project_get_distinct_filenames(
        self, tmp_path, monkeypatch
    ):
        bioinfocore = tmp_path / "bioinfocore"
        bioinfocore.mkdir()
        config = make_config(
            {
                "Paths": {
                    "baseData": str(tmp_path),
                    "groupData": str(tmp_path),
                    "bioinfoCoreDir": str(bioinfocore),
                }
            }
        )
        tuples = [["lib1", "sample1", "protocol", False]]

        for libraryType, org_label in (
            ("ChIP-Seq", "hg38"),
            ("stranded mRNA-Seq", "mm10"),
        ):
            outputDir = tmp_path / "out" / f"{libraryType}_{org_label}"
            (outputDir / "multiQC").mkdir(parents=True)
            (outputDir / "multiQC" / "multiqc_report.html").write_text(
                libraryType + org_label
            )
            monkeypatch.setattr(
                PushButton, "createPath", lambda *a, _o=outputDir, **k: str(_o)
            )
            monkeypatch.setattr(PushButton, "linkFiles", lambda *a, **k: False)
            PushButton.relinkFiles(
                config, "group", "proj", org_label, libraryType, tuples
            )

        assert sorted(p.name for p in bioinfocore.iterdir()) == [
            "Analysisproj_ChIP-Seq_hg38_multiqc.html",
            "Analysisproj_strandedmRNA-Seq_mm10_multiqc.html",
        ]


class TestManagedSubprocessCallSites:
    def test_no_check_call_remains(self):
        source = Path(PushButton.__file__).read_text()
        assert "subprocess.check_call" not in source

    def test_all_twelve_sites_converted(self):
        source = Path(PushButton.__file__).read_text()
        assert source.count("runManagedSubprocess(") == 12

    def test_the_two_cwd_sites_keep_their_cwd(self):
        source = Path(PushButton.__file__).read_text()
        assert source.count("runManagedSubprocess(CMD, cwd=outputDir)") == 1
        assert source.count('runManagedSubprocess(" ".join(CMD), cwd=outputDir)') == 1

    def test_rna_uses_the_managed_runner(self, tmp_path, monkeypatch):
        config = make_config()
        calls = []
        monkeypatch.setattr(PushButton, "createPath", lambda *a, **k: str(tmp_path))
        monkeypatch.setattr(PushButton, "linkFiles", lambda *a, **k: True)
        monkeypatch.setattr(PushButton, "removeLinkFiles", lambda *a, **k: None)
        monkeypatch.setattr(PushButton, "relinkFiles", lambda *a, **k: None)
        monkeypatch.setattr(PushButton, "tidyUpABit", lambda *a, **k: None)
        monkeypatch.setattr(
            PushButton,
            "runManagedSubprocess",
            lambda cmd, **kw: calls.append((cmd, kw)) or 0,
        )
        _outputDir, rv, _samba = PushButton.RNA(
            config,
            "grp",
            "1_Doe_Smith",
            ("mouse", "GRCm38", "/yaml/GRCm38.yaml"),
            "stranded mRNA-Seq",
            [["18L001", "sampleA", "TruSeq", False]],
        )
        assert rv == 0
        assert len(calls) == 1
        assert "mRNAseq" in calls[0][0]
        assert calls[0][1] == {}

    def test_rna_still_returns_rv_1_when_the_command_fails(self, tmp_path, monkeypatch):
        config = make_config()

        def boom(cmd, **kw):
            raise subprocess.CalledProcessError(1, cmd)

        monkeypatch.setattr(PushButton, "createPath", lambda *a, **k: str(tmp_path))
        monkeypatch.setattr(PushButton, "linkFiles", lambda *a, **k: True)
        monkeypatch.setattr(PushButton, "runManagedSubprocess", boom)
        _outputDir, rv, _samba = PushButton.RNA(
            config,
            "grp",
            "1_Doe_Smith",
            ("mouse", "GRCm38", "/yaml/GRCm38.yaml"),
            "stranded mRNA-Seq",
            [["18L001", "sampleA", "TruSeq", False]],
        )
        assert rv == 1


from BRB import jobtrack


def make_item(tmp_path, pipeline="RNA"):
    return PushButton.WorkItem(
        project="1_Doe_Smith",
        group="smith",
        pipeline=pipeline,
        organism=("mouse", "GRCm38", "/yaml/GRCm38.yaml"),
        libraryType="stranded mRNA-Seq",
        tuples=[["18L001", "sampleA", "TruSeq", False]],
    )


def _dead_pid_kill(dead):
    def _kill(pid, sig):
        if pid == dead:
            raise ProcessLookupError

    return _kill


class TestRunOneGroupMarkerGate:
    @pytest.fixture(autouse=True)
    def _stub(self, tmp_path, monkeypatch):
        self.dispatched = []
        monkeypatch.setattr(PushButton, "createPath", lambda *a, **k: str(tmp_path))

        def fakeRNA(config, group, project, organism, libraryType, tuples):
            self.dispatched.append(project)
            return str(tmp_path), 0, False

        monkeypatch.setattr(PushButton, "RNA", fakeRNA)
        monkeypatch.setattr(
            PushButton.BRB.ET,
            "phoneHome",
            lambda *a, **k: [
                "1_Doe_Smith",
                "mouse",
                "stranded mRNA-Seq",
                "RNA",
                "OK",
                "updated",
            ],
        )
        self.tmp_path = tmp_path

    def test_absent_marker_dispatches_and_cleans_up(self):
        reg = jobtrack.JobRegistry()
        row = PushButton.runOneGroup(make_config(), make_item(self.tmp_path), reg)
        assert self.dispatched == ["1_Doe_Smith"]
        assert row[0][4] == "OK"
        assert jobtrack.markerState(self.tmp_path)[0] == jobtrack.MARKER_ABSENT
        assert reg.active_groups() == []

    def test_marker_exists_and_handle_bound_while_the_pipeline_runs(self, monkeypatch):
        seen = {}

        def peek(config, group, project, organism, libraryType, tuples):
            seen["state"] = jobtrack.markerState(self.tmp_path)[0]
            seen["bound"] = jobtrack.currentHandle().outputDir
            return str(self.tmp_path), 0, False

        monkeypatch.setattr(PushButton, "RNA", peek)
        PushButton.runOneGroup(
            make_config(), make_item(self.tmp_path), jobtrack.JobRegistry()
        )
        assert seen["state"] == jobtrack.MARKER_LIVE
        assert seen["bound"] == str(self.tmp_path)

    def test_live_marker_skips_without_dispatching(self):
        jobtrack.writeMarker(self.tmp_path)
        row = PushButton.runOneGroup(
            make_config(), make_item(self.tmp_path), jobtrack.JobRegistry()
        )
        assert row == [
            [
                "1_Doe_Smith",
                "mouse",
                "stranded mRNA-Seq",
                "RNA",
                "SKIPPED (owned by live PID)",
                "not updated",
                False,
                0,
            ]
        ]
        assert self.dispatched == []
        assert jobtrack.markerState(self.tmp_path)[0] == jobtrack.MARKER_LIVE

    def test_abandoned_marker_is_removed_and_group_dispatched(self, monkeypatch):
        jobtrack.writeMarker(self.tmp_path, pid=999999, job_ids=["5"])
        monkeypatch.setattr(jobtrack.os, "kill", _dead_pid_kill(999999))
        row = PushButton.runOneGroup(
            make_config(), make_item(self.tmp_path), jobtrack.JobRegistry()
        )
        assert self.dispatched == ["1_Doe_Smith"]
        assert row[0][4] == "OK"

    def test_corrupt_marker_skips_without_dispatching(self):
        (self.tmp_path / "running.pid").write_text('{"pid": 12')
        row = PushButton.runOneGroup(
            make_config(), make_item(self.tmp_path), jobtrack.JobRegistry()
        )
        assert row[0][4] == "SKIPPED (unreadable marker)"
        assert self.dispatched == []
        assert (self.tmp_path / "running.pid").exists()

    def test_aborted_registry_skips_without_dispatching(self):
        reg = jobtrack.JobRegistry()
        reg.abort()
        row = PushButton.runOneGroup(make_config(), make_item(self.tmp_path), reg)
        assert row[0][4] == "SKIPPED (flowcell teardown in progress)"
        assert self.dispatched == []

    def test_marker_is_cleared_when_the_pipeline_raises(self, monkeypatch):
        def boom(*a, **k):
            raise RuntimeError("driver died")

        monkeypatch.setattr(PushButton, "RNA", boom)
        reg = jobtrack.JobRegistry()
        with pytest.raises(RuntimeError):
            PushButton.runOneGroup(make_config(), make_item(self.tmp_path), reg)
        assert jobtrack.markerState(self.tmp_path)[0] == jobtrack.MARKER_ABSENT
        assert reg.active_groups() == []

    def test_failure_retries_once_then_reports_failed(self, monkeypatch):
        attempts = []

        def failing(config, group, project, organism, libraryType, tuples):
            attempts.append(1)
            return str(self.tmp_path), 1, False

        monkeypatch.setattr(PushButton, "RNA", failing)
        row = PushButton.runOneGroup(
            make_config(), make_item(self.tmp_path), jobtrack.JobRegistry()
        )
        assert len(attempts) == 2
        assert row[0][4] == "FAILED"
        assert row[0][7] == 1

    def test_abort_between_attempts_suppresses_the_retry(self, monkeypatch):
        reg = jobtrack.JobRegistry()
        attempts = []

        def failThenAbort(config, group, project, organism, libraryType, tuples):
            attempts.append(1)
            reg.abort()
            return str(self.tmp_path), 1, False

        monkeypatch.setattr(PushButton, "RNA", failThenAbort)
        row = PushButton.runOneGroup(make_config(), make_item(self.tmp_path), reg)
        assert len(attempts) == 1
        assert row[0][4] == "SKIPPED (flowcell teardown in progress)"
        assert reg.active_groups() == []


class TestRunFlowcellRegistryAndCancellation:
    """
    Covers Task 10's actual scope: `runFlowcell` accepting an injectable
    `registry`, threading it into every `runOneGroup` call, and calling
    `BRB.jobtrack.cancelAllGroups(registry)` exactly once -- before the pool
    is shut down -- the moment the first worker exception is observed.

    `runFlowcell`'s crash path (the break-on-first-exception loop, the
    second drain pass collecting other already-done failures, wrapping into
    `GroupDispatchError`, `shutdown(wait=False, cancel_futures=True)`) is
    Phase 1 / Task 7's already-reviewed behaviour and is exercised by
    tests/test_PushButton_dispatch.py::TestRunFlowcell. These tests only
    check the new cancellation wiring on top of it.
    """

    def _parkourDict(self, name="1_Doe_Smith"):
        return {
            name: {
                "18L001": [
                    "sampleA",
                    "stranded mRNA-Seq",
                    "TruSeq",
                    ("mouse", "GRCm38", "/yaml/GRCm38.yaml"),
                    "single",
                    10,
                ]
            }
        }

    @pytest.fixture(autouse=True)
    def _projectDirExists(self, monkeypatch):
        # runFlowcell's own os.path.exists check (does the project directory
        # exist on this lane) is irrelevant to this test class -- always say
        # yes so every project reaches GetResults.
        monkeypatch.setattr(os.path, "exists", lambda p: True)

    def test_explicit_registry_is_used_not_a_fresh_one(self, tmp_path, monkeypatch):
        reg = jobtrack.JobRegistry()
        cancelled = []
        monkeypatch.setattr(
            PushButton.BRB.jobtrack,
            "cancelAllGroups",
            lambda registry, **kw: cancelled.append(registry),
        )
        item = make_item(tmp_path)
        monkeypatch.setattr(PushButton, "GetResults", lambda *a, **k: ([item], []))
        monkeypatch.setattr(
            PushButton,
            "runOneGroup",
            lambda config, it, registry: (_ for _ in ()).throw(RuntimeError("boom")),
        )
        with pytest.raises(PushButton.GroupDispatchError):
            PushButton.runFlowcell(make_config(), self._parkourDict(), reg)
        assert cancelled == [reg]

    def test_registry_is_threaded_into_every_runonegroup_call(
        self, tmp_path, monkeypatch
    ):
        seen = []
        items = [make_item(tmp_path), make_item(tmp_path)]
        monkeypatch.setattr(PushButton, "GetResults", lambda *a, **k: (items, []))
        monkeypatch.setattr(
            PushButton,
            "runOneGroup",
            lambda config, it, registry: seen.append(registry) or ["row"],
        )
        reg = jobtrack.JobRegistry()
        PushButton.runFlowcell(make_config(), self._parkourDict(), reg)
        assert seen == [reg, reg]

    def test_default_registry_is_constructed_and_still_threaded(
        self, tmp_path, monkeypatch
    ):
        seen = []
        item = make_item(tmp_path)
        monkeypatch.setattr(PushButton, "GetResults", lambda *a, **k: ([item], []))
        monkeypatch.setattr(
            PushButton,
            "runOneGroup",
            lambda config, it, registry: seen.append(registry) or ["row"],
        )
        PushButton.runFlowcell(make_config(), self._parkourDict())
        assert len(seen) == 1
        assert isinstance(seen[0], jobtrack.JobRegistry)

    def test_worker_exception_cancels_once_and_raises_groupdispatcherror(
        self, tmp_path, monkeypatch
    ):
        reg = jobtrack.JobRegistry()
        cancelled = []
        monkeypatch.setattr(
            PushButton.BRB.jobtrack,
            "cancelAllGroups",
            lambda registry, **kw: cancelled.append(registry),
        )
        items = [make_item(tmp_path) for _ in range(4)]
        monkeypatch.setattr(PushButton, "GetResults", lambda *a, **k: (items, []))
        monkeypatch.setattr(PushButton, "POOL_SIZE", 4)

        def boom(config, item, registry):
            raise RuntimeError("driver died")

        monkeypatch.setattr(PushButton, "runOneGroup", boom)
        with pytest.raises(PushButton.GroupDispatchError) as excinfo:
            PushButton.runFlowcell(make_config(), self._parkourDict(), reg)
        assert cancelled == [reg]
        assert isinstance(excinfo.value.__cause__, RuntimeError)
        assert len(excinfo.value.failures) >= 1

    def test_every_captured_worker_exception_is_logged_and_in_failures(
        self, tmp_path, monkeypatch, caplog
    ):
        # Not asserting a specific failure count: how many sibling futures
        # have already reached done() at the instant the first exception
        # surfaces is scheduler-dependent (see the equivalent caveat in
        # test_PushButton_dispatch.py::TestRunFlowcell). What must hold
        # regardless is that every failure the drain pass collected is both
        # logged at CRITICAL and present in GroupDispatchError.failures.
        monkeypatch.setattr(
            PushButton.BRB.jobtrack, "cancelAllGroups", lambda registry, **kw: None
        )
        items = [make_item(tmp_path), make_item(tmp_path)]
        monkeypatch.setattr(PushButton, "GetResults", lambda *a, **k: (items, []))
        monkeypatch.setattr(PushButton, "POOL_SIZE", 2)
        errors = iter([RuntimeError("first"), ValueError("second")])

        def boom(config, item, registry):
            raise next(errors)

        monkeypatch.setattr(PushButton, "runOneGroup", boom)
        with (
            caplog.at_level("CRITICAL"),
            pytest.raises(PushButton.GroupDispatchError) as excinfo,
        ):
            PushButton.runFlowcell(make_config(), self._parkourDict())

        criticals = [
            r.getMessage() for r in caplog.records if r.levelname == "CRITICAL"
        ]
        failures = excinfo.value.failures
        assert len(failures) >= 1
        assert len(criticals) == len(failures)
        for item, exc in failures:
            assert any(repr(exc) in m for m in criticals)
            assert item.project in str(excinfo.value)

    def test_cancel_happens_before_pool_shutdown(self, tmp_path, monkeypatch):
        order = []
        monkeypatch.setattr(
            PushButton.BRB.jobtrack,
            "cancelAllGroups",
            lambda registry, **kw: order.append("cancel"),
        )
        realExecutor = PushButton.ThreadPoolExecutor

        class SpyExecutor(realExecutor):
            def shutdown(self, wait=True, cancel_futures=False):
                order.append("shutdown")
                return super().shutdown(wait=wait, cancel_futures=cancel_futures)

        monkeypatch.setattr(PushButton, "ThreadPoolExecutor", SpyExecutor)
        item = make_item(tmp_path)
        monkeypatch.setattr(PushButton, "GetResults", lambda *a, **k: ([item], []))
        monkeypatch.setattr(
            PushButton,
            "runOneGroup",
            lambda *a, **k: (_ for _ in ()).throw(RuntimeError("x")),
        )
        with pytest.raises(PushButton.GroupDispatchError):
            PushButton.runFlowcell(make_config(), self._parkourDict())
        assert order == ["cancel", "shutdown"]

    def test_pool_is_shut_down_without_waiting_on_crash(self, tmp_path, monkeypatch):
        recorded = {}
        realExecutor = PushButton.ThreadPoolExecutor

        class SpyExecutor(realExecutor):
            def shutdown(self, wait=True, cancel_futures=False):
                recorded["wait"] = wait
                recorded["cancel_futures"] = cancel_futures
                return super().shutdown(wait=False, cancel_futures=True)

        monkeypatch.setattr(PushButton, "ThreadPoolExecutor", SpyExecutor)
        monkeypatch.setattr(
            PushButton.BRB.jobtrack, "cancelAllGroups", lambda registry, **kw: None
        )
        item = make_item(tmp_path)
        monkeypatch.setattr(PushButton, "GetResults", lambda *a, **k: ([item], []))
        monkeypatch.setattr(
            PushButton,
            "runOneGroup",
            lambda *a, **k: (_ for _ in ()).throw(RuntimeError("x")),
        )
        with pytest.raises(PushButton.GroupDispatchError):
            PushButton.runFlowcell(make_config(), self._parkourDict())
        assert recorded == {"wait": False, "cancel_futures": True}

    def test_clean_run_does_not_cancel(self, tmp_path, monkeypatch):
        cancelled = []
        monkeypatch.setattr(
            PushButton.BRB.jobtrack,
            "cancelAllGroups",
            lambda registry, **kw: cancelled.append(registry),
        )
        item = make_item(tmp_path)
        monkeypatch.setattr(PushButton, "GetResults", lambda *a, **k: ([item], []))
        monkeypatch.setattr(PushButton, "runOneGroup", lambda *a, **k: ["ok"])
        msg = PushButton.runFlowcell(make_config(), self._parkourDict())
        assert msg == ["ok"]
        assert cancelled == []


class TestPoolSize:
    def test_default_is_the_module_constant(self):
        assert PushButton.poolSize(make_config()) == PushButton.POOL_SIZE
        assert PushButton.POOL_SIZE == 2

    def test_ini_value_wins(self):
        config = make_config({"Options": {"poolSize": "4"}})
        assert PushButton.poolSize(config) == 4

    def test_per_sequencer_run_halves_the_configured_value(self):
        config = make_config({"Options": {"poolSize": "4"}})
        assert PushButton.poolSize(config, sequencer="illumina") == 2

    def test_per_sequencer_never_drops_below_one(self):
        config = make_config({"Options": {"poolSize": "1"}})
        assert PushButton.poolSize(config, sequencer="aviti") == 1

    def test_garbage_falls_back_to_the_constant(self):
        config = make_config({"Options": {"poolSize": "lots"}})
        assert PushButton.poolSize(config) == PushButton.POOL_SIZE

    def test_zero_or_negative_falls_back_to_the_constant(self):
        assert (
            PushButton.poolSize(make_config({"Options": {"poolSize": "0"}}))
            == PushButton.POOL_SIZE
        )
        assert (
            PushButton.poolSize(make_config({"Options": {"poolSize": "-3"}}))
            == PushButton.POOL_SIZE
        )

    def test_runflowcell_honours_maxworkers(self, tmp_path, monkeypatch):
        widths = []
        realExecutor = PushButton.ThreadPoolExecutor
        monkeypatch.setattr(
            PushButton,
            "ThreadPoolExecutor",
            lambda max_workers: (
                widths.append(max_workers) or realExecutor(max_workers=max_workers)
            ),
        )
        monkeypatch.setattr(
            PushButton, "GetResults", lambda *a, **k: ([make_item(tmp_path)], [])
        )
        monkeypatch.setattr(os.path, "exists", lambda p: True)
        monkeypatch.setattr(PushButton, "runOneGroup", lambda *a, **k: ["ok"])
        PushButton.runFlowcell(make_config(), {"1_Doe_Smith": {}}, None, maxWorkers=3)
        assert widths == [3]
