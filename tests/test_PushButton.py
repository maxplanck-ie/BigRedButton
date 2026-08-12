import configparser
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
        with patch("BRB.PushButton.subprocess.check_call") as mock_call:
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
        with patch("BRB.PushButton.subprocess.check_call") as mock_call:
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
        with patch("BRB.PushButton.subprocess.check_call"):
            deliverExternalRELACS(config, str(outputDir), "99_Foo_Bar")
        assert (outputDir / "external_delivery.done").exists()

    def test_second_call_is_a_noop(self, tmp_path):
        outputDir = self._make_output_dir(tmp_path)
        (outputDir / "external_delivery.done").touch()
        config = make_config()
        with patch("BRB.PushButton.subprocess.check_call") as mock_call:
            deliverExternalRELACS(config, str(outputDir), "99_Foo_Bar")
        mock_call.assert_not_called()

    def test_failure_is_logged_not_raised(self, tmp_path):
        outputDir = self._make_output_dir(tmp_path)
        config = make_config()
        with patch("BRB.PushButton.subprocess.check_call", side_effect=OSError("boom")):
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
    @patch("BRB.PushButton.subprocess.check_call")
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
    @patch("BRB.PushButton.subprocess.check_call")
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
