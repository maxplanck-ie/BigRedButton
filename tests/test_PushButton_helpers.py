import configparser
import os

from BRB import PushButton


def make_config(overrides=None):
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
    for section, kv in (overrides or {}).items():
        if section not in config:
            config[section] = {}
        config[section].update(kv)
    return config


class TestCreatePath:
    def test_external_project_roots_under_basedata(self, tmp_path):
        config = make_config({"Paths": {"baseData": str(tmp_path)}})
        tuples = [["lib1", "s1", "WGS", True]]

        outputDir = PushButton.createPath(config, "grp", "Proj", "h38", "WGS", tuples)

        assert outputDir == str(
            tmp_path
            / "20260101_AV999999_1234567_lanes_1_2"
            / "Analysis_Proj"
            / "WGS_h38"
        )
        assert os.path.isdir(outputDir)

    def test_internal_project_roots_under_groupdata_seqdir(self, tmp_path, monkeypatch):
        config = make_config({"Paths": {"groupData": str(tmp_path)}})
        (tmp_path / "grp").mkdir()
        tuples = [["lib1", "s1", "WGS", False]]

        outputDir = PushButton.createPath(config, "grp", "Proj", "h38", "WGS", tuples)

        assert outputDir == str(
            tmp_path
            / "grp"
            / "sequencing_data"
            / "20260101_AV999999_1234567_lanes_1_2"
            / "Analysis_Proj"
            / "WGS_h38"
        )
        assert os.path.isdir(outputDir)


class TestLinkFiles:
    def test_external_project_links_both_mates_and_reports_paired_end(self, tmp_path):
        config = make_config({"Paths": {"baseData": str(tmp_path)}})
        runID = config.get("Options", "runID")
        sampleDir = tmp_path / runID / "Project_Proj" / "Sample_s1"
        sampleDir.mkdir(parents=True)
        (sampleDir / "lib1_R1.fastq.gz").write_text("r1")
        (sampleDir / "lib1_R2.fastq.gz").write_text("r2")
        odir = tmp_path / "out"
        odir.mkdir()
        tuples = [["s1", "lib1", "WGS", True]]

        PE = PushButton.linkFiles(config, "grp", "Proj", str(odir), tuples)

        assert PE is True
        assert (odir / "lib1_R1.fastq.gz").is_symlink()
        assert (odir / "lib1_R2.fastq.gz").is_symlink()

    def test_internal_project_links_from_groupdata_seqdir(self, tmp_path):
        config = make_config({"Paths": {"groupData": str(tmp_path)}})
        (tmp_path / "grp").mkdir()
        runID = config.get("Options", "runID")
        sampleDir = (
            tmp_path / "grp" / "sequencing_data" / runID / "Project_Proj" / "Sample_s1"
        )
        sampleDir.mkdir(parents=True)
        (sampleDir / "lib1_R1.fastq.gz").write_text("r1")
        odir = tmp_path / "out"
        odir.mkdir()
        tuples = [["s1", "lib1", "WGS", False]]

        PE = PushButton.linkFiles(config, "grp", "Proj", str(odir), tuples)

        assert PE is False  # only R1 present
        assert (odir / "lib1_R1.fastq.gz").is_symlink()

    def test_existing_symlink_is_not_recreated(self, tmp_path):
        config = make_config({"Paths": {"baseData": str(tmp_path)}})
        runID = config.get("Options", "runID")
        sampleDir = tmp_path / runID / "Project_Proj" / "Sample_s1"
        sampleDir.mkdir(parents=True)
        target = sampleDir / "lib1_R1.fastq.gz"
        target.write_text("r1")
        odir = tmp_path / "out"
        odir.mkdir()
        existingLink = odir / "lib1_R1.fastq.gz"
        existingLink.symlink_to(target)
        tuples = [["s1", "lib1", "WGS", True]]

        # Must not raise FileExistsError.
        PushButton.linkFiles(config, "grp", "Proj", str(odir), tuples)
        assert existingLink.is_symlink()


class TestRemoveLinkFiles:
    def test_removes_symlinks_from_originalfastq_and_top_level(self, tmp_path):
        d = tmp_path / "outdir"
        d.mkdir()
        elsewhere = tmp_path / "elsewhere"
        elsewhere.mkdir()
        target = elsewhere / "real.fastq.gz"
        target.write_text("data")
        originalFastq = d / "originalFASTQ"
        originalFastq.mkdir()
        (originalFastq / "lib1_R1.fastq.gz").symlink_to(target)
        (d / "lib1_R2.fastq.gz").symlink_to(target)

        PushButton.removeLinkFiles(str(d))

        assert not (originalFastq / "lib1_R1.fastq.gz").exists()
        assert not (d / "lib1_R2.fastq.gz").exists()
        assert target.exists()  # the real file is untouched

    def test_no_symlinks_present_does_not_raise(self, tmp_path):
        PushButton.removeLinkFiles(str(tmp_path))  # must not raise


class TestRelinkFiles:
    def test_copies_multiqc_report_when_present(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        (outputDir / "multiQC").mkdir(parents=True)
        (outputDir / "multiQC" / "multiqc_report.html").write_text("<html></html>")
        bioinfoCoreDir = tmp_path / "bioinfocore"
        bioinfoCoreDir.mkdir()
        config = make_config({"Paths": {"bioinfoCoreDir": str(bioinfoCoreDir)}})
        monkeypatch.setattr(PushButton, "createPath", lambda *a, **k: str(outputDir))
        monkeypatch.setattr(PushButton, "linkFiles", lambda *a, **k: None)
        tuples = [["s1", "lib1", "WGS", False]]

        PushButton.relinkFiles(config, "grp", "Proj", "h38", "WGS", tuples)

        expected = bioinfoCoreDir / "AnalysisProj_WGS_h38_multiqc.html"
        assert expected.exists()

    def test_missing_multiqc_report_is_logged_and_skipped(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        outputDir.mkdir()
        bioinfoCoreDir = tmp_path / "bioinfocore"
        bioinfoCoreDir.mkdir()
        config = make_config({"Paths": {"bioinfoCoreDir": str(bioinfoCoreDir)}})
        monkeypatch.setattr(PushButton, "createPath", lambda *a, **k: str(outputDir))
        monkeypatch.setattr(PushButton, "linkFiles", lambda *a, **k: None)
        tuples = [["s1", "lib1", "WGS", False]]

        PushButton.relinkFiles(config, "grp", "Proj", "h38", "WGS", tuples)

        assert list(bioinfoCoreDir.iterdir()) == []

    def test_external_project_skips_multiqc_copy_entirely(self, tmp_path, monkeypatch):
        outputDir = tmp_path / "out"
        (outputDir / "multiQC").mkdir(parents=True)
        (outputDir / "multiQC" / "multiqc_report.html").write_text("<html></html>")
        bioinfoCoreDir = tmp_path / "bioinfocore"
        bioinfoCoreDir.mkdir()
        config = make_config({"Paths": {"bioinfoCoreDir": str(bioinfoCoreDir)}})
        monkeypatch.setattr(PushButton, "createPath", lambda *a, **k: str(outputDir))
        monkeypatch.setattr(PushButton, "linkFiles", lambda *a, **k: None)
        tuples = [["s1", "lib1", "WGS", True]]

        PushButton.relinkFiles(config, "grp", "Proj", "h38", "WGS", tuples)

        assert list(bioinfoCoreDir.iterdir()) == []


class TestGetSambaPath:
    def test_illumina_uses_two_digit_year_prefix(self):
        current_year, year_postfix = PushButton.getsambaPath(
            "210608_A00931", "Illumina"
        )
        assert current_year == "2021"
        assert str(year_postfix) == "Sequence_Quality_2021/Illumina_2021"

    def test_aviti_known_machine_maps_to_facility_name(self):
        current_year, year_postfix = PushButton.getsambaPath(
            "2026_AV251009_lanes_1", "Aviti"
        )
        assert current_year == "2026"
        assert str(year_postfix) == "Sequence_Quality_2026/AVITI24_2026"

    def test_aviti_unknown_machine_falls_back_to_serial(self, caplog):
        _current_year, year_postfix = PushButton.getsambaPath(
            "2026_AVUNKNOWN_lanes_1", "Aviti"
        )
        assert str(year_postfix) == "Sequence_Quality_2026/AVUNKNOWN_2026"


class TestCopyCellRanger:
    def test_ignore_true_skips_entirely(self, tmp_path):
        config = make_config()
        PushButton.copyCellRanger(config, str(tmp_path), ignore=True)
        # nothing to assert on filesystem; must simply not raise / not glob

    def test_copies_web_summary_to_seqfac_and_bioinfocore(self, tmp_path):
        seqFacDir = tmp_path / "seqfac"
        bioinfoCoreDir = tmp_path / "bioinfocore"
        seqFacDir.mkdir()
        bioinfoCoreDir.mkdir()
        config = make_config(
            {
                "Paths": {
                    "seqFacDir": str(seqFacDir),
                    "bioinfoCoreDir": str(bioinfoCoreDir),
                }
            }
        )
        # d = .../<lane_dir>/Analysis_Proj/RNA-Seq_human ; lane_dir is
        # Path(d).parents[1].stem
        laneDir = "210608_A00931_0309_HCCMWDRXY"
        d = tmp_path / laneDir / "Analysis_Proj" / "RNA-Seq_human"
        sampleOuts = d / "sample1" / "outs"
        sampleOuts.mkdir(parents=True)
        (sampleOuts / "web_summary.html").write_text("<html></html>")

        PushButton.copyCellRanger(config, str(d), ignore=False)

        seqfacFiles = list((seqFacDir).rglob("*.html"))
        bioinfoFiles = list(bioinfoCoreDir.glob("*.html"))
        assert len(seqfacFiles) == 1
        assert len(bioinfoFiles) == 1


class TestTidyUpABit:
    def test_removes_clusters_logs_and_snakemake_dirs(self, tmp_path):
        (tmp_path / "clusters_logs").mkdir()
        (tmp_path / "clusters_logs" / "log.txt").write_text("x")
        (tmp_path / ".snakemake").mkdir()

        PushButton.tidyUpABit(str(tmp_path))

        assert not (tmp_path / "clusters_logs").exists()
        assert not (tmp_path / ".snakemake").exists()

    def test_missing_files_do_not_raise(self, tmp_path):
        PushButton.tidyUpABit(str(tmp_path))  # nothing exists -- must not raise

    def test_removes_config_yaml_and_multiqc_logs(self, tmp_path):
        (tmp_path / "config.yaml").write_text("x")
        multiqcDir = tmp_path / "multiQC" / "multiqc_data"
        multiqcDir.mkdir(parents=True)
        (multiqcDir / "multiqc.log").write_text("x")

        PushButton.tidyUpABit(str(tmp_path))

        assert not (tmp_path / "config.yaml").exists()
        assert not (multiqcDir / "multiqc.log").exists()


class TestRemoveDone:
    def test_removes_existing_marker(self, tmp_path):
        (tmp_path / "analysis.done").touch()
        PushButton.removeDone(str(tmp_path))
        assert not (tmp_path / "analysis.done").exists()

    def test_missing_marker_does_not_raise(self, tmp_path):
        PushButton.removeDone(str(tmp_path))
