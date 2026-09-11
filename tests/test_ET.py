import configparser
import gzip
import threading
import time

import BRB.ET


def parkourConfig(tmp_path):
    config = configparser.ConfigParser()
    config["Paths"] = {"baseData": str(tmp_path)}
    config["Options"] = {"runID": "210608_A00931_0309_BHCCMWDRXY"}
    config["Parkour"] = {
        "ResultsURL": "https://parkour.example.org/api/results/",
        "user": "u",
        "password": "p",
        "cert": "",
    }
    return config


class TestSendToParkourSerialisation:
    def test_concurrent_posts_are_serialised(self, tmp_path, monkeypatch):
        config = parkourConfig(tmp_path)
        state = {"inFlight": 0, "maxInFlight": 0, "posts": 0}
        guard = threading.Lock()

        def fakePost(url, auth=None, data=None, verify=None):
            with guard:
                state["inFlight"] += 1
                state["posts"] += 1
                state["maxInFlight"] = max(state["maxInFlight"], state["inFlight"])
            time.sleep(0.05)
            with guard:
                state["inFlight"] -= 1
            return "RESPONSE"

        monkeypatch.setattr(BRB.ET.requests, "post", fakePost)

        threads = [
            threading.Thread(
                target=BRB.ET.sendToParkour, args=(config, [{"barcode": str(i)}])
            )
            for i in range(4)
        ]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=10)

        assert state["posts"] == 4
        assert state["maxInFlight"] == 1

    def test_return_value_unchanged(self, tmp_path, monkeypatch):
        config = parkourConfig(tmp_path)
        monkeypatch.setattr(BRB.ET.requests, "post", lambda *a, **k: "RESPONSE")
        assert BRB.ET.sendToParkour(config, [{"barcode": "1"}]) == "RESPONSE"


class TestGetOffSpeciesRate:
    def _writeReport(self, d, lines):
        (d / "sample.rep").write_text("\n".join(lines) + "\n")

    def test_known_label_matches_its_kraken_group(self, tmp_path):
        self._writeReport(
            tmp_path,
            ["97.00\t123\t45\tG1\t0\tflygrp", "1.00\t1\t1\tG1\t0\tother"],
        )
        rate = BRB.ET.getOffSpeciesRate(str(tmp_path), "drosophila")
        assert rate == 1 - 97.00 / 100

    def test_unknown_label_returns_zero_and_does_not_raise(self, tmp_path):
        self._writeReport(tmp_path, ["97.00\t123\t45\tG1\t0\tflygrp"])
        assert BRB.ET.getOffSpeciesRate(str(tmp_path), "ricefish") == 0
        assert BRB.ET.getOffSpeciesRate(str(tmp_path), None) == 0

    def test_known_label_whose_group_is_absent_from_report_returns_zero(self, tmp_path):
        # regression: org_map[org_label] never matches a line in the report,
        # `off` used to stay unbound and raise UnboundLocalError instead of
        # returning a sane default.
        self._writeReport(tmp_path, ["100.00\t123\t45\tG1\t0\tsomeothergrp"])
        assert BRB.ET.getOffSpeciesRate(str(tmp_path), "drosophila") == 0


class TestGetNReads:
    def test_duplicate_file_present_computes_optical_dupe_fraction(self, tmp_path):
        (tmp_path / "sample.duplicate.txt").write_text("10 100\n")
        total, optFrac = BRB.ET.getNReads(str(tmp_path))
        assert total == 90
        assert optFrac == 10.0

    def test_duplicate_file_with_zero_total_returns_minus_one_fraction(self, tmp_path):
        (tmp_path / "sample.duplicate.txt").write_text("0 0\n")
        total, optFrac = BRB.ET.getNReads(str(tmp_path))
        assert total == 0
        assert optFrac == -1.0

    def test_falls_back_to_counting_fastq_lines_when_no_duplicate_file(self, tmp_path):
        # getNReads strips a leading "FASTQC_" from the dir name when no
        # *.duplicate.txt exists, and looks for *_R1.fastq.gz there instead.
        fastqDir = tmp_path / "Sample_S1"
        fastqDir.mkdir()
        with gzip.open(fastqDir / "S1_R1.fastq.gz", "wt") as fh:
            fh.write("line\n" * 40)  # 10 reads x 4 lines/read

        d = str(tmp_path / "FASTQC_Sample_S1")
        total, optFrac = BRB.ET.getNReads(d)

        assert total == 10
        assert optFrac == 0.0


class TestGetBaseStatistics:
    def _buildSampleDir(
        self, config, project, sample, dupeText="10 100", repLines=None
    ):
        baseData = config.get("Paths", "baseData")
        runID = config.get("Options", "runID")
        sampleDir = f"{baseData}/{runID}/FASTQC_Project_{project}/Sample_{sample}"
        import os

        os.makedirs(sampleDir, exist_ok=True)
        with open(f"{sampleDir}/{sample}_R1_fastqc.zip", "w") as fh:
            fh.write("")
        with open(f"{sampleDir}/{sample}.duplicate.txt", "w") as fh:
            fh.write(dupeText)
        with open(f"{sampleDir}/{sample}.rep", "w") as fh:
            fh.write("\n".join(repLines or ["100.00\t1\t1\tG1\t0\tnone"]) + "\n")
        return sampleDir

    def test_parses_known_sample_into_base_dict(self, tmp_path):
        config = parkourConfig(tmp_path)
        self._buildSampleDir(config, "ProjectX", "S1")
        outputDir = (
            f"{tmp_path}/{config.get('Options', 'runID')}/Analysis_ProjectX/DNA_human"
        )

        baseDict, s2l = BRB.ET.getBaseStatistics(config, outputDir, ["S1"])

        assert baseDict["S1"][0] == "S1"  # sampleName
        assert baseDict["S1"][1] == 90  # nReads (100 - 10 dupes)
        assert s2l["S1"] == "S1"

    def test_sample_without_fastqc_zip_is_skipped(self, tmp_path):
        config = parkourConfig(tmp_path)
        import os

        sampleDir = (
            f"{tmp_path}/{config.get('Options', 'runID')}"
            "/FASTQC_Project_ProjectX/Sample_S1"
        )
        os.makedirs(sampleDir, exist_ok=True)  # no *_R1_fastqc.zip written
        outputDir = (
            f"{tmp_path}/{config.get('Options', 'runID')}/Analysis_ProjectX/DNA_human"
        )

        baseDict, s2l = BRB.ET.getBaseStatistics(config, outputDir, ["S1"])

        assert baseDict == {}
        assert s2l == {}


class TestSendToParkourFlowcellId:
    def test_illumina_run_extracts_fcid_from_fourth_field(self, tmp_path, monkeypatch):
        config = parkourConfig(tmp_path)  # runID 210608_A00931_0309_BHCCMWDRXY
        seen = {}
        monkeypatch.setattr(
            BRB.ET.requests,
            "post",
            lambda url, auth=None, data=None, verify=None: seen.update(data) or "OK",
        )
        BRB.ET.sendToParkour(config, [{"barcode": "1"}])
        assert seen["flowcell_id"] == "HCCMWDRXY"

    def test_aviti_run_extracts_fcid_from_run_manifest(self, tmp_path, monkeypatch):
        config = parkourConfig(tmp_path)
        config.set("Options", "runID", "20260101_AV999999_1234-5-6")
        (tmp_path / "20260101_AV999999_1234-5-6").mkdir()
        (tmp_path / "20260101_AV999999_1234-5-6" / "RunManifest.csv").write_text("")
        seen = {}
        monkeypatch.setattr(
            BRB.ET.requests,
            "post",
            lambda url, auth=None, data=None, verify=None: seen.update(data) or "OK",
        )
        BRB.ET.sendToParkour(config, [{"barcode": "1"}])
        assert seen["flowcell_id"] == "6"

    def test_illumina_run_with_dash_in_fcid_keeps_only_the_suffix(
        self, tmp_path, monkeypatch
    ):
        config = parkourConfig(tmp_path)
        config.set("Options", "runID", "210608_A00931_0309_BAAA-HCCMWDRXY")
        seen = {}
        monkeypatch.setattr(
            BRB.ET.requests,
            "post",
            lambda url, auth=None, data=None, verify=None: seen.update(data) or "OK",
        )
        BRB.ET.sendToParkour(config, [{"barcode": "1"}])
        assert seen["flowcell_id"] == "HCCMWDRXY"


class TestDNA:
    def test_relacs_output_returns_base_fields_only(self, tmp_path):
        config = parkourConfig(tmp_path)
        outputDir = tmp_path / "DNA_human"
        outputDir.mkdir()
        (outputDir / "RELACS_sampleSheet.txt").write_text("L1\tbarcode\n")
        baseDict = {"L1": ["S1", 90, 0.05, 10.0]}

        m = BRB.ET.DNA(config, str(outputDir), baseDict, {"S1": "L1"})

        assert m == [
            {
                "barcode": "L1",
                "reads_pf_sequenced": 90,
                "confident_reads": 0.05,
                "optical_duplicates": 10.0,
            }
        ]

    def test_non_relacs_output_adds_mapping_dup_and_insert_size(self, tmp_path):
        import os

        config = parkourConfig(tmp_path)
        outputDir = tmp_path / "DNA_human"
        bowtieDir = outputDir / "Bowtie2"
        multiqcDir = outputDir / "multiQC" / "multiqc_data"
        deeptoolsDir = outputDir / "deepTools_qc" / "bamPEFragmentSize"
        for d in (bowtieDir, multiqcDir, deeptoolsDir):
            os.makedirs(d, exist_ok=True)

        (bowtieDir / "S1.Bowtie2_summary.txt").write_text(
            "some header\n95.50% overall alignment rate\n"
        )
        (multiqcDir / "multiqc_samtools_flagstat.txt").write_text(
            "Sample\ttotal_passed\tduplicates_passed\nS1\t100\t10\n"
        )
        (deeptoolsDir / "fragmentSize.metric.tsv").write_text(
            "Unnamed: 0\tFrag. Len. Median\nfiltered_bam/S1.filtered.bam\t250\n"
        )

        baseDict = {"L1": ["S1", 90, 0.05, 10.0]}
        m = BRB.ET.DNA(config, str(outputDir), baseDict, {"S1": "L1"})

        assert m == [
            {
                "barcode": "L1",
                "reads_pf_sequenced": 90,
                "confident_reads": 0.05,
                "optical_duplicates": 10.0,
                "dupped_reads": 10.0,
                "mapped_reads": 95.50,
                "insert_size": 250,
            }
        ]

    def test_sample_missing_from_dup_df_defaults_dup_rate_to_zero(self, tmp_path):
        import os

        config = parkourConfig(tmp_path)
        outputDir = tmp_path / "DNA_human"
        bowtieDir = outputDir / "Bowtie2"
        multiqcDir = outputDir / "multiQC" / "multiqc_data"
        deeptoolsDir = outputDir / "deepTools_qc" / "bamPEFragmentSize"
        for d in (bowtieDir, multiqcDir, deeptoolsDir):
            os.makedirs(d, exist_ok=True)

        (bowtieDir / "S1.Bowtie2_summary.txt").write_text(
            "some header\n95.50% overall alignment rate\n"
        )
        (multiqcDir / "multiqc_samtools_flagstat.txt").write_text(
            "Sample\ttotal_passed\tduplicates_passed\nOTHER\t100\t10\n"
        )
        (deeptoolsDir / "fragmentSize.metric.tsv").write_text(
            "Unnamed: 0\tFrag. Len. Median\nfiltered_bam/S1.filtered.bam\t250\n"
        )

        baseDict = {"L1": ["S1", 90, 0.05, 10.0]}
        m = BRB.ET.DNA(config, str(outputDir), baseDict, {"S1": "L1"})

        assert m[0]["dupped_reads"] == 0


class TestRNA:
    def test_parses_star_multiqc_and_featurecounts_into_message(self, tmp_path):
        import os

        config = parkourConfig(tmp_path)
        outputDir = tmp_path / "RNA_human"
        starDir = outputDir / "STAR" / "S1"
        multiqcDir = outputDir / "multiQC" / "multiqc_data"
        for d in (starDir, multiqcDir):
            os.makedirs(d, exist_ok=True)

        (starDir / "S1.Log.final.out").write_text(
            "Uniquely mapped reads %\t85.00%\n"
            "% of reads mapped to multiple loci\t5.00%\n"
        )
        (multiqcDir / "multiqc_samtools_flagstat.txt").write_text(
            "Sample\ttotal_passed\tduplicates_passed\nS1\t100\t10\n"
        )
        (multiqcDir / "multiqc_featurecounts.txt").write_text(
            "Sample\tTotal\tAssigned\nS1.filtered\t100\t80\n"
        )

        baseDict = {"L1": ["S1", 90, 0.05, 10.0]}
        m = BRB.ET.RNA(config, str(outputDir), baseDict, {"S1": "L1"})

        assert m == [
            {
                "barcode": "L1",
                "reads_pf_sequenced": 90,
                "confident_reads": 0.05,
                "optical_duplicates": 10.0,
                "mapped_reads": 90.0,
                "uniq_mapped": 85.0,
                "multi_mapped": 5.0,
                "dupped_reads": 10.0,
                "assigned_reads": 80.0,
            }
        ]


class TestPhoneHome:
    def test_dna_pipeline_dispatches_to_dna_module_and_sends_to_parkour(
        self, tmp_path, monkeypatch
    ):
        config = parkourConfig(tmp_path)
        monkeypatch.setattr(
            BRB.ET,
            "getBaseStatistics",
            lambda *a, **k: ({"L1": ["S1", 1, 2, 3]}, {"S1": "L1"}),
        )
        monkeypatch.setattr(BRB.ET, "DNA", lambda *a, **k: [{"barcode": "L1"}])
        sent = {}

        def fakeSend(c, m):
            sent["msg"] = m
            return "OK"

        monkeypatch.setattr(BRB.ET, "sendToParkour", fakeSend)

        result = BRB.ET.phoneHome(
            config,
            "outdir",
            "DNA",
            [("S1", "libtype")],
            "human",
            "h38",
            "ProjectX",
            "WGS",
        )

        assert sent["msg"] == [{"barcode": "L1"}]
        assert result == ["ProjectX", "human", "WGS", "DNA", "success", "OK"]

    def test_rna_pipeline_dispatches_to_rna_module(self, tmp_path, monkeypatch):
        config = parkourConfig(tmp_path)
        monkeypatch.setattr(
            BRB.ET,
            "getBaseStatistics",
            lambda *a, **k: ({"L1": ["S1", 1, 2, 3]}, {"S1": "L1"}),
        )
        monkeypatch.setattr(BRB.ET, "RNA", lambda *a, **k: [{"barcode": "L1"}])
        monkeypatch.setattr(BRB.ET, "sendToParkour", lambda c, m: m)

        result = BRB.ET.phoneHome(
            config, "outdir", "RNA", [("S1", "libtype")], "human", "h38", "P", "mRNA"
        )

        assert result == ["P", "human", "mRNA", "RNA", "success", [{"barcode": "L1"}]]

    def test_other_pipeline_uses_base_fields_only(self, tmp_path, monkeypatch):
        config = parkourConfig(tmp_path)
        monkeypatch.setattr(
            BRB.ET,
            "getBaseStatistics",
            lambda *a, **k: ({"L1": ["S1", 1, 2, 3]}, {"S1": "L1"}),
        )
        monkeypatch.setattr(BRB.ET, "sendToParkour", lambda c, m: m)

        result = BRB.ET.phoneHome(
            config, "outdir", "ATAC", [("S1", "libtype")], "human", "h38", "P", "ATAC"
        )

        assert result[5] == [
            {
                "barcode": "L1",
                "reads_pf_sequenced": 1,
                "confident_reads": 2,
                "optical_duplicates": 3,
            }
        ]


class TestTelegraphHome:
    def test_builds_matrix_and_returns_skip_report(self, tmp_path, monkeypatch):
        config = parkourConfig(tmp_path)
        config["Paths"]["groupData"] = str(tmp_path / "groupdata")
        import os

        piDir = tmp_path / "groupdata" / "PIName"
        os.makedirs(piDir / "sequencing_data", exist_ok=True)

        monkeypatch.setattr(
            BRB.ET,
            "getBaseStatistics",
            lambda *a, **k: ({"L1": ["S1", 1, 2, 3]}, {"S1": "L1"}),
        )
        monkeypatch.setattr(BRB.ET, "sendToParkour", lambda c, m: m)

        skipList = [("S1", "S1", "WGS")]
        result = BRB.ET.telegraphHome(
            config, "PIName", "ProjectX", skipList, organism="human", org_label="h38"
        )

        assert result[0] == "ProjectX"
        assert result[1] == "human"
        assert result[2] == "WGS"
        assert result[3] is None
        assert result[4] is None
        assert result[5] == [
            {
                "barcode": "L1",
                "reads_pf_sequenced": 1,
                "confident_reads": 2,
                "optical_duplicates": 3,
            }
        ]
        assert result[6] is False
