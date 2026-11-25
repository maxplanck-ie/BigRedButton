import configparser
import gzip
import math
import os
import threading
import time
from pathlib import Path

import pytest
import requests

import BRB.ET


def fakeResponse(status=200, reason="OK", body=""):
    res = requests.Response()
    res.status_code = status
    res.reason = reason
    res._content = body.encode()
    return res


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

        def fakePost(url, auth=None, data=None, verify=None, timeout=None):
            with guard:
                state["inFlight"] += 1
                state["posts"] += 1
                state["maxInFlight"] = max(state["maxInFlight"], state["inFlight"])
            time.sleep(0.05)
            with guard:
                state["inFlight"] -= 1
            return fakeResponse()

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
        monkeypatch.setattr(
            BRB.ET.requests, "post", lambda *a, **k: fakeResponse(body="RESPONSE")
        )
        assert BRB.ET.sendToParkour(config, [{"barcode": "1"}]).text == "RESPONSE"


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

    def test_unknown_label_is_nan_and_does_not_raise(self, tmp_path):
        self._writeReport(tmp_path, ["97.00\t123\t45\tG1\t0\tflygrp"])
        assert math.isnan(BRB.ET.getOffSpeciesRate(str(tmp_path), "ricefish"))
        assert math.isnan(BRB.ET.getOffSpeciesRate(str(tmp_path), None))

    def test_known_label_whose_group_is_absent_from_report_is_nan(self, tmp_path):
        # regression: org_map[org_label] never matches a line in the report,
        # `off` used to stay unbound and raise UnboundLocalError. 0 would claim
        # a perfect library, so the rate is unknown (NaN) instead.
        self._writeReport(tmp_path, ["100.00\t123\t45\tG1\t0\tsomeothergrp"])
        assert math.isnan(BRB.ET.getOffSpeciesRate(str(tmp_path), "drosophila"))

    def test_missing_kraken_report_is_nan_and_does_not_raise(self, tmp_path):
        assert math.isnan(BRB.ET.getOffSpeciesRate(str(tmp_path), "drosophila"))


class TestGetNReads:
    def test_duplicate_file_present_computes_optical_dupe_fraction(self, tmp_path):
        (tmp_path / "sample.duplicate.txt").write_text("10 100\n")
        total, optFrac = BRB.ET.getNReads(str(tmp_path))
        assert total == 90
        assert optFrac == 10.0

    def test_duplicate_file_with_zero_total_has_no_fraction(self, tmp_path):
        (tmp_path / "sample.duplicate.txt").write_text("0 0\n")
        total, optFrac = BRB.ET.getNReads(str(tmp_path))
        assert total == 0
        assert math.isnan(optFrac)

    def test_metrics_file_computes_optical_dupe_fraction(self, tmp_path):
        # splitFastq now writes <sample>.metrics instead of <sample>.duplicate.txt
        (tmp_path / "sample.metrics").write_text("144129\t14488041\n")
        total, optFrac = BRB.ET.getNReads(str(tmp_path))
        assert total == 14488041 - 144129
        assert optFrac == pytest.approx(100.0 * 144129 / 14488041)

    def test_metrics_file_left_in_unprefixed_project_dir_is_found(self, tmp_path):
        # runs from before dissectBCL v1.3.5 never had .metrics moved to FASTQC_*
        sampleDir = tmp_path / "Sample_S1"
        sampleDir.mkdir()
        (sampleDir / "S1.metrics").write_text("10 100\n")
        total, optFrac = BRB.ET.getNReads(str(tmp_path / "FASTQC_Sample_S1"))
        assert total == 90
        assert optFrac == 10.0

    def test_no_duplicate_file_counts_fastq_reads_and_rate_is_unknown(self, tmp_path):
        # getNReads strips a leading "FASTQC_" from the dir name when no
        # *.duplicate.txt exists, and looks for *_R1.fastq.gz there instead.
        # No duplicate marking means the rate is unknown, not 0 %.
        fastqDir = tmp_path / "Sample_S1"
        fastqDir.mkdir()
        with gzip.open(fastqDir / "S1_R1.fastq.gz", "wt") as fh:
            fh.write("line\n" * 40)  # 10 reads x 4 lines/read

        d = str(tmp_path / "FASTQC_Sample_S1")
        total, optFrac = BRB.ET.getNReads(d)

        assert total == 10
        assert math.isnan(optFrac)


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
            lambda url, auth=None, data=None, verify=None, timeout=None: (
                seen.update(data) or fakeResponse()
            ),
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
            lambda url, auth=None, data=None, verify=None, timeout=None: (
                seen.update(data) or fakeResponse()
            ),
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
            lambda url, auth=None, data=None, verify=None, timeout=None: (
                seen.update(data) or fakeResponse()
            ),
        )
        BRB.ET.sendToParkour(config, [{"barcode": "1"}])
        assert seen["flowcell_id"] == "HCCMWDRXY"


class TestDNA:
    def test_relacs_output_without_demuxed_sample_files_returns_base_fields_only(
        self, tmp_path
    ):
        config = parkourConfig(tmp_path)
        outputDir = tmp_path / "DNA_human"
        outputDir.mkdir()
        (outputDir / "RELACS_sampleSheet.txt").write_text("Sample_L1\tAAAA\tS1_Pol2\n")
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

    def test_relacs_output_sends_stats_per_demuxed_sample(self, tmp_path):
        import os

        config = parkourConfig(tmp_path)
        outputDir = tmp_path / "DNA_human"
        for d in (
            "Bowtie2",
            "multiQC/multiqc_data",
            "deepTools_qc/bamPEFragmentSize",
            "MACS2_QC",
        ):
            os.makedirs(outputDir / d)
        (outputDir / "RELACS_sampleSheet.txt").write_text(
            "Sample_L1\tAAAA\tUninf_Pol2\n"
            "Sample_L1\tCCCC\tCovid_Pol2\n"
            "Sample_L1\tdefault\tunused\n"
            "Sample_L2\tAAAA\tUninf_Input\n"
        )
        (outputDir / "Bowtie2" / "Uninf_Pol2.Bowtie2_summary.txt").write_text(
            "some header\n95.50% overall alignment rate\n"
        )
        (outputDir / "Bowtie2" / "Covid_Pol2.Bowtie2_summary.txt").write_text(
            "some header\n80.00% overall alignment rate\n"
        )
        # Covid_Pol2 has no flagstat row, Uninf_Pol2 no fragment row
        (outputDir / "multiQC/multiqc_data/multiqc_samtools_flagstat.txt").write_text(
            "Sample\ttotal_passed\tduplicates_passed\tread1_passed\n"
            "Uninf_Pol2.markdup\t100\t10\t60\n"
        )
        (
            outputDir / "deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv"
        ).write_text(
            "\tFrag. Len. Min.\tFrag. Len. 1st. Qu.\tFrag. Len. Mean\t"
            "Frag. Len. Median\tFrag. Len. 3rd Qu.\tFrag. Len. 10%\t"
            "Frag. Len. 90%\n"
            "filtered_bam/Covid_Pol2.filtered.bam\t1\t100\t200\t210\t300\t50\t400\n"
        )
        (outputDir / "MACS2_QC" / "Uninf_Pol2.filtered.BAM_peaks.qc.txt").write_text(
            "peak_count\tFRiP\tpeak_genome_coverage\n12\t0.25\t0.01\n"
        )
        baseDict = {"L1": ["S1", 90, 0.05, 10.0], "L2": ["S2", 50, 0.01, 1.0]}

        m = BRB.ET.DNA(config, str(outputDir), baseDict, {"S1": "L1", "S2": "L2"})

        assert m == [
            {
                "barcode": "L1",
                "reads_pf_sequenced": 90,
                "confident_reads": 0.05,
                "optical_duplicates": 10.0,
                "dupped_reads": 10.0,
                "fragments_total": 60,
                # Covid_Pol2 has no flagstat row: no weight, so it is left out
                "mapped_reads": 95.5,
                "frip": 0.25,
                "demuxed_samples": {
                    "Uninf_Pol2": {
                        "mapped_reads": 95.5,
                        "dupped_reads": 10.0,
                        "fragments_total": 60,
                        "peak_count": 12,
                        "frip": 0.25,
                        "peak_genome_coverage": 0.01,
                    },
                    "Covid_Pol2": {
                        "mapped_reads": 80.0,
                        "fragment_q1": 100,
                        "fragment_median": 210,
                        "fragment_q3": 300,
                        "fragment_p10": 50,
                        "fragment_p90": 400,
                        "insert_size": 210,
                    },
                },
            },
            {
                "barcode": "L2",
                "reads_pf_sequenced": 50,
                "confident_reads": 0.01,
                "optical_duplicates": 1.0,
            },
        ]

    def test_relacs_peaks_read_from_macs2_dir(self, tmp_path):
        import os

        config = parkourConfig(tmp_path)
        outputDir = tmp_path / "DNA_human"
        for d in ("Bowtie2", "multiQC/multiqc_data", "MACS2"):
            os.makedirs(outputDir / d)
        (outputDir / "RELACS_sampleSheet.txt").write_text(
            "Sample_L1\tAAAA\tUninf_Pol2\n"
        )
        (outputDir / "Bowtie2" / "Uninf_Pol2.Bowtie2_summary.txt").write_text(
            "some header\n95.50% overall alignment rate\n"
        )
        (outputDir / "MACS2" / "Uninf_Pol2.filtered.BAM_peaks.qc.txt").write_text(
            "peak_count\tFRiP\tpeak_genome_coverage\n12\t0.25\t0.01\n"
        )
        baseDict = {"L1": ["S1", 90, 0.05, 10.0]}

        m = BRB.ET.DNA(config, str(outputDir), baseDict, {"S1": "L1"})

        assert m[0]["demuxed_samples"] == {
            "Uninf_Pol2": {
                "mapped_reads": 95.5,
                "peak_count": 12,
                "frip": 0.25,
                "peak_genome_coverage": 0.01,
            }
        }

    def test_relacs_four_column_sheet_with_spaces_in_the_label(self, tmp_path):
        import os

        config = parkourConfig(tmp_path)
        outputDir = tmp_path / "DNA_mouse"
        for d in ("Bowtie2", "multiQC/multiqc_data"):
            os.makedirs(outputDir / d)
        # sample, bc position, barcode, label; files carry the label with _ for spaces
        (outputDir / "RELACS_sampleSheet.txt").write_text(
            "Sample_L1\tA1-21\tTACCGATG\tPellet 1: Rep1 INPUT\n"
        )
        (outputDir / "Bowtie2" / "Pellet_1:_Rep1_INPUT.Bowtie2_summary.txt").write_text(
            "some header\n95.50% overall alignment rate\n"
        )
        (outputDir / "multiQC/multiqc_data/multiqc_samtools_flagstat.txt").write_text(
            "Sample\ttotal_passed\tduplicates_passed\tread1_passed\n"
            "Pellet_1:_Rep1_INPUT\t100\t10\t60\n"
        )
        baseDict = {"L1": ["S1", 90, 0.05, 10.0]}

        m = BRB.ET.DNA(config, str(outputDir), baseDict, {"S1": "L1"})

        assert m[0]["demuxed_samples"] == {
            "Pellet_1:_Rep1_INPUT": {
                "mapped_reads": 95.5,
                "dupped_reads": 10.0,
                "fragments_total": 60,
            }
        }

    def test_relacs_demuxed_sample_with_a_broken_flagstat_row_keeps_the_others(
        self, tmp_path
    ):
        import os

        config = parkourConfig(tmp_path)
        outputDir = tmp_path / "DNA_mouse"
        for d in ("Bowtie2", "multiQC/multiqc_data"):
            os.makedirs(outputDir / d)
        (outputDir / "RELACS_sampleSheet.txt").write_text(
            "Sample_L1\tA1-21\tTACCGATG\tD1\nSample_L1\tA1-22\tTACCGATC\tD2\n"
        )
        for name in ("D1", "D2"):
            (outputDir / "Bowtie2" / f"{name}.Bowtie2_summary.txt").write_text(
                "some header\n90.00% overall alignment rate\n"
            )
        (outputDir / "multiQC/multiqc_data/multiqc_samtools_flagstat.txt").write_text(
            "Sample\ttotal_passed\tduplicates_passed\tread1_passed\n"
            "D1\t100\t10\t60\n"
            "D2\tnot_a_number\t10\t40\n"
        )
        baseDict = {"L1": ["S1", 90, 0.05, 10.0]}

        m = BRB.ET.DNA(config, str(outputDir), baseDict, {"S1": "L1"})

        assert m[0]["demuxed_samples"] == {
            "D1": {"mapped_reads": 90.0, "dupped_reads": 10.0, "fragments_total": 60},
            "D2": {"mapped_reads": 90.0, "fragments_total": 40},
        }

    def test_relacs_library_values_are_weighted_means_of_the_demuxed_samples(self):
        samples = {
            "A": {
                "dupped_reads": 10.0,
                "mapped_reads": 90.0,
                "fragments_total": 10,
                "fragment_median": 100,
                "frip": 0.1,
                "peak_count": 5,
            },
            "B": {
                "dupped_reads": 30.0,
                "mapped_reads": 60.0,
                "fragments_total": 30,
                "fragment_median": 200,
                "frip": 0.3,
                "peak_count": 7,
            },
            "C": {"dupped_reads": float("nan"), "fragments_total": 5},  # no weight
        }
        weights = {"A": 100, "B": 300}

        out = BRB.ET._relacsLibrary(samples, weights)

        assert out == {
            "fragments_total": 45,  # exact sum
            "dupped_reads": 25.0,  # (100*10 + 300*30) / 400
            "mapped_reads": 67.5,  # (100*90 + 300*60) / 400
            "fragment_median": 175.0,  # (10*100 + 30*200) / 40
            "frip": 0.25,  # (10*0.1 + 30*0.3) / 40
        }  # no peak_count: the samples' peak sets overlap

    def test_dropNaN_cleans_nested_demuxed_samples(self):
        m = [{"barcode": "L1", "demuxed_samples": {"S": {"a": float("nan"), "b": 1}}}]
        assert BRB.ET.dropNaN(m) == [
            {"barcode": "L1", "demuxed_samples": {"S": {"b": 1}}}
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

    def test_sample_missing_from_run_is_skipped(self, tmp_path):
        import os

        config = parkourConfig(tmp_path)
        outputDir = tmp_path / "DNA_human"
        multiqcDir = outputDir / "multiQC" / "multiqc_data"
        fragDir = outputDir / "deepTools_qc" / "bamPEFragmentSize"
        for d in (outputDir / "Bowtie2", multiqcDir, fragDir):
            os.makedirs(d, exist_ok=True)
        for s in ("S1", "S9"):
            (outputDir / "Bowtie2" / f"{s}.Bowtie2_summary.txt").write_text(
                "some header\n95.50% overall alignment rate\n"
            )
        (multiqcDir / "multiqc_samtools_flagstat.txt").write_text(
            "Sample\ttotal_passed\tduplicates_passed\nS1\t100\t10\n"
        )
        (fragDir / "fragmentSize.metric.tsv").write_text(
            "Unnamed: 0\tFrag. Len. Median\nfiltered_bam/S1.filtered.bam\t250\n"
        )

        baseDict = {"L1": ["S1", 90, 0.05, 10.0]}
        m = BRB.ET.DNA(config, str(outputDir), baseDict, {"S1": "L1"})

        assert len(m) == 1
        assert m[0]["barcode"] == "L1"
        assert m[0]["mapped_reads"] == 95.50

    def test_sample_missing_from_dup_df_leaves_dup_rate_unknown(self, tmp_path):
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

        assert math.isnan(m[0]["dupped_reads"])

    def test_sample_without_fragment_row_keeps_its_other_values(self, tmp_path):
        import os

        config = parkourConfig(tmp_path)
        outputDir = tmp_path / "DNA_human"
        multiqcDir = outputDir / "multiQC" / "multiqc_data"
        fragDir = outputDir / "deepTools_qc" / "bamPEFragmentSize"
        for d in (outputDir / "Bowtie2", multiqcDir, fragDir):
            os.makedirs(d, exist_ok=True)
        for s in ("S1", "S2"):
            (outputDir / "Bowtie2" / f"{s}.Bowtie2_summary.txt").write_text(
                "some header\n95.50% overall alignment rate\n"
            )
        (multiqcDir / "multiqc_samtools_flagstat.txt").write_text(
            "Sample\ttotal_passed\tduplicates_passed\nS1\t100\t10\nS2\t100\t20\n"
        )
        (fragDir / "fragmentSize.metric.tsv").write_text(
            "Unnamed: 0\tFrag. Len. Median\nfiltered_bam/S1.filtered.bam\t250\n"
        )

        baseDict = {"L1": ["S1", 90, 0.05, 10.0], "L2": ["S2", 80, 0.04, 8.0]}
        m = {
            e["barcode"]: e
            for e in BRB.ET.DNA(
                config, str(outputDir), baseDict, {"S1": "L1", "S2": "L2"}
            )
        }

        assert m["L1"]["insert_size"] == 250
        assert m["L2"]["mapped_reads"] == 95.50
        assert m["L2"]["dupped_reads"] == 20.0
        assert "insert_size" not in m["L2"]

    def test_run_without_bowtie2_output_sends_base_fields_only(self, tmp_path):
        config = parkourConfig(tmp_path)
        baseDict = {"L1": ["S1", 90, 0.05, 10.0]}
        m = BRB.ET.DNA(config, str(tmp_path / "DNA_human"), baseDict, {"S1": "L1"})

        assert m == [
            {
                "barcode": "L1",
                "reads_pf_sequenced": 90,
                "confident_reads": 0.05,
                "optical_duplicates": 10.0,
            }
        ]


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

    def test_sample_missing_from_run_is_skipped(self, tmp_path):
        import os

        config = parkourConfig(tmp_path)
        outputDir = tmp_path / "RNA_human"
        multiqcDir = outputDir / "multiQC" / "multiqc_data"
        os.makedirs(multiqcDir, exist_ok=True)
        for s in ("S1", "S9"):
            os.makedirs(outputDir / "STAR" / s, exist_ok=True)
            (outputDir / "STAR" / s / f"{s}.Log.final.out").write_text(
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

        assert len(m) == 1
        assert m[0]["barcode"] == "L1"
        assert m[0]["uniq_mapped"] == 85.0


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
            return fakeResponse(body="OK")

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
        assert result[:5] == ["ProjectX", "human", "WGS", "DNA", "success"]
        assert result[5].text == "OK"

    def test_rna_pipeline_dispatches_to_rna_module(self, tmp_path, monkeypatch):
        config = parkourConfig(tmp_path)
        monkeypatch.setattr(
            BRB.ET,
            "getBaseStatistics",
            lambda *a, **k: ({"L1": ["S1", 1, 2, 3]}, {"S1": "L1"}),
        )
        monkeypatch.setattr(BRB.ET, "RNA", lambda *a, **k: [{"barcode": "L1"}])
        sent = {}
        monkeypatch.setattr(
            BRB.ET, "sendToParkour", lambda c, m: sent.update(msg=m) or fakeResponse()
        )

        result = BRB.ET.phoneHome(
            config, "outdir", "RNA", [("S1", "libtype")], "human", "h38", "P", "mRNA"
        )

        assert sent["msg"] == [{"barcode": "L1"}]
        assert result[:5] == ["P", "human", "mRNA", "RNA", "success"]

    def test_other_pipeline_uses_base_fields_only(self, tmp_path, monkeypatch):
        config = parkourConfig(tmp_path)
        monkeypatch.setattr(
            BRB.ET,
            "getBaseStatistics",
            lambda *a, **k: ({"L1": ["S1", 1, 2, 3]}, {"S1": "L1"}),
        )
        sent = {}
        monkeypatch.setattr(
            BRB.ET, "sendToParkour", lambda c, m: sent.update(msg=m) or fakeResponse()
        )

        result = BRB.ET.phoneHome(
            config, "outdir", "HiC", [("S1", "libtype")], "human", "h38", "P", "HiC"
        )

        assert sent["msg"] == [
            {
                "barcode": "L1",
                "reads_pf_sequenced": 1,
                "confident_reads": 2,
                "optical_duplicates": 3,
            }
        ]
        assert result[4] == "success"


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
        sent = {}
        monkeypatch.setattr(
            BRB.ET, "sendToParkour", lambda c, m: sent.update(msg=m) or fakeResponse()
        )

        skipList = [("S1", "S1", "WGS")]
        result = BRB.ET.telegraphHome(
            config, "PIName", "ProjectX", skipList, organism="human", org_label="h38"
        )

        assert result[0] == "ProjectX"
        assert result[1] == "human"
        assert result[2] == "WGS"
        assert result[3] is None
        assert result[4] is None
        assert sent["msg"] == [
            {
                "barcode": "L1",
                "reads_pf_sequenced": 1,
                "confident_reads": 2,
                "optical_duplicates": 3,
            }
        ]
        assert result[6] is False


class TestPeakQC:
    def test_maps_macs2_qc_columns(self, tmp_path):
        qc = tmp_path / "MACS2_QC"
        qc.mkdir()
        (qc / "s1.filtered.BAM_peaks.qc.txt").write_text(
            "peak_count\tFRiP\tpeak_genome_coverage\n3364\t0.198\t0.0123\n"
        )
        assert BRB.ET.peakQC(str(tmp_path), "s1") == {
            "peak_count": 3364,
            "frip": 0.198,
            "peak_genome_coverage": 0.0123,
        }

    def test_peak_file_with_only_a_header_adds_no_keys(self, tmp_path):
        qc = tmp_path / "MACS2_QC"
        qc.mkdir()
        (qc / "s1.filtered.BAM_peaks.qc.txt").write_text(
            "peak_count\tFRiP\tpeak_genome_coverage\n"
        )
        assert BRB.ET.peakQC(str(tmp_path), "s1") == {}

    def test_peak_file_with_a_missing_column_keeps_the_other_keys(self, tmp_path):
        qc = tmp_path / "MACS2_QC"
        qc.mkdir()
        (qc / "s1.filtered.BAM_peaks.qc.txt").write_text(
            "peak_count\tFRiP\n3364\t0.198\n"
        )
        assert BRB.ET.peakQC(str(tmp_path), "s1") == {
            "peak_count": 3364,
            "frip": 0.198,
        }

    def test_run_without_peak_calling_adds_no_keys(self, tmp_path):
        assert BRB.ET.peakQC(str(tmp_path), "s1") == {}

    def test_chip_run_reads_macs2_dir(self, tmp_path):
        macs = tmp_path / "MACS2"
        macs.mkdir()
        (macs / "s1.filtered.BAM_peaks.qc.txt").write_text(
            "peak_count\tFRiP\tpeak_genome_coverage\n110371\t0.125\t0.0004\n"
        )
        assert BRB.ET.peakQC(str(tmp_path), "s1") == {
            "peak_count": 110371,
            "frip": 0.125,
            "peak_genome_coverage": 0.0004,
        }

    def test_macs2_qc_wins_when_both_dirs_exist(self, tmp_path):
        for d, count in (("MACS2_QC", 1), ("MACS2", 2)):
            (tmp_path / d).mkdir()
            (tmp_path / d / "s1.filtered.BAM_peaks.qc.txt").write_text(
                f"peak_count\tFRiP\tpeak_genome_coverage\n{count}\t0.1\t0.1\n"
            )
        assert BRB.ET.peakQC(str(tmp_path), "s1")["peak_count"] == 1

    def test_chip_chromatin_sends_peak_and_fragment_keys(self, tmp_path):
        config = parkourConfig(tmp_path)
        outputDir = tmp_path / "ChIP_mm10"
        for d in (
            "Bowtie2",
            "multiQC/multiqc_data",
            "deepTools_qc/bamPEFragmentSize",
            "MACS2",
        ):
            os.makedirs(outputDir / d)
        (outputDir / "Bowtie2" / "S1.Bowtie2_summary.txt").write_text(
            "some header\n95.50% overall alignment rate\n"
        )
        (outputDir / "multiQC/multiqc_data/multiqc_samtools_flagstat.txt").write_text(
            "Sample\ttotal_passed\tduplicates_passed\tread1_passed\nS1\t100\t10\t60\n"
        )
        (
            outputDir / "deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv"
        ).write_text(
            "\tFrag. Len. Min.\tFrag. Len. 1st. Qu.\tFrag. Len. Mean\t"
            "Frag. Len. Median\tFrag. Len. 3rd Qu.\tFrag. Len. 10%\t"
            "Frag. Len. 90%\n"
            "filtered_bam/S1.filtered.bam\t1\t100\t200\t210\t300\t50\t400\n"
        )
        (outputDir / "MACS2" / "S1.filtered.BAM_peaks.qc.txt").write_text(
            "peak_count\tFRiP\tpeak_genome_coverage\n12\t0.25\t0.01\n"
        )
        baseDict = {"L1": ["S1", 90, 0.05, 10.0]}

        m = BRB.ET.chromatin(config, str(outputDir), baseDict, {"S1": "L1"})

        assert m == [
            {
                "barcode": "L1",
                "reads_pf_sequenced": 90,
                "confident_reads": 0.05,
                "optical_duplicates": 10.0,
                "dupped_reads": 10.0,
                "mapped_reads": 95.5,
                "insert_size": 210,
                "fragment_q1": 100,
                "fragment_median": 210,
                "fragment_q3": 300,
                "fragment_p10": 50,
                "fragment_p90": 400,
                "fragments_total": 60,
                "peak_count": 12,
                "frip": 0.25,
                "peak_genome_coverage": 0.01,
            }
        ]


class TestFragmentQC:
    def test_reads_quantiles_and_total_pairs(self, tmp_path):
        frag = tmp_path / "deepTools_qc" / "bamPEFragmentSize"
        frag.mkdir(parents=True)
        (frag / "fragmentSize.metric.tsv").write_text(
            "\tFrag. Len. 1st. Qu.\tFrag. Len. Median\tFrag. Len. 3rd Qu."
            "\tFrag. Len. 10%\tFrag. Len. 90%\n"
            "filtered_bam/s1.filtered.bam\t57.0\t85.0\t151.0\t44.0\t246.0\n"
        )
        mq = tmp_path / "multiQC" / "multiqc_data"
        mq.mkdir(parents=True)
        (mq / "multiqc_samtools_flagstat.txt").write_text(
            "Sample\ttotal_passed\tread1_passed\ns1\t16595486\t8297743\n"
        )
        assert BRB.ET.fragmentQC(str(tmp_path), "s1") == {
            "fragment_q1": 57,
            "fragment_median": 85,
            "fragment_q3": 151,
            "fragment_p10": 44,
            "fragment_p90": 246,
            "fragments_total": 8297743,
        }

    def test_missing_sources_add_no_keys(self, tmp_path):
        assert BRB.ET.fragmentQC(str(tmp_path), "s1") == {}


HIC_BASE = {"25L1": ["hicA", 1000, 1000, 0.5], "25L2": ["hicB", 2000, 2000, 0.5]}
HIC_SAMPLES = {"hicA": "25L1", "hicB": "25L2"}


class TestHiC:
    def writeMultiqc(self, tmp_path, rows):
        mq = tmp_path / "multiQC" / "multiqc_data"
        mq.mkdir(parents=True)
        (mq / "multiqc_hicexplorer.txt").write_text(
            "Sample\tPairs mappable, unique and high quality_pct\n" + rows
        )

    def test_sends_10kb_pair_rate_per_sample(self, tmp_path):
        self.writeMultiqc(
            tmp_path,
            "hicA_bs10kb.cool\t87.7\nhicA_bs50kb.cool\t12.3\nhicB_bs10kb.cool\t61.0\n",
        )
        m = BRB.ET.HiC(None, str(tmp_path), HIC_BASE, HIC_SAMPLES)
        assert {e["barcode"]: e["uniq_mapped"] for e in m} == {
            "25L1": 87.7,
            "25L2": 61.0,
        }
        assert m[0]["reads_pf_sequenced"] == 1000

    def test_library_without_a_row_keeps_base_stats_only(self, tmp_path):
        self.writeMultiqc(tmp_path, "hicA_bs10kb.cool\t87.7\n")
        m = BRB.ET.HiC(None, str(tmp_path), HIC_BASE, HIC_SAMPLES)
        assert "uniq_mapped" not in m[1]
        assert m[1]["optical_duplicates"] == 0.5

    def test_run_without_multiqc_sends_base_stats(self, tmp_path):
        m = BRB.ET.HiC(None, str(tmp_path), HIC_BASE, HIC_SAMPLES)
        assert [e["barcode"] for e in m] == ["25L1", "25L2"]
        assert all("uniq_mapped" not in e for e in m)

    def test_10kb_row_beats_an_earlier_other_bin_size(self, tmp_path):
        self.writeMultiqc(tmp_path, "hicA_bs100kb.cool\t12.0\nhicA_bs10kb.cool\t87.7\n")
        m = BRB.ET.HiC(None, str(tmp_path), HIC_BASE, HIC_SAMPLES)
        assert m[0]["uniq_mapped"] == 87.7

    def test_first_row_used_when_no_10kb_row(self, tmp_path):
        self.writeMultiqc(tmp_path, "hicA_bs50kb.cool\t12.3\nhicA_bs100kb.cool\t9.0\n")
        m = BRB.ET.HiC(None, str(tmp_path), HIC_BASE, HIC_SAMPLES)
        assert m[0]["uniq_mapped"] == 12.3

    def test_shared_barcode_uses_rate_of_the_sample_in_base_stats(self, tmp_path):
        self.writeMultiqc(tmp_path, "hicA_bs10kb.cool\t87.7\nhicB_bs10kb.cool\t61.0\n")
        shared = {"hicA": "25L1", "hicB": "25L1"}
        for sample, rate in (("hicA", 87.7), ("hicB", 61.0)):
            base = {"25L1": [sample, 1000, 1000, 0.5]}
            m = BRB.ET.HiC(None, str(tmp_path), base, shared)
            assert m[0]["uniq_mapped"] == rate


def _fastqcLayout(tmp_path, stats, reads=5):
    import gzip

    qc = tmp_path / "FASTQC_Project_1" / "Sample_B1"
    qc.mkdir(parents=True)
    for name, text in stats.items():
        (qc / name).write_text(text)
    raw = tmp_path / "Project_1" / "Sample_B1"
    raw.mkdir(parents=True)
    with gzip.open(raw / "s1_R1.fastq.gz", "wt") as f:
        f.write("@r\nA\n+\nI\n" * reads)
    return qc


@pytest.mark.parametrize("text", ["", "10.0 100.0", "opt total", "5 3", "1 2 3"])
def test_getNReads_unusable_stats_file_gives_read_count_and_unknown_rate(
    tmp_path, text
):
    qc = _fastqcLayout(tmp_path, {"s1.duplicate.txt": text})
    total, opt = BRB.ET.getNReads(str(qc))
    assert total == 5
    assert math.isnan(opt)


def test_getNReads_skips_empty_file_and_uses_a_valid_sibling(tmp_path):
    qc = _fastqcLayout(tmp_path, {"a.metrics": "", "b.metrics": "10 100\n"})
    assert BRB.ET.getNReads(str(qc)) == (90, 10.0)


def test_getNReads_unusable_stats_and_no_fastq_gives_unknown_not_crash(tmp_path):
    qc = tmp_path / "FASTQC_Project_1" / "Sample_B1"
    qc.mkdir(parents=True)
    (qc / "s1.metrics").write_text("")
    total, opt = BRB.ET.getNReads(str(qc))
    assert math.isnan(total) and math.isnan(opt)


def test_dropNaN_leaves_out_undefined_rates_only():
    m = [
        {
            "barcode": "B",
            "dupped_reads": float("nan"),
            "optical_duplicates": 0.0,
            "mapped_reads": 3.5,
        }
    ]
    assert BRB.ET.dropNaN(m) == [
        {"barcode": "B", "optical_duplicates": 0.0, "mapped_reads": 3.5}
    ]


def test_phoneHome_sends_base_stats_when_the_parser_raises(tmp_path, monkeypatch):
    config = parkourConfig(tmp_path)
    monkeypatch.setattr(
        BRB.ET,
        "getBaseStatistics",
        lambda *a, **k: ({"L1": ["S1", 10, 8, 0.5]}, {"S1": "L1"}),
    )

    def boom(*a, **k):
        raise IndexError("empty peaks.qc")

    monkeypatch.setattr(BRB.ET, "chromatin", boom)
    sent = {}
    monkeypatch.setattr(
        BRB.ET, "sendToParkour", lambda c, m: sent.update(msg=m) or fakeResponse()
    )
    res = BRB.ET.phoneHome(
        config,
        "outdir",
        "ATAC",
        [("S1", "libtype")],
        "human",
        "h38",
        "ProjectX",
        "ATAC",
    )
    assert sent["msg"] == [
        {
            "barcode": "L1",
            "reads_pf_sequenced": 10,
            "confident_reads": 8,
            "optical_duplicates": 0.5,
        }
    ]
    assert res[4] == "success"


def test_telegraphHome_looks_up_the_analysis_prefixed_project_dir(
    tmp_path, monkeypatch
):
    config = parkourConfig(tmp_path)
    config["Paths"]["groupData"] = str(tmp_path / "groupdata")
    os.makedirs(tmp_path / "groupdata" / "PIName" / "sequencing_data", exist_ok=True)
    got = {}

    def fakeBase(cfg, outputDir, samples, label):
        got["dir"] = Path(outputDir)
        return {"L1": ["S1", 1, 2, 3]}, {"S1": "L1"}

    monkeypatch.setattr(BRB.ET, "getBaseStatistics", fakeBase)
    monkeypatch.setattr(BRB.ET, "sendToParkour", lambda c, m: fakeResponse())
    BRB.ET.telegraphHome(
        config, "PIName", "ProjectX", [("S1", "S1", "WGS")], "human", "h38"
    )
    assert got["dir"].parent.name == "Analysis_ProjectX"


def _answerPost(monkeypatch, answer):
    def fakePost(*a, **k):
        if isinstance(answer, Exception):
            raise answer
        return answer

    monkeypatch.setattr(BRB.ET.requests, "post", fakePost)


class TestParkourOutcomeStatus:
    @pytest.mark.parametrize(
        "answer, status",
        [
            (fakeResponse(200), "success"),
            (fakeResponse(500, "Internal Server Error"), "FAILED"),
            (requests.ConnectionError("parkour down"), "FAILED"),
        ],
    )
    def test_phoneHome_row_reports_parkour_outcome(
        self, tmp_path, monkeypatch, answer, status
    ):
        config = parkourConfig(tmp_path)
        monkeypatch.setattr(
            BRB.ET,
            "getBaseStatistics",
            lambda *a, **k: ({"L1": ["S1", 1, 2, 3]}, {"S1": "L1"}),
        )
        _answerPost(monkeypatch, answer)

        result = BRB.ET.phoneHome(
            config, "outdir", "WGS", [("S1", "libtype")], "human", "h38", "P", "WGS"
        )

        assert result[4] == status

    @pytest.mark.parametrize(
        "answer, status",
        [
            (fakeResponse(200), None),
            (fakeResponse(500, "Internal Server Error"), "FAILED"),
            (requests.ConnectionError("parkour down"), "FAILED"),
        ],
    )
    def test_telegraphHome_row_reports_parkour_outcome(
        self, tmp_path, monkeypatch, answer, status
    ):
        config = parkourConfig(tmp_path)
        config["Paths"]["groupData"] = str(tmp_path / "groupdata")
        os.makedirs(
            tmp_path / "groupdata" / "PIName" / "sequencing_data", exist_ok=True
        )
        monkeypatch.setattr(
            BRB.ET,
            "getBaseStatistics",
            lambda *a, **k: ({"L1": ["S1", 1, 2, 3]}, {"S1": "L1"}),
        )
        _answerPost(monkeypatch, answer)

        result = BRB.ET.telegraphHome(
            config, "PIName", "ProjectX", [("S1", "S1", "WGS")], "human", "h38"
        )

        assert result[4] == status


class TestStarMapping:
    def test_missing_lines_are_nan_not_zero(self, tmp_path):
        f = tmp_path / "S1.Log.final.out"
        f.write_text("Uniquely mapped reads %\t80.00%\n")
        mapped, uniq, multi = BRB.ET._starMapping(str(f))
        assert uniq == 80.0
        assert math.isnan(multi)
        assert math.isnan(mapped)

    def test_empty_log_gives_all_nan(self, tmp_path):
        f = tmp_path / "S1.Log.final.out"
        f.write_text("")
        assert all(math.isnan(v) for v in BRB.ET._starMapping(str(f)))


class TestRNAMissingSources:
    def _writeStar(self, outputDir, sample, uniq):
        d = outputDir / "STAR" / sample
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{sample}.Log.final.out").write_text(
            f"Uniquely mapped reads %\t{uniq}%\n"
            "% of reads mapped to multiple loci\t5.00%\n"
        )

    def test_library_missing_a_source_keeps_base_stats_others_keep_rna(
        self, tmp_path, caplog
    ):
        outputDir = tmp_path / "RNA_human"
        self._writeStar(outputDir, "S.1", "85.00")
        self._writeStar(outputDir, "S2", "70.00")
        mq = outputDir / "multiQC" / "multiqc_data"
        mq.mkdir(parents=True)
        (mq / "multiqc_samtools_flagstat.txt").write_text(
            "Sample\ttotal_passed\tduplicates_passed\nS.1\t100\t10\nS2\t100\t20\n"
        )
        (mq / "multiqc_featurecounts.txt").write_text(
            "Sample\tTotal\tAssigned\nS.1.filtered\t100\t80\n"
        )
        baseDict = {
            "L1": ["S.1", 90, 0.05, 10.0],
            "L2": ["S2", 80, 0.04, 20.0],
            "L3": ["S3", 70, 0.03, 30.0],
        }

        m = BRB.ET.RNA(
            None, str(outputDir), baseDict, {"S.1": "L1", "S2": "L2", "S3": "L3"}
        )

        byLib = {e["barcode"]: e for e in m}
        assert byLib["L1"]["mapped_reads"] == 90.0
        assert byLib["L1"]["uniq_mapped"] == 85.0
        assert byLib["L1"]["dupped_reads"] == 10.0
        assert byLib["L1"]["assigned_reads"] == 80.0
        assert "mapped_reads" not in byLib["L2"]
        assert byLib["L2"]["reads_pf_sequenced"] == 80
        assert "mapped_reads" not in byLib["L3"]
        assert byLib["L3"]["optical_duplicates"] == 30.0
        assert "L2" in caplog.text and "multiQC featureCounts" in caplog.text
        assert "L3" in caplog.text and "STAR Log.final.out" in caplog.text

    def test_missing_multiqc_tables_send_base_stats_for_every_library(self, tmp_path):
        outputDir = tmp_path / "RNA_human"
        self._writeStar(outputDir, "S1", "85.00")

        m = BRB.ET.RNA(
            None, str(outputDir), {"L1": ["S1", 90, 0.05, 10.0]}, {"S1": "L1"}
        )

        assert m == [
            {
                "barcode": "L1",
                "reads_pf_sequenced": 90,
                "confident_reads": 0.05,
                "optical_duplicates": 10.0,
            }
        ]


PAIRS_STATS = (
    "total\t1000\ntotal_unmapped\t200\ntotal_single_sided_mapped\t100\n"
    "total_mapped\t700\ntotal_dups\t140\ntotal_nodups\t560\ncis\t500\n"
    "pair_types/MM\t10\n"
)


PAIRS_BASE = {"25L1": ["mcA", 1000, 1000, 0.5], "25L2": ["mcB", 2000, 2000, 0.5]}
PAIRS_SAMPLES = {"mcA": "25L1", "mcB": "25L2"}


class TestMakePairs:
    def writeStats(self, tmp_path, sample, text):
        pairs = tmp_path / "pairs"
        pairs.mkdir(exist_ok=True)
        (pairs / f"{sample}.pairs.dedup.stats").write_text(text)

    def test_sends_mapped_and_dup_rate_per_sample(self, tmp_path):
        self.writeStats(tmp_path, "mcA", PAIRS_STATS)
        m = BRB.ET.makePairs(None, str(tmp_path), PAIRS_BASE, PAIRS_SAMPLES)
        assert m[0]["mapped_reads"] == 70.0
        assert m[0]["dupped_reads"] == 20.0
        assert m[0]["reads_pf_sequenced"] == 1000

    def test_library_without_stats_keeps_base_stats_only(self, tmp_path):
        self.writeStats(tmp_path, "mcA", PAIRS_STATS)
        m = BRB.ET.makePairs(None, str(tmp_path), PAIRS_BASE, PAIRS_SAMPLES)
        assert "mapped_reads" not in m[1] and "dupped_reads" not in m[1]
        assert m[1]["optical_duplicates"] == 0.5

    def test_unusable_stats_file_is_not_sent_as_zero(self, tmp_path):
        self.writeStats(tmp_path, "mcA", "total\t0\n")
        m = BRB.ET.makePairs(None, str(tmp_path), PAIRS_BASE, PAIRS_SAMPLES)
        assert "mapped_reads" not in m[0] and "dupped_reads" not in m[0]

    def test_no_mapped_pairs_leaves_dup_rate_out_of_payload(self, tmp_path):
        self.writeStats(tmp_path, "mcA", "total\t10\ntotal_mapped\t0\ntotal_dups\t0\n")
        m = BRB.ET.dropNaN(
            BRB.ET.makePairs(None, str(tmp_path), PAIRS_BASE, PAIRS_SAMPLES)
        )
        assert m[0]["mapped_reads"] == 0.0
        assert "dupped_reads" not in m[0]
