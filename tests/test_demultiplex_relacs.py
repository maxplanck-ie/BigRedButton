import argparse
import io

import pytest

import BRB.demultiplex_relacs as dr


def _args(umiLength=0, buffer=1):
    return argparse.Namespace(umiLength=umiLength, buffer=buffer)


class TestParseArgs:
    def test_parses_required_positional_arguments(self):
        args = dr.parseArgs(["sampleTable.txt", "outdir"])
        assert args.sampleTable == "sampleTable.txt"
        assert args.output == "outdir"
        assert args.numThreads == 1
        assert args.buffer == 1
        assert args.umiLength == 0

    def test_parses_optional_flags(self):
        args = dr.parseArgs(
            ["-p", "4", "-b", "2", "--umiLength", "8", "table.txt", "out"]
        )
        assert args.numThreads == 4
        assert args.buffer == 2
        assert args.umiLength == 8

    def test_missing_positional_arguments_exits(self):
        with pytest.raises(SystemExit):
            dr.parseArgs([])


class TestCheckDuplicatedLabels:
    def test_unique_labels_pass_without_exit(self, capsys):
        data = {
            "SampleA": {
                "ACTACT": ["Sample1", ""],
                "TGACTG": ["Sample2", ""],
            }
        }
        dr.checkDuplicatedLabels(data)  # must not raise / exit
        assert "OK" in capsys.readouterr().out

    def test_duplicated_labels_across_samples_exits(self, capsys):
        data = {
            "SampleA": {"ACTACT": ["Sample1", ""]},
            "SampleB": {"TGACTG": ["Sample1", ""]},
        }
        with pytest.raises(SystemExit) as excinfo:
            dr.checkDuplicatedLabels(data)
        assert excinfo.value.code == 1
        assert "duplicated labels" in capsys.readouterr().out


class TestReadSampleTable:
    def test_three_column_table_parses_barcode_and_label(self, tmp_path):
        table = tmp_path / "samples.txt"
        table.write_text(
            "Sample_Lib_C838_11\tACTACT\tSample1\n"
            "Sample_Lib_C838_11\tTGACTG\tSample2\n"
            "Sample_Lib_C838_11\tdefault\tSample3\n"
        )
        d, bcLen = dr.readSampleTable(str(table))

        assert bcLen == 6
        assert d["Sample_Lib_C838_11"]["ACTACT"] == ["Sample1", ""]
        assert d["Sample_Lib_C838_11"]["default"] == ["Sample3", ""]

    def test_four_column_table_includes_bc_position(self, tmp_path):
        table = tmp_path / "samples.txt"
        table.write_text("SampleA\tA1-21\tACTACT\tSample1\n")
        d, bcLen = dr.readSampleTable(str(table))

        assert d["SampleA"]["ACTACT"] == ["Sample1", "A1-21"]
        assert bcLen == 6

    def test_label_with_spaces_is_sanitized(self, tmp_path):
        table = tmp_path / "samples.txt"
        table.write_text("SampleA\tACTACT\tSample One\n")
        d, _bcLen = dr.readSampleTable(str(table))
        assert d["SampleA"]["ACTACT"][0] == "Sample_One"

    def test_lines_with_too_few_columns_are_skipped(self, tmp_path):
        table = tmp_path / "samples.txt"
        table.write_text("garbage line\nSampleA\tACTACT\tSample1\n")
        d, bcLen = dr.readSampleTable(str(table))
        assert list(d.keys()) == ["SampleA"]
        assert bcLen == 6

    def test_empty_samplesheet_exits(self, tmp_path, capsys):
        table = tmp_path / "empty.txt"
        table.write_text("")
        with pytest.raises(SystemExit) as excinfo:
            dr.readSampleTable(str(table))
        assert excinfo.value.code == 1
        assert "Samplesheet is empty" in capsys.readouterr().out

    def test_duplicated_labels_in_samplesheet_exits(self, tmp_path):
        table = tmp_path / "samples.txt"
        table.write_text("SampleA\tACTACT\tSample1\nSampleB\tTGACTG\tSample1\n")
        with pytest.raises(SystemExit):
            dr.readSampleTable(str(table))


class TestMatchSample:
    def _oDict(self):
        return {"ACTACT": ["Sample1", ""], "TGACTG": ["Sample2", ""]}

    def test_exact_match_single_end(self):
        seq = "ACTACTGGGG\n"
        bc, isDefault = dr.matchSample(seq, None, self._oDict(), 6, 0)
        assert bc == "ACTACT"
        assert isDefault is True

    def test_exact_match_paired_end_agrees(self):
        seq1 = "ACTACTGGGG\n"
        seq2 = "ACTACTCCCC\n"
        bc, isDefault = dr.matchSample(seq1, seq2, self._oDict(), 6, 0)
        assert bc == "ACTACT"
        assert isDefault is True

    def test_one_mismatch_is_tolerated(self):
        # AGTACT differs from ACTACT by one base
        seq = "AGTACTGGGG\n"
        bc, isDefault = dr.matchSample(seq, None, self._oDict(), 6, 0)
        assert bc == "ACTACT"
        assert isDefault is True

    def test_unrecognized_barcode_falls_back_to_default(self):
        seq = "GGGGGGGGGG\n"
        bc, isDefault = dr.matchSample(seq, None, self._oDict(), 6, 0)
        assert bc == "default"
        assert isDefault is False

    def test_umi_length_offsets_the_barcode_window(self):
        # 4-base UMI precedes the barcode
        seq = "NNNNACTACTGGGG\n"
        bc, isDefault = dr.matchSample(seq, None, self._oDict(), 6, 4)
        assert bc == "ACTACT"
        assert isDefault is True


class TestWriteRead:
    def test_trims_barcode_and_tags_read_name(self):
        of = io.BytesIO()
        lineList = ["@read1 1:N:0\n", "ACTACTGGGGAAAA\n", "+\n", "IIIIIIIIIIIIII\n"]
        args = _args(umiLength=0, buffer=1)

        rname, bc = dr.writeRead(lineList, of, "ACTACT", 6, args, doTrim=True)

        assert rname == "@read1_ACTACT 1:N:0\n"
        assert bc == "ACTACT"
        written = of.getvalue().decode()
        assert written.startswith("@read1_ACTACT 1:N:0\n")
        assert "GGGAAAA\n" in written  # barcode (6) + 1 buffer base trimmed

    def test_umi_is_appended_to_read_name(self):
        of = io.BytesIO()
        lineList = ["@read1 1:N:0\n", "NNNNACTACTGGGG\n", "+\n", "IIIIIIIIIIIIII\n"]
        args = _args(umiLength=4, buffer=1)

        rname, _bc = dr.writeRead(lineList, of, "ACTACT", 6, args, doTrim=True)

        assert rname == "@read1_ACTACT_NNNN 1:N:0\n"

    def test_no_trim_writes_read_unmodified(self):
        of = io.BytesIO()
        lineList = ["@read1 1:N:0\n", "ACTACTGGGG\n", "+\n", "IIIIIIIIII\n"]
        args = _args()

        rname, _bc = dr.writeRead(lineList, of, "default", 6, args, doTrim=False)

        assert rname == "@read1 1:N:0\n"
        written = of.getvalue().decode()
        assert "ACTACTGGGG\n" in written


class TestWriteRead2:
    def test_renumbers_mate_to_read_two_and_trims(self):
        of = io.BytesIO()
        lineList = ["@read1 1:N:0\n", "ACTACTGGGG\n", "+\n", "IIIIIIIIII\n"]
        args = _args(umiLength=0, buffer=1)

        dr.writeRead2(lineList, of, 6, args, doTrim=True)

        written = of.getvalue().decode()
        assert written.startswith("@read1 2:N:0\n")
        assert "GGG\n" in written


class TestWritePaired:
    def test_writes_both_mates_with_shared_read_name(self):
        of1, of2 = io.BytesIO(), io.BytesIO()
        read1 = ["@read1 1:N:0\n", "ACTACTGGGG\n", "+\n", "IIIIIIIIII\n"]
        read2 = ["@read1 2:N:0\n", "ACTACTCCCC\n", "+\n", "IIIIIIIIII\n"]
        args = _args(umiLength=0, buffer=1)

        bc = dr.writePaired(read1, read2, (of1, of2), "ACTACT", 6, args, doTrim=True)

        assert bc == "ACTACT"
        w1 = of1.getvalue().decode()
        w2 = of2.getvalue().decode()
        assert w1.startswith("@read1_ACTACT 1:N:0\n")
        assert w2.startswith("@read1_ACTACT 1:N:0\n")
        assert "GGG\n" in w1
        assert "CCC\n" in w2

    def test_umi_from_both_mates_is_concatenated_into_read_name(self):
        of1, of2 = io.BytesIO(), io.BytesIO()
        read1 = ["@read1 1:N:0\n", "AAAAACTACTGGGG\n", "+\n", "IIIIIIIIIIIIII\n"]
        read2 = ["@read1 2:N:0\n", "TTTTACTACTCCCC\n", "+\n", "IIIIIIIIIIIIII\n"]
        args = _args(umiLength=4, buffer=1)

        bc = dr.writePaired(read1, read2, (of1, of2), "ACTACT", 6, args, doTrim=True)

        assert bc == "ACTACT"
        w1 = of1.getvalue().decode()
        assert w1.startswith("@read1_ACTACT_AAAATTTT 1:N:0\n")
