import argparse
import gzip
import io
from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

import BRB.demultiplex_relacs as dr


def _args(umiLength=0, buffer=1):
    return argparse.Namespace(umiLength=umiLength, buffer=buffer)


def _gzFastq(path, records):
    with gzip.open(path, "wb") as fh:
        for name, seq, qual in records:
            fh.write(f"{name}\n{seq}\n+\n{qual}\n".encode())


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

    def test_one_mismatch_paired_end_agrees(self):
        # AGTACT differs from ACTACT (the key) by one base, and mate 2 agrees
        # with mate 1's own (mismatched) barcode -- exercises the paired-end
        # branch of the 1-mismatch fallback loop, not just the single-end one.
        seq1 = "AGTACTGGGG\n"
        seq2 = "AGTACTCCCC\n"
        bc, isDefault = dr.matchSample(seq1, seq2, self._oDict(), 6, 0)
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
    def test_writes_both_mates_with_shared_name_and_own_pair_flag(self):
        of1, of2 = io.BytesIO(), io.BytesIO()
        read1 = ["@read1 1:N:0\n", "ACTACTGGGG\n", "+\n", "IIIIIIIIII\n"]
        read2 = ["@read1 2:N:0\n", "ACTACTCCCC\n", "+\n", "IIIIIIIIII\n"]
        args = _args(umiLength=0, buffer=1)

        bc = dr.writePaired(read1, read2, (of1, of2), "ACTACT", 6, args, doTrim=True)

        assert bc == "ACTACT"
        w1 = of1.getvalue().decode()
        w2 = of2.getvalue().decode()
        # Same barcode-suffixed name, differing only in the pair-number field.
        assert w1.startswith("@read1_ACTACT 1:N:0\n")
        assert w2.startswith("@read1_ACTACT 2:N:0\n")
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


class TestProcessSingle:
    """
    Regression test for a bug where processSingle called writeRead with the
    raw oDict list value (as wired up by wrapper: oDict[k] = [stdin]) instead
    of unwrapping it, crashing with AttributeError on every read -- the
    single-end demux path was unconditionally broken.
    """

    def test_routes_matched_and_unmatched_reads_to_their_files(self, tmp_path):
        r1 = tmp_path / "s_R1.fastq.gz"
        _gzFastq(
            r1,
            [
                ("@read1 1:N:0:1", "ACGTNTTTTTTTT", "IIIIIIIIIIIII"),
                ("@read2 1:N:0:1", "GGGGNCCCCCCCC", "IIIIIIIIIIIII"),
            ],
        )
        matched = io.BytesIO()
        default = io.BytesIO()
        sDict = {"ACGT": [matched], "default": [default]}

        dr.processSingle(_args(), sDict, 4, str(r1))

        matchedText = matched.getvalue().decode()
        assert matchedText == "@read1_ACGT 1:N:0:1\nTTTTTTTT\n+\nIIIIIIII\n"
        defaultText = default.getvalue().decode()
        assert defaultText == "@read2 1:N:0:1\nGGGGNCCCCCCCC\n+\nIIIIIIIIIIIII\n"


class TestProcessPaired:
    def test_routes_matched_and_unmatched_pairs_and_reports_bc_occurance(
        self, tmp_path, monkeypatch
    ):
        r1 = tmp_path / "s_R1.fastq.gz"
        r2 = tmp_path / "s_R2.fastq.gz"
        _gzFastq(
            r1,
            [
                ("@read1 1:N:0:1", "ACGTNTTTTTTTT", "IIIIIIIIIIIII"),
                ("@read2 1:N:0:1", "GGGGNCCCCCCCC", "IIIIIIIIIIIII"),
            ],
        )
        _gzFastq(
            r2,
            [
                ("@read1 2:N:0:1", "ACGTNAAAAAAAA", "IIIIIIIIIIIII"),
                ("@read2 2:N:0:1", "GGGGNTTTTTTTT", "IIIIIIIIIIIII"),
            ],
        )
        matched1, matched2 = io.BytesIO(), io.BytesIO()
        default1, default2 = io.BytesIO(), io.BytesIO()
        oDict = {"ACGT": [matched1, matched2], "default": [default1, default2]}
        bc_dict = {}
        seenPlot = {}
        monkeypatch.setattr(
            dr,
            "plot_bc_occurance",
            lambda read1, bcd, falseBc, output, oriDict: seenPlot.update(
                read1=read1, bc_dict=dict(bcd), false_bc=falseBc
            ),
        )

        args = _args()
        args.output = str(tmp_path)
        dr.processPaired(args, oDict, 4, str(r1), str(r2), bc_dict, oDict)

        m1 = matched1.getvalue().decode()
        assert m1 == "@read1_ACGT 1:N:0:1\nTTTTTTTT\n+\nIIIIIIII\n"
        assert (
            matched2.getvalue().decode()
            == "@read1_ACGT 2:N:0:1\nAAAAAAAA\n+\nIIIIIIII\n"
        )
        assert default1.getvalue().decode() == (
            "@read2 1:N:0:1\nGGGGNCCCCCCCC\n+\nIIIIIIIIIIIII\n"
        )
        assert bc_dict == {"ACGT": 1}
        assert seenPlot == {"read1": str(r1), "bc_dict": {"ACGT": 1}, "false_bc": 1}

    def test_repeated_barcode_increments_its_bc_dict_count(self, tmp_path, monkeypatch):
        r1 = tmp_path / "s_R1.fastq.gz"
        r2 = tmp_path / "s_R2.fastq.gz"
        _gzFastq(
            r1,
            [
                ("@read1 1:N:0:1", "ACGTNTTTTTTTT", "IIIIIIIIIIIII"),
                ("@read2 1:N:0:1", "ACGTNCCCCCCCC", "IIIIIIIIIIIII"),
            ],
        )
        _gzFastq(
            r2,
            [
                ("@read1 2:N:0:1", "ACGTNAAAAAAAA", "IIIIIIIIIIIII"),
                ("@read2 2:N:0:1", "ACGTNGGGGGGGG", "IIIIIIIIIIIII"),
            ],
        )
        matched1, matched2 = io.BytesIO(), io.BytesIO()
        oDict = {"ACGT": [matched1, matched2]}
        bc_dict = {}
        monkeypatch.setattr(dr, "plot_bc_occurance", lambda *a, **k: None)

        args = _args()
        args.output = str(tmp_path)
        dr.processPaired(args, oDict, 4, str(r1), str(r2), bc_dict, oDict)

        assert bc_dict == {"ACGT": 2}


class TestPlotBcOccurance:
    def test_writes_a_png_named_after_the_r1_sample(self, tmp_path):
        sDict = {"ACTACT": ["Sample1", ""], "TGACTG": ["Sample2", ""]}
        bc_dict = {"ACTACT": 5, "TGACTG": 3}

        dr.plot_bc_occurance("run_R1.fastq.gz", bc_dict, 2, str(tmp_path), sDict)

        assert (tmp_path / "run_fig.png").exists()

    def test_uses_the_bc_pos_label_when_set(self, tmp_path):
        sDict = {"ACTACT": ["Sample1", "A1-21"], "TGACTG": ["Sample2", ""]}
        bc_dict = {"ACTACT": 5, "TGACTG": 3}

        dr.plot_bc_occurance("run_R1.fastq.gz", bc_dict, 2, str(tmp_path), sDict)

        assert (tmp_path / "run_fig.png").exists()


class TestWrapper:
    """
    `wrapper` builds `{args.output}/{d}` by plain string formatting (not
    os.path.join), so passing an absolute `d` would nest the whole absolute
    path under args.output. Its real caller (main) always passes `d` as a
    samplesheet-relative directory name, so these tests match that by
    chdir-ing into tmp_path and using a bare relative dirname throughout.

    Also regression-tests two bugs found while writing these tests: the
    output gzip subprocesses were never closed/waited on (so their output
    was never flushed deterministically -- production "worked" only because
    the whole worker process eventually exited, forcing the OS to close the
    fds), and the auto-added "default" row for single-end samples set
    v = "unknown" (a bare string) instead of ["unknown", ""], so v[0]
    indexed into the string and produced a file named "u_R1.fastq.gz".
    """

    def test_single_end_sample_with_no_default_row_gets_one_added(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        d = Path("Sample_lib1")
        d.mkdir()
        _gzFastq(
            d / "sample_R1.fastq.gz",
            [
                ("@read1 1:N:0:1", "ACGTNTTTTTTTT", "IIIIIIIIIIIII"),
                ("@read2 1:N:0:1", "GGGGNCCCCCCCC", "IIIIIIIIIIIII"),
            ],
        )
        args = _args()
        args.output = "out"
        sDict = {"ACGT": ["Sample1", ""]}
        bc_dict = {}

        dr.wrapper((str(d), args, sDict, 4, bc_dict))

        matched = gzip.open(Path("out") / d / "Sample1_R1.fastq.gz").read()
        assert matched.decode() == "@read1_ACGT 1:N:0:1\nTTTTTTTT\n+\nIIIIIIII\n"
        default = gzip.open(Path("out") / d / "unknown_R1.fastq.gz").read()
        assert default.decode() == "@read2 1:N:0:1\nGGGGNCCCCCCCC\n+\nIIIIIIIIIIIII\n"

    def test_paired_end_sample_with_no_default_row_gets_one_added(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        d = Path("Sample_lib1")
        d.mkdir()
        _gzFastq(
            d / "sample_R1.fastq.gz",
            [("@read1 1:N:0:1", "ACGTNTTTTTTTT", "IIIIIIIIIIIII")],
        )
        _gzFastq(
            d / "sample_R2.fastq.gz",
            [("@read1 2:N:0:1", "ACGTNAAAAAAAA", "IIIIIIIIIIIII")],
        )
        args = _args()
        args.output = "out"
        sDict = {"ACGT": ["Sample1", ""]}
        bc_dict = {}
        monkeypatch.setattr(dr, "plot_bc_occurance", lambda *a, **k: None)

        dr.wrapper((str(d), args, sDict, 4, bc_dict))

        matched1 = gzip.open(Path("out") / d / "Sample1_R1.fastq.gz").read()
        assert matched1.decode() == "@read1_ACGT 1:N:0:1\nTTTTTTTT\n+\nIIIIIIII\n"
        matched2 = gzip.open(Path("out") / d / "Sample1_R2.fastq.gz").read()
        assert matched2.decode() == "@read1_ACGT 2:N:0:1\nAAAAAAAA\n+\nIIIIIIII\n"
        assert (Path("out") / d / "unknown_R2.fastq.gz").exists()

    def test_warns_and_uses_first_when_multiple_r1_files_found(
        self, tmp_path, monkeypatch, capsys
    ):
        monkeypatch.chdir(tmp_path)
        d = Path("Sample_lib1")
        d.mkdir()
        _gzFastq(
            d / "a_R1.fastq.gz",
            [("@read1 1:N:0:1", "ACGTNTTTTTTTT", "IIIIIIIIIIIII")],
        )
        _gzFastq(
            d / "b_R1.fastq.gz",
            [("@read1 1:N:0:1", "ACGTNAAAAAAAA", "IIIIIIIIIIIII")],
        )
        args = _args()
        args.output = "out"
        sDict = {"ACGT": ["Sample1", ""]}

        dr.wrapper((str(d), args, sDict, 4, {}))

        assert "more than 1 sample found" in capsys.readouterr().out
        assert (Path("out") / d / "Sample1_R1.fastq.gz").exists()


class TestMain:
    def test_end_to_end_single_end_demux(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        sampleDir = Path("Sample_lib1")
        sampleDir.mkdir()
        _gzFastq(
            sampleDir / "sample_R1.fastq.gz",
            [
                ("@read1 1:N:0:1", "ACGTNTTTTTTTT", "IIIIIIIIIIIII"),
                ("@read2 1:N:0:1", "GGGGNCCCCCCCC", "IIIIIIIIIIIII"),
            ],
        )
        table = Path("samples.txt")
        table.write_text("Sample_lib1\tACGT\tSample1\n")
        Path("out").mkdir()

        dr.main([str(table), "out", "-p", "1"])

        matched = gzip.open(Path("out") / sampleDir / "Sample1_R1.fastq.gz").read()
        assert matched.decode() == "@read1_ACGT 1:N:0:1\nTTTTTTTT\n+\nIIIIIIII\n"
        default = gzip.open(Path("out") / sampleDir / "unknown_R1.fastq.gz").read()
        assert default.decode() == "@read2 1:N:0:1\nGGGGNCCCCCCCC\n+\nIIIIIIIIIIIII\n"
