from BRB.demultiplex_relacs import edit_distance, matchSample


def test_edit_distance_exact_and_one_mismatch():
    assert edit_distance("ACTACT", "ACTACT", max_dist=1) == 0
    assert edit_distance("ACTACT", "ACTGCT", max_dist=1) == 1
    assert edit_distance("ACTACT", "TTTTTT", max_dist=1) > 1


def test_match_sample_allows_one_mismatch():
    sample_sheet = {"ACTACT": ["Sample1", ""]}

    # exact barcode match
    assert matchSample("ACTACTNN", None, sample_sheet, bcLen=6, umiLength=0) == (
        "ACTACT",
        True,
    )

    # one mismatch should still match
    assert matchSample("ACTGCTNN", None, sample_sheet, bcLen=6, umiLength=0) == (
        "ACTACT",
        True,
    )

    # two mismatches should fall through to default
    assert matchSample("TTTTTTNN", None, sample_sheet, bcLen=6, umiLength=0) == (
        "default",
        False,
    )


def test_match_sample_with_paired_read_check():
    sample_sheet = {"ACTACT": ["Sample1", ""]}
    result = matchSample("ACTACTNN", "ACTGCTNN", sample_sheet, bcLen=6, umiLength=0)
    assert result == ("ACTACT", True)
