import pytest

from BRB.editdistance import edit_distance, eval as editdistance_eval


@pytest.mark.parametrize(
    "a,b,expected",
    [
        ("", "", 0),
        ("abc", "abc", 0),
        ("abc", "ab", 1),
        ("kitten", "sitting", 3),
        ("GATTACA", "GCATGCU", 4),
    ],
)
def test_edit_distance(a, b, expected):
    assert edit_distance(a, b) == expected
    assert editdistance_eval(a, b) == expected


def test_edit_distance_is_symmetric():
    assert edit_distance("barcode", "barc0de") == edit_distance("barc0de", "barcode")


def test_edit_distance_max_dist_short_circuits():
    assert edit_distance("abc", "xyz", max_dist=1) == 2
    assert edit_distance("abc", "ab", max_dist=1) == 1
