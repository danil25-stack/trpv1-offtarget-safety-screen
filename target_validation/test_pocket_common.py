"""Regression test for the zero-padding bug described in
pocket_common.filter_common_positions' docstring and
esm_pocket_similarity.md: a paralog's gapped pocket position must be
dropped from the comparison entirely, never zero-padded in. No torch/
fair-esm/network needed -- run with: pytest test_pocket_common.py
"""
from pocket_common import filter_common_positions


def _row(trpv1_pos, paralog_pos):
    return {"trpv1_pos": str(trpv1_pos), "paralog_pos": str(paralog_pos) if paralog_pos else ""}


def test_no_gaps_keeps_all_positions():
    rows = [_row(511, 469), _row(512, 470), _row(550, 508), _row(570, 528)]
    trpv1_pos, paralog_pos = filter_common_positions(rows)
    assert trpv1_pos == [511, 512, 550, 570]
    assert paralog_pos == [469, 470, 508, 528]


def test_gap_is_excluded_not_zero_padded():
    """The exact TRPV6/S512 case: one row has no aligned paralog position."""
    rows = [_row(511, 465), _row(512, None), _row(550, 503), _row(570, 523)]
    trpv1_pos, paralog_pos = filter_common_positions(rows)

    # the gapped TRPV1 position (512) must not appear on either side
    assert 512 not in trpv1_pos
    assert len(trpv1_pos) == len(paralog_pos) == 3
    assert trpv1_pos == [511, 550, 570]
    assert paralog_pos == [465, 503, 523]

    # both lists are the same length as each other by construction --
    # this is what makes it safe to zip them into equal-length vectors
    # without ever needing a placeholder for the missing position
    assert len(trpv1_pos) == len(paralog_pos)


def test_output_sorted_by_trpv1_position_regardless_of_input_order():
    rows = [_row(570, 523), _row(511, 465), _row(550, 503)]
    trpv1_pos, _ = filter_common_positions(rows)
    assert trpv1_pos == sorted(trpv1_pos)


def test_all_gaps_returns_empty_lists():
    rows = [_row(511, None), _row(512, None)]
    trpv1_pos, paralog_pos = filter_common_positions(rows)
    assert trpv1_pos == []
    assert paralog_pos == []
