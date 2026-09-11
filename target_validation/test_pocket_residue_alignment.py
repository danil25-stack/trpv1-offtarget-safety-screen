"""Unit tests for align_and_map using small synthetic sequences -- no
network fetch, no real TRPV1/paralog data needed. Requires biopython
(already a base project dependency); run with:
    pytest test_pocket_residue_alignment.py
"""
from pocket_residue_alignment import align_and_map


def test_identical_sequence_maps_positions_to_themselves():
    seq = "MAAVLGWTNMLYYGRK"
    #      1234567890123456
    mapping = align_and_map(seq, seq, positions_of_interest={7, 8})
    assert mapping == {7: 7, 8: 8}


def test_substitution_without_indel_keeps_position_number():
    trpv1_like = "MAAVLGWTNMLYYGRK"
    #                    ^ position 8 = T
    paralog_like = "MAAVLGWLNMLYYGRK"
    #                      ^ position 8 = L (substituted, no length change)
    mapping = align_and_map(trpv1_like, paralog_like, positions_of_interest={7, 8})
    assert mapping[7] == 7  # unaffected flanking residue
    assert mapping[8] == 8  # substituted but still aligned 1:1
    assert paralog_like[mapping[8] - 1] == "L"


def test_upstream_deletion_shifts_downstream_position():
    trpv1_like = "MAAVLGWTNMLYYGRK"
    #              1234567890123456
    # remove 2 residues upstream of position 8 -> everything from there
    # on should shift left by 2 in the paralog's own numbering
    paralog_like = "MAVLGWTNMLYYGRK"  # "AA" -> "A": one residue removed
    mapping = align_and_map(trpv1_like, paralog_like, positions_of_interest={8})
    assert mapping[8] is not None
    assert paralog_like[mapping[8] - 1] == "T"


def test_gap_at_pocket_position_maps_to_none():
    """If the paralog is missing the residue aligned to a pocket
    position entirely (a true indel at that exact site), the mapping
    must be None there -- this is what pocket_common.filter_common_positions
    is later required to drop rather than zero-pad."""
    trpv1_like = "MAAVLGWTNMLYYGRK"
    #                    ^ position 8
    # delete exactly the residue that aligns to position 8
    paralog_like = "MAAVLGWNMLYYGRK"
    mapping = align_and_map(trpv1_like, paralog_like, positions_of_interest={7, 8, 9})
    assert mapping[7] is not None
    assert mapping[9] is not None
    assert mapping[8] is None


def test_default_positions_of_interest_is_the_real_pocket_residues():
    from pocket_residue_alignment import POCKET_RESIDUES

    assert POCKET_RESIDUES == {511: "Y", 512: "S", 550: "T", 570: "E"}
