"""Pure, dependency-light helpers shared by the pocket-residue-alignment
and ESM-2 embedding-comparison scripts.

Deliberately has no torch/fair-esm/requests/biopython imports so it (and
its regression tests) run in a plain Python environment -- no GPU, no
model download, no network -- unlike the two scripts that use it.
"""


def filter_common_positions(paralog_rows: list) -> tuple:
    """Given a paralog's rows from pocket_residue_alignment.csv (each a
    dict with at least 'trpv1_pos' and 'paralog_pos', the latter possibly
    an empty string for a gap), return (trpv1_positions, paralog_positions)
    -- two same-length int lists containing *only* the positions aligned
    (non-gap) in this paralog, sorted by trpv1_pos.

    This is the fix for a real bug: an earlier version zero-padded a
    gapped position into the compared embedding vector instead of
    dropping it, which mechanically depressed cosine similarity for any
    paralog with a gap (caught via TRPV6's S512 gap -- see
    esm_pocket_similarity.md). A gap must be *excluded* from both sides,
    never zero-filled -- that's the behavior this function exists to lock
    in, and what test_pocket_common.py's regression test checks directly.
    """
    rows_sorted = sorted(paralog_rows, key=lambda r: int(r["trpv1_pos"]))
    present = [r for r in rows_sorted if r["paralog_pos"]]
    trpv1_positions = [int(r["trpv1_pos"]) for r in present]
    paralog_positions = [int(r["paralog_pos"]) for r in present]
    return trpv1_positions, paralog_positions
