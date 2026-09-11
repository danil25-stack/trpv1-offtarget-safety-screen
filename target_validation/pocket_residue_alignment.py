"""Map TRPV1's four vanilloid-pocket-lining residues (Y511, S512, T550,
E570 -- confirmed against UniProt Q8NER1 canonical numbering, matching the
residues used throughout md_residence_time/STATUS.md and
OFFTARGET_SAFETY_WRITEUP.md, Section 2.3/3.2) onto the aligned position in
each TRPV2-6 paralog via pairwise global alignment.

This is the prerequisite for a local (pocket-residue-only) ESM-2 embedding
comparison -- comparing residue numbers directly across paralogs is wrong
whenever there's an indel between TRPV1 and a paralog upstream of the
pocket, which a single global alignment resolves.

Output: target_validation/pocket_residue_alignment.csv (one row per
paralog per pocket residue: TRPV1 position/residue, aligned paralog
position/residue, a +/-5-residue window in each sequence for manual
sanity-checking, and whether the aligned residue actually matches the
physicochemical class of the TRPV1 residue).
"""
import csv
from pathlib import Path

import requests
from Bio import Align
from Bio.Align import substitution_matrices

HERE = Path(__file__).resolve().parent
FASTA_DIR = HERE / "data_pocket_esm" / "fasta"
FASTA_DIR.mkdir(parents=True, exist_ok=True)

UNIPROT = {
    "TRPV1": "Q8NER1",
    "TRPV2": "Q9Y5S1",
    "TRPV3": "Q8NET8",
    "TRPV4": "Q9HBA0",
    "TRPV5": "Q9NQA5",
    "TRPV6": "Q9H1D0",
}

# TRPV1 (Q8NER1) vanilloid-pocket-lining residues, 1-indexed, confirmed
# against the canonical sequence in this script's header comment.
POCKET_RESIDUES = {511: "Y", 512: "S", 550: "T", 570: "E"}

AA_CLASS = {
    **{a: "aromatic" for a in "FWY"},
    **{a: "polar" for a in "STNQCH"},
    **{a: "acidic" for a in "DE"},
    **{a: "basic" for a in "KRH"},
    **{a: "aliphatic" for a in "AVLIMG"},
    "P": "proline",
}


def fetch_sequence(accession: str) -> str:
    cache = FASTA_DIR / f"{accession}.fasta"
    if not cache.exists():
        r = requests.get(f"https://rest.uniprot.org/uniprotkb/{accession}.fasta", timeout=30)
        r.raise_for_status()
        cache.write_text(r.text)
    lines = cache.read_text().splitlines()
    return "".join(l.strip() for l in lines if not l.startswith(">"))


def align_and_map(trpv1_seq: str, other_seq: str, positions_of_interest=None) -> dict:
    """Global pairwise alignment; returns {trpv1_1indexed_pos: other_1indexed_pos_or_None}
    for each position in `positions_of_interest` (defaults to the real
    vanilloid-pocket positions, POCKET_RESIDUES) -- a value of None means
    that TRPV1 position aligns to a gap in `other_seq`.

    `positions_of_interest` is a parameter (not hardcoded to
    POCKET_RESIDUES) specifically so this function is unit-testable with
    small synthetic sequences/positions -- see
    test_pocket_residue_alignment.py.
    """
    if positions_of_interest is None:
        positions_of_interest = set(POCKET_RESIDUES)
    else:
        positions_of_interest = set(positions_of_interest)

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "global"
    alignment = aligner.align(trpv1_seq, other_seq)[0]
    aligned_trpv1, aligned_other = alignment[0], alignment[1]

    mapping = {}
    trpv1_pos = 0  # 1-indexed position in the *ungapped* TRPV1 sequence
    other_pos = 0
    for a_res, b_res in zip(aligned_trpv1, aligned_other):
        if a_res != "-":
            trpv1_pos += 1
        if b_res != "-":
            other_pos += 1
        if a_res != "-" and trpv1_pos in positions_of_interest:
            mapping[trpv1_pos] = other_pos if b_res != "-" else None
    return mapping


def window(seq: str, pos_1indexed, half: int = 5) -> str:
    if pos_1indexed is None:
        return "(gap)"
    i = pos_1indexed - 1
    lo, hi = max(0, i - half), min(len(seq), i + half + 1)
    marked = seq[lo:i] + "[" + seq[i] + "]" + seq[i + 1 : hi]
    return marked


def main():
    seqs = {name: fetch_sequence(acc) for name, acc in UNIPROT.items()}
    trpv1_seq = seqs["TRPV1"]

    for pos, res in POCKET_RESIDUES.items():
        assert trpv1_seq[pos - 1] == res, f"TRPV1 position {pos} is {trpv1_seq[pos-1]!r}, expected {res!r}"

    rows = []
    for name in ["TRPV2", "TRPV3", "TRPV4", "TRPV5", "TRPV6"]:
        mapping = align_and_map(trpv1_seq, seqs[name])
        for trpv1_pos, trpv1_res in POCKET_RESIDUES.items():
            other_pos = mapping.get(trpv1_pos)
            other_res = seqs[name][other_pos - 1] if other_pos else None
            same_class = (
                other_res is not None
                and AA_CLASS.get(trpv1_res) == AA_CLASS.get(other_res)
            )
            rows.append(
                {
                    "paralog": name,
                    "trpv1_pos": trpv1_pos,
                    "trpv1_res": trpv1_res,
                    "paralog_pos": other_pos or "",
                    "paralog_res": other_res or "",
                    "same_physicochem_class": same_class,
                    "trpv1_window": window(trpv1_seq, trpv1_pos),
                    "paralog_window": window(seqs[name], other_pos),
                }
            )

    out_csv = HERE / "pocket_residue_alignment.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {out_csv}")
    for r in rows:
        print(
            f"{r['paralog']:6s} TRPV1 {r['trpv1_res']}{r['trpv1_pos']} -> "
            f"{r['paralog_res'] or '-'}{r['paralog_pos'] or ''}  "
            f"same_class={r['same_physicochem_class']}  {r['paralog_window']}"
        )


if __name__ == "__main__":
    main()
