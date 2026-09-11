"""Local (pocket-residue-only) ESM-2 embedding comparison of TRPV1 vs its
five paralogs (TRPV2-6), as a quantitative complement to the qualitative
"vanilloid pocket is conserved across TRPV1-3, TRPV4 diverges by point
substitution" claim in OFFTARGET_SAFETY_WRITEUP.md Section 2.3.

Whole-sequence mean-pooled ESM-2 similarity would just re-derive the same
whole-protein identity ranking already in Table 1 (that's the thing the
writeup already flagged as a poor proxy for pocket-level risk) -- so
instead this extracts per-residue ESM-2 representations only at the four
vanilloid-pocket-lining positions (Y511/S512/T550/E570 in TRPV1, mapped to
each paralog's aligned position by pocket_residue_alignment.py) and
compares those local vectors directly.

Run on a GPU instance (pocket_residue_alignment.py's alignment step is
cheap/local; this script's ESM-2 forward pass is the actual compute and
runs on vast.ai per project convention). Requires: pip install fair-esm
torch.

Usage: python esm_pocket_embedding_similarity.py
Inputs: pocket_residue_alignment.csv (from pocket_residue_alignment.py)
        data_pocket_esm/fasta/*.fasta (same cache used by that script)
Output: esm_pocket_similarity.csv
"""
import csv
from pathlib import Path

import torch
import esm

HERE = Path(__file__).resolve().parent
FASTA_DIR = HERE / "data_pocket_esm" / "fasta"
ALIGNMENT_CSV = HERE / "pocket_residue_alignment.csv"
OUT_CSV = HERE / "esm_pocket_similarity.csv"

UNIPROT = {
    "TRPV1": "Q8NER1",
    "TRPV2": "Q9Y5S1",
    "TRPV3": "Q8NET8",
    "TRPV4": "Q9HBA0",
    "TRPV5": "Q9NQA5",
    "TRPV6": "Q9H1D0",
}

# Window half-width (residues) averaged around each pocket position, to
# smooth over single-residue embedding noise and any off-by-one alignment
# uncertainty -- 2 gives a 5-residue window, small enough to stay local.
WINDOW_HALF = 2
MODEL_NAME = "esm2_t33_650M_UR50D"


def load_sequence(name: str) -> str:
    fasta = FASTA_DIR / f"{UNIPROT[name]}.fasta"
    lines = fasta.read_text().splitlines()
    return "".join(l.strip() for l in lines if not l.startswith(">"))


def load_model():
    model, alphabet = getattr(esm.pretrained, MODEL_NAME)()
    model.eval()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(device)
    return model, alphabet, device


def embed_sequence(model, alphabet, device, seq: str) -> torch.Tensor:
    """Returns per-residue representation, shape (len(seq), embed_dim),
    from the model's final layer, batch/BOS/EOS tokens stripped."""
    batch_converter = alphabet.get_batch_converter()
    _, _, tokens = batch_converter([("query", seq)])
    tokens = tokens.to(device)
    with torch.no_grad():
        out = model(tokens, repr_layers=[model.num_layers], return_contacts=False)
    reps = out["representations"][model.num_layers][0]
    # strip BOS (position 0) and EOS (last position)
    return reps[1 : 1 + len(seq)].cpu()


def local_pocket_vector(reps: torch.Tensor, positions: list) -> torch.Tensor:
    """Mean-pool a window around each 1-indexed position, then concatenate
    across the 4 pocket positions into one vector."""
    vecs = []
    for pos in positions:
        if pos is None:
            vecs.append(torch.zeros(reps.shape[1]))
            continue
        i = pos - 1
        lo, hi = max(0, i - WINDOW_HALF), min(reps.shape[0], i + WINDOW_HALF + 1)
        vecs.append(reps[lo:hi].mean(dim=0))
    return torch.cat(vecs)


def main():
    with open(ALIGNMENT_CSV) as f:
        alignment_rows = list(csv.DictReader(f))

    print(f"loading {MODEL_NAME} ...")
    model, alphabet, device = load_model()
    print(f"using device={device}")

    seqs = {name: load_sequence(name) for name in UNIPROT}
    reps = {}
    for name, seq in seqs.items():
        print(f"embedding {name} ({len(seq)} aa) ...")
        reps[name] = embed_sequence(model, alphabet, device, seq)

    # whole-protein mean-pooled baseline, for explicit comparison against
    # the local pocket-residue result (expected to be uninformative/
    # redundant with Table 1's identity ranking -- see module docstring)
    trpv1_global = reps["TRPV1"].mean(dim=0)

    rows = []
    for name in ["TRPV2", "TRPV3", "TRPV4", "TRPV5", "TRPV6"]:
        paralog_rows = [r for r in alignment_rows if r["paralog"] == name]
        paralog_rows.sort(key=lambda r: int(r["trpv1_pos"]))

        # Compare only positions aligned (non-gap) in *this* paralog. A
        # gap position zero-padded into the concatenated vector would
        # inject a whole zero-block into one operand and mechanically
        # depress cosine similarity regardless of true local similarity
        # at the other positions -- an artifact, not a biological signal
        # (caught when TRPV6, which gaps at S512, came out a clear
        # outlier -- confirmed here it isn't just that).
        present = [r for r in paralog_rows if r["paralog_pos"]]
        trpv1_positions_present = [int(r["trpv1_pos"]) for r in present]
        paralog_positions_present = [int(r["paralog_pos"]) for r in present]

        trpv1_vec = local_pocket_vector(reps["TRPV1"], trpv1_positions_present)
        paralog_vec = local_pocket_vector(reps[name], paralog_positions_present)
        paralog_global = reps[name].mean(dim=0)

        local_cos = torch.nn.functional.cosine_similarity(trpv1_vec, paralog_vec, dim=0).item()
        global_cos = torch.nn.functional.cosine_similarity(trpv1_global, paralog_global, dim=0).item()

        rows.append(
            {
                "paralog": name,
                "local_pocket_cosine_sim": round(local_cos, 4),
                "global_meanpool_cosine_sim": round(global_cos, 4),
                "n_positions_compared": len(present),
                "n_gap_positions": len(paralog_rows) - len(present),
            }
        )
        print(
            f"{name}: local_pocket_cos={local_cos:.4f}  global_meanpool_cos={global_cos:.4f}  "
            f"(n_compared={len(present)}, n_gap={len(paralog_rows) - len(present)})"
        )

    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {OUT_CSV}")


if __name__ == "__main__":
    main()
