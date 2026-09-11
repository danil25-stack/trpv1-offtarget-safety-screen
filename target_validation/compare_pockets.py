"""Parse fpocket info.txt outputs and compare every detected pocket's
descriptor fingerprint against the TRPV1 vanilloid pocket reference.

No docking/ligand poses involved -- pure cavity geometry+chemistry comparison.
"""
import re
import csv
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data" / "processed" / "skin_offtarget_screen"
STRUCT_DIR = DATA_DIR / "fpocket_out"
TRPV1_OUT_DIR = DATA_DIR / "fpocket_out_trpv1_reference" / "TRPV1_AD_out"
TRPV1_INFO = TRPV1_OUT_DIR / "TRPV1_AD_info.txt"
OUT_CSV = DATA_DIR / "pocket_comparison_results.csv"

# Descriptors used for the fingerprint: chemically meaningful, scale-comparable
# across proteins (excludes raw counts like "Number of Alpha Spheres" that
# scale with pocket size redundantly with Volume).
FIELDS = [
    "Druggability Score",
    "Volume",
    "Hydrophobicity score",
    "Polarity score",
    "Charge score",
    "Proportion of polar atoms",
    "Apolar alpha sphere proportion",
    "Mean local hydrophobic density",
]

TRPV1_VANILLOID_RESIDUES = {("A", 511), ("A", 512), ("A", 550), ("A", 570)}


def parse_info_file(path):
    text = "\n" + path.read_text()
    blocks = re.split(r"\nPocket (\d+) :\n", text)
    pockets = {}
    for i in range(1, len(blocks), 2):
        pnum = int(blocks[i])
        body = blocks[i + 1]
        vals = {}
        for line in body.strip().split("\n"):
            if ":" not in line:
                continue
            key, val = line.split(":", 1)
            key = key.strip().lstrip("\t").strip()
            val = val.strip()
            try:
                vals[key] = float(val)
            except ValueError:
                pass
        pockets[pnum] = vals
    return pockets


def pocket_residues(atm_pdb_path):
    residues = set()
    for line in atm_pdb_path.read_text().splitlines():
        if line.startswith("ATOM"):
            chain = line[21].strip()
            resnum = int(line[22:26])
            residues.add((chain, resnum))
    return residues


def vector(pocket_vals):
    return np.array([pocket_vals.get(f, np.nan) for f in FIELDS])


def main():
    trpv1_pockets = parse_info_file(TRPV1_INFO)
    trpv1_atm_dir = TRPV1_OUT_DIR / "pockets"
    ref_pocket_num = None
    for pnum in trpv1_pockets:
        atm = trpv1_atm_dir / f"pocket{pnum}_atm.pdb"
        if not atm.exists():
            continue
        residues = pocket_residues(atm)
        if TRPV1_VANILLOID_RESIDUES <= residues:
            ref_pocket_num = pnum
            break
    assert ref_pocket_num is not None, "could not relocate vanilloid pocket"
    ref_vals = trpv1_pockets[ref_pocket_num]
    print(f"TRPV1 reference: pocket {ref_pocket_num}, druggability={ref_vals['Druggability Score']:.3f}, volume={ref_vals['Volume']:.1f}")
    ref_vec = vector(ref_vals)

    # collect all pockets (all targets + TRPV1 itself) to compute a shared
    # normalization (z-score) so no single descriptor's raw scale dominates.
    all_records = []  # (gene, pocket_num, vals)
    for struct_dir in sorted(STRUCT_DIR.glob("*_out")):
        gene = struct_dir.name[: -len("_out")]
        info_path = struct_dir / f"{gene}_info.txt"
        if not info_path.exists():
            print(f"  WARNING: no info file for {gene}")
            continue
        pockets = parse_info_file(info_path)
        for pnum, vals in pockets.items():
            all_records.append((gene, pnum, vals))
    all_records.append(("TRPV1", ref_pocket_num, ref_vals))

    matrix = np.array([vector(vals) for _, _, vals in all_records])
    mean = np.nanmean(matrix, axis=0)
    std = np.nanstd(matrix, axis=0)
    std[std == 0] = 1.0
    norm_matrix = (matrix - mean) / std
    norm_matrix = np.nan_to_num(norm_matrix, nan=0.0)

    ref_idx = len(all_records) - 1
    ref_norm = norm_matrix[ref_idx]

    results = []
    for i, (gene, pnum, vals) in enumerate(all_records):
        if gene == "TRPV1":
            continue
        dist = np.linalg.norm(norm_matrix[i] - ref_norm)
        cos_sim = np.dot(norm_matrix[i], ref_norm) / (np.linalg.norm(norm_matrix[i]) * np.linalg.norm(ref_norm) + 1e-9)
        results.append({
            "gene": gene, "pocket_num": pnum, "euclidean_dist": dist, "cosine_sim": cos_sim,
            **{f: vals.get(f) for f in FIELDS},
        })

    # best (closest) pocket ANYWHERE on the structure per gene -- exploratory,
    # inflated by proteins with many detected pockets (multiple-comparisons).
    best_per_gene = {}
    for r in results:
        g = r["gene"]
        if g not in best_per_gene or r["euclidean_dist"] < best_per_gene[g]["euclidean_dist"]:
            best_per_gene[g] = r
    ranked_any = sorted(best_per_gene.values(), key=lambda r: r["euclidean_dist"])

    # top-ranked pocket per gene by fpocket's own internal Score (pocket #1,
    # i.e. what fpocket itself considers the single most likely real cavity --
    # this is what was validated against TRPV1's known vanilloid pocket, so
    # it's the fairer one-shot-per-target comparison, not multiple-testing inflated.
    top1_per_gene = {r["gene"]: r for r in results if r["pocket_num"] == 1}
    ranked_top1 = sorted(top1_per_gene.values(), key=lambda r: r["euclidean_dist"])

    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["gene", "pocket_num", "euclidean_dist", "cosine_sim"] + FIELDS)
        writer.writeheader()
        writer.writerows(ranked_any)

    top1_csv = OUT_CSV.with_name("pocket_comparison_top1_only.csv")
    with open(top1_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["gene", "pocket_num", "euclidean_dist", "cosine_sim"] + FIELDS)
        writer.writeheader()
        writer.writerows(ranked_top1)

    print(f"\n=== BEST POCKET ANYWHERE ON STRUCTURE (exploratory, multiple-comparisons-inflated) ===")
    print(f"{'Gene':10s} {'best pocket':12s} {'eucl.dist':10s} {'cos.sim':8s} {'drug.score':10s} {'volume':8s}")
    for r in ranked_any:
        print(f"{r['gene']:10s} #{r['pocket_num']:<11d} {r['euclidean_dist']:<10.3f} {r['cosine_sim']:<8.3f} {r['Druggability Score']:<10.3f} {r['Volume']:<8.1f}")

    print(f"\n=== TOP FPOCKET-RANKED POCKET ONLY (fairer one-shot comparison) ===")
    print(f"{'Gene':10s} {'eucl.dist':10s} {'cos.sim':8s} {'drug.score':10s} {'volume':8s}")
    for r in ranked_top1:
        print(f"{r['gene']:10s} {r['euclidean_dist']:<10.3f} {r['cosine_sim']:<8.3f} {r['Druggability Score']:<10.3f} {r['Volume']:<8.1f}")

    print(f"\nWrote {OUT_CSV} and {top1_csv}")


if __name__ == "__main__":
    main()
