"""Redock the 7 TRPV1 validation ligands (md_residence_time/ligands.py) into
one off-target receptor -- the actual cross-docking validation for
target_validation/skin_offtarget_pocket_comparison.md's fpocket-based screen.

Expects receptor.pdbqt + box.txt (Vina box config) alongside this script.
receptor.pdbqt: `mk_prepare_receptor.py --read_pdb <target>.pdb -o receptor -p
--box_center ... --box_size ... -v` (box from target_validation/
box_from_pocket.py run on the target's fpocket pocket1_vert.pqr).

Run one instance of this per off-target (compute is cheap enough to just
split across a vast.ai instance per target rather than looping in-process --
see project convention: all non-trivial compute runs on vast.ai, not
locally).
"""
import re
import csv
from pathlib import Path

from ligand_prep import smiles_to_sdf, sdf_to_pdbqt
from run_vina import dock, read_box
from ligands import LIGANDS

HERE = Path(__file__).resolve().parent
LIGAND_DIR = HERE / "ligands_prepped"
LIGAND_DIR.mkdir(exist_ok=True)
OUT_DIR = HERE / "docked"
OUT_DIR.mkdir(exist_ok=True)


def best_score(pdbqt_out: Path) -> float:
    for line in pdbqt_out.read_text().splitlines():
        m = re.match(r"REMARK VINA RESULT:\s*(-?\d+\.?\d*)", line)
        if m:
            return float(m.group(1))
    raise RuntimeError(f"no VINA RESULT found in {pdbqt_out}")


def main():
    receptor = HERE / "receptor.pdbqt"
    box = read_box(HERE / "box.txt")

    results = []
    for name, (smiles, ec50, cid) in LIGANDS.items():
        sdf = LIGAND_DIR / f"{name}.sdf"
        pdbqt = LIGAND_DIR / f"{name}.pdbqt"
        if not pdbqt.exists():
            print(f"preparing ligand {name}...")
            smiles_to_sdf(smiles, name, sdf)
            sdf_to_pdbqt(sdf, pdbqt)

        out_pdbqt = OUT_DIR / f"{name}_docked.pdbqt"
        print(f"docking {name}...")
        try:
            dock(receptor, pdbqt, box, out_pdbqt, exhaustiveness=32, num_modes=9)
            score = best_score(out_pdbqt)
        except Exception as e:
            print(f"  FAILED: {e}")
            score = None
        results.append({"ligand": name, "ec50_nM": ec50, "vina_score": score})
        print(f"  {name}: {score}")

    out_csv = HERE / "redocking_results.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["ligand", "ec50_nM", "vina_score"])
        writer.writeheader()
        writer.writerows(results)
    print(f"\nWrote {out_csv}")
    for r in results:
        print(r)


if __name__ == "__main__":
    main()
