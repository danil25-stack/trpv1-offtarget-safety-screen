"""Dock a property-matched active/inactive ChEMBL ligand set into RARG or
MMP3 to check whether Vina score actually discriminates known binders from
non-binders on these off-target folds -- the missing ground-truth check for
skin_offtarget_pocket_comparison.md / OFFTARGET_SAFETY_ARTICLE_EN.md's
cross-docking section (Section 3.3), which so far only compares absolute
scores *across* receptors with no within-receptor discrimination check.

Runs entirely inside the phd-docking conda env (vina, meeko, rdkit,
biopython, openbabel -- see docker/environment-docking.yml). Intended to run
on a vast.ai instance, one process per target (RARG and MMP3 in parallel on
two instances), per project convention.

Box parameters and receptor structures are the *same* ones already validated
in the article (OFFTARGET_SAFETY_ARTICLE_EN.md Section 2.4/3.3) -- reusing
them keeps this run directly comparable to the existing cross-docking table,
not a new, differently-parameterized experiment.

Usage:
    python run_offtarget_validation.py RARG ligands.csv out_dir/
    python run_offtarget_validation.py MMP3 ligands.csv out_dir/
"""
import csv
import re
import sys
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ligand_prep import smiles_to_sdf, sdf_to_pdbqt
from run_vina import dock

from Bio.PDB import MMCIFParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa

# Same PDB structures, chain/hetero selection, and boxes already validated in
# OFFTARGET_SAFETY_ARTICLE_EN.md Section 2.4/3.3 (redocking RMSD 1.42 A for
# RARG/E9T, 2.10 A for MMP3/L04) -- reused as-is, not rederived here.
TARGETS = {
    "RARG": {
        "pdb_id": "6FX0",
        "box": {"center_x": 12.91, "center_y": -8.71, "center_z": -8.53,
                "size_x": 22.3, "size_y": 19.7, "size_z": 19.8},
        "retain_metals": False,
    },
    "MMP3": {
        "pdb_id": "1HFS",
        "box": {"center_x": 35.34, "center_y": 36.10, "center_z": 24.54,
                "size_x": 25.7, "size_y": 19.6, "size_z": 21.4},
        "retain_metals": True,
    },
}

EXHAUSTIVENESS = 32
NUM_MODES = 9
# These instances are fractional-CPU GPU rentals -- nproc reports the host's
# full core count, but the container's cgroup CFS quota is much smaller
# (observed 1.65 cores on both instances used for this run). Running 8
# concurrent Vina jobs (each defaulting to nproc threads) against that quota
# thrashes instead of parallelizing -- 20+ min with zero completions. Match
# concurrency to the real quota instead: 2 processes x 1 thread each.
MAX_WORKERS = 2
VINA_CPU = 1


class ProteinPlusMetals(Select):
    """Protein-only, optionally keeping ZN/CA hetero ions (metalloenzyme
    catalytic/structural metals) -- everything else hetero (waters, the
    co-crystallized ligand, cryoprotectant) is dropped."""

    def __init__(self, retain_metals):
        self.retain_metals = retain_metals

    def accept_residue(self, residue):
        if is_aa(residue, standard=True):
            return True
        if self.retain_metals and residue.get_resname() in ("ZN", "CA"):
            return True
        return False


def fetch_cif(pdb_id, out_path):
    if out_path.exists():
        return out_path
    subprocess.run(
        ["curl", "-sf", "-o", str(out_path), f"https://files.rcsb.org/download/{pdb_id}.cif"],
        check=True,
    )
    return out_path


def prep_receptor(target_name, work_dir):
    cfg = TARGETS[target_name]
    cif_path = fetch_cif(cfg["pdb_id"], work_dir / f"{cfg['pdb_id']}.cif")
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(cfg["pdb_id"], str(cif_path))

    receptor_pdb = work_dir / "receptor_raw.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(receptor_pdb), ProteinPlusMetals(cfg["retain_metals"]))

    receptor_base = work_dir / "receptor"
    box = cfg["box"]
    cmd = [
        "mk_prepare_receptor.py", "--read_pdb", str(receptor_pdb),
        "-o", str(receptor_base), "-p", "-v",
        "--box_center", str(box["center_x"]), str(box["center_y"]), str(box["center_z"]),
        "--box_size", str(box["size_x"]), str(box["size_y"]), str(box["size_z"]),
        "--allow_bad_res", "--default_altloc", "A",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"mk_prepare_receptor.py failed:\n{result.stdout}\n{result.stderr}")
    receptor_pdbqt = work_dir / "receptor.pdbqt"
    if not receptor_pdbqt.exists():
        # meeko names the output <o>.pdbqt -- handle either convention
        candidates = list(work_dir.glob("receptor*.pdbqt"))
        if not candidates:
            raise RuntimeError(f"no receptor pdbqt produced:\n{result.stdout}")
        receptor_pdbqt = candidates[0]
    return receptor_pdbqt, structure, cfg


def best_score(pdbqt_out: Path) -> float:
    for line in pdbqt_out.read_text().splitlines():
        m = re.match(r"REMARK VINA RESULT:\s*(-?\d+\.?\d*)", line)
        if m:
            return float(m.group(1))
    raise RuntimeError(f"no VINA RESULT found in {pdbqt_out}")


def dock_one(args):
    mid, smiles, receptor_pdbqt, box, ligand_dir, docked_dir = args
    sdf = ligand_dir / f"{mid}.sdf"
    pdbqt = ligand_dir / f"{mid}.pdbqt"
    try:
        if not pdbqt.exists():
            smiles_to_sdf(smiles, mid, sdf)
            sdf_to_pdbqt(sdf, pdbqt)
        out_pdbqt = docked_dir / f"{mid}_docked.pdbqt"
        dock(receptor_pdbqt, pdbqt, box, out_pdbqt, exhaustiveness=EXHAUSTIVENESS, num_modes=NUM_MODES, cpu=VINA_CPU)
        score = best_score(out_pdbqt)
        return mid, score, None
    except Exception as e:
        return mid, None, str(e)


def main():
    target_name, ligand_csv, out_dir = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    out_dir.mkdir(parents=True, exist_ok=True)
    ligand_dir = out_dir / "ligands_prepped"
    ligand_dir.mkdir(exist_ok=True)
    docked_dir = out_dir / "docked"
    docked_dir.mkdir(exist_ok=True)

    print(f"[{target_name}] preparing receptor...")
    receptor_pdbqt, structure, cfg = prep_receptor(target_name, out_dir)
    print(f"[{target_name}] receptor ready: {receptor_pdbqt}")

    with open(ligand_csv) as f:
        rows = list(csv.DictReader(f))
    print(f"[{target_name}] {len(rows)} ligands to dock")

    jobs = [
        (r["molecule_chembl_id"], r["smiles"], receptor_pdbqt, cfg["box"], ligand_dir, docked_dir)
        for r in rows
    ]
    results = {}
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futures = {ex.submit(dock_one, j): j[0] for j in jobs}
        done = 0
        for fut in as_completed(futures):
            mid, score, err = fut.result()
            done += 1
            if err:
                print(f"[{target_name}] {done}/{len(jobs)} {mid}: FAILED ({err})")
            else:
                print(f"[{target_name}] {done}/{len(jobs)} {mid}: {score}")
            results[mid] = (score, err)

    out_csv = out_dir / f"{target_name}_docking_validation_results.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "molecule_chembl_id", "smiles", "label", "standard_value_nM", "vina_score", "error",
        ])
        writer.writeheader()
        for r in rows:
            mid = r["molecule_chembl_id"]
            score, err = results.get(mid, (None, "not run"))
            writer.writerow({
                "molecule_chembl_id": mid, "smiles": r["smiles"], "label": r["label"],
                "standard_value_nM": r["standard_value_nM"], "vina_score": score, "error": err or "",
            })
    n_ok = sum(1 for s, e in results.values() if e is None)
    print(f"[{target_name}] done: {n_ok}/{len(jobs)} succeeded -> {out_csv}")


if __name__ == "__main__":
    main()
