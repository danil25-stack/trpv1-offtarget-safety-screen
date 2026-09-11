"""Re-run the MMP3 actives/inactives docking validation
(run_offtarget_validation.py) with AD4Zn scoring instead of plain Vina.

Motivation: the plain-Vina MMP3 run (data/processed/skin_offtarget_screen/
MMP3_docking_validation_results.csv) scored ROC-AUC = 0.503 against known
ChEMBL actives/inactives -- indistinguishable from random. Standard Vina has
no metal-specific potential; retaining the catalytic Zn2+/Ca2+ ions as inert
rigid-body atoms (what OFFTARGET_SAFETY_ARTICLE_EN.md's "corrected" MMP3
receptor did) isn't the same as actually scoring zinc coordination chemistry.
AD4Zn (Santos-Martins et al. 2014, 10.1021/ci500209e) adds a directional
tetrahedral-coordination potential via zinc pseudo-atoms; Vina >=1.2 can use
these precomputed AD4 affinity maps directly via --scoring ad4 --maps,
without needing the AutoDock4 GA search engine itself.

Pipeline (see docking_zinc.html in the Vina docs for the reference procedure
this follows, minus the MolKit/AutoDockTools-dependent gpf writer, replaced
here by prepare_gpf4zn_lite.py):
  1. mk_prepare_receptor.py -> receptor.pdbqt (protein + Zn2+/Ca2+, same
     selection as the plain-Vina MMP3 run)
  2. zinc_pseudo.py -> receptor_TZ.pdbqt (adds tetrahedral zinc pseudo-atoms)
  3. prepare ligand pdbqts for the full 60-ligand set (same RDKit/Meeko
     pipeline as before)
  4. prepare_gpf4zn_lite.py -> protein_tz.gpf (grid covering the union of
     atom types across all 60 ligands, so one map set serves the whole batch)
  5. autogrid4 -> affinity maps (one-time cost, not per-ligand)
  6. vina --scoring ad4 --maps protein_tz --ligand X.pdbqt, once per ligand
     (parallelized across MAX_WORKERS processes)

Usage:
    python run_ad4zn_validation.py MMP3_docking_validation_set.csv ad4zn_run/
"""
import csv
import re
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ligand_prep import smiles_to_sdf
from run_offtarget_validation import ProteinPlusMetals, fetch_cif

from Bio.PDB import MMCIFParser, PDBIO

PDB_ID = "1HFS"
BOX = {"center_x": 35.34, "center_y": 36.10, "center_z": 24.54,
       "size_x": 25.7, "size_y": 19.6, "size_z": 21.4}
MAX_WORKERS = int(__import__("os").environ.get("AD4ZN_MAX_WORKERS", "2"))
VINA_CPU = __import__("os").environ.get("AD4ZN_VINA_CPU")  # None = vina's own default (nproc-detected)
HERE = Path(__file__).resolve().parent


def run(cmd, **kw):
    print(f"$ {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if result.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed:\n{result.stdout}\n{result.stderr}")
    return result


def prep_receptor_tz(work_dir):
    cif_path = fetch_cif(PDB_ID, work_dir / f"{PDB_ID}.cif")
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(PDB_ID, str(cif_path))
    receptor_pdb = work_dir / "receptor_raw.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(receptor_pdb), ProteinPlusMetals(retain_metals=True))

    receptor_base = work_dir / "receptor"
    run(["mk_prepare_receptor.py", "--read_pdb", str(receptor_pdb),
         "-o", str(receptor_base), "-p", "--allow_bad_res"])
    receptor_pdbqt = work_dir / "receptor.pdbqt"

    receptor_tz = work_dir / "receptor_TZ.pdbqt"
    run(["python3", str(HERE / "zinc_pseudo.py"), "-r", str(receptor_pdbqt), "-o", str(receptor_tz)])
    return receptor_tz


def sdf_to_pdbqt_rigid(sdf_path, out_path):
    # AD4Zn.dat only defines the classic AutoDock4 atom-type set -- it
    # predates Vina 1.2's macrocycle-flexibility feature, which emits
    # synthetic ring-glue atom types (e.g. "CG0"/"G0") that autogrid4
    # rejects outright ("unknown ligand atom type"). Keep macrocycles rigid
    # for this pipeline so every ligand only uses classic AD4 atom types.
    result = subprocess.run(
        ["mk_prepare_ligand.py", "-i", str(sdf_path), "-o", str(out_path), "--rigid_macrocycles"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"mk_prepare_ligand.py failed:\n{result.stdout}\n{result.stderr}")
    return out_path


def prep_ligands(rows, work_dir):
    ligand_dir = work_dir / "ligands_prepped"
    ligand_dir.mkdir(exist_ok=True)
    for r in rows:
        mid = r["molecule_chembl_id"]
        pdbqt = ligand_dir / f"{mid}.pdbqt"
        if pdbqt.exists():
            continue
        sdf = ligand_dir / f"{mid}.sdf"
        try:
            smiles_to_sdf(r["smiles"], mid, sdf)
            sdf_to_pdbqt_rigid(sdf, pdbqt)
        except Exception as e:
            print(f"  ligand prep FAILED for {mid}: {e}")
    return ligand_dir


def build_maps(receptor_tz, ligand_dir, work_dir):
    # Both this script and autogrid4 are run with cwd=work_dir below, and
    # the .gpf format embeds bare relative filenames (autogrid4 convention:
    # run from the directory containing every referenced file) -- so every
    # path handed to prepare_gpf4zn_lite.py must already be relative to
    # work_dir, not to this script's original cwd, or autogrid4 resolves a
    # doubled-up path (work_dir/work_dir/...) and fails to find the files.
    receptor_rel = receptor_tz.relative_to(work_dir)
    ligand_glob_rel = (ligand_dir / "*.pdbqt").relative_to(work_dir)
    run([
        "python3", str(HERE / "prepare_gpf4zn_lite.py"),
        "--receptor", str(receptor_rel),
        "--ligand-glob", str(ligand_glob_rel),
        "--center", str(BOX["center_x"]), str(BOX["center_y"]), str(BOX["center_z"]),
        "--size", str(BOX["size_x"]), str(BOX["size_y"]), str(BOX["size_z"]),
        "--name", "protein_tz", "--out", "protein_tz.gpf",
    ], cwd=work_dir)
    # AD4Zn.dat must be in the working directory (referenced by relative
    # path inside the .gpf, per autogrid4 convention)
    ad4zn = HERE / "AD4Zn.dat"
    (work_dir / "AD4Zn.dat").write_text(ad4zn.read_text())
    run(["autogrid4", "-p", "protein_tz.gpf", "-l", "protein_tz.glg"], cwd=work_dir)


def best_score(pdbqt_out: Path) -> float:
    for line in pdbqt_out.read_text().splitlines():
        m = re.match(r"REMARK VINA RESULT:\s*(-?\d+\.?\d*)", line)
        if m:
            return float(m.group(1))
    raise RuntimeError(f"no VINA RESULT found in {pdbqt_out}")


def dock_one(args):
    mid, ligand_pdbqt, work_dir, docked_dir = args
    out_pdbqt = docked_dir / f"{mid}_docked.pdbqt"
    try:
        # vina runs with cwd=work_dir (matches --maps protein_tz, a bare
        # relative name written by build_maps) -- ligand/out paths must be
        # relative to work_dir too, or they double up with the work_dir
        # prefix already baked into ligand_pdbqt/out_pdbqt (see build_maps).
        cmd = [
            "vina", "--ligand", str(ligand_pdbqt.relative_to(work_dir)), "--maps", "protein_tz",
            "--scoring", "ad4", "--exhaustiveness", "32", "--num_modes", "9",
            "--out", str(out_pdbqt.relative_to(work_dir)),
        ]
        if VINA_CPU:
            cmd += ["--cpu", VINA_CPU]
        run(cmd, cwd=work_dir)
        return mid, best_score(out_pdbqt), None
    except Exception as e:
        return mid, None, str(e)


def main():
    ligand_csv, out_dir = Path(sys.argv[1]), Path(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)
    docked_dir = out_dir / "docked"
    docked_dir.mkdir(exist_ok=True)

    with open(ligand_csv) as f:
        rows = list(csv.DictReader(f))
    print(f"{len(rows)} ligands")

    print("preparing receptor + zinc pseudo-atoms...")
    receptor_tz = prep_receptor_tz(out_dir)

    print("preparing ligands...")
    ligand_dir = prep_ligands(rows, out_dir)

    print("building AD4Zn affinity maps (one-time)...")
    build_maps(receptor_tz, ligand_dir, out_dir)

    jobs = []
    for r in rows:
        mid = r["molecule_chembl_id"]
        lp = ligand_dir / f"{mid}.pdbqt"
        if lp.exists():
            jobs.append((mid, lp, out_dir, docked_dir))
    print(f"docking {len(jobs)} ligands (ad4 scoring, maps precomputed)...")

    results = {}
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futures = {ex.submit(dock_one, j): j[0] for j in jobs}
        done = 0
        for fut in as_completed(futures):
            mid, score, err = fut.result()
            done += 1
            print(f"{done}/{len(jobs)} {mid}: {score if err is None else 'FAILED ' + err}")
            results[mid] = (score, err)

    out_csv = out_dir / "MMP3_ad4zn_docking_validation_results.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "molecule_chembl_id", "smiles", "label", "standard_value_nM", "vina_score", "error",
        ])
        writer.writeheader()
        for r in rows:
            mid = r["molecule_chembl_id"]
            score, err = results.get(mid, (None, "ligand prep failed"))
            writer.writerow({
                "molecule_chembl_id": mid, "smiles": r["smiles"], "label": r["label"],
                "standard_value_nM": r["standard_value_nM"], "vina_score": score, "error": err or "",
            })
    n_ok = sum(1 for s, e in results.values() if e is None)
    print(f"done: {n_ok}/{len(jobs)} succeeded -> {out_csv}")


if __name__ == "__main__":
    main()
