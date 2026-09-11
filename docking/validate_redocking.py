"""Redocking validation: dock a ligand back into the pocket it was
crystallized in, then RMSD the best pose against the crystal coordinates.
< ~2 A on heavy atoms is the usual bar for "the docking setup reproduces
reality" before trusting it on unknown compounds.
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS


def to_sdf_via_obabel(in_path: Path, out_sdf: Path) -> Path:
    # obabel does distance-based bond perception, which sidesteps RDKit's
    # strict fixed-column PDB parser -- that parser breaks on 5-character
    # extended CCD residue codes (e.g. A1C8K), which overflow the classic
    # 3-character resName column.
    result = subprocess.run(["obabel", str(in_path), "-O", str(out_sdf)], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"obabel failed:\n{result.stdout}\n{result.stderr}")
    return out_sdf


def best_pose_sdf(vina_out_pdbqt: Path, out_sdf: Path) -> Path:
    text = vina_out_pdbqt.read_text()
    first_model = text.split("ENDMDL")[0] + "ENDMDL\n"
    tmp_pdbqt = out_sdf.with_suffix(".best.pdbqt")
    tmp_pdbqt.write_text(first_model)
    return to_sdf_via_obabel(tmp_pdbqt, out_sdf)


def rmsd_to_reference(smiles: str, reference_pdb: Path, docked_pdbqt: Path, work_dir: Path) -> float:
    template = Chem.MolFromSmiles(smiles)

    # Reference (crystal) pose: obabel's bond guess is clean here (no extra
    # explicit polar Hs), so template-based bond-order reassignment works
    # directly and gives a mol with unambiguous, correctly-ordered bonds.
    ref_sdf = to_sdf_via_obabel(reference_pdb, work_dir / "reference_ligand.sdf")
    ref_raw = Chem.MolFromMolFile(str(ref_sdf), sanitize=False, removeHs=True)
    ref_mol = AllChem.AssignBondOrdersFromTemplate(template, ref_raw)
    Chem.SanitizeMol(ref_mol)

    # Docked pose: PDBQT keeps polar hydrogens explicit and its alkyl tail's
    # local symmetry (terminal isopropyl) makes AssignBondOrdersFromTemplate's
    # matcher pick an inconsistent mapping (atoms come out as radicals after
    # sanitization). Bond *orders* don't matter for an atom correspondence,
    # so match on topology alone via FMCS (bondCompare=CompareAny) instead --
    # robust to the same ambiguity that trips up strict bond-order matching.
    docked_sdf = best_pose_sdf(docked_pdbqt, work_dir / "best_pose.sdf")
    docked_raw = Chem.MolFromMolFile(str(docked_sdf), sanitize=False, removeHs=True)
    Chem.SanitizeMol(
        docked_raw,
        sanitizeOps=Chem.SANITIZE_FINDRADICALS | Chem.SANITIZE_SETAROMATICITY
        | Chem.SANITIZE_ADJUSTHS | Chem.SANITIZE_SYMMRINGS,
    )

    mcs = rdFMCS.FindMCS(
        [ref_mol, docked_raw],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        ringMatchesRingOnly=False, completeRingsOnly=False, timeout=30,
    )
    if mcs.numAtoms != ref_mol.GetNumAtoms():
        raise RuntimeError(
            f"MCS only covers {mcs.numAtoms}/{ref_mol.GetNumAtoms()} atoms -- "
            "docked pose topology doesn't match the reference ligand."
        )
    patt = Chem.MolFromSmarts(mcs.smartsString)
    ref_match = ref_mol.GetSubstructMatch(patt)
    docked_match = docked_raw.GetSubstructMatch(patt)

    ref_conf = ref_mol.GetConformer()
    docked_conf = docked_raw.GetConformer()
    sq_dists = [
        np.sum((np.array(ref_conf.GetAtomPosition(a)) - np.array(docked_conf.GetAtomPosition(b))) ** 2)
        for a, b in zip(ref_match, docked_match)
    ]
    return float(np.sqrt(np.mean(sq_dists)))


if __name__ == "__main__":
    smiles, ref_pdb, vina_out, work_dir = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])
    rmsd = rmsd_to_reference(smiles, ref_pdb, vina_out, work_dir)
    print(f"Redocking RMSD (best pose vs crystal, heavy atoms): {rmsd:.2f} A")
