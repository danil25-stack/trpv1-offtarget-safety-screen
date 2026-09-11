"""Build a 3D ligand structure from SMILES and prepare it as PDBQT for Vina.

Docking should start from an independently-generated conformer, not the
crystal pose, so redocking against a known structure is a meaningful
validation (compare the docked pose back to the crystal pose via RMSD).
"""
import subprocess
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


def smiles_to_sdf(smiles: str, name: str, out_path: Path, seed: int = 42) -> Path:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    mol = Chem.AddHs(mol)
    mol.SetProp("_Name", name)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError(f"3D embedding failed for {name} ({smiles})")
    AllChem.MMFFOptimizeMolecule(mol)
    writer = Chem.SDWriter(str(out_path))
    writer.write(mol)
    writer.close()
    return out_path


def sdf_to_pdbqt(sdf_path: Path, out_path: Path) -> Path:
    result = subprocess.run(
        ["mk_prepare_ligand.py", "-i", str(sdf_path), "-o", str(out_path)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"mk_prepare_ligand.py failed:\n{result.stdout}\n{result.stderr}")
    return out_path


if __name__ == "__main__":
    smiles, name, out_dir = sys.argv[1], sys.argv[2], Path(sys.argv[3])
    out_dir.mkdir(parents=True, exist_ok=True)
    sdf = smiles_to_sdf(smiles, name, out_dir / f"{name}.sdf")
    pdbqt = sdf_to_pdbqt(sdf, out_dir / f"{name}.pdbqt")
    print(f"Wrote {pdbqt}")
