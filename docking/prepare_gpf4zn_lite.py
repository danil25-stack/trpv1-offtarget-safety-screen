"""Write an AutoGrid4 .gpf file for AD4Zn (tetrahedral-zinc-pseudoatom)
scoring, without the legacy MGLTools/AutoDockTools dependency that the
official prepare_gpf4zn.py script requires (it imports MolKit and
AutoDockTools.GridParameters, which aren't packaged for modern Python/conda
and normally require installing the full ADFRsuite).

The .gpf format itself is plain text with no computation beyond: box
geometry, the set of atom types actually present in the receptor and ligand
pdbqt files, and a fixed block of AD4Zn pairwise parameters that are
force-field constants, not receptor-specific -- copied verbatim from the
worked example in the AutoDock-Vina repo (ccsb-scripps/AutoDock-Vina,
example/docking_with_zinc_metalloproteins/solution/protein_tz.gpf), which is
the authoritative source for this force field
(https://autodock-vina.readthedocs.io/en/latest/docking_zinc.html).

Usage:
    python prepare_gpf4zn_lite.py --receptor receptor_tz.pdbqt \\
        --ligand-types-from ligands_prepped/*.pdbqt \\
        --center CX CY CZ --size SX SY SZ \\
        --name protein_tz --out protein_tz.gpf
"""
import argparse
import glob
from pathlib import Path

SPACING = 0.375  # AD4 default grid spacing (A)

# AD4Zn's tetrahedral-zinc pairwise parameters -- force-field constants from
# Santos-Martins et al. 2014 (10.1021/ci500209e), not receptor- or
# ligand-specific. Copied verbatim from the official worked example.
NBP_R_EPS_BLOCK = """nbp_r_eps 0.25 23.2135 12 6 NA TZ
nbp_r_eps 2.1   3.8453 12 6 OA Zn
nbp_r_eps 2.25  7.5914 12 6 SA Zn
nbp_r_eps 1.0   0.0    12 6 HD Zn
nbp_r_eps 2.0   0.0060 12 6 NA Zn
nbp_r_eps 2.0   0.2966 12 6  N Zn"""


def read_pdbqt_atom_types(path):
    types = set()
    for line in Path(path).read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            types.add(line.split()[-1])
    return types


def npts_for(size):
    # AutoGrid requires an even number of grid points per axis.
    n = int(size / SPACING) + 1
    return n + 1 if n % 2 else n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", required=True, help="receptor_tz.pdbqt (post zinc_pseudo.py)")
    ap.add_argument("--ligand-glob", required=True, help="glob pattern matching all ligand pdbqt files to be docked against these maps")
    ap.add_argument("--center", nargs=3, type=float, required=True)
    ap.add_argument("--size", nargs=3, type=float, required=True)
    ap.add_argument("--name", required=True, help="basename for map files (e.g. protein_tz)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    receptor_types = sorted(read_pdbqt_atom_types(args.receptor))
    ligand_files = sorted(glob.glob(args.ligand_glob))
    if not ligand_files:
        raise RuntimeError(f"no ligand files matched {args.ligand_glob!r}")
    ligand_types = set()
    for lf in ligand_files:
        ligand_types |= read_pdbqt_atom_types(lf)
    ligand_types = sorted(ligand_types)
    print(f"receptor_types: {receptor_types}")
    print(f"ligand_types ({len(ligand_files)} ligands scanned): {ligand_types}")

    npts = [npts_for(s) for s in args.size]
    cx, cy, cz = args.center
    name = args.name

    lines = [
        f"npts {npts[0]} {npts[1]} {npts[2]}                        # num.grid points in xyz",
        "parameter_file AD4Zn.dat             # force field default parameter file",
        f"gridfld {name}.maps.fld          # grid_data_file",
        f"spacing {SPACING}                        # spacing(A)",
        f"receptor_types {' '.join(receptor_types)} # receptor atom types",
        f"ligand_types {' '.join(ligand_types)}       # ligand atom types",
        f"receptor {args.receptor}            # macromolecule",
        f"gridcenter {cx} {cy} {cz}                 # xyz-coordinates or auto",
        "smooth 0.5                           # store minimum energy w/in rad(A)",
    ]
    for lt in ligand_types:
        lines.append(f"map {name}.{lt}.map                 # atom-specific affinity map")
    lines.append(f"elecmap {name}.e.map             # electrostatic potential map")
    lines.append(f"dsolvmap {name}.d.map              # desolvation potential map")
    lines.append("dielectric -0.1465                   # <0, AD4 distance-dep.diel;>0, constant")
    lines.append(NBP_R_EPS_BLOCK)

    Path(args.out).write_text("\n".join(lines) + "\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
