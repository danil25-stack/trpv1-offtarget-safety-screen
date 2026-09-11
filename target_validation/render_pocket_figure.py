"""Render TRPV1's vanilloid pocket and RARG's retinoid pocket side by side
from their fpocket detections, for the article figure requested to
illustrate "pocket similarity" visually.

Honesty note: TRPV1 (ion channel) and RARG (nuclear receptor) share no
fold homology, so there is no valid 3D frame to structurally *superimpose*
them in -- the article's pocket-similarity claim is a comparison of an
8-descriptor fpocket fingerprint (volume, druggability, hydrophobicity,
etc.), not a spatial/geometric alignment (see the "purpose-built
pocket-comparison methods... would provide a more rigorous... comparison"
limitation already stated in the article). This renders each detected
pocket independently, side by side at matched relative scale, rather than
overlaid -- an overlay would visually claim an alignment that was never
computed or validated.

Pocket shape is approximated as the convex hull of fpocket's alpha-sphere
centers (pocketN_vert.pqr) -- a simplification of fpocket's own alpha-shape
(it will not show concavities), disclosed as such. Lining residues
(pocketN_atm.pdb) are rendered as backbone-atom scatter.

Inputs (already computed by the group's existing fpocket pipeline,
compare_pockets.py / skin_offtarget_pocket_comparison.md):
    data/processed/skin_offtarget_screen/fpocket_out_trpv1_reference/TRPV1_AD_out/
    data/processed/skin_offtarget_screen/fpocket_out/RARG_out/

Usage:
    python render_pocket_figure.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull

DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "processed" / "skin_offtarget_screen"
OUT_PNG = Path(__file__).resolve().parent / "figures" / "pocket_shapes.png"

TARGETS = {
    "TRPV1 (vanilloid pocket)": {
        "pocket_dir": DATA_DIR / "fpocket_out_trpv1_reference" / "TRPV1_AD_out" / "pockets",
        "pocket_num": 1,
        "druggability": 0.783,
        "volume": 758,
        "highlight_resnames_seqids": {("TYR", 511), ("SER", 512), ("THR", 550), ("GLU", 570)},
        "color": "#0f6b5e",
    },
    "RARG (retinoid pocket)": {
        "pocket_dir": DATA_DIR / "fpocket_out" / "RARG_out" / "pockets",
        "pocket_num": 1,
        "druggability": 0.982,
        "volume": 1032,
        "highlight_resnames_seqids": set(),
        "color": "#9a4a2e",
    },
}


def parse_pqr_coords(path: Path):
    coords = []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        parts = line.split()
        x, y, z = float(parts[-5]), float(parts[-4]), float(parts[-3])
        coords.append((x, y, z))
    return np.array(coords)


def parse_atm_pdb(path: Path):
    """Return (coords, resname, resid) for each CA/backbone-representative
    atom of residues lining the pocket."""
    atoms = []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        resname = line[17:20].strip()
        try:
            resid = int(line[22:26])
        except ValueError:
            continue
        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        atoms.append((np.array([x, y, z]), resname, resid))
    return atoms


def plot_pocket(ax, name, cfg):
    vert_path = cfg["pocket_dir"] / f"pocket{cfg['pocket_num']}_vert.pqr"
    atm_path = cfg["pocket_dir"] / f"pocket{cfg['pocket_num']}_atm.pdb"
    alpha_spheres = parse_pqr_coords(vert_path)
    lining_atoms = parse_atm_pdb(atm_path)

    centroid = alpha_spheres.mean(axis=0)
    alpha_c = alpha_spheres - centroid

    hull = ConvexHull(alpha_c)
    poly = Poly3DCollection(alpha_c[hull.simplices], alpha=0.28, facecolor=cfg["color"], edgecolor=cfg["color"], linewidths=0.2)
    ax.add_collection3d(poly)
    ax.scatter(*alpha_c.T, s=4, color=cfg["color"], alpha=0.5, depthshade=True, label="alpha spheres (fpocket)")

    lining_xyz = np.array([a[0] for a in lining_atoms]) - centroid
    highlighted = cfg["highlight_resnames_seqids"]
    hi_mask = np.array([(rn, rid) in highlighted for _, rn, rid in lining_atoms])
    if hi_mask.any():
        ax.scatter(*lining_xyz[hi_mask].T, s=60, color="#d4a017", edgecolor="black", linewidths=0.5,
                   depthshade=False, label="literature-validated residues")
        for (xyz, rn, rid), hi in zip(lining_atoms, hi_mask):
            if hi:
                ax.text(*(xyz - centroid), f"{rn}{rid}", fontsize=7, color="#7a5c00")
    ax.scatter(*lining_xyz[~hi_mask].T, s=14, color="#444444", alpha=0.6, depthshade=True, label="pocket-lining Cα")

    span = np.abs(alpha_c).max() * 1.15
    ax.set_xlim(-span, span); ax.set_ylim(-span, span); ax.set_zlim(-span, span)
    ax.set_box_aspect([1, 1, 1])
    ax.set_title(f"{name}\ndruggability {cfg['druggability']:.3f}  ·  volume {cfg['volume']} Å³  ·  "
                 f"{len(alpha_spheres)} alpha spheres", fontsize=10)
    ax.set_axis_off()


def main():
    fig = plt.figure(figsize=(11, 5.5), dpi=200)
    for i, (name, cfg) in enumerate(TARGETS.items()):
        ax = fig.add_subplot(1, 2, i + 1, projection="3d")
        plot_pocket(ax, name, cfg)
        if i == 0:
            ax.legend(loc="upper left", fontsize=7, framealpha=0.9)

    fig.suptitle(
        "Detected pocket shape (fpocket alpha-sphere convex hull) — shown\n"
        "independently, not superimposed: TRPV1 and RARG share no fold homology,\n"
        "so no structural alignment between them was computed or is implied here",
        fontsize=9, color="#565f68"
    )
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    OUT_PNG.parent.mkdir(exist_ok=True)
    fig.savefig(OUT_PNG)
    print(f"wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
