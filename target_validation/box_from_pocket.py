"""Compute a Vina docking box (center + size) from an fpocket pocket*_vert.pqr
file (the alpha-sphere/voronoi-vertex cloud fpocket detected for that pocket).
"""
import sys
from pathlib import Path

PADDING = 6.0  # angstrom, added to the raw pocket extent on each side


def parse_pqr_coords(path: Path):
    coords = []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        parts = line.split()
        x, y, z = float(parts[-5]), float(parts[-4]), float(parts[-3])
        coords.append((x, y, z))
    return coords


def box_from_coords(coords):
    xs, ys, zs = zip(*coords)
    cx, cy, cz = sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)
    sx = (max(xs) - min(xs)) + PADDING
    sy = (max(ys) - min(ys)) + PADDING
    sz = (max(zs) - min(zs)) + PADDING
    return cx, cy, cz, sx, sy, sz


if __name__ == "__main__":
    pqr_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])
    coords = parse_pqr_coords(pqr_path)
    cx, cy, cz, sx, sy, sz = box_from_coords(coords)
    text = (
        f"center_x = {cx:.3f}\ncenter_y = {cy:.3f}\ncenter_z = {cz:.3f}\n"
        f"size_x = {sx:.3f}\nsize_y = {sy:.3f}\nsize_z = {sz:.3f}\n"
    )
    out_path.write_text(text)
    print(f"{pqr_path.parent.parent.name}: center=({cx:.2f},{cy:.2f},{cz:.2f}) size=({sx:.1f},{sy:.1f},{sz:.1f})")
