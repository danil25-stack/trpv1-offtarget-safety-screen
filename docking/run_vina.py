"""Run AutoDock Vina for one ligand against a prepared receptor + box."""
import subprocess
import sys
from pathlib import Path


def read_box(box_txt: Path) -> dict:
    """Parse the Vina-style box file meeko writes (receptor.box.txt)."""
    values = {}
    for line in box_txt.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        key, val = line.split("=")
        values[key.strip()] = float(val.strip())
    return values


def dock(receptor_pdbqt: Path, ligand_pdbqt: Path, box: dict, out_pdbqt: Path,
         exhaustiveness: int = 8, num_modes: int = 9, cpu: int | None = None) -> str:
    cmd = [
        "vina",
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", str(box["center_x"]), "--center_y", str(box["center_y"]), "--center_z", str(box["center_z"]),
        "--size_x", str(box["size_x"]), "--size_y", str(box["size_y"]), "--size_z", str(box["size_z"]),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--out", str(out_pdbqt),
    ]
    if cpu is not None:
        # Vina defaults to using every CPU it *detects* (nproc), which
        # ignores a cgroup CFS quota -- on a fractional-CPU cloud instance
        # (e.g. a GPU rental with a 1-2 core quota but nproc reporting the
        # full host count) that means every concurrent docking job spawns
        # far more threads than the container can actually run at once,
        # thrashing instead of parallelizing. Pin explicitly to the real quota.
        cmd += ["--cpu", str(cpu)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"vina failed:\n{result.stdout}\n{result.stderr}")
    return result.stdout


if __name__ == "__main__":
    receptor, ligand, box_txt, out_pdbqt = (Path(a) for a in sys.argv[1:5])
    box = read_box(box_txt)
    log = dock(receptor, ligand, box, out_pdbqt)
    print(log)
