"""Render the candidate-selection funnel (Figure 1 of the manuscript).

Every stage count below is stated in the article's Section 3.2 and traces
back to the screening scripts in this directory:

    604  skin tissue-specific genes          fetch_skin_targets.py
    118  with an Open Targets target class   fetch_skin_targets.py
     26  with confirmed pocket evidence      compare_pockets.py
      3  flagged by both rankings + biology  compare_pockets.py
      2  carried into cross-docking          (MMP3, RARG)

This figure is pure bookkeeping -- it needs no data files, so it can be
regenerated anywhere:

    python render_pipeline_funnel.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

OUT_PNG = Path(__file__).resolve().parent / "figures" / "pipeline_funnel.png"

# (count, label, note) from widest to narrowest.
STAGES = [
    (604, "Skin tissue-specific genes", "Human Protein Atlas expression filter"),
    (118, "Druggable target class", "Open Targets annotation (486 excluded)"),
    (26, "Confirmed pocket evidence", "highest-confidence tier, carried to fpocket"),
    (3, "Flagged candidates", "top of both rankings: MMP3, NTF4, BLMH"),
    (2, "Cross-docked", "MMP3 and RARG"),
]

# Sequential ramp, dark at the top of the funnel to light at the tip; the
# count is also printed on every band, so colour never carries meaning alone.
COLORS = ["#1b4f4a", "#2f7d6e", "#5aa392", "#93c4b4", "#c9e0d6"]

BAND_H = 1.0          # height of each band in data units
GAP = 0.28            # vertical gap between bands
MIN_W = 0.22          # width of the narrowest band, as a fraction of the widest


def band_width(count: int, widest: int) -> float:
    """Width for a stage, log-scaled so 2 of 604 stays visible."""
    import math

    frac = math.log10(count + 1) / math.log10(widest + 1)
    return MIN_W + (1.0 - MIN_W) * frac


def main() -> None:
    widest = STAGES[0][0]
    fig, ax = plt.subplots(figsize=(7.2, 5.4), dpi=300)

    top = 0.0
    for i, (count, label, note) in enumerate(STAGES):
        w = band_width(count, widest)
        w_next = band_width(STAGES[i + 1][0], widest) if i + 1 < len(STAGES) else w
        bottom = top - BAND_H

        # Trapezoid tapering towards the next stage's width.
        ax.add_patch(
            Polygon(
                [(-w / 2, top), (w / 2, top), (w_next / 2, bottom), (-w_next / 2, bottom)],
                closed=True,
                facecolor=COLORS[i],
                edgecolor="white",
                linewidth=1.5,
            )
        )

        # n inside the band; label and note outside, so nothing can collide.
        ax.text(
            0, top - BAND_H / 2, f"n = {count:,}".replace(",", " "),
            ha="center", va="center", fontsize=13, fontweight="bold",
            color="white" if i < 3 else "#17352f",
        )
        ax.text(
            w / 2 + 0.04, top - BAND_H / 2 + 0.13, label,
            ha="left", va="center", fontsize=11, color="#17352f",
        )
        ax.text(
            w / 2 + 0.04, top - BAND_H / 2 - 0.17, note,
            ha="left", va="center", fontsize=8.5, color="#5c6b68",
        )

        top = bottom - GAP

    ax.set_xlim(-0.62, 1.5)
    ax.set_ylim(top + GAP - 0.25, 0.3)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(OUT_PNG, bbox_inches="tight", facecolor="white")
    print(f"wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
