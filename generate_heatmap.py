#!/usr/bin/env python3
"""Generate a GS2 catalytic-domain variant heatmap from the ranking CSV."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("/tmp/matplotlib-cache")))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

DEFAULT_INPUT = Path("gs2_catalytic_domain_snp_ranking.csv")
DEFAULT_OUTPUT = Path("gs2_efficiency_heatmap.png")
AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def ordered_amino_acids(frame: pd.DataFrame) -> list[str]:
    observed = set(frame["mut_aa"]).union(frame["wt_aa"])
    ordered = [aa for aa in AA_ORDER if aa in observed]
    extras = sorted(observed.difference(AA_ORDER))
    return ordered + extras


def build_heatmap_data(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    amino_acids = ordered_amino_acids(frame)
    positions = sorted(frame["aa_pos_1based"].unique())

    heatmap = frame.pivot_table(
        index="mut_aa",
        columns="aa_pos_1based",
        values="delta_log_likelihood",
        aggfunc="mean",
    )
    heatmap = heatmap.reindex(index=amino_acids, columns=positions)

    wildtype = (
        frame[["aa_pos_1based", "wt_aa"]]
        .drop_duplicates()
        .sort_values("aa_pos_1based")
        .set_index("aa_pos_1based")
    )
    return heatmap, wildtype


def draw_heatmap(heatmap: pd.DataFrame, wildtype: pd.DataFrame, output_path: Path) -> None:
    sns.set_theme(style="white")
    fig_width = max(12, len(heatmap.columns) * 0.18)
    fig, ax = plt.subplots(figsize=(fig_width, 6.5))

    sns.heatmap(
        heatmap,
        ax=ax,
        cmap="RdBu_r",
        center=0,
        linewidths=0.25,
        linecolor="#f2f2f2",
        cbar_kws={"label": "Delta log-likelihood"},
    )

    x_positions = list(range(len(heatmap.columns)))
    tick_labels = [str(pos) if pos % 10 == 0 else "" for pos in heatmap.columns]
    ax.set_xticks([x + 0.5 for x in x_positions])
    ax.set_xticklabels(tick_labels, rotation=0)
    ax.set_xlabel("GS2 Residue Position")
    ax.set_ylabel("Mutant Amino Acid")
    ax.set_title("GS2 Catalytic-Domain Variant Effect Heatmap")

    aa_to_row = {aa: idx for idx, aa in enumerate(heatmap.index)}
    marker_x = []
    marker_y = []
    for col_index, position in enumerate(heatmap.columns):
        wt_aa = wildtype.loc[position, "wt_aa"]
        if wt_aa in aa_to_row:
            marker_x.append(col_index + 0.5)
            marker_y.append(aa_to_row[wt_aa] + 0.5)

    ax.scatter(
        marker_x,
        marker_y,
        s=18,
        facecolors="none",
        edgecolors="black",
        linewidths=0.6,
        zorder=3,
    )

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    frame = pd.read_csv(args.input)
    heatmap, wildtype = build_heatmap_data(frame)
    draw_heatmap(heatmap, wildtype, args.output)
    print(f"Saved heatmap to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
