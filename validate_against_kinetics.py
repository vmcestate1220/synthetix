#!/usr/bin/env python3
"""Anchor Synthetix Δ-LL scores against experimentally-measured mutant effects.

The pipeline's ``delta_log_likelihood`` is an evolutionary-plausibility score
from Evo 2. Before trusting the top of the ranking, we need evidence that it
correlates with actually-measured functional quantities on GS2 mutants.

This script joins:

  * the Synthetix ranking CSV (produced by ``predict_efficiency.py``)
  * a user-curated CSV of known mutants with measured effects

and reports Pearson + Spearman correlations between measured effect and
Δ-LL. The correlation floor below which the ranking should not be trusted
is a judgment call, but r < ~0.3 on a dataset of 10+ mutants is a strong
signal that the scoring proxy is not measuring catalytic fitness.

See ``known_mutants_template.csv`` for the expected input schema.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--ranking",
        type=Path,
        default=Path("gs2_catalytic_domain_snp_ranking.csv"),
        help="Synthetix ranking CSV with delta_log_likelihood column.",
    )
    parser.add_argument(
        "--known",
        type=Path,
        required=True,
        help=(
            "User-curated CSV of known mutants. Required columns: "
            "mutation (e.g. D125N), measured_effect (float; higher = more "
            "favorable by convention). Optional: effect_type, citation, notes."
        ),
    )
    parser.add_argument(
        "--plot",
        type=Path,
        default=None,
        help="Optional scatter-plot output path (PNG).",
    )
    parser.add_argument(
        "--min-observations",
        type=int,
        default=5,
        help="Refuse to report correlation below this many matched rows.",
    )
    return parser.parse_args()


def split_mutation(label: str) -> tuple[str, int, str]:
    label = label.strip().upper()
    if len(label) < 3:
        raise ValueError(f"Mutation label too short: {label!r}")
    wt_aa = label[0]
    mut_aa = label[-1]
    position = int(label[1:-1])
    return wt_aa, position, mut_aa


def aggregate_ranking_by_mutation(ranking: pd.DataFrame) -> pd.DataFrame:
    """Collapse multiple SNPs that produce the same AA change to one row.

    A missense mutation can arise from up to three distinct SNPs; we take
    the best (maximum) Δ-LL across them, matching the convention that
    ``predict_efficiency.py`` sorts by (rank 1 = max delta).
    """
    missense = ranking[ranking["effect"] == "missense"].copy()
    missense["mutation"] = (
        missense["wt_aa"] + missense["aa_pos_1based"].astype(str) + missense["mut_aa"]
    )
    grouped = (
        missense.groupby("mutation", as_index=False)
        .agg(delta_log_likelihood=("delta_log_likelihood", "max"))
    )
    return grouped


def pearson(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    cov = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    var_x = sum((x - mean_x) ** 2 for x in xs)
    var_y = sum((y - mean_y) ** 2 for y in ys)
    denom = math.sqrt(var_x * var_y)
    return cov / denom if denom > 0 else float("nan")


def spearman(xs: list[float], ys: list[float]) -> float:
    def ranks(values: list[float]) -> list[float]:
        indexed = sorted(range(len(values)), key=lambda i: values[i])
        out = [0.0] * len(values)
        i = 0
        while i < len(values):
            j = i
            while j + 1 < len(values) and values[indexed[j + 1]] == values[indexed[i]]:
                j += 1
            avg_rank = (i + j) / 2 + 1
            for k in range(i, j + 1):
                out[indexed[k]] = avg_rank
            i = j + 1
        return out

    return pearson(ranks(xs), ranks(ys))


def main() -> int:
    args = parse_args()
    ranking = pd.read_csv(args.ranking)
    known = pd.read_csv(args.known)

    required = {"mutation", "measured_effect"}
    missing = required.difference(known.columns)
    if missing:
        print(f"Error: --known CSV is missing required columns: {sorted(missing)}")
        return 2

    # Validate mutation labels early so bad rows surface before the join.
    for label in known["mutation"]:
        split_mutation(str(label))

    aggregated = aggregate_ranking_by_mutation(ranking)
    merged = known.merge(aggregated, on="mutation", how="inner")

    unmatched = set(known["mutation"]) - set(merged["mutation"])
    if unmatched:
        print(
            f"Warning: {len(unmatched)} mutation(s) in the known set have no "
            f"Synthetix score and were dropped: {sorted(unmatched)}"
        )

    if len(merged) < args.min_observations:
        print(
            f"Not enough matched observations: {len(merged)} < "
            f"{args.min_observations}. Refusing to report a correlation that "
            "would be noise-dominated. Populate more rows in --known and retry."
        )
        return 3

    xs = merged["measured_effect"].astype(float).tolist()
    ys = merged["delta_log_likelihood"].astype(float).tolist()

    r_pearson = pearson(xs, ys)
    r_spearman = spearman(xs, ys)

    print(f"Matched {len(merged)} mutations.")
    print(f"Pearson  r = {r_pearson:+.3f}")
    print(f"Spearman ρ = {r_spearman:+.3f}")
    print()
    print(merged[["mutation", "measured_effect", "delta_log_likelihood"]].to_string(index=False))

    if args.plot is not None:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.scatter(xs, ys, s=40)
        for _, row in merged.iterrows():
            ax.annotate(row["mutation"], (row["measured_effect"], row["delta_log_likelihood"]), fontsize=8)
        ax.set_xlabel("Measured effect (from --known)")
        ax.set_ylabel("Synthetix delta_log_likelihood")
        ax.set_title(
            f"Evo 2 Δ-LL vs measured (Pearson r={r_pearson:+.2f}, Spearman ρ={r_spearman:+.2f})"
        )
        ax.axhline(0, color="#888", linewidth=0.5)
        ax.axvline(0, color="#888", linewidth=0.5)
        fig.tight_layout()
        fig.savefig(args.plot, dpi=200)
        print(f"Scatter saved to {args.plot}")

    # Exit non-zero if correlation is too weak to act on — makes it usable
    # as a CI gate between scoring and downstream structural triage.
    if r_spearman < 0.3:
        print(
            "\nSpearman ρ < 0.3: the ranking is not reliably aligned with "
            "measured function on this benchmark. Do not triage structural "
            "audits off the top of the Synthetix CSV without orthogonal "
            "evidence (stability models, MD, or wet-lab).",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
