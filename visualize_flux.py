#!/usr/bin/env python3
"""Boilerplate VPython scene for GS2 mutation score visualization."""

from __future__ import annotations

import csv
from pathlib import Path

from vpython import canvas, color, cylinder, label, rate, sphere, vector

DEFAULT_SCORES = Path("gs2_catalytic_domain_snp_ranking.csv")


def load_scores(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle))


def make_scene() -> None:
    scene = canvas(
        title="Project Synthetix: GS2 variant-energy landscape",
        width=1200,
        height=720,
        background=vector(0.06, 0.08, 0.10),
        center=vector(0, 2, 0),
    )

    stem = cylinder(pos=vector(0, -4, 0), axis=vector(0, 8, 0), radius=0.35, color=color.green)
    node = sphere(pos=vector(0, 1.2, 0), radius=1.2, color=vector(0.38, 0.75, 0.32), opacity=0.7)
    root = cylinder(pos=vector(0, -4, 0), axis=vector(0, -2.5, 0), radius=0.2, color=vector(0.45, 0.30, 0.12))

    label(
        pos=node.pos + vector(0, 1.8, 0),
        text="GS2 catalytic-domain energy peaks",
        box=False,
        height=18,
        color=color.white,
    )

    scores = load_scores(DEFAULT_SCORES)
    top_scores = scores[:12]
    peak_objects = []
    pulse_steps = []

    for index, row in enumerate(top_scores):
        delta = float(row["delta_log_likelihood"])
        peak_height = max(0.25, min(4.0, 0.25 + max(delta, 0.0)))
        x_offset = -3.0 + index * 0.55
        peak = cylinder(
            pos=vector(x_offset, 1.2, 0),
            axis=vector(0, peak_height, 0),
            radius=0.10,
            color=vector(1.0, 0.45, 0.10),
            opacity=0.85,
        )
        peak_objects.append(peak)
        pulse_steps.append((index % 3) + 1)

    if not top_scores:
        label(
            pos=vector(0, -1.5, 0),
            text="Run predict_efficiency.py first to generate ranking data.",
            box=False,
            height=14,
            color=color.white,
        )

    while True:
        rate(30)
        for peak, pulse_step in zip(peak_objects, pulse_steps, strict=True):
            pulse = 0.02 * pulse_step
            peak.radius = 0.10 + pulse


if __name__ == "__main__":
    make_scene()
