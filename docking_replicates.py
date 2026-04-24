#!/usr/bin/env python3
"""Run AutoDock Vina multiple times with distinct seeds and summarize results.

AutoDock Vina is stochastic — a single run gives you a point estimate of the
best-mode binding affinity, with no sense of its variance. This wrapper
runs Vina N times with distinct ``--seed`` values, parses each output, and
reports the mean and standard deviation of the top-mode affinity across
replicates, plus the RMSD of each replicate's best pose against the first.

A wide distribution (std > ~0.5 kcal/mol across 5 runs, or pose RMSD >
~2.0 Å) is a signal that ``exhaustiveness`` is too low for the pocket, or
that the search box is under-constraining the ligand.
"""

from __future__ import annotations

import argparse
import math
import re
import statistics
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--vina",
        default="./vina",
        help="Path to the Vina binary (default: ./vina from repo root).",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("docking_config.txt"),
        help="Vina config file (receptor/ligand/box).",
    )
    parser.add_argument(
        "--replicates",
        type=int,
        default=5,
        help="Number of replicate runs (default: 5).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("vina_replicates"),
        help="Directory for per-replicate outputs.",
    )
    parser.add_argument(
        "--exhaustiveness",
        type=int,
        default=None,
        help="Override exhaustiveness from the config file.",
    )
    return parser.parse_args()


MODE_LINE = re.compile(
    r"^\s*1\s+(-?\d+\.\d+)"  # mode 1's first numeric column is the affinity
)


def parse_best_affinity(output_pdbqt: Path) -> float:
    """Read the REMARK VINA RESULT from the top docked mode."""
    with output_pdbqt.open() as handle:
        for line in handle:
            if line.startswith("REMARK VINA RESULT:"):
                # REMARK VINA RESULT:   -7.839  0.000  0.000
                parts = line.split()
                # parts: ['REMARK', 'VINA', 'RESULT:', affinity, rmsd_lb, rmsd_ub]
                return float(parts[3])
    raise ValueError(f"No REMARK VINA RESULT found in {output_pdbqt}")


def parse_first_model_coords(output_pdbqt: Path) -> list[tuple[float, float, float]]:
    """Extract heavy-atom XYZ for the first docked MODEL as a baseline pose."""
    coords: list[tuple[float, float, float]] = []
    in_first_model = False
    with output_pdbqt.open() as handle:
        for line in handle:
            if line.startswith("MODEL"):
                if in_first_model:
                    break
                in_first_model = True
                continue
            if line.startswith("ENDMDL"):
                if in_first_model:
                    break
            if in_first_model and (line.startswith("ATOM") or line.startswith("HETATM")):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append((x, y, z))
    return coords


def pose_rmsd(a: list[tuple[float, float, float]], b: list[tuple[float, float, float]]) -> float:
    if len(a) != len(b) or not a:
        return float("nan")
    acc = 0.0
    for (ax, ay, az), (bx, by, bz) in zip(a, b):
        acc += (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2
    return math.sqrt(acc / len(a))


def run_vina(
    vina: str,
    config: Path,
    out_path: Path,
    seed: int,
    exhaustiveness: int | None,
) -> None:
    cmd = [vina, "--config", str(config), "--out", str(out_path), "--seed", str(seed)]
    if exhaustiveness is not None:
        cmd += ["--exhaustiveness", str(exhaustiveness)]
    subprocess.run(cmd, check=True, capture_output=True)


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    affinities: list[float] = []
    first_pose: list[tuple[float, float, float]] | None = None
    rmsds: list[float] = []

    for seed in range(1, args.replicates + 1):
        out_path = args.out_dir / f"out_seed{seed}.pdbqt"
        run_vina(args.vina, args.config, out_path, seed, args.exhaustiveness)
        affinity = parse_best_affinity(out_path)
        affinities.append(affinity)
        pose = parse_first_model_coords(out_path)
        if first_pose is None:
            first_pose = pose
            rmsds.append(0.0)
        else:
            rmsds.append(pose_rmsd(first_pose, pose))
        print(f"seed={seed}  affinity={affinity:+.3f} kcal/mol  rmsd_vs_seed1={rmsds[-1]:.2f} Å")

    mean = statistics.mean(affinities)
    std = statistics.stdev(affinities) if len(affinities) > 1 else 0.0
    print()
    print(f"Affinity over {args.replicates} replicates: {mean:+.3f} ± {std:.3f} kcal/mol")
    print(f"Best-pose RMSD vs seed=1: mean {statistics.mean(rmsds):.2f} Å, max {max(rmsds):.2f} Å")
    if std > 0.5:
        print(
            "WARNING: affinity std > 0.5 kcal/mol — the pocket is under-sampled. "
            "Consider raising --exhaustiveness or narrowing the search box."
        )
    if max(rmsds) > 2.0:
        print(
            "WARNING: best-pose RMSD > 2.0 Å between replicates — Vina is "
            "finding structurally distinct poses. Review the search box."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
