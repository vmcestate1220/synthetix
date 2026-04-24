#!/usr/bin/env python3
"""Build structural mutants of GS2 with real side-chain rebuilding.

Uses PDBFixer's ``applyMutations`` to rewrite the residue identity and
rebuild the side chain from a rotamer library, then adds hydrogens at the
chloroplast-stroma pH (8.0) so the output is directly consumable by the
docking pipeline.

Supersedes an earlier version that only rewrote the three-letter residue
label while leaving all side-chain atoms in their wild-type positions.
Any PDB produced by that older version was wild-type geometry in disguise
and should be regenerated here before use.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from pdbfixer import PDBFixer
from openmm.app import PDBFile

# (residue_number, wt_three_letter, mut_three_letter, label)
MUTATIONS: list[tuple[int, str, str, str]] = [
    (112, "ASN", "TYR", "N112Y"),
    (117, "SER", "CYS", "S117C"),
    (170, "GLU", "ALA", "E170A"),
    (125, "ASP", "ASN", "D125N"),
    (111, "TRP", "GLY", "W111G"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--template",
        type=Path,
        default=Path("AF-Q43127-F1.pdb"),
        help="Wild-type template PDB (default: AlphaFold Q43127).",
    )
    parser.add_argument(
        "--chain",
        default="A",
        help="Chain identifier on which to apply mutations (default: A).",
    )
    parser.add_argument(
        "--ph",
        type=float,
        default=8.0,
        help="pH for protonation; 8.0 simulates chloroplast stroma.",
    )
    parser.add_argument(
        "--output-prefix",
        default="gs2_model_",
        help="Output filename prefix; full name is <prefix><label>.pdb.",
    )
    return parser.parse_args()


def build_mutant(
    template: Path,
    chain: str,
    residue_number: int,
    wt_three_letter: str,
    mut_three_letter: str,
    label: str,
    ph: float,
    output_path: Path,
) -> None:
    fixer = PDBFixer(filename=str(template))

    # PDBFixer's mutation spec is "WT-RESNUM-MUT" using three-letter codes.
    spec = f"{wt_three_letter}-{residue_number}-{mut_three_letter}"
    fixer.applyMutations([spec], chain)

    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=ph)

    with output_path.open("w") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle)
    print(f"Wrote {output_path} ({label}: {spec} on chain {chain}, pH {ph})")


def main() -> int:
    args = parse_args()
    if not args.template.exists():
        print(f"Error: template {args.template} not found.")
        return 1

    for residue_number, wt, mut, label in MUTATIONS:
        output_path = Path(f"{args.output_prefix}{label}.pdb")
        build_mutant(
            template=args.template,
            chain=args.chain,
            residue_number=residue_number,
            wt_three_letter=wt,
            mut_three_letter=mut,
            label=label,
            ph=args.ph,
            output_path=output_path,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
