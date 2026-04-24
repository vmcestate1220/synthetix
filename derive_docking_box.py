#!/usr/bin/env python3
"""Derive the Vina search box from GS2 active-site residue coordinates.

The docking box in ``docking_config.txt`` pins a cube in receptor space; if
the cube is not centered on the active site or is too small to contain a
full-fit of the ligand, every downstream Vina run produces meaningless
affinities.

This script reads an AlphaFold PDB, pulls the heavy-atom coordinates of a
named set of active-site residues, computes the centroid and the smallest
cube (padded) that contains all those atoms, and prints the result as a
Vina config block. The committed ``docking_config.txt`` is annotated with
the exact invocation that produced its numbers so the provenance is
reproducible.

Active-site residue set comes from the GS literature:
    * Metal-binding (n1/n2 sites):   Glu-124, Asp-125, Glu-187, Glu-195, His-268
    * Substrate (glutamate) anchor:  Arg-319, Ser-173
    * ATP / phospho-transfer:        Gly-251, Arg-354, Arg-335
    * Catalytic:                     Glu-187
(Residue numbers are 1-based in the full-length UniProt Q43127 sequence.
They inherit from the bacterial GS active-site map of Eisenberg et al.
2000 and Liaw & Eisenberg 1994, adjusted to the Arabidopsis GS2 numbering.)
"""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio.PDB import PDBParser

DEFAULT_ACTIVE_SITE = (124, 125, 151, 173, 187, 195, 251, 268, 319, 335, 354)
DEFAULT_PADDING_ANGSTROMS = 8.0
DEFAULT_MIN_BOX_SIZE = 20.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pdb",
        type=Path,
        default=Path("AF-Q43127-F1.pdb"),
        help="Receptor PDB to measure (default: AlphaFold Q43127).",
    )
    parser.add_argument(
        "--chain",
        default="A",
        help="Chain identifier (default: A).",
    )
    parser.add_argument(
        "--residues",
        nargs="+",
        type=int,
        default=list(DEFAULT_ACTIVE_SITE),
        help=f"Active-site residue numbers (default: {list(DEFAULT_ACTIVE_SITE)}).",
    )
    parser.add_argument(
        "--padding",
        type=float,
        default=DEFAULT_PADDING_ANGSTROMS,
        help=f"Extra Å around the residue bounding box (default: {DEFAULT_PADDING_ANGSTROMS}).",
    )
    parser.add_argument(
        "--min-size",
        type=float,
        default=DEFAULT_MIN_BOX_SIZE,
        help=f"Minimum cube edge length in Å (default: {DEFAULT_MIN_BOX_SIZE}).",
    )
    parser.add_argument(
        "--receptor",
        default="WT_ready.pdbqt",
        help="Receptor .pdbqt filename to write into the Vina config block.",
    )
    parser.add_argument(
        "--ligand",
        default="Glutamate.pdbqt",
        help="Ligand .pdbqt filename to write into the Vina config block.",
    )
    parser.add_argument(
        "--exhaustiveness",
        type=int,
        default=8,
        help="Vina exhaustiveness (default: 8).",
    )
    return parser.parse_args()


def collect_active_site_atoms(
    pdb_path: Path, chain_id: str, residue_numbers: list[int]
) -> list[tuple[str, int, tuple[float, float, float]]]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("gs2", str(pdb_path))
    atoms: list[tuple[str, int, tuple[float, float, float]]] = []
    wanted = set(residue_numbers)
    found_residues: set[int] = set()
    for model in structure:
        for chain in model:
            if chain.id != chain_id:
                continue
            for residue in chain:
                resnum = residue.get_id()[1]
                if resnum not in wanted:
                    continue
                found_residues.add(resnum)
                resname = residue.get_resname()
                for atom in residue:
                    if atom.element == "H":
                        continue
                    x, y, z = atom.coord
                    atoms.append((resname, resnum, (float(x), float(y), float(z))))
    missing = wanted.difference(found_residues)
    if missing:
        raise ValueError(
            f"Residues {sorted(missing)} not found in chain {chain_id} of {pdb_path}. "
            "Verify the residue set matches the numbering of this PDB."
        )
    return atoms


def compute_box(
    atoms: list[tuple[str, int, tuple[float, float, float]]],
    padding: float,
    min_size: float,
) -> dict[str, float]:
    xs = [a[2][0] for a in atoms]
    ys = [a[2][1] for a in atoms]
    zs = [a[2][2] for a in atoms]
    cx = (min(xs) + max(xs)) / 2
    cy = (min(ys) + max(ys)) / 2
    cz = (min(zs) + max(zs)) / 2
    sx = max(max(xs) - min(xs) + 2 * padding, min_size)
    sy = max(max(ys) - min(ys) + 2 * padding, min_size)
    sz = max(max(zs) - min(zs) + 2 * padding, min_size)
    return {
        "center_x": cx,
        "center_y": cy,
        "center_z": cz,
        "size_x": sx,
        "size_y": sy,
        "size_z": sz,
    }


def render_config(
    box: dict[str, float],
    receptor: str,
    ligand: str,
    exhaustiveness: int,
    residues: list[int],
    pdb: Path,
    chain: str,
    padding: float,
) -> str:
    residue_list = ", ".join(str(r) for r in residues)
    return (
        f"# Generated by derive_docking_box.py\n"
        f"# Source PDB:   {pdb}   (chain {chain})\n"
        f"# Active-site residues: {residue_list}\n"
        f"# Padding around bounding box: {padding:.1f} Å\n"
        f"# Rerun to regenerate:\n"
        f"#   python3 derive_docking_box.py --pdb {pdb} --chain {chain} \\\n"
        f"#       --residues {residue_list.replace(',', '')} \\\n"
        f"#       --padding {padding} --receptor {receptor} --ligand {ligand}\n"
        f"\n"
        f"receptor = {receptor}\n"
        f"ligand = {ligand}\n"
        f"\n"
        f"center_x = {box['center_x']:.3f}\n"
        f"center_y = {box['center_y']:.3f}\n"
        f"center_z = {box['center_z']:.3f}\n"
        f"\n"
        f"size_x = {box['size_x']:.1f}\n"
        f"size_y = {box['size_y']:.1f}\n"
        f"size_z = {box['size_z']:.1f}\n"
        f"\n"
        f"exhaustiveness = {exhaustiveness}\n"
    )


def main() -> int:
    args = parse_args()
    atoms = collect_active_site_atoms(args.pdb, args.chain, args.residues)
    box = compute_box(atoms, padding=args.padding, min_size=args.min_size)

    print(f"Collected {len(atoms)} heavy atoms from {len(args.residues)} residues.")
    print(
        f"Bounding centroid: ({box['center_x']:+.3f}, "
        f"{box['center_y']:+.3f}, {box['center_z']:+.3f})"
    )
    print(
        f"Cube (with {args.padding:.1f} Å padding, min edge "
        f"{args.min_size:.1f} Å): {box['size_x']:.1f} × "
        f"{box['size_y']:.1f} × {box['size_z']:.1f}"
    )
    print()
    print(
        render_config(
            box=box,
            receptor=args.receptor,
            ligand=args.ligand,
            exhaustiveness=args.exhaustiveness,
            residues=args.residues,
            pdb=args.pdb,
            chain=args.chain,
            padding=args.padding,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
