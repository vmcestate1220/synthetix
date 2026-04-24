#!/usr/bin/env python3
"""Map mutation residue numbers from a GS homolog onto Arabidopsis GS2.

The bulk of the published GS mutagenesis record lives on *E. coli*, *Salmonella*,
or plant GS1 (cytosolic) — not on *Arabidopsis* GS2 directly. To use that data
to validate a Synthetix scan (which is numbered against UniProt Q43127), we
need to translate residue indices across homologs via a pairwise alignment.

This script wraps Biopython's ``PairwiseAligner`` (BLOSUM62, global) and:

    1. Aligns a source FASTA (e.g. E. coli GS GlnA or maize GS1a) against
       the GS2 reference protein sequence.
    2. For each source mutation (e.g. ``E297A``), reports the aligned GS2
       residue, the GS2 residue's WT identity, and the local alignment
       identity in a ±10 window (a quick sanity check that the mapping is
       in a well-conserved region).
    3. Emits a CSV row ready to paste into ``known_mutants.csv``.

Usage::

    python3 align_to_gs2.py \\
        --source EcoliGlnA.fasta \\
        --mutations E297A D50A \\
        --measured-effects -4.6 -3.1 \\
        --source-label EcoliGS

Confidence guidance: if the local identity at the mapped position is below
~50% or the source WT residue differs from the GS2 WT residue, do NOT use
the mapping for validation — the position is outside the conserved core.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Iterable

from Bio import Align, SeqIO

# Hardcoded from predict_efficiency.py so this script is usable even when
# the user runs it before generating any scan.
GS2_REFERENCE_PROTEIN = (
    "MAQILAASPTCQMRVPKHSSVIASSSKLWSSVVLKQKKQSNNKVRGFRVLALQSDNSTVNRVETLL"
    "NLDTKPYSDRIIAEYIWIGGSGIDLRSKSRTIEKPVEDPSELPKWNYDGSSTGQAPGEDSEVILYP"
    "QAIFRDPFRGGNNILVICDTWTPAGEPIPTNKRAKAAEIFSNKKVSGEVPWFGIEQEYTLLQQNVK"
    "WPLGWPVGAFPGPQGPYYCGVGADKIWGRDISDAHYKACLYAGINISGTNGEVMPGQWEFQVGPSV"
    "GIDAGDHVWCARYLLERITEQAGVVLTLDPKPIEGDWNGAGCHTNYSTKSMREEGGFEVIKKAILN"
    "LSLRHKEHISAYGEGNERRLTGKHETASIDQFSWGVANRGCSIRVGRDTEAKGKGYLEDRRPASNM"
    "DPYIVTSLLAETTLLWEPTLEAEALAAQKLSLNV"
)

MUTATION_PATTERN = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source", type=Path, required=True,
        help="FASTA with the homolog protein sequence (e.g. E. coli GlnA).",
    )
    parser.add_argument(
        "--target", type=Path, default=None,
        help="Target GS2 FASTA. If omitted, the hardcoded Q43127 reference is used.",
    )
    parser.add_argument(
        "--source-label", required=True,
        help="Short tag recorded in CSV rows (e.g. EcoliGS, PhaseolusGS1).",
    )
    parser.add_argument(
        "--mutations", nargs="+", required=True,
        help="Source-numbering mutations (e.g. E297A D50A).",
    )
    parser.add_argument(
        "--measured-effects", nargs="+", type=float, required=True,
        help="Measured effect per mutation (same order as --mutations).",
    )
    parser.add_argument(
        "--effect-type", default="log2_relative_kcat",
        help="Units/convention for the measured effect (consistent across rows).",
    )
    parser.add_argument(
        "--citation", default="",
        help="Primary citation string for the measurement source.",
    )
    parser.add_argument(
        "--identity-window", type=int, default=10,
        help="Residue window around each mapped position for local-identity sanity check.",
    )
    parser.add_argument(
        "--out-csv", type=Path, default=None,
        help="Optional output CSV path (ready for known_mutants.csv).",
    )
    return parser.parse_args()


def load_protein(path: Path) -> str:
    record = next(SeqIO.parse(str(path), "fasta"), None)
    if record is None:
        raise ValueError(f"No records in {path}")
    return str(record.seq).upper().rstrip("*")


def align(source: str, target: str) -> Align.Alignment:
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(source, target)
    return alignments[0]


def build_source_to_target_map(alignment: Align.Alignment) -> dict[int, int]:
    """Return a 1-based residue index map from source to target.

    Gaps map to ``None`` implicitly (omitted from the returned dict).
    """
    aligned_source, aligned_target = str(alignment[0]), str(alignment[1])
    mapping: dict[int, int] = {}
    src_idx = 0
    tgt_idx = 0
    for s, t in zip(aligned_source, aligned_target):
        if s != "-":
            src_idx += 1
        if t != "-":
            tgt_idx += 1
        if s != "-" and t != "-":
            mapping[src_idx] = tgt_idx
    return mapping


def local_identity(
    alignment: Align.Alignment, source_residue_1based: int, window: int
) -> float:
    aligned_source, aligned_target = str(alignment[0]), str(alignment[1])
    # Find the alignment column that holds this source residue.
    src_idx = 0
    col = None
    for i, s in enumerate(aligned_source):
        if s != "-":
            src_idx += 1
        if src_idx == source_residue_1based:
            col = i
            break
    if col is None:
        return 0.0
    lo = max(0, col - window)
    hi = min(len(aligned_source), col + window + 1)
    paired = matched = 0
    for s, t in zip(aligned_source[lo:hi], aligned_target[lo:hi]):
        if s != "-" and t != "-":
            paired += 1
            if s == t:
                matched += 1
    return matched / paired if paired else 0.0


def parse_mutation(label: str) -> tuple[str, int, str]:
    m = MUTATION_PATTERN.match(label.strip().upper())
    if not m:
        raise ValueError(f"Mutation must look like 'E297A', got {label!r}")
    return m.group(1), int(m.group(2)), m.group(3)


def main() -> int:
    args = parse_args()
    if len(args.mutations) != len(args.measured_effects):
        print("--mutations and --measured-effects must have equal length.", file=sys.stderr)
        return 2

    source = load_protein(args.source)
    if args.target is not None:
        target = load_protein(args.target)
    else:
        target = GS2_REFERENCE_PROTEIN
    alignment = align(source, target)
    mapping = build_source_to_target_map(alignment)

    header = ["mutation", "measured_effect", "effect_type", "citation", "notes"]
    rows: list[list[str]] = []
    print(f"# Alignment score: {alignment.score}")
    print(f"# Source length: {len(source)}   Target length: {len(target)}")
    print(f"# {'src':^8} {'tgt':^8} {'src_WT':^7} {'tgt_WT':^7} {'mut':^5} {'local_id':^9} {'effect':^10}")

    for label, effect in zip(args.mutations, args.measured_effects):
        src_wt, src_pos, new_aa = parse_mutation(label)
        if src_pos > len(source):
            print(f"# WARNING: position {src_pos} exceeds source length {len(source)}", file=sys.stderr)
            continue
        actual_src_wt = source[src_pos - 1]
        if actual_src_wt != src_wt:
            print(
                f"# WARNING: mutation {label} claims WT={src_wt} but source[{src_pos}]={actual_src_wt}",
                file=sys.stderr,
            )
        tgt_pos = mapping.get(src_pos)
        if tgt_pos is None:
            print(f"# {label}: source residue aligns to a gap — no GS2 equivalent")
            continue
        tgt_wt = target[tgt_pos - 1]
        ident = local_identity(alignment, src_pos, args.identity_window)
        gs2_mutation = f"{tgt_wt}{tgt_pos}{new_aa}"
        notes = (
            f"mapped from {args.source_label} {label}; "
            f"local identity {ident:.0%} in ±{args.identity_window} window"
            + ("; WT mismatch — verify alignment" if tgt_wt != src_wt else "")
        )
        print(
            f"# {src_pos:>6}    {tgt_pos:>6}    {actual_src_wt:^7} {tgt_wt:^7} {new_aa:^5} {ident*100:>6.1f}%    {effect:+.3f}"
        )
        rows.append([gs2_mutation, f"{effect:+.4f}", args.effect_type, args.citation, notes])

    print()
    print(",".join(header))
    for r in rows:
        print(",".join(r))

    if args.out_csv is not None:
        import csv
        write_header = not args.out_csv.exists() or args.out_csv.stat().st_size == 0
        with args.out_csv.open("a", newline="") as handle:
            writer = csv.writer(handle)
            if write_header:
                writer.writerow(header)
            for r in rows:
                writer.writerow(r)
        print(f"\nAppended {len(rows)} rows to {args.out_csv}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
