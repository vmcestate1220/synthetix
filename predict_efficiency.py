#!/usr/bin/env python3
"""Rank GS2 catalytic-domain SNPs via direct Evo 2 logit analysis.

Scoring method
--------------
For each candidate SNP at DNA position ``i`` (0-based), the script prompts
Evo 2 40B with the *left-context only* (bases ``0..i-1``) and requests one
continuation token with ``enable_logits=True``. The logits at position ``i``
score every possible next byte; we read the A/C/G/T entries directly and
report

    delta_log_likelihood = logit(alt) - logit(ref)

which equals ``log P(alt | left) - log P(ref | left)`` under any shared
softmax normalization (the partition function cancels in the difference).
No wild-type baseline call is needed: reference and alternate are scored from
the same logit row, so per-position variance cancels exactly.

Interpretation caveat
---------------------
Evo 2 is a genomic foundation model trained on observed nucleotide sequences.
``delta_log_likelihood`` measures the **evolutionary plausibility** of the
alternate base in its local context — not catalytic improvement, not thermal
stability, not docking affinity. A positive Δ-LL at a conserved active-site
residue is most often a signal that the proxy is poorly calibrated for that
position, not that the mutation enhances the enzyme. See
``validate_against_kinetics.py`` and README § Methodology for how to anchor
these scores against experimentally-measured functional data before trusting
the ranking.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import requests
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm

WILD_TYPE_PROTEIN = (
    "MAQILAASPTCQMRVPKHSSVIASSSKLWSSVVLKQKKQSNNKVRGFRVLALQSDNSTVNRVETLL"
    "NLDTKPYSDRIIAEYIWIGGSGIDLRSKSRTIEKPVEDPSELPKWNYDGSSTGQAPGEDSEVILYP"
    "QAIFRDPFRGGNNILVICDTWTPAGEPIPTNKRAKAAEIFSNKKVSGEVPWFGIEQEYTLLQQNVK"
    "WPLGWPVGAFPGPQGPYYCGVGADKIWGRDISDAHYKACLYAGINISGTNGEVMPGQWEFQVGPSV"
    "GIDAGDHVWCARYLLERITEQAGVVLTLDPKPIEGDWNGAGCHTNYSTKSMREEGGFEVIKKAILN"
    "LSLRHKEHISAYGEGNERRLTGKHETASIDQFSWGVANRGCSIRVGRDTEAKGKGYLEDRRPASNM"
    "DPYIVTSLLAETTLLWEPTLEAEALAAQKLSLNV"
)

CATALYTIC_DOMAIN_START_AA = 111
CATALYTIC_DOMAIN_END_AA = 430
DEFAULT_API_URL = "https://health.api.nvidia.com/v1/biology/arc/evo2-40b/generate"
DNA_ALPHABET = ("A", "C", "G", "T")
# Evo 2 tokenizes DNA at the byte level; these indices correspond to ASCII
# codes for the uppercase bases (A=65, C=67, G=71, T=84) in the returned
# logits row.
DNA_TO_LOGIT_INDEX = {"A": 65, "C": 67, "G": 71, "T": 84}

ARABIDOPSIS_PREFERRED_CODONS = {
    "A": "GCT", "C": "TGC", "D": "GAT", "E": "GAA", "F": "TTC",
    "G": "GGA", "H": "CAC", "I": "ATC", "K": "AAG", "L": "TTG",
    "M": "ATG", "N": "AAC", "P": "CCA", "Q": "CAG", "R": "AGA",
    "S": "AGC", "T": "ACC", "V": "GTG", "W": "TGG", "Y": "TAC",
}


@dataclass(frozen=True)
class VariantScore:
    rank: int
    dna_pos_1based: int
    ref_nt: str
    alt_nt: str
    aa_pos_1based: int
    wt_codon: str
    mut_codon: str
    wt_aa: str
    mut_aa: str
    effect: str
    ref_logit: float
    alt_logit: float
    ref_probability: float
    alt_probability: float
    delta_log_likelihood: float
    elapsed_ms: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cds-fasta",
        type=Path,
        default=None,
        help=(
            "Path to a FASTA containing the real GS2 CDS (e.g. TAIR "
            "AT5G35630.1 or NCBI NM_122499). When omitted, the script "
            "falls back to a synthetic Arabidopsis-preferred-codon "
            "back-translation and prints a warning."
        ),
    )
    parser.add_argument(
        "--protein-fasta",
        type=Path,
        default=None,
        help="Optional protein FASTA used to validate the supplied CDS.",
    )
    parser.add_argument("--protein-sequence", default=WILD_TYPE_PROTEIN)
    parser.add_argument("--api-url", default=DEFAULT_API_URL)
    parser.add_argument(
        "--concurrency",
        type=int,
        default=8,
        help=(
            "Number of API requests submitted in parallel per round. Each "
            "round flushes to the output CSV before the next begins, so this "
            "also serves as the resume/checkpoint granularity. Tune down if "
            "you hit rate limits (HTTP 429) or up if the API tolerates it."
        ),
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=None,
        help="Deprecated alias for --concurrency. Prefer --concurrency.",
    )
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--top-k", type=int, default=50)
    parser.add_argument("--max-positions", type=int, default=None)
    parser.add_argument(
        "--include-synonymous",
        action="store_true",
        help="Include synonymous SNPs. By default they are excluded.",
    )
    parser.add_argument(
        "--residues",
        nargs="+",
        type=int,
        help="Specific residues (1-based) to scan. If not provided, scans the whole catalytic domain.",
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("gs2_catalytic_domain_snp_ranking.csv"),
    )
    return parser.parse_args()


def load_fasta_record(path: Path) -> str:
    record = next(SeqIO.parse(str(path), "fasta"), None)
    if record is None:
        raise ValueError(f"No records found in FASTA: {path}")
    return str(record.seq).upper().replace("U", "T")


def back_translate(protein_sequence: str) -> str:
    try:
        return "".join(ARABIDOPSIS_PREFERRED_CODONS[aa] for aa in protein_sequence)
    except KeyError as exc:
        raise ValueError(f"Unsupported residue in protein sequence: {exc.args[0]}") from exc


def validate_translation(dna_sequence: str, protein_sequence: str) -> None:
    # Trim to a multiple of 3 and drop a trailing stop codon if present.
    usable = dna_sequence[: len(dna_sequence) - (len(dna_sequence) % 3)]
    translated = str(Seq(usable).translate(to_stop=False)).rstrip("*")
    if not translated.startswith(protein_sequence):
        raise ValueError(
            "CDS does not translate to the expected protein. "
            f"Expected prefix: {protein_sequence[:30]}...; got: {translated[:30]}..."
        )


def load_dna_sequence(args: argparse.Namespace, protein_sequence: str) -> str:
    if args.cds_fasta is not None:
        dna = load_fasta_record(args.cds_fasta)
        validate_translation(dna, protein_sequence)
        return dna
    print(
        "WARNING: --cds-fasta not provided. Using a synthetic "
        "Arabidopsis-preferred-codon back-translation. This is out-of-distribution "
        "for Evo 2 and will bias scores — supply the real AT5G35630 CDS via "
        "--cds-fasta for trustworthy results.",
        file=sys.stderr,
    )
    dna = back_translate(protein_sequence)
    validate_translation(dna, protein_sequence)
    return dna


def iter_domain_positions(
    dna_sequence: str, residue_list: list[int] | None = None
) -> list[int]:
    if residue_list:
        indices: list[int] = []
        for residue in residue_list:
            start_nt = (residue - 1) * 3
            indices.extend([start_nt, start_nt + 1, start_nt + 2])
    else:
        start_nt = (CATALYTIC_DOMAIN_START_AA - 1) * 3
        end_nt = CATALYTIC_DOMAIN_END_AA * 3
        indices = list(range(start_nt, end_nt))
    return [i for i in indices if 0 <= i < len(dna_sequence)]


def classify_variant(
    dna_sequence: str, dna_index: int, alt_nt: str, protein_sequence: str
) -> dict[str, object]:
    codon_index = dna_index // 3
    codon_start = codon_index * 3
    wt_codon = dna_sequence[codon_start:codon_start + 3]
    mutated = list(dna_sequence)
    mutated[dna_index] = alt_nt
    mut_codon = "".join(mutated[codon_start:codon_start + 3])
    wt_aa = protein_sequence[codon_index] if codon_index < len(protein_sequence) else "?"
    mut_aa = str(Seq(mut_codon).translate())
    if mut_aa == "*":
        effect = "nonsense"
    elif mut_aa == wt_aa:
        effect = "synonymous"
    else:
        effect = "missense"
    return {
        "dna_index": dna_index,
        "alt_nt": alt_nt,
        "codon_index": codon_index,
        "wt_codon": wt_codon,
        "mut_codon": mut_codon,
        "wt_aa": wt_aa,
        "mut_aa": mut_aa,
        "effect": effect,
    }


def variant_key(dna_pos_1based: int, alt_nt: str) -> str:
    return f"{dna_pos_1based}:{alt_nt}"


def load_existing_scores(output_path: Path) -> list[VariantScore]:
    if not output_path.exists():
        return []

    scores: list[VariantScore] = []
    with output_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        expected = set(VariantScore.__dataclass_fields__.keys())
        actual = set(reader.fieldnames or [])
        if not expected.issubset(actual):
            missing = expected.difference(actual)
            raise ValueError(
                f"Existing CSV {output_path} is missing required columns: "
                f"{sorted(missing)}. This file was produced by an older "
                "version of predict_efficiency.py with incompatible scoring; "
                "delete it and re-run to regenerate under the current schema."
            )
        for row in reader:
            scores.append(
                VariantScore(
                    rank=int(row["rank"]),
                    dna_pos_1based=int(row["dna_pos_1based"]),
                    ref_nt=row["ref_nt"],
                    alt_nt=row["alt_nt"],
                    aa_pos_1based=int(row["aa_pos_1based"]),
                    wt_codon=row["wt_codon"],
                    mut_codon=row["mut_codon"],
                    wt_aa=row["wt_aa"],
                    mut_aa=row["mut_aa"],
                    effect=row["effect"],
                    ref_logit=float(row["ref_logit"]),
                    alt_logit=float(row["alt_logit"]),
                    ref_probability=float(row["ref_probability"]),
                    alt_probability=float(row["alt_probability"]),
                    delta_log_likelihood=float(row["delta_log_likelihood"]),
                    elapsed_ms=int(row["elapsed_ms"]),
                )
            )
    return scores


def append_scores(output_path: Path, scores: list[VariantScore]) -> None:
    if not scores:
        return
    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not output_path.exists() or output_path.stat().st_size == 0
    with output_path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(VariantScore.__dataclass_fields__.keys()))
        if write_header:
            writer.writeheader()
        for score in scores:
            writer.writerow(score.__dict__)


def softmax_over_bases(logits: dict[str, float]) -> dict[str, float]:
    max_logit = max(logits.values())
    exps = {b: math.exp(l - max_logit) for b, l in logits.items()}
    total = sum(exps.values())
    return {b: e / total for b, e in exps.items()}


class Evo2GenerateClient:
    def __init__(self, api_url: str, timeout: int, dry_run: bool) -> None:
        self.api_url = api_url
        self.timeout = timeout
        self.dry_run = dry_run
        self.api_key = os.getenv("NVCF_RUN_KEY")
        if not dry_run and not self.api_key:
            raise RuntimeError("NVCF_RUN_KEY is not set.")

    def build_payload(self, left_context: str) -> dict[str, object]:
        return {
            "sequence": left_context,
            "num_tokens": 1,
            "top_k": 1,
            "temperature": 1.0,
            "enable_logits": True,
            "enable_sampled_probs": False,
        }

    def generate(self, left_context: str) -> dict[str, object]:
        payload = self.build_payload(left_context)
        if self.dry_run:
            import json
            print(json.dumps({**payload, "sequence": f"<{len(left_context)} nt>"}, indent=2))
            return {"logits": [[0.0] * (max(DNA_TO_LOGIT_INDEX.values()) + 1)], "elapsed_ms": 0}

        last_exc: Exception | None = None
        for attempt in range(3):
            try:
                response = requests.post(
                    self.api_url,
                    headers={
                        "Authorization": f"Bearer {self.api_key}",
                        "Accept": "application/json",
                        "Content-Type": "application/json",
                    },
                    json=payload,
                    timeout=self.timeout,
                )
                response.raise_for_status()
                return response.json()
            except (requests.exceptions.HTTPError, requests.exceptions.ConnectionError) as exc:
                last_exc = exc
                if attempt == 2:
                    break
                import time
                time.sleep(2 ** attempt)
        assert last_exc is not None
        raise last_exc


def extract_base_logits(response: dict[str, object]) -> dict[str, float]:
    """Read the A/C/G/T logits from the first (and only) generated position.

    Raises on any missing / malformed field so a bad API response never
    silently degrades into a plausible-looking Δ-LL.
    """
    logits = response.get("logits")
    if not isinstance(logits, list) or not logits:
        raise ValueError(f"Generate response has no 'logits' field: {response!r}")
    row = logits[0]
    if not isinstance(row, list) or len(row) <= max(DNA_TO_LOGIT_INDEX.values()):
        raise ValueError(
            f"Logits row is too short to contain DNA byte indices: len={len(row) if isinstance(row, list) else type(row).__name__}"
        )
    try:
        return {base: float(row[idx]) for base, idx in DNA_TO_LOGIT_INDEX.items()}
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Non-numeric logit value in response: {exc}") from exc


def _emit_for_position(
    dna_sequence: str,
    protein_sequence: str,
    dna_index: int,
    include_synonymous: bool,
    completed: set[str],
) -> list[dict[str, object]]:
    """Classify the three alt variants at one DNA position and return the
    ones that need to be scored (not nonsense, not skipped synonymous, not
    already in the resume set)."""
    ref_nt = dna_sequence[dna_index]
    emit: list[dict[str, object]] = []
    for alt in DNA_ALPHABET:
        if alt == ref_nt:
            continue
        info = classify_variant(dna_sequence, dna_index, alt, protein_sequence)
        if info["effect"] == "nonsense":
            continue
        if info["effect"] == "synonymous" and not include_synonymous:
            continue
        if variant_key(dna_index + 1, str(info["alt_nt"])) in completed:
            continue
        emit.append(info)
    return emit


def _scores_from_response(
    dna_sequence: str,
    dna_index: int,
    emit: list[dict[str, object]],
    response: dict[str, object],
) -> list[VariantScore]:
    ref_nt = dna_sequence[dna_index]
    base_logits = extract_base_logits(response)
    base_probs = softmax_over_bases(base_logits)
    elapsed = int(response.get("elapsed_ms", 0) or 0)
    ref_logit = base_logits[ref_nt]
    ref_prob = base_probs[ref_nt]

    out: list[VariantScore] = []
    for info in emit:
        alt = str(info["alt_nt"])
        out.append(
            VariantScore(
                rank=0,
                dna_pos_1based=dna_index + 1,
                ref_nt=ref_nt,
                alt_nt=alt,
                aa_pos_1based=int(info["codon_index"]) + 1,
                wt_codon=str(info["wt_codon"]),
                mut_codon=str(info["mut_codon"]),
                wt_aa=str(info["wt_aa"]),
                mut_aa=str(info["mut_aa"]),
                effect=str(info["effect"]),
                ref_logit=ref_logit,
                alt_logit=base_logits[alt],
                ref_probability=ref_prob,
                alt_probability=base_probs[alt],
                delta_log_likelihood=base_logits[alt] - ref_logit,
                elapsed_ms=elapsed,
            )
        )
    return out


def score_positions(
    client: Evo2GenerateClient,
    dna_sequence: str,
    protein_sequence: str,
    positions: list[int],
    include_synonymous: bool,
    concurrency: int,
    output_path: Path,
) -> tuple[list[VariantScore], int]:
    if concurrency < 1:
        raise ValueError(f"concurrency must be >= 1, got {concurrency}")

    existing = load_existing_scores(output_path)
    filtered_existing = [
        s for s in existing
        if s.effect != "nonsense" and (include_synonymous or s.effect != "synonymous")
    ]
    completed = {variant_key(s.dna_pos_1based, s.alt_nt) for s in filtered_existing}

    # Build a work list of (position, emit_infos) skipping positions where
    # every alt is already scored or filtered out. Doing this up front lets
    # the progress bar reflect real work and keeps the thread pool busy.
    work: list[tuple[int, list[dict[str, object]]]] = []
    for dna_index in positions:
        emit = _emit_for_position(
            dna_sequence, protein_sequence, dna_index, include_synonymous, completed
        )
        if emit:
            work.append((dna_index, emit))

    new_results: list[VariantScore] = []
    progress = tqdm(total=len(work), desc="Scoring positions")

    # ThreadPoolExecutor — Evo 2 calls are IO-bound HTTP, the GIL is
    # released around requests.post, and the endpoint handles its own
    # concurrency. If the API starts returning HTTP 429s the existing
    # retry-with-backoff in Evo2GenerateClient.generate handles it; drop
    # --concurrency further if the backoff becomes the bottleneck.
    with ThreadPoolExecutor(max_workers=concurrency) as executor:
        for chunk_start in range(0, len(work), concurrency):
            chunk = work[chunk_start:chunk_start + concurrency]
            futures = [
                executor.submit(client.generate, dna_sequence[:dna_index])
                for dna_index, _ in chunk
            ]
            chunk_scores: list[VariantScore] = []
            for (dna_index, emit), future in zip(chunk, futures):
                response = future.result()
                chunk_scores.extend(
                    _scores_from_response(dna_sequence, dna_index, emit, response)
                )
                progress.update(1)

            append_scores(output_path, chunk_scores)
            new_results.extend(chunk_scores)
            for score in chunk_scores:
                completed.add(variant_key(score.dna_pos_1based, score.alt_nt))

    progress.close()

    combined = filtered_existing + new_results
    combined.sort(key=lambda item: item.delta_log_likelihood, reverse=True)
    ranked: list[VariantScore] = []
    for index, score in enumerate(combined, start=1):
        payload = dict(score.__dict__)
        payload["rank"] = index
        ranked.append(VariantScore(**payload))
    return ranked, len(new_results)


def write_csv(output_path: Path, scores: list[VariantScore]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(VariantScore.__dataclass_fields__.keys()))
        writer.writeheader()
        for score in scores:
            writer.writerow(score.__dict__)


def main() -> int:
    args = parse_args()
    if args.protein_fasta is not None:
        protein_sequence = load_fasta_record(args.protein_fasta).strip("*")
    else:
        protein_sequence = args.protein_sequence.strip().upper()
    dna_sequence = load_dna_sequence(args, protein_sequence)

    client = Evo2GenerateClient(
        api_url=args.api_url,
        timeout=args.timeout,
        dry_run=args.dry_run,
    )
    positions = iter_domain_positions(dna_sequence, residue_list=args.residues)
    if args.max_positions is not None:
        positions = positions[: args.max_positions]

    concurrency = args.batch_size if args.batch_size is not None else args.concurrency
    if args.batch_size is not None:
        print(
            "WARNING: --batch-size is deprecated; it now aliases --concurrency.",
            file=sys.stderr,
        )

    ranked, new_count = score_positions(
        client=client,
        dna_sequence=dna_sequence,
        protein_sequence=protein_sequence,
        positions=positions,
        include_synonymous=args.include_synonymous,
        concurrency=concurrency,
        output_path=args.output,
    )
    write_csv(args.output, ranked)

    domain_label = (
        f"residues {args.residues}"
        if args.residues
        else f"GS2 catalytic-domain residues {CATALYTIC_DOMAIN_START_AA}-{CATALYTIC_DOMAIN_END_AA}"
    )
    print(f"Scored {len(ranked)} SNPs across {domain_label}.")
    print(f"New variants evaluated in this run: {new_count}")
    print(f"Ranked variants written to {args.output}")
    print(
        "Reminder: delta_log_likelihood reflects Evo 2's evolutionary prior, "
        "not catalytic improvement. Run validate_against_kinetics.py against "
        "published mutant kinetic data before acting on the top of this list."
    )
    for score in ranked[: args.top_k]:
        print(
            f"{score.rank:>3}  nt {score.dna_pos_1based}:{score.ref_nt}>{score.alt_nt}  "
            f"aa {score.aa_pos_1based}:{score.wt_aa}>{score.mut_aa}  "
            f"{score.effect:<10}  delta={score.delta_log_likelihood:+.4f}  "
            f"p(alt)={score.alt_probability:.4f}  p(ref)={score.ref_probability:.4f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
