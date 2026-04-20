#!/usr/bin/env python3
"""Rank GS2 catalytic-domain SNPs with the NVIDIA-hosted Evo 2 generate API.

This workflow uses the hosted Evo 2 40B generate endpoint rather than local
CUDA inference. Because the `generate` API returns continuation statistics for
newly generated tokens, this script scores each SNP by a one-token continuation
proxy rather than a full teacher-forced sequence log-likelihood.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import requests
from Bio.Seq import Seq

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
DNA_TO_LOGIT_INDEX = {"A": 65, "C": 67, "G": 71, "T": 84}

ARABIDOPSIS_PREFERRED_CODONS = {
    "A": "GCT",
    "C": "TGC",
    "D": "GAT",
    "E": "GAA",
    "F": "TTC",
    "G": "GGA",
    "H": "CAC",
    "I": "ATC",
    "K": "AAG",
    "L": "TTG",
    "M": "ATG",
    "N": "AAC",
    "P": "CCA",
    "Q": "CAG",
    "R": "AGA",
    "S": "AGC",
    "T": "ACC",
    "V": "GTG",
    "W": "TGG",
    "Y": "TAC",
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
    prompt_sequence: str
    generated_base: str
    sampled_probability: float
    continuation_log_likelihood: float
    wildtype_continuation_log_likelihood: float
    delta_log_likelihood: float
    elapsed_ms: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protein-sequence", default=WILD_TYPE_PROTEIN)
    parser.add_argument("--api-url", default=DEFAULT_API_URL)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--top-k", type=int, default=50)
    parser.add_argument("--max-variants", type=int, default=None)
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


def back_translate(protein_sequence: str) -> str:
    try:
        return "".join(ARABIDOPSIS_PREFERRED_CODONS[aa] for aa in protein_sequence)
    except KeyError as exc:
        raise ValueError(f"Unsupported residue in protein sequence: {exc.args[0]}") from exc


def validate_translation(dna_sequence: str, protein_sequence: str) -> None:
    translated = str(Seq(dna_sequence).translate(to_stop=False))
    if translated != protein_sequence:
        raise ValueError("Back-translation check failed: translated DNA does not match protein.")


def iter_domain_snps(dna_sequence: str, residue_list: list[int] | None = None) -> Iterable[tuple[int, str]]:
    if residue_list:
        indices_to_scan = []
        for residue in residue_list:
            start_nt = (residue - 1) * 3
            indices_to_scan.extend([start_nt, start_nt + 1, start_nt + 2])
    else:
        start_nt = (CATALYTIC_DOMAIN_START_AA - 1) * 3
        end_nt = CATALYTIC_DOMAIN_END_AA * 3
        indices_to_scan = range(start_nt, end_nt)

    for index in indices_to_scan:
        if index < 0 or index >= len(dna_sequence):
            continue
        ref_nt = dna_sequence[index]
        for alt_nt in DNA_ALPHABET:
            if alt_nt != ref_nt:
                yield index, alt_nt


def build_variant_records(
    dna_sequence: str,
    protein_sequence: str,
    max_variants: int | None,
    include_synonymous: bool,
    residue_list: list[int] | None = None,
) -> list[dict[str, object]]:
    variants: list[dict[str, object]] = []
    for dna_index, alt_nt in iter_domain_snps(dna_sequence, residue_list=residue_list):
        mutated = list(dna_sequence)
        mutated[dna_index] = alt_nt
        mutated_sequence = "".join(mutated)
        codon_index = dna_index // 3
        codon_start = codon_index * 3
        wt_codon = dna_sequence[codon_start:codon_start + 3]
        mut_codon = mutated_sequence[codon_start:codon_start + 3]
        wt_aa = protein_sequence[codon_index]
        mut_aa = str(Seq(mut_codon).translate())
        effect = "synonymous" if mut_aa == wt_aa else ("nonsense" if mut_aa == "*" else "missense")
        if effect == "nonsense":
            continue
        if effect == "synonymous" and not include_synonymous:
            continue

        variants.append(
            {
                "sequence": mutated_sequence,
                "dna_index": dna_index,
                "alt_nt": alt_nt,
                "codon_index": codon_index,
                "wt_codon": wt_codon,
                "mut_codon": mut_codon,
                "wt_aa": wt_aa,
                "mut_aa": mut_aa,
                "effect": effect,
            }
        )
        if max_variants is not None and len(variants) >= max_variants:
            break
    return variants


def chunked(items: list[dict[str, object]], size: int) -> Iterable[list[dict[str, object]]]:
    for start in range(0, len(items), size):
        yield items[start:start + size]


def safe_log(probability: float) -> float:
    clamped = max(min(probability, 1.0), 1e-12)
    return math.log(clamped)


def variant_key(dna_pos_1based: int, alt_nt: str) -> str:
    return f"{dna_pos_1based}:{alt_nt}"


def load_existing_scores(output_path: Path) -> list[VariantScore]:
    if not output_path.exists():
        return []

    scores: list[VariantScore] = []
    with output_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
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
                    prompt_sequence=row["prompt_sequence"],
                    generated_base=row["generated_base"],
                    sampled_probability=float(row["sampled_probability"]),
                    continuation_log_likelihood=float(row["continuation_log_likelihood"]),
                    wildtype_continuation_log_likelihood=float(row["wildtype_continuation_log_likelihood"]),
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


class Evo2GenerateClient:
    def __init__(self, api_url: str, timeout: int, dry_run: bool) -> None:
        self.api_url = api_url
        self.timeout = timeout
        self.dry_run = dry_run
        self.api_key = os.getenv("NVCF_RUN_KEY")
        if not dry_run and not self.api_key:
            raise RuntimeError("NVCF_RUN_KEY is not set.")

    def build_payload(self, sequence: str) -> dict[str, object]:
        return {
            "sequence": sequence,
            "num_tokens": 1,
            "top_k": 1,
            "temperature": 1.0,
            "enable_logits": True,
            "enable_sampled_probs": True,
        }

    def generate(self, sequence: str) -> dict[str, object]:
        payload = self.build_payload(sequence)
        if self.dry_run:
            print(json.dumps(payload, indent=2))
            raise SystemExit(0)

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
        data = response.json()
        if "sequence" not in data:
            raise ValueError(f"Generate response missing 'sequence': {data}")
        return data


def parse_generate_score(prompt_sequence: str, response: dict[str, object]) -> dict[str, object]:
    generated = str(response["sequence"])
    generated_base = generated[0] if generated else ""

    sampled_probs = response.get("sampled_probs") or []
    sampled_probability = float(sampled_probs[0]) if sampled_probs else 0.0

    logits = response.get("logits") or []
    target_token_logit = None
    if logits and isinstance(logits, list) and logits[0]:
        row = logits[0]
        if len(row) > max(DNA_TO_LOGIT_INDEX.values()):
            target_token_logit = float(row[DNA_TO_LOGIT_INDEX.get(generated_base, 65)])

    return {
        "generated_base": generated_base,
        "sampled_probability": sampled_probability,
        "continuation_log_likelihood": safe_log(sampled_probability),
        "target_token_logit": target_token_logit,
        "elapsed_ms": int(response.get("elapsed_ms", 0)),
    }


def score_variants(
    client: Evo2GenerateClient,
    dna_sequence: str,
    protein_sequence: str,
    batch_size: int,
    max_variants: int | None,
    include_synonymous: bool,
    output_path: Path,
    residue_list: list[int] | None = None,
) -> tuple[list[VariantScore], int]:
    variants = build_variant_records(
        dna_sequence=dna_sequence,
        protein_sequence=protein_sequence,
        max_variants=max_variants,
        include_synonymous=include_synonymous,
        residue_list=residue_list,
    )
    existing_scores = [
        score
        for score in load_existing_scores(output_path)
        if score.effect != "nonsense" and (include_synonymous or score.effect != "synonymous")
    ]
    completed_keys = {
        variant_key(score.dna_pos_1based, score.alt_nt)
        for score in existing_scores
    }
    pending_variants = [
        variant
        for variant in variants
        if variant_key(int(variant["dna_index"]) + 1, str(variant["alt_nt"])) not in completed_keys
    ]

    if existing_scores:
        wildtype_log_likelihood = existing_scores[0].wildtype_continuation_log_likelihood
    else:
        wildtype_response = client.generate(dna_sequence)
        wildtype_score = parse_generate_score(dna_sequence, wildtype_response)
        wildtype_log_likelihood = float(wildtype_score["continuation_log_likelihood"])

    new_results: list[VariantScore] = []
    for batch in chunked(pending_variants, batch_size):
        batch_results: list[VariantScore] = []
        for variant in batch:
            mutant_response = client.generate(str(variant["sequence"]))
            mutant_score = parse_generate_score(str(variant["sequence"]), mutant_response)

            batch_results.append(
                VariantScore(
                    rank=0,
                    dna_pos_1based=int(variant["dna_index"]) + 1,
                    ref_nt=dna_sequence[int(variant["dna_index"])],
                    alt_nt=str(variant["alt_nt"]),
                    aa_pos_1based=int(variant["codon_index"]) + 1,
                    wt_codon=str(variant["wt_codon"]),
                    mut_codon=str(variant["mut_codon"]),
                    wt_aa=str(variant["wt_aa"]),
                    mut_aa=str(variant["mut_aa"]),
                    effect=str(variant["effect"]),
                    prompt_sequence=str(variant["sequence"]),
                    generated_base=str(mutant_score["generated_base"]),
                    sampled_probability=float(mutant_score["sampled_probability"]),
                    continuation_log_likelihood=float(mutant_score["continuation_log_likelihood"]),
                    wildtype_continuation_log_likelihood=wildtype_log_likelihood,
                    delta_log_likelihood=float(mutant_score["continuation_log_likelihood"])
                    - wildtype_log_likelihood,
                    elapsed_ms=int(mutant_score["elapsed_ms"]),
                )
            )

        append_scores(output_path, batch_results)
        new_results.extend(batch_results)
        for score in batch_results:
            completed_keys.add(variant_key(score.dna_pos_1based, score.alt_nt))

    results = existing_scores + new_results
    results.sort(key=lambda item: item.delta_log_likelihood, reverse=True)
    ranked_results: list[VariantScore] = []
    for index, score in enumerate(results, start=1):
        payload = dict(score.__dict__)
        payload["rank"] = index
        ranked_results.append(VariantScore(**payload))
    return ranked_results, len(pending_variants)


def write_csv(output_path: Path, scores: list[VariantScore]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(VariantScore.__dataclass_fields__.keys()))
        writer.writeheader()
        for score in scores:
            writer.writerow(score.__dict__)


def main() -> int:
    args = parse_args()
    protein_sequence = args.protein_sequence.strip().upper()
    dna_sequence = back_translate(protein_sequence)
    validate_translation(dna_sequence, protein_sequence)

    client = Evo2GenerateClient(
        api_url=args.api_url,
        timeout=args.timeout,
        dry_run=args.dry_run,
    )
    scores = score_variants(
        client=client,
        dna_sequence=dna_sequence,
        protein_sequence=protein_sequence,
        batch_size=args.batch_size,
        max_variants=args.max_variants,
        include_synonymous=args.include_synonymous,
        output_path=args.output,
        residue_list=args.residues,
    )
    ranked_scores, pending_count = scores
    write_csv(args.output, ranked_scores)

    domain_label = (
        f"residues {args.residues}"
        if args.residues
        else f"GS2 catalytic-domain residues {CATALYTIC_DOMAIN_START_AA}-{CATALYTIC_DOMAIN_END_AA}"
    )
    print(f"Scored {len(ranked_scores)} SNPs across {domain_label}.")
    print(f"New variants evaluated in this run: {pending_count}")
    print(f"Ranked variants written to {args.output}")
    for score in ranked_scores[: args.top_k]:
        print(
            f"{score.rank:>3}  nt {score.dna_pos_1based}:{score.ref_nt}>{score.alt_nt}  "
            f"aa {score.aa_pos_1based}:{score.wt_aa}>{score.mut_aa}  "
            f"{score.effect:<10}  delta={score.delta_log_likelihood:.4f}  "
            f"p={score.sampled_probability:.6f}  next={score.generated_base or '?'}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
