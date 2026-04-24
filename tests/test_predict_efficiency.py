"""Unit tests for predict_efficiency core helpers.

These tests exercise the code paths that don't hit the Evo 2 API — variant
enumeration, classification, logit parsing, CSV round-trip — so they run
offline in <1s and catch schema/offset regressions before any API call.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import pytest

from predict_efficiency import (
    ARABIDOPSIS_PREFERRED_CODONS,
    CATALYTIC_DOMAIN_END_AA,
    CATALYTIC_DOMAIN_START_AA,
    DNA_ALPHABET,
    DNA_TO_LOGIT_INDEX,
    VariantScore,
    append_scores,
    back_translate,
    classify_variant,
    extract_base_logits,
    iter_domain_positions,
    load_existing_scores,
    softmax_over_bases,
    validate_translation,
    variant_key,
)


# Small synthetic protein/DNA fixture so the tests run fast and don't depend
# on the 430-aa GS2 constant.
TEST_PROTEIN = "MKWA"  # ATG AAG TGG GCT under Arabidopsis-preferred codons
TEST_DNA = "ATGAAGTGGGCT"


def test_back_translate_round_trips():
    dna = back_translate(TEST_PROTEIN)
    assert dna == TEST_DNA
    validate_translation(dna, TEST_PROTEIN)


def test_back_translate_rejects_unknown_residue():
    with pytest.raises(ValueError, match="Unsupported residue"):
        back_translate("MKX")


def test_validate_translation_rejects_mismatch():
    with pytest.raises(ValueError, match="does not translate"):
        validate_translation("AAAAAA", TEST_PROTEIN)


def test_validate_translation_tolerates_trailing_stop():
    dna_with_stop = TEST_DNA + "TAA"
    validate_translation(dna_with_stop, TEST_PROTEIN)  # should not raise


def test_iter_domain_positions_residue_list_maps_correctly():
    # Residue 2 (1-based) → DNA indices 3, 4, 5
    positions = iter_domain_positions(TEST_DNA, residue_list=[2])
    assert positions == [3, 4, 5]


def test_iter_domain_positions_full_range_uses_catalytic_window():
    # Synthesize a sequence longer than the catalytic window so the range
    # is fully materialized.
    dna = "N" * (CATALYTIC_DOMAIN_END_AA * 3 + 9)
    positions = iter_domain_positions(dna, residue_list=None)
    expected_start = (CATALYTIC_DOMAIN_START_AA - 1) * 3
    expected_end = CATALYTIC_DOMAIN_END_AA * 3
    assert positions[0] == expected_start
    assert positions[-1] == expected_end - 1
    assert len(positions) == expected_end - expected_start


def test_iter_domain_positions_filters_out_of_range():
    # Residue 100 is past the end of a 4-aa sequence; no positions should
    # survive the in-range filter.
    assert iter_domain_positions(TEST_DNA, residue_list=[100]) == []


def test_classify_variant_detects_synonymous_missense_nonsense():
    # TGG (W) → TAG (*) is nonsense: mutate position 7 (1-based 8) G→A
    info = classify_variant(TEST_DNA, dna_index=7, alt_nt="A", protein_sequence=TEST_PROTEIN)
    assert info["effect"] == "nonsense"
    assert info["wt_aa"] == "W"
    assert info["mut_aa"] == "*"

    # AAG (K) → AAA (K) at codon 2 position 2 is synonymous
    info = classify_variant(TEST_DNA, dna_index=5, alt_nt="A", protein_sequence=TEST_PROTEIN)
    assert info["effect"] == "synonymous"

    # ATG (M) → GTG (V) is missense
    info = classify_variant(TEST_DNA, dna_index=0, alt_nt="G", protein_sequence=TEST_PROTEIN)
    assert info["effect"] == "missense"
    assert info["wt_aa"] == "M"
    assert info["mut_aa"] == "V"


def test_variant_key_is_stable():
    assert variant_key(42, "T") == "42:T"


def test_softmax_over_bases_sums_to_one_and_is_numerically_stable():
    logits = {"A": 1000.0, "C": 999.0, "G": 998.0, "T": 997.0}
    probs = softmax_over_bases(logits)
    assert math.isclose(sum(probs.values()), 1.0, abs_tol=1e-9)
    # A should dominate.
    assert probs["A"] > probs["C"] > probs["G"] > probs["T"]


def test_extract_base_logits_reads_correct_indices():
    row = [0.0] * (max(DNA_TO_LOGIT_INDEX.values()) + 1)
    row[DNA_TO_LOGIT_INDEX["A"]] = 1.5
    row[DNA_TO_LOGIT_INDEX["C"]] = -0.5
    row[DNA_TO_LOGIT_INDEX["G"]] = 2.0
    row[DNA_TO_LOGIT_INDEX["T"]] = 0.0
    result = extract_base_logits({"logits": [row]})
    assert result == {"A": 1.5, "C": -0.5, "G": 2.0, "T": 0.0}


def test_extract_base_logits_raises_on_missing_field():
    with pytest.raises(ValueError, match="no 'logits' field"):
        extract_base_logits({})


def test_extract_base_logits_raises_on_short_row():
    with pytest.raises(ValueError, match="too short"):
        extract_base_logits({"logits": [[0.0] * 10]})


def test_extract_base_logits_raises_on_non_numeric():
    row = [0.0] * (max(DNA_TO_LOGIT_INDEX.values()) + 1)
    row[DNA_TO_LOGIT_INDEX["A"]] = "not a float"  # type: ignore[assignment]
    with pytest.raises(ValueError, match="Non-numeric"):
        extract_base_logits({"logits": [row]})


def _sample_score(**overrides) -> VariantScore:
    defaults = dict(
        rank=1,
        dna_pos_1based=100,
        ref_nt="A",
        alt_nt="C",
        aa_pos_1based=34,
        wt_codon="GAT",
        mut_codon="GCT",
        wt_aa="D",
        mut_aa="A",
        effect="missense",
        ref_logit=1.23,
        alt_logit=0.45,
        ref_probability=0.6,
        alt_probability=0.3,
        delta_log_likelihood=-0.78,
        elapsed_ms=42,
    )
    defaults.update(overrides)
    return VariantScore(**defaults)


def test_csv_round_trip_preserves_values(tmp_path: Path):
    output = tmp_path / "scores.csv"
    originals = [_sample_score(), _sample_score(alt_nt="G", delta_log_likelihood=1.1)]
    append_scores(output, originals)
    reloaded = load_existing_scores(output)
    assert reloaded == originals


def test_load_existing_scores_rejects_legacy_schema(tmp_path: Path):
    # Simulate the pre-rewrite CSV that lacked the new ref_logit/alt_logit
    # columns. Loading it should fail loudly rather than silently dropping
    # fields.
    legacy = tmp_path / "legacy.csv"
    with legacy.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["rank", "dna_pos_1based", "delta_log_likelihood"])
        writer.writerow([1, 100, 0.5])
    with pytest.raises(ValueError, match="missing required columns"):
        load_existing_scores(legacy)


def test_preferred_codons_cover_all_standard_amino_acids():
    assert set(ARABIDOPSIS_PREFERRED_CODONS.keys()) == set("ACDEFGHIKLMNPQRSTVWY")


def test_dna_alphabet_and_logit_index_are_consistent():
    assert set(DNA_ALPHABET) == set(DNA_TO_LOGIT_INDEX.keys())
