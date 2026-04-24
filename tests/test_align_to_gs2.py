"""Unit tests for align_to_gs2 — exercise the mapping + identity logic on
small synthetic alignments so the behaviour at gaps is specified."""

from __future__ import annotations

from dataclasses import dataclass

import pytest

from align_to_gs2 import (
    build_source_to_target_map,
    local_identity,
    parse_mutation,
)


@dataclass
class FakeAlignment:
    """Minimal stand-in for Bio.Align.Alignment for unit testing.

    PairwiseAligner returns an Alignment object whose ``str(aln[0])`` and
    ``str(aln[1])`` give the aligned source/target strings. We mimic that
    with a subscript-able tuple.
    """
    aligned_source: str
    aligned_target: str

    def __getitem__(self, i):
        return [self.aligned_source, self.aligned_target][i]


def test_parse_mutation():
    assert parse_mutation("D125N") == ("D", 125, "N")
    assert parse_mutation(" e187a ") == ("E", 187, "A")


def test_parse_mutation_rejects_bad_format():
    with pytest.raises(ValueError):
        parse_mutation("D125")
    with pytest.raises(ValueError):
        parse_mutation("125N")
    with pytest.raises(ValueError):
        parse_mutation("XYZ")


def test_mapping_identity_case():
    aln = FakeAlignment("ABCDE", "ABCDE")
    assert build_source_to_target_map(aln) == {1: 1, 2: 2, 3: 3, 4: 4, 5: 5}


def test_mapping_with_gaps_in_target():
    # Source has 5 residues, target has 3 — positions 2 and 4 of source
    # don't have a target equivalent.
    aln = FakeAlignment("ABCDE", "A-C-E")
    assert build_source_to_target_map(aln) == {1: 1, 3: 2, 5: 3}


def test_mapping_with_gaps_in_source():
    aln = FakeAlignment("A-C-E", "ABCDE")
    assert build_source_to_target_map(aln) == {1: 1, 2: 3, 3: 5}


def test_local_identity_perfect():
    aln = FakeAlignment("MAQILAAS", "MAQILAAS")
    # Any source position should report 100% identity.
    assert local_identity(aln, source_residue_1based=3, window=2) == pytest.approx(1.0)


def test_local_identity_partial():
    # Window of ±1 around position 3: positions 2,3,4 → source QIL vs target RIL
    # identity = 2/3
    aln = FakeAlignment("MQILAAS", "MRILAAS")
    assert local_identity(aln, source_residue_1based=3, window=1) == pytest.approx(2 / 3)


def test_local_identity_ignores_gap_columns():
    # Source position 3 is 'C'; target aligned column at that row is '-'.
    # That column is not counted as "paired" so it shouldn't drag identity.
    aln = FakeAlignment("AB CDE".replace(" ", ""), "AB-DE")
    aln = FakeAlignment("ABCDE", "AB-DE")
    # Around position 3 with window 1: columns for source-idx 2,3,4.
    # Column 2 = B/B match, column 3 = C/- gap (not paired), column 4 = D/D match.
    # paired = 2 (B,D), matched = 2 → identity 100%.
    assert local_identity(aln, source_residue_1based=3, window=1) == pytest.approx(1.0)
