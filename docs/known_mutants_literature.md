# Known GS mutants from the literature

This document collects published GS mutagenesis data that can feed
`validate_against_kinetics.py`. It exists because the vast majority of the
experimental literature is on bacterial GS (*E. coli*, *Salmonella*) or plant
GS1 (cytosolic) — **not** on *Arabidopsis* GS2 (chloroplastic) directly, which
is what Synthetix is scoring. To use any of these entries for validation you
must first map the residue number from the source homolog onto the GS2
numbering with `align_to_gs2.py` and then manually review the mapping.

## Honesty caveat

This review was conducted from article abstracts and PubMed/PMC summaries.
Several papers report quantitative kinetic tables only in the full text,
which this reviewer did not have license to reproduce. **Do not paste rows
into `known_mutants.csv` without opening the primary source and verifying the
number** — an abstract-only citation chain is not durable enough to gate a
scoring decision.

The entries below are organized by source organism and flag how confidently
a number can be extracted without the full text.

---

## Direct *Arabidopsis* GS2 data

**None found** in a targeted literature search (April 2026). Published work
on *Arabidopsis* GS2 / GLN2 focuses on knockout phenotypes (`gln2-1`,
`gln2-2`; Igarashi et al. 2022; Potel et al. 2019), dual organelle targeting
(Taira et al. 2004), and activation by ACR11 (Osanai et al. 2017). None of
these reports deposit point-mutant kinetic tables suitable for Δ-LL
validation.

This is the central gap that makes validation hard.

---

## *Arabidopsis* GS1 (cytosolic; homologous but not the Synthetix target)

### Ishiyama, K. et al. (2006) *J. Biol. Chem.* 281, 29287–29296. PMID 16338958.

Mutagenesis of GLN1;3 (low-affinity cytosolic isoform).

| Source mutation | Reported effect (qualitative; abstract only) |
|-----------------|-----------------------------------------------|
| K49Q            | Decreased Km for ammonium; increased catalytic efficiency |
| A174S           | Decreased Km for ammonium; increased catalytic efficiency |

The abstract states the effect direction; the magnitudes live in full-text
tables. Do not record quantitative numbers here until someone extracts them
from the primary source.

Residue mapping: GLN1;3 has no chloroplast transit peptide, so its
residue 49 is near the start of the mature enzyme. GS2's mature enzyme
begins after a ~58-residue transit peptide (position ≈59 of Q43127),
with the catalytic domain conventionally starting at residue 111. Run
`align_to_gs2.py` on the GLN1;3 sequence (UniProt Q9LVI3) to resolve the
GS2 equivalents of K49 and A174.

---

## *Phaseolus vulgaris* GS1 α-polypeptide

### Clemente M.T. & Márquez A.J. (2000?) *Plant Mol. Biol.* (doi:10.1023/A:1006257323624)

"Site-directed mutagenesis of Glu-297 from the α-polypeptide of Phaseolus
vulgaris glutamine synthetase alters kinetic and structural properties and
confers resistance to L-methionine sulfoximine."

| Source mutation | Reported effect |
|-----------------|-----------------|
| E297A           | **Km for ammonium (biosynthetic) increased 100-fold**; other parameters "not greatly altered"; complements *E. coli glnA* lesion |

The 100-fold Km_NH4 shift is the single quantitative number extractable
from the abstract. On a log₂ scale it would be `log2(100) ≈ 6.6` for the
Km increase (i.e. catalytic efficiency loss of roughly `-6.6` in the
`log2_relative_kcat_over_km_NH4` convention).

Residue mapping: the α-polypeptide is GS1. Phaseolus GS1 residue 297 is
homologous to *E. coli* GlnA residue 327 (authors note this explicitly).
Run `align_to_gs2.py` against GS2 to get the Q43127 equivalent.

---

## *E. coli* GlnA (bacterial GS)

### Alibhai, M. & Villafranca, J.J. (1994) *Biochemistry* 33, 682–686. PMID 7904829.

Mutagenesis of two key active-site residues.

| Source mutation | Reported effect (qualitative; abstract only) |
|-----------------|-----------------------------------------------|
| D50A            | Increased Km for NH4+; destabilizes both ground and transition states for phosphoryl transfer |
| D50E            | Active with Mn²⁺ but very low activity with Mg²⁺ (physiological metal); all three kcat/Km substantially altered |
| E327A           | Decreased kcat/Km for NH4+ |

Again: numbers are in the full text, not the abstract. Do not record
quantitative rows without pulling the primary source.

Residue mapping: *E. coli* GlnA (UniProt P0A9C5) has no transit peptide.
D50 and E327 are canonical active-site residues conserved across all
GSI-β enzymes; both map into the GS2 catalytic domain. Use
`align_to_gs2.py`.

---

## *Providencia vermicola* GS

### Mishra, S.K. et al. (2018) *Sci. Rep.* 8, 15549. PMID 30353099 / PMC 6199252.

**Only paper in this review with a full quantitative kinetic table accessible
without paywall.**

Wild-type vs S54A (from Table 1 of the paper):

| Substrate       | WT Km (mM)          | WT kcat (s⁻¹)   | S54A Km (mM)        | S54A kcat (s⁻¹) |
|-----------------|---------------------|-----------------|---------------------|-----------------|
| Hydroxylamine   | 15.7 ± 1.1          | 17.0 ± 0.6      | 4.6 ± 0.3           | 19.6 ± 0.5      |
| ADP-Na₂         | (25.2 ± 1.5)×10⁻⁵   | 9.14 ± 0.12     | (10.2 ± 5.6)×10⁻⁶   | 15.2 ± 0.2      |
| L-Glutamine     | 32.6 ± 1.7          | 30.5 ± 1.0      | 23.8 ± 1.2          | 41.9 ± 1.1      |

- S54A kcat/Km for hydroxylamine: **3.1× WT** (transferase activity).
- S54A kcat/Km for glutamine: **2.9× WT** (glutaminase activity).
- Also improved thermostability: residual activity +10.7× at 0 °C, +3.8× at 10 and 50 °C vs WT.

Caveat: the paper measures the γ-glutamyl transferase and glutaminase
activities, not the biosynthetic (forward) reaction most relevant to in vivo
nitrogen assimilation. Using these kcat/Km ratios as a stand-in for
biosynthetic improvement is a non-trivial assumption that should be flagged
in `notes`.

Residue mapping: PveGS S54 is homologous to *E. coli* GlnA S53 (Mishra
et al. provide an alignment in their Fig. 1). From there, `align_to_gs2.py`
will place it on GS2.

---

## *Salmonella typhimurium* GS (the Eisenberg-lab structural mutants)

### Liaw, S.-H. & Eisenberg, D. (1994) *Biochemistry* 33, 675–681.
### Pinkofsky, H.B. et al. (1993) *Biochemistry* 32, 10853–10858. PMID 8102250.

These are the crystallographic mechanism papers that define the active-site
architecture GS2 inherits. They describe structures of enzyme-substrate
complexes and the H269N mutant (oxidative modification study). The H269N
paper does not tabulate kcat/Km in the abstract — again, full text required.

### Eisenberg, D. et al. (2000) *Biochim. Biophys. Acta* 1477, 122–145.

Review of GS structure-function; recommended as the source for the
**residue identities of the metal-binding site (n1/n2) and glutamate /
ammonium subsites**, which are what `derive_docking_box.py` points its
active-site residue list at.

---

## *Zea mays* GS1a (plant crystal structure)

### Unno, H. et al. (2006) *J. Biol. Chem.* 281, 29287–29296. PMID 16829528.

First plant GS crystal structure (PDB 2D3A/2D3B/2D3C). Not a mutagenesis
paper per se, but establishes the plant-specific active-site architecture
and the 10-mer oligomerization. Use the structure to sanity-check the
active-site residue numbers in `derive_docking_box.py` for plant GS —
they're taken from Eisenberg et al. (bacterial-numbered) and mapped
across a conserved alignment.

---

## Workflow for populating `known_mutants.csv`

1. Open the primary source (not the abstract).
2. Extract the source mutation label + WT-normalized quantitative effect (`log2(kcat/kcat_WT)` or `-ΔΔG` — pick one convention and hold it constant).
3. Obtain the homolog's FASTA (UniProt: `P0A9C5` for *E. coli*, `Q9LVI3` for *Arabidopsis* GLN1;3, etc.).
4. Run `align_to_gs2.py --source homolog.fasta --source-label <tag> --mutations <list> --measured-effects <list> --citation <ref>`.
5. Inspect the printed local-identity column. If identity < ~50% or the source WT doesn't match GS2 WT, treat the mapping as unreliable.
6. Append the surviving rows to `known_mutants.csv`.

A correlation computed on <5 high-confidence rows is noise-dominated; aim
for ≥10 before drawing conclusions about the Δ-LL ranking's calibration.
