# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

Synthetix ranks single-nucleotide variants in the catalytic domain (residues 111–430) of *Arabidopsis thaliana* Glutamine Synthetase 2 (GS2; UniProt Q43127, TAIR AT5G35630). Scoring is delegated to the **NVIDIA-hosted Evo 2 40B** `generate` endpoint — there is no local model inference. Top-ranked variants are then promoted to structural audits via AutoDock Vina.

Lunar food security is the long-range motivation but **not** a current technical constraint of the pipeline. The pH 8.0 in `synthetix_receptor_prep.py` simulates the alkaline chloroplast stroma (not anything lunar). If someone asks to wire lunar-specific physics in (pressure, gas exchange, cosmic-radiation damage), treat that as a real design change — don't pretend the current code already captures it.

Companion context files: `AGENTS.md` is obsolete ("no source yet" — written before the codebase existed); `GEMINI.md` overlaps the README but is pre-v0.2 (references the old sampled-probability scoring). Prefer this file and the README.

## Environment

- Python 3.11, all work runs inside the in-tree venv at `synthetix/` (git-ignored; the venv happens to share the project name, don't confuse it with a package directory).
- Activate with `source synthetix/bin/activate` before running anything.
- Install with `pip install -e '.[structure,viz3d,dev]'` for the full stack. `pyproject.toml` is the source of truth for dependencies.
- `predict_efficiency.py` requires `NVCF_RUN_KEY` exported in the shell; every other script runs offline.
- Tests: `pytest` (offline, <1s).

## Common commands

```bash
# Targeted scan with the real CDS
python3 predict_efficiency.py --cds-fasta AT5G35630.fasta --residues 124 125 127 173 187 --output results.csv

# Tune parallelism if the endpoint tolerates it (default --concurrency=8). Drop if you hit HTTP 429.
python3 predict_efficiency.py --cds-fasta AT5G35630.fasta --concurrency 4 --output results.csv

# Full catalytic-domain scan with the real CDS
python3 predict_efficiency.py --cds-fasta AT5G35630.fasta --output gs2_catalytic_domain_snp_ranking.csv

# Fallback (synthetic back-translation) — emits warning, exploratory only
python3 predict_efficiency.py --residues 125 --output scratch.csv

# Inspect payload without calling the API
python3 predict_efficiency.py --dry-run --max-positions 2

# Validate Δ-LL against published mutant kinetic data (CI gate: exits non-zero if ρ < 0.3)
python3 validate_against_kinetics.py --ranking results.csv --known known_mutants.csv --plot validation.png

# Map a homolog mutation onto GS2 numbering (homolog→GS2 alignment + local-identity sanity check)
python3 align_to_gs2.py --source EcoliGlnA.fasta --source-label EcoliGS \
    --mutations D50A E327A --measured-effects -2.1 -3.4 \
    --citation "Alibhai & Villafranca 1994" --out-csv known_mutants.csv

# Derive the Vina search box from active-site residue coordinates
python3 derive_docking_box.py --pdb AF-Q43127-F1.pdb --residues 125 --padding 0 --min-size 20

# Heatmap from a ranking CSV
python3 generate_heatmap.py --input gs2_catalytic_domain_snp_ranking.csv --output gs2_efficiency_heatmap.png

# Structural audit track
python3 generate_mutants.py                    # PDBFixer applyMutations → real side-chain repacking
python3 synthetix_receptor_prep.py             # PDBFixer: gs2_model_*.pdb → *_ready.pdb, pH 8.0
python3 synthetix_ligand_prep.py               # RDKit SMILES → *_ligand.pdb
python3 synthetix_pdbqt_conv.py                # Meeko: *.pdb → *.pdbqt
python3 docking_replicates.py --replicates 5   # ./vina × N seeds, mean ± std, pose RMSD

# Tests
pytest
```

## Architecture

Two independent tracks that both feed the top-ranked variants:

### Track 1 — Sequence-level scoring (`predict_efficiency.py`)

1. **Load the CDS.** Real DNA from `--cds-fasta` (preferred) or synthetic back-translation via `ARABIDOPSIS_PREFERRED_CODONS` (fallback, warns). `validate_translation` checks the DNA translates to the expected protein prefix.
2. **Enumerate DNA positions** across residues 111–430 (or the `--residues` subset). One API call per *position*, not per variant — the returned logit row scores all four bases at once.
3. **Score via direct logit read.** For DNA position *i*, prompt Evo 2 with `sequence=dna[:i]`, `num_tokens=1`, `enable_logits=True`. `extract_base_logits` pulls A/C/G/T entries from `logits[0]` at the `DNA_TO_LOGIT_INDEX` (ASCII-code) positions. Δ-LL = `logit(alt) − logit(ref)` — no softmax normalization needed since the partition function cancels in the difference. No wild-type baseline call. Missing/malformed API fields raise loudly.
4. **Persist incrementally.** API calls are submitted in parallel via `ThreadPoolExecutor` with `max_workers=--concurrency` (default 8). Each round waits for all futures to return, flushes the chunk to the output CSV via `append_scores`, then starts the next round. `load_existing_scores` enables resume. Final ranks are recomputed over the union and the CSV is rewritten at the end. `--batch-size` is a deprecated alias for `--concurrency`.

Key constants: `CATALYTIC_DOMAIN_START_AA`/`_END_AA` (1-based, inclusive start / inclusive end), `DNA_TO_LOGIT_INDEX` (ASCII codes 65/67/71/84 for A/C/G/T).

**Methodology reality check.** Evo 2's prior is evolutionary plausibility, not catalytic fitness. A positive Δ-LL at a conserved active-site residue (D125, E187, S173) is much more likely to indicate a miscalibrated proxy than a real rate-enhancing mutation. Always pair scoring with `validate_against_kinetics.py` before acting on ranks. If someone proposes removing that caveat from the README, push back.

**Validation-data gap.** As of April 2026 there is essentially no direct Arabidopsis GS2 point-mutant kinetic data in the literature — most GS mutagenesis is on E. coli GlnA or plant GS1. `docs/known_mutants_literature.md` catalogs the available sources with abstract-level summaries; numbers in full-text tables weren't extracted during the triage (requires paper access). To build a validation set, use `align_to_gs2.py` to map homolog mutations onto GS2 numbering and verify the local alignment identity. Discard mappings with <50% local identity.

### Track 2 — Structural audit

- `generate_mutants.py` uses PDBFixer's `applyMutations` to actually rebuild mutant side chains on `AF-Q43127-F1.pdb`, then protonates at pH 8.0. **Important**: earlier pre-v0.2 versions of this script only rewrote the residue label and left side-chain atoms untouched; any `gs2_model_*.pdb` produced before the rewrite is WT geometry in disguise and must be regenerated before use.
- `synthetix_receptor_prep.py` repairs gaps and re-protonates → `*_ready.pdb`.
- `synthetix_ligand_prep.py` embeds ATP/Glutamate from SMILES via RDKit ETKDG.
- `synthetix_pdbqt_conv.py` runs Meeko → AutoDock `.pdbqt`.
- `docking_config.txt` pins the search box (center ≈ 19.6, −15.9, 16.3; 20 Å cube) centered on Asp-125 (metal-binding Asp; Eisenberg et al. 2000). The file is annotated with the `derive_docking_box.py` invocation that reproduces its numbers; use that script to regenerate with a different residue set (e.g. the N-site Gly-251/Arg-354 cluster for ATP docking).
- `docking_replicates.py` wraps `./vina` across N seeds and reports mean ± std affinity plus pose RMSD — single-run Vina numbers are noise.

### Visualization

- `generate_heatmap.py` pivots the ranking CSV into a position × mutant-AA matrix (`aggfunc="mean"` over Δ-LL), draws with seaborn `RdBu_r` centered at 0, overlays open circles on the WT residue per column. Headless (`Agg`, `/tmp/matplotlib-cache`).
- `visualize_flux.py` reads the top 12 rows and spawns a VPython canvas with an infinite pulse loop — does not exit on its own.

## CSV schema (v0.2)

The on-disk format is locked to the `VariantScore` dataclass. Columns: `rank`, `dna_pos_1based`, `ref_nt`, `alt_nt`, `aa_pos_1based`, `wt_codon`, `mut_codon`, `wt_aa`, `mut_aa`, `effect`, `ref_logit`, `alt_logit`, `ref_probability`, `alt_probability`, `delta_log_likelihood`, `elapsed_ms`. Positions are always 1-based. Rank 1 = highest Δ-LL.

When adding columns, update the dataclass, `load_existing_scores`, and the CSV writer together — `DictWriter` uses `VariantScore.__dataclass_fields__.keys()` as the source of truth. `load_existing_scores` already fails loudly on the pre-v0.2 schema (which had `sampled_probability` / `continuation_log_likelihood` / `wildtype_continuation_log_likelihood`), so legacy CSVs must be deleted and regenerated.
