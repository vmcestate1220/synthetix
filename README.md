<p align="center">
  <img src="GS2_AlphaFold_Structure.png" alt="AlphaFold structure of chloroplastic Glutamine Synthetase 2 (GS2) from Arabidopsis thaliana (UniProt Q43127), colored by pLDDT" width="720"/>
  <br/>
  <sub><em>AlphaFold structure of the GS2 catalytic domain (residues 111–430; <a href="https://alphafold.ebi.ac.uk/entry/Q43127">UniProt Q43127</a>).<br/>Colored by pLDDT: dark blue ≥90, light blue 70–90.</em></sub>
</p>

# Synthetix: GS2 Catalytic-Domain Variant Analysis Pipeline

[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/release/python-3110/)
[![Model: NVIDIA Evo 2 40B](https://img.shields.io/badge/model-NVIDIA%20Evo%202%2040B-76B900.svg)](https://build.nvidia.com/arc/evo2-40b)
[![Biopython](https://img.shields.io/badge/Biopython-1.8x-informational.svg)](https://biopython.org/)
[![Seaborn](https://img.shields.io/badge/Seaborn-plots-9cf.svg)](https://seaborn.pydata.org/)
[![Matplotlib](https://img.shields.io/badge/Matplotlib-visualization-11557c.svg)](https://matplotlib.org/)
[![VPython](https://img.shields.io/badge/VPython-3D-orange.svg)](https://vpython.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](#license)
[![Status](https://img.shields.io/badge/status-active%20research-success.svg)](#)

---

## Contents

- [Overview](#overview)
- [Setup](#setup)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Dataset](#dataset)
- [Visualization](#visualization)
- [References](#references)
- [Author](#author)

---

## Overview

Synthetix is a research pipeline for ranking and visualizing the functional impact of Single Nucleotide Polymorphisms (SNPs) within the catalytic domain of chloroplastic **Glutamine Synthetase 2 (GS2)** from *Arabidopsis thaliana* (UniProt **Q43127**). Its goal is to identify non-synonymous substitutions that may plausibly **improve catalytic efficiency, substrate affinity, or fold stability** — and to narrow a combinatorially large mutation landscape down to a prioritized shortlist of candidates for structural audit and wet-lab validation.

GS2 catalyzes the ATP-dependent condensation of glutamate and ammonia into glutamine, the entry point for nitrogen assimilation in plant chloroplasts. Even small improvements in its catalytic throughput have implications for **nitrogen-use efficiency** in crops, and the enzyme's active site is densely constrained by metal-binding, substrate-anchoring, and transition-state residues — making it a rigorous testbed for computational variant-effect prediction.

### Role of Evo 2

Synthetix uses **NVIDIA Evo 2 40B** — an autoregressive genomic foundation model trained on ~9.3 trillion nucleotides spanning all domains of life — as the scoring engine behind its variant predictions. Evo 2 models DNA directly at single-nucleotide resolution, which lets Synthetix evaluate SNPs *in their native codon and genomic context* rather than after translation to protein. This preserves signals (codon usage bias, local sequence composition, mutational neighborhood) that are lost in protein-only language models.

Concretely, the pipeline:

1. Back-translates the wild-type GS2 protein into a DNA sequence using *Arabidopsis*-preferred codons.
2. Generates each candidate single-nucleotide variant across the catalytic domain (residues 111–430).
3. Prompts Evo 2 with the sequence context preceding the variant position and reads out the model's **continuation log-likelihood** for the reference vs. alternate base.
4. Ranks variants by the **delta log-likelihood** (Δ-LL) relative to wild type: positive Δ-LL ≈ the variant fits the evolutionary/functional prior learned by Evo 2 better than the WT base; negative Δ-LL ≈ the opposite.

Because inference runs against the hosted NVIDIA API, no local GPU is required, and the full catalytic-domain scan is reproducible from a single command.

### Project Aims

1. **Variant discovery.** Produce a ranked, reproducible list of catalytic-domain SNPs most likely to modulate GS2 function.
2. **Active-site focus.** Targeted scans over functionally critical residues — metal-binding (Glu-124, Asp-125, Asp-151), substrate affinity (Ser-173), and the catalytic core (Glu-187) — to triage leads in the regions where fitness effects are most interpretable.
3. **Structural validation.** Generate mutant PDB models (e.g., `gs2_model_D125N.pdb`) and run docking audits against ATP and glutamate to cross-check the sequence-level predictions with geometry-level plausibility.
4. **Communicable results.** Deliver publication-quality heatmaps and 3D visualizations so that leads can be inspected, compared, and handed off for experimental follow-up.

### Tech Stack

- **Model:** NVIDIA Evo 2 40B (hosted genomic foundation model)
- **Language:** Python 3.11
- **Libraries:** Biopython, Pandas, Seaborn, Matplotlib, VPython, Requests
- **Docking:** AutoDock Vina with `.pdbqt`-converted receptor and ligand inputs

---

## Setup

The project requires Python 3.11. Activate the local virtual environment and install dependencies:

```bash
source synthetix/bin/activate
pip install -r requirements.txt  # If requirements.txt is provided
```

Set your NVIDIA API key as an environment variable:

```bash
export NVCF_RUN_KEY="your_nvidia_api_key_here"
```

---

## Quick Start

### 1. Environment Setup
The project requires Python 3.11. Activate the local virtual environment:

```bash
source synthetix/bin/activate
pip install -r requirements.txt  # If requirements.txt is provided
```

### 2. API Configuration
Set your NVIDIA API key as an environment variable:

```bash
export NVCF_RUN_KEY="your_nvidia_api_key_here"
```

### 3. Run a Scan
Predict the efficiency of specific residues (e.g., active site residues 124, 125, 127):

```bash
python3 predict_efficiency.py --residues 124 125 127 --output results.csv
```

### 4. Visualize Results
Generate a publication-quality heatmap of the mutation landscape:

```bash
python3 generate_heatmap.py --input results.csv --output gs2_heatmap.png
```

Or launch an interactive 3D visualization:

```bash
python3 visualize_flux.py
```

---

## Usage

### Project Structure

| File | Description |
| :--- | :--- |
| `predict_efficiency.py` | Core scoring engine using the Evo 2 API. |
| `generate_heatmap.py` | Heatmap generation using Seaborn and Matplotlib. |
| `generate_mutants.py` | Utility to construct SNP-perturbed sequence sets. |
| `visualize_flux.py` | Interactive 3D visualization using VPython. |
| `synthetix_ligand_prep.py` | Prepares ligand geometries for downstream docking. |
| `synthetix_receptor_prep.py` | Prepares GS2 receptor structures for docking. |
| `synthetix_pdbqt_conv.py` | Converts PDB structures to AutoDock `.pdbqt` format. |
| `veron_baseline_audit.py` | Baseline audit utility for the Veron variant sub-project. |
| `AGENTS.md` | Developer-specific guidelines and environment notes. |
| `GEMINI.md` | Instructional context for AI-assisted development. |

### Scoring Variants

`predict_efficiency.py` is the main entry point. It back-translates the GS2 protein sequence using *Arabidopsis*-preferred codons, generates all single-nucleotide variants within the catalytic domain (or an explicit residue list), and scores each one against the Evo 2 40B `generate` endpoint. Results are ranked by delta log-likelihood relative to the wild-type continuation.

Key flags:

- `--residues` — limit the scan to specific 1-based residues
- `--include-synonymous` — include synonymous SNPs (excluded by default)
- `--max-variants` — cap the number of variants evaluated
- `--batch-size` — API batching granularity
- `--dry-run` — print the request payload without calling the API
- `--output` — CSV path for ranked results

### Generating Mutant Structures

`generate_mutants.py` produces perturbed PDB models (e.g., `gs2_model_D125N.pdb`, `gs2_model_N112Y.pdb`) used for docking audits and figure generation.

### Docking Preparation

Use `synthetix_receptor_prep.py` and `synthetix_ligand_prep.py` to prepare AutoDock-compatible inputs. `synthetix_pdbqt_conv.py` converts between formats; `docking_config.txt` provides the reference configuration.

---

## Dataset

Ranked variant scores are serialized to CSV and can be regenerated with `predict_efficiency.py`. Columns match the `VariantScore` dataclass in `predict_efficiency.py:66`.

| Artifact | Description |
| :--- | :--- |
| `gs2_catalytic_domain_snp_ranking.csv` | Full catalytic-domain (residues 111–430) SNP ranking with delta log-likelihoods. |
| `gs2_active_site_scan.csv` | Targeted scan over the active-site residues. |
| `test_api.csv` | Small fixture used for API smoke tests. |
| `AF-Q43127-F1.pdb` / `AF-GS2.pdb` | AlphaFold GS2 reference structures. |
| `gs2_wildtype.pdb` | Wild-type GS2 model. |
| `gs2_model_{D125N,E170A,N112Y,S117C,W111G}.pdb` | Top-ranked mutant structural models. |
| `ATP.pdbqt`, `Glutamate.pdbqt` | Prepared ligands for docking. |

Key score columns include `dna_pos_1based`, `aa_pos_1based`, `wt_aa`, `mut_aa`, `effect`, `sampled_probability`, `continuation_log_likelihood`, `wildtype_continuation_log_likelihood`, and `delta_log_likelihood`. Variants are sorted so that rank 1 is the most favorable (highest delta) relative to wild type.

---

## Visualization

Several rendering paths are available:

- **Variant-effect heatmaps** — `generate_heatmap.py` builds a position × amino-acid matrix of mean delta log-likelihoods and overlays the wild-type residue per column. Produces `gs2_efficiency_heatmap.png`, `gs2_active_site_heatmap.png`, and `gs2_full_landscape_heatmap.png`.
- **Interactive 3D flux view** — `visualize_flux.py` launches a VPython scene for inspecting substrate flow / residue geometry.
- **Structural audits** — static renderings such as `D125N_Audit_Synthetix.png`, `N112Y_Audit.png`, and `Synthetix_Structural_Portfolio.png` summarize docking and mutational outcomes for selected variants.

Example:

```bash
python3 generate_heatmap.py \
  --input gs2_catalytic_domain_snp_ranking.csv \
  --output gs2_efficiency_heatmap.png
```

---

## References

- Nguyen, E. *et al.* **Evo 2: Genome modeling and design across all domains of life.** Arc Institute / NVIDIA, 2024. <https://build.nvidia.com/arc/evo2-40b>
- Jumper, J. *et al.* **Highly accurate protein structure prediction with AlphaFold.** *Nature* 596, 583–589 (2021).
- Eisenberg, D., Gill, H. S., Pfluegl, G. M. U. & Rotstein, S. H. **Structure–function relationships of glutamine synthetases.** *Biochim. Biophys. Acta* 1477, 122–145 (2000).
- Trott, O. & Olson, A. J. **AutoDock Vina: Improving the speed and accuracy of docking with a new scoring function.** *J. Comput. Chem.* 31, 455–461 (2010).
- Cock, P. J. A. *et al.* **Biopython: freely available Python tools for computational molecular biology and bioinformatics.** *Bioinformatics* 25, 1422–1423 (2009).
- UniProt entry **Q43127** — Glutamine synthetase leaf isozyme, chloroplastic (*Arabidopsis thaliana*). <https://www.uniprot.org/uniprotkb/Q43127>

---

## Author

**Christopher Cocchiaraley**
Maintainer of the Synthetix pipeline.

- GitHub: [@vmcestate1220](https://github.com/vmcestate1220)

### License

Released under the **MIT License**. See [`LICENSE`](LICENSE) for the full text.
