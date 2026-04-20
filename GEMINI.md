# Project Synthetix

Synthetix is a scientific analysis pipeline designed to rank and visualize the effects of Single Nucleotide Polymorphisms (SNPs) within the catalytic domain of chloroplastic Glutamine Synthetase 2 (GS2). It leverages the **NVIDIA Evo 2 40B** genomic foundation model to predict the evolutionary fitness and efficiency impact of missense mutations.

## Architecture & Workflow

The project follows a Research -> Score -> Visualize workflow:
1.  **Prediction:** `predict_efficiency.py` performs back-translation of the GS2 protein sequence (Arabidopsis preference) and queries the NVIDIA hosted Evo 2 API to obtain continuation log-likelihoods for variants.
2.  **Analysis:** Results are saved in a structured CSV (`gs2_catalytic_domain_snp_ranking.csv`) including rank, delta log-likelihood, and sampled probabilities.
3.  **Visualization:**
    *   **Heatmaps:** `generate_heatmap.py` produces a mutation landscape heatmap (`gs2_efficiency_heatmap.png`) using Seaborn/Matplotlib.
    *   **3D Scene:** `visualize_flux.py` creates a VPython-based interactive visualization of the top-ranked variants.

## Development Setup

### Environment
The project uses a Python 3.11 virtual environment located in the `synthetix/` directory (ignored by Git).

```bash
# Activate the environment
source synthetix/bin/activate

# Verify setup
python --version
pip list
```

### Key Commands

- **Score Variants:**
  ```bash
  python3 predict_efficiency.py --residues <list of residues> --output results.csv
  ```
  *Requires `NVCF_RUN_KEY` environment variable for API access.*

- **Generate Heatmap:**
  ```bash
  python3 generate_heatmap.py --input results.csv --output heatmap.png
  ```

- **Run 3D Visualization:**
  ```bash
  python3 visualize_flux.py
  ```

## Project Structure

- `predict_efficiency.py`: Main scoring script using Evo 2.
- `generate_heatmap.py`: Heatmap generation tool.
- `visualize_flux.py`: VPython visualization script.
- `AGENTS.md`: Detailed developer guidelines and environment management notes.
- `gs2_catalytic_domain_snp_ranking.csv`: Default output/input for scoring data.
- `.gitignore`: Configured to ignore the `synthetix/` venv and model/cache artifacts.

## Conventions

- **Python Version:** 3.11+
- **Style:** PEP 8 with 4-space indentation.
- **Dependencies:** Core dependencies include `requests`, `biopython`, `pandas`, `seaborn`, `matplotlib`, and `vpython`.
