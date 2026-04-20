# Synthetix: GS2 Catalytic-Domain Variant Analysis Pipeline

Synthetix is a genomic analysis toolset for ranking and visualizing the effects of Single Nucleotide Polymorphisms (SNPs) within the catalytic domain of chloroplastic **Glutamine Synthetase 2 (GS2)**. 

Using the **NVIDIA Evo 2 40B** genomic foundation model, Synthetix predicts the evolutionary fitness and enzymatic efficiency impact of missense mutations by calculating continuation log-likelihoods.

## 🚀 Quick Start

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

## 📂 Project Structure

| File | Description |
| :--- | :--- |
| `predict_efficiency.py` | Core scoring engine using the Evo 2 API. |
| `generate_heatmap.py` | Heatmap generation using Seaborn and Matplotlib. |
| `visualize_flux.py` | Interactive 3D visualization using VPython. |
| `AGENTS.md` | Developer-specific guidelines and environment notes. |
| `GEMINI.md` | Instructional context for AI-assisted development. |

## 🧬 Scientific Context

The catalytic domain of GS2 (residues 111-430) is responsible for nitrogen assimilation in plants. Synthetix targets key functional residues including:
- **Metal Binding:** Glu-124, Asp-125, Asp-151.
- **Substrate Affinity:** Ser-173.
- **Catalytic Core:** Glu-187.

## 🛠 Tech Stack

- **Model:** NVIDIA Evo 2 40B (Genomic Foundation Model)
- **Language:** Python 3.11
- **Libraries:** Biopython, Pandas, Seaborn, Matplotlib, VPython, Requests

## 📄 License
*Specify license here (e.g., MIT, Apache 2.0)*
