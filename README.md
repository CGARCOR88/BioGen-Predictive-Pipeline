# 🧬 BioGen Predictive Pipeline

![Tests](https://github.com/tu_usuario/BioGen-Predictive-Pipeline/actions/workflows/tests.yml/badge.svg)
![Python](https://img.shields.io/badge/python-3.11-blue)
![License](https://img.shields.io/badge/license-MIT-green)

Bioinformatics pipeline for genomic analysis of FASTA sequences: DNA/protein metrics calculation, NCBI homology search, statistical modeling with linear regression, and anomaly detection.

---

## Description

EThis project analyzes human DNA sequences (FASTA format) by applying a 4-phase pipeline:

1. **Genomic Processing** — GC%, length, DNA→protein translation, molecular weight, aromaticity.
2. **Online Homology** (optional) — BLAST search against the NCBI `nt` database and reference downloading.
3. **AI Model** — Linear regression (length → GC%) with anomaly detection (>2σ).
4. **Visualization** — Correlation heatmap generated with seaborn.

The included sample genes are: **TP53**, **BRCA1**, **APOE**, and **KCNJ1**.

---

## Project Structure

```
text
BioGen-Predictive-Pipeline/
├── data/
│   ├── raw/                  # Input FASTAs (.fasta / .fa)
│   ├── processed/            # features_genes.csv (generated)
│   └── reference/            # References downloaded from NCBI (online mode)
├── graficos/                 # heatmap_correlacion.png (generated)
├── logs/                     # pipeline.log (generated)
├── resultados/               # Results CSV and PNG (generated)
├── src/
│   ├── __init__.py           # CLI entry point
│   ├── descarga_sec.py       # Bulk download of FASTAs from NCBI
│   └── fetch_tools.py        # Reference FASTA downloading from NCBI
│   ├── homology_search.py    # BLAST + GenBank details extraction
│   ├── model_ia.py           # Linear regression + anomaly detector
│   ├── procesado_sec.py      # DNA/Protein metrics orchestrator
│   └── visualizer.py         # Correlation heatmap with seaborn
├── tests/
│   └── test_pipeline.py      # 29 unit tests (pytest)
├── .env.example              # Environment variables template
├── requirements.txt
├── main.py                   # Orchestrator: DNA + protein metrics
└── setup.py                  # Folder initialization script

```

---

## Installation

```bash
git clone [https://github.com/tu_usuario/BioGen-Predictive-Pipeline.git](https://github.com/tu_usuario/BioGen-Predictive-Pipeline.git)
cd BioGen-Predictive-Pipeline

python -m venv venv
# Windows:
venv\Scripts\activate
# Linux/Mac:
source venv/bin/activate

pip install -r requirements.txt
```

### Configure NCBI Credentials
Copy `.env.example` to `.env` and add your email:

```bash
cp .env.example .env
```

```env
ENTREZ_EMAIL=tu_email@ejemplo.com
```

> An email address is required for NCBI API calls. Without it, the `--online` mode will fail.

---

## Usage

### Offline Mode (Local analysis, no NCBI)

```bash
python scripts/main_pipeline.py
```

### Online Mode (BLAST + reference download)

```bash
python scripts/main_pipeline.py --online
```

---

## Tests

```bash
pytest tests/ -v
```
Expected output: **29 passed**.

---

## Tech Stack

| Library | Usage |
|---|---|
| BioPython ≥1.81 | FASTA parsing, BLAST, Entrez, ProteinAnalysis |
| pandas ≥2.0 | DataFrames and CSV export |
| scikit-learn ≥1.3 | Linear regression |
| seaborn / matplotlib | Correlation heatmap |
| python-dotenv | Secure credentials management |
| pytest | Unit testing |

---

## Dataset

The sequences included in `data/raw/` are public RefSeq records from NCBI:

| Gene | Accession | Function |
|---|---|---|
| TP53 | NM_000546 | Supresor tumoral |
| BRCA1 | NM_007294 | Reparación de ADN |
| APOE | NM_000041 | Metabolismo lipídico |
| KCNJ1 | NM_001301717 | Canal de potasio |

---

## Author
Carlos Garcia Corona
