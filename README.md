# 🧬 BioGen Predictive Pipeline

![Tests](https://github.com/tu_usuario/BioGen-Predictive-Pipeline/actions/workflows/tests.yml/badge.svg)
![Python](https://img.shields.io/badge/python-3.11-blue)
![License](https://img.shields.io/badge/license-MIT-green)

Pipeline bioinformático para el análisis genómico de secuencias FASTA: cálculo de métricas de ADN/proteína, búsqueda de homología en NCBI, modelado estadístico con regresión lineal y detección de anomalías.

---

## Descripción

Este proyecto analiza secuencias de ADN humano (formato FASTA) aplicando un pipeline de 4 fases:

1. **Procesamiento genómico** — GC%, longitud, traducción ADN→proteína, peso molecular, aromaticidad
2. **Homología online** (opcional) — BLAST contra la base de datos `nt` del NCBI y descarga de referencias
3. **Modelo de IA** — Regresión lineal (longitud → GC%) con detección de anomalías (>2σ)
4. **Visualización** — Heatmap de correlación generado con seaborn

Los genes de ejemplo incluidos son: **TP53**, **BRCA1**, **APOE** y **KCNJ1**.

---

## Estructura del proyecto

```
BioGen-Predictive-Pipeline/
├── data/
│   ├── raw/                  # FASTAs de entrada (.fasta / .fa)
│   ├── processed/            # features_genes.csv (generado)
│   └── reference/            # Referencias descargadas de NCBI (modo online)
├── graficos/                 # heatmap_correlacion.png (generado)
├── logs/                     # pipeline.log (generado)
├── resultados/               # CSV y PNG de resultados (generados)
├── src/
│   ├── __init__.py           # Punto de entrada CLI
│   ├── descarga_sec.py       # Descarga masiva de FASTAs desde NCBI
│   └── fetch_tools.py        # Descarga de referencias FASTA desde NCBI
│   ├── homology_search.py    # BLAST + extracción de detalles GenBank
│   ├── model_ia.py           # Regresión lineal + detector de anomalías
│   ├── procesado_sec.py      # Orquestador de métricas ADN/Proteína
│   └── visualizer.py         # Heatmap de correlación con seaborn
├── tests/
│   └── test_pipeline.py      # 29 tests unitarios (pytest)
├── .env.example              # Plantilla de variables de entorno
├── requirements.txt
├── main.py                   # Orquestador: métricas ADN + proteína
└── setup.py                  # Script de inicialización de carpetas

```

---

## Instalación

```bash
git clone https://github.com/tu_usuario/BioGen-Predictive-Pipeline.git
cd BioGen-Predictive-Pipeline

python -m venv venv
# Windows:
venv\Scripts\activate
# Linux/Mac:
source venv/bin/activate

pip install -r requirements.txt
```

### Configurar credenciales NCBI

Copia `.env.example` a `.env` y añade tu email:

```bash
cp .env.example .env
```

```env
ENTREZ_EMAIL=tu_email@ejemplo.com
```

> El email es obligatorio para las llamadas a la API de NCBI. Sin él, el modo `--online` fallará.

---

## Uso

### Modo offline (análisis local, sin NCBI)

```bash
python scripts/main_pipeline.py
```

### Modo online (BLAST + descarga de referencias)

```bash
python scripts/main_pipeline.py --online
```

---

## Tests

```bash
pytest tests/ -v
```

Resultado esperado: **29 passed**.

---

## Tecnologías

| Librería | Uso |
|---|---|
| BioPython ≥1.81 | Parsing FASTA, BLAST, Entrez, ProteinAnalysis |
| pandas ≥2.0 | DataFrames y exportación CSV |
| scikit-learn ≥1.3 | Regresión lineal |
| seaborn / matplotlib | Heatmap de correlación |
| python-dotenv | Gestión segura de credenciales |
| pytest | Tests unitarios |

---

## Dataset

Las secuencias incluidas en `data/raw/` son registros RefSeq públicos del NCBI:

| Gen | Accesión | Función |
|---|---|---|
| TP53 | NM_000546 | Supresor tumoral |
| BRCA1 | NM_007294 | Reparación de ADN |
| APOE | NM_000041 | Metabolismo lipídico |
| KCNJ1 | NM_001301717 | Canal de potasio |

---

## Autor
Carlos Garcia Corona
