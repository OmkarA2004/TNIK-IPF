# TNIK INSIGHT
### Machine Learning-Guided TNIK Inhibitor Discovery for Idiopathic Pulmonary Fibrosis (IPF)

**TNIK INSIGHT** is a computational research project focused on exploring potential TNIK inhibitors for Idiopathic Pulmonary Fibrosis (IPF) using machine learning, cheminformatics, molecular docking, molecular dynamics, and explainable AI.

## Project Overview

Transcriptional intermediary factor 1 (TNIK) is investigated as a potential therapeutic target in IPF. This project combines computational drug-discovery approaches to investigate compounds and predict their potential TNIK inhibitory activity.

The project aims to:

- Develop machine-learning models for predicting TNIK inhibitory activity.
- Analyse chemical compounds using cheminformatics.
- Investigate potential TNIK-binding compounds through molecular docking.
- Explore protein–ligand stability using molecular dynamics simulations.
- Interpret model predictions using explainable AI.
- Develop a web-based TNIK activity prediction platform.

## Computational Workflow

1. **Literature Review** – Research on TNIK, IPF, and potential therapeutic strategies.
2. **Data Collection** – Collection and organisation of chemical and biological datasets.
3. **Data Processing** – Preparation of compounds and molecular descriptors.
4. **Machine Learning** – Development and evaluation of predictive models.
5. **Molecular Docking** – Computational investigation of protein–ligand interactions.
6. **Molecular Dynamics** – Investigation of selected protein–ligand systems.
7. **Explainable AI** – Interpretation of machine-learning predictions.
8. **Analogue Exploration** – Investigation of structurally related compounds.
9. **Web Application** – Development of the TNIK INSIGHT prediction interface.

## Machine Learning Performance

The Random Forest baseline model was evaluated using a held-out test set.

| Evaluation metric | Test result |
|---|---:|
| R² | 0.4914 |
| RMSE | 0.7976 pIC50 |
| MAE | 0.5917 pIC50 |
| Median absolute error | 0.3857 |
| Maximum absolute error | 2.3396 |

**Dataset split:**
- Total compounds: 453
- Training compounds: 362 (79.9%)
- Test compounds: 91 (20.1%)

These values represent the reported evaluation results for the Random Forest model. They should not be interpreted as experimental confirmation of inhibitor activity or clinical efficacy.

## Repository Structure

```text
TNIK-INSIGHT/
│
├── 01_literature/
├── 02_raw_data/
├── 03_processed_data/
├── 04_code/
├── 05_models/
├── 06_docking/
├── 07_explainability/
├── 08_analogues/
├── 09_results/
├── 10_report/
├── 11_archive/
├── 12_project_documentation/
│
└── 12_TNIK_Predict_Website/
```

## Web Application

The project includes a web application in:

`12_TNIK_Predict_Website/`

It is intended to provide an interface for TNIK activity prediction using the developed computational models.

## Technologies and Methods

The project involves computational methods and tools across:

- Python
- Machine Learning
- Cheminformatics
- Molecular Descriptors
- Random Forest
- Molecular Docking
- Molecular Dynamics
- Explainable Artificial Intelligence
- Web Application Development

## Research Significance

TNIK INSIGHT explores an integrated computational workflow for investigating potential TNIK inhibitors in the context of IPF.

The results are computational predictions and require further experimental validation before any therapeutic conclusions can be made.

## Important Note

Large molecular-dynamics trajectory files and selected large log files have been excluded from this GitHub repository to keep the repository manageable. Original files are retained separately in the local project directory.

## Disclaimer

This repository is intended for academic and research purposes. Computational predictions do not establish biological activity, safety, or clinical efficacy.

## Author

**Omkar Aravattu**

GitHub: [@OmkarA2004](https://github.com/OmkarA2004)

---

*TNIK INSIGHT – Computational exploration of potential TNIK inhibitors for IPF.*
