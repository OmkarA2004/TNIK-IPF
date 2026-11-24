# Materials and Methods

## 1. Study design and computational workflow

A computational drug-repurposing workflow was developed to investigate potential interactions between approved drugs and TRAF2- and NCK-interacting kinase (TNIK) in the context of idiopathic pulmonary fibrosis (IPF).

The workflow comprised TNIK bioactivity retrieval and curation, molecular fingerprint generation, Random Forest model development, model validation, explainability analysis, applicability-domain assessment, molecular docking, protein–ligand interaction analysis and preliminary molecular-dynamics equilibration.

All computational findings were treated as hypothesis-generating evidence rather than experimental confirmation of TNIK inhibition or therapeutic efficacy.

## 2. TNIK bioactivity data retrieval

TNIK bioactivity records were retrieved from the ChEMBL database using the human TNIK target identifier CHEMBL4527.

The initial collection contained 2,177 activity records covering multiple endpoints, including IC50, Kd, association and dissociation rates, inhibition measurements and thermal-shift measurements.

Activity endpoints were preserved as distinct measurement types rather than pooled as interchangeable quantities. The initial machine-learning dataset was restricted to IC50 measurements reported in nM with an exact equality relation.

Records with censored relations, such as greater-than or less-than measurements, were not included in the exact-IC50 regression dataset. Missing numerical values were also excluded from pIC50 calculation.

The resulting dataset contained 554 exact-IC50 records representing 453 unique compounds.

## 3. pIC50 calculation and dataset preparation

IC50 measurements reported in nM were converted into pIC50 values using:

pIC50 = 9 − log10(IC50 in nM)

Compound identifiers and assay metadata were retained to support traceability. The curated dataset included information on assay identifiers, assay confidence, target organism and activity records.

Repeated measurements were identified during compound-level dataset preparation. The final unique-compound modelling dataset contained 453 compounds.

The final cleaned activity data were stored under:

03_processed_data/activity_clean/

## 4. Molecular fingerprint generation

Molecular structures were represented using Morgan circular fingerprints generated with RDKit.

Fingerprints were calculated using a radius of 2 and a bit length of 2,048. These fingerprints encode local molecular environments and were used as numerical input features for machine-learning modelling.

The unique-compound fingerprint dataset contained 453 compounds and 2,048 fingerprint features, together with associated compound and activity information.

## 5. Random Forest model development

A Random Forest regression model was developed using the Morgan fingerprint features as input variables and pIC50 as the continuous prediction target.

The dataset was divided into:
- Training set: 362 compounds.
- Independent test set: 91 compounds.

The training and test partitions contained no overlapping compounds.

Model performance was evaluated using:
- Coefficient of determination (R²).
- Root mean squared error (RMSE).
- Mean absolute error (MAE).

Five-fold cross-validation was also performed to evaluate performance across multiple data partitions.

Additional evaluation included scaffold-based validation to investigate generalisation to structurally distinct chemical series.

Y-randomization was performed using 100 permutations to investigate whether the observed model performance depended on the original activity labels.

## 6. Model explainability

SHAP analysis was used to investigate the contributions of molecular fingerprint features to Random Forest predictions.

Feature-importance and summary outputs were generated to describe the influence of fingerprint features on model behaviour.

Because Morgan fingerprint bits represent encoded molecular environments, their contributions were interpreted as model-level explanations rather than direct evidence of specific biological mechanisms.

## 7. Applicability-domain assessment

A nearest-neighbour molecular similarity approach was used to assess whether test compounds were represented within the chemical space of the training dataset.

Tanimoto similarity was used to quantify molecular fingerprint similarity. A threshold of 0.40 was used in the documented applicability-domain assessment.

Test compounds were classified as inside or outside the applicability domain according to the selected similarity criterion.

Prediction errors were compared between the two groups using MAE. The small number of outside-domain compounds was considered when interpreting this comparison.

## 8. Protein structure preparation and docking

The TNIK kinase-domain crystal structure PDB 6RA7 was used for the structure-based investigation.

The co-crystallized reference ligand JWK was used to evaluate the docking protocol through redocking. The redocked pose was compared with the crystallographic reference pose using heavy-atom root mean square deviation (RMSD).

AutoDock Vina was used for molecular docking and virtual screening. Meeko was used in the receptor and ligand preparation workflow to generate docking-compatible molecular files.

An approved-drug library obtained from e-Drug3D contained 2,197 compounds and was used for virtual screening.

Four compounds were selected for integrated structural investigation:
- Ponatinib.
- Tepotinib.
- Entrectinib.
- Nilotinib.

Docking scores were treated as approximate scoring-function outputs rather than experimentally measured binding affinities.

Repeated-seed docking was used to examine score consistency for selected compounds.

## 9. Protein–ligand interaction analysis

Predicted ligand poses were examined using protein–ligand interaction analysis and three-dimensional molecular visualisation.

PLIP was used to identify and summarise predicted protein–ligand interactions, including hydrogen-bond and hydrophobic interaction assignments.

PyMOL was used to inspect and visualise representative docking poses and binding-pocket regions.

Additional geometric analysis considered nearby residues, potential hydrogen-bond geometry candidates, carbon contacts and distances between ligand and protein atoms.

These descriptors were interpreted as computational observations. Potential contacts were not assumed to represent experimentally confirmed interactions.

## 10. Experimental TNIK affinity record assessment

Available ChEMBL TNIK Kd records for the shortlisted drugs were examined separately from docking results.

The original activity relation, numerical limit, assay identifier, year and source document were retained when reporting these measurements.

Censored measurements were preserved using their original inequality signs rather than treated as exact affinity values.

The assay context of each measurement was considered when interpreting its relevance to direct TNIK inhibition.

## 11. Molecular-dynamics system preparation and equilibration

Selected TNIK–ligand systems were prepared for molecular-dynamics analysis.

Protein parameters were assigned using Amber protein force-field parameters, while ligand parameters were generated using GAFF2. Explicit TIP3P water was used for solvation.

Ligand protonation microstates were considered at nominal pH 7.4. The prepared systems included two singly protonated Ponatinib alternatives (N32-H and N37-H) and a singly protonated Tepotinib model.

AmberTools and tleap were used for system preparation and topology generation. GROMACS was used for energy minimisation and molecular-dynamics equilibration and trajectory analysis.

The documented restrained equilibration workflow included:
- 500 ps NVT equilibration.
- 1 ns NPT equilibration.
- A 2 fs integration timestep.
- Heavy-atom position restraints.

Trajectory analyses included backbone RMSD, ligand RMSD, hydrogen-bond analysis, ligand–protein contacts and other structural descriptors.

Because these simulations were short and restrained, their outputs were interpreted as preliminary equilibration observations rather than evidence of long-term stability under unrestrained production conditions.

## 12. Software and computational environment

The computational workflow used the following software and tools:

| Software/tool | Application |
|---|---|
| Python | Data processing and computational analysis |
| pandas and NumPy | Data manipulation and numerical operations |
| RDKit | Molecular fingerprint generation |
| scikit-learn | Random Forest modelling and validation |
| SHAP | Model explainability |
| AutoDock Vina | Molecular docking and virtual screening |
| Meeko | Receptor and ligand preparation |
| Open Babel | Molecular conversion and protonation preparation |
| AmberTools / tleap | Molecular-dynamics system preparation |
| GAFF2 | Ligand force-field parameters |
| GROMACS | Molecular-dynamics equilibration and analysis |
| PLIP | Protein–ligand interaction analysis |
| PyMOL | Three-dimensional molecular visualisation |
| Matplotlib | Scientific plotting |

Exact software versions and computational settings should be reported from the archived environment and run configuration files.

## 13. Reproducibility and data management

Raw ChEMBL and approved-drug data were retained separately from processed datasets.

Processed activity data, molecular fingerprints, model files, train/test partitions, docking outputs, interaction analyses, molecular-dynamics outputs and final evidence tables were organised within the project directory structure.

Analysis outputs were retained in their corresponding project folders to support traceability between input data, computational procedures and reported results.
