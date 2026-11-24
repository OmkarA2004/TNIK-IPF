# Computational Prioritisation of Approved Drugs for TNIK in the Context of Idiopathic Pulmonary Fibrosis

# Abstract

**Background:** Tumour necrosis factor receptor-associated factor 2 and NCK-interacting kinase (TNIK) has emerged as a potential therapeutic target in idiopathic pulmonary fibrosis (IPF) (Ren et al., 2025; Moua et al., 2024). Computational drug repurposing may help prioritise approved compounds for further investigation.

**Objective:** This study aimed to develop and evaluate a TNIK-focused machine-learning workflow and integrate its findings with molecular docking, interaction analysis and preliminary molecular dynamics simulations to prioritise approved drugs for experimental follow-up.

**Methods:** TNIK bioactivity data were retrieved from ChEMBL and curated to obtain exact IC50 measurements. A dataset of 453 unique compounds was represented using 2,048-bit Morgan fingerprints. A Random Forest regression model (Breiman, 2001) was developed to predict pIC50 values and evaluated using a held-out test set, cross-validation and scaffold-based validation. Model interpretability and applicability-domain analyses were performed. An approved-drug library from e-Drug3D was screened against the TNIK crystal structure (PDB: 6RA7) using AutoDock Vina. Selected compounds underwent interaction analysis and preliminary restrained molecular dynamics equilibration.

**Results:** The Random Forest model achieved a test-set R² of 0.4914, RMSE of 0.7976 and MAE of 0.5917. Five-fold random cross-validation produced a mean R² of 0.5439 ± 0.0931, while scaffold-based validation showed weaker performance, indicating limitations in generalisation to structurally distinct compounds. Ponatinib, tepotinib, entrectinib and nilotinib were prioritised for further examination. Their integrated docking scores ranged from −10.523 to −11.149 kcal/mol. ChEMBL records included a reported TNIK Kd of 385 nM for ponatinib; however, this measurement originated from a Kinobeads cell-lysate competition assay and does not establish inhibition in a purified-TNIK enzyme assay. Molecular dynamics analyses were limited to short, restrained equilibration simulations.

**Conclusion:** This study establishes an integrated computational workflow for TNIK-focused approved-drug prioritisation. Ponatinib warrants further experimental investigation based on the combined computational and reported binding evidence. Nevertheless, docking predictions, model outputs and restrained simulations do not establish TNIK inhibition, cellular efficacy or therapeutic benefit in IPF. Biochemical and cellular validation is required.


---


# Introduction

## 1. Idiopathic pulmonary fibrosis

Idiopathic pulmonary fibrosis (IPF) is a chronic, progressive interstitial lung disease characterized by progressive pulmonary scarring, excessive extracellular matrix deposition and deterioration of lung function. The disease is associated with substantial morbidity and mortality, and its underlying biological mechanisms involve interactions between alveolar epithelial injury, abnormal tissue repair, fibroblast activation and persistent extracellular matrix accumulation.

Although the understanding of IPF pathogenesis has improved, treatment remains challenging. Nintedanib and pirfenidone are established antifibrotic therapies that can slow the decline in lung function in some patients, but they do not reliably reverse established fibrosis. Their limitations in efficacy and tolerability highlight the continuing need to investigate additional molecular targets and therapeutic strategies.

## 2. Molecular signalling pathways involved in pulmonary fibrosis

Pulmonary fibrosis is associated with dysregulation of several interconnected signalling pathways. Transforming growth factor-beta (TGF-β), Wnt/β-catenin, Hippo-YAP/TAZ, c-Jun N-terminal kinase (JNK) and nuclear factor-kappa B (NF-κB) signalling have all been implicated in processes relevant to fibrotic disease.

These pathways influence cellular responses such as fibroblast activation, myofibroblast differentiation, cell migration, inflammation and extracellular matrix production. Their interactions contribute to the complex biological processes involved in the initiation and progression of pulmonary fibrosis.

The involvement of multiple signalling pathways provides a rationale for investigating regulatory proteins that may influence more than one profibrotic cellular process.

## 3. TNIK as a potential target in fibrosis

TRAF2- and NCK-interacting kinase (TNIK) is a serine/threonine protein kinase associated with several intracellular signalling pathways, including Wnt-related signalling. TNIK has been linked to pathways relevant to fibrosis, including Wnt, TGF-β, Hippo-YAP/TAZ, JNK and NF-κB signalling.

A study published online in 2024 and in the 2025 journal volume by the developers of a small-molecule TNIK inhibitor identified TNIK as a potential antifibrotic target and reported preclinical activity of INS018_055, subsequently named rentosertib (Ren et al., 2025). The study described antifibrotic effects across experimental models and reported safety, tolerability and pharmacokinetic findings from phase I evaluation.

Subsequently, a randomized, placebo-controlled phase 2a study of rentosertib in patients with IPF was published in 2025. The study reported a dose-dependent pattern in lung-function outcomes over 12 weeks, including a mean forced vital capacity change of +98.4 ml in the 60 mg once-daily group compared with −20.3 ml in the placebo group. The authors concluded that the findings warranted further investigation in larger and longer clinical studies (Xu et al., 2025).

These developments strengthen the rationale for studying TNIK in the context of IPF. However, the clinical development of a specific TNIK inhibitor does not establish that other drugs can engage TNIK or reproduce its biological effects. Independent evaluation of potential compounds remains important.

## 4. Drug repurposing and computational approaches

Drug repurposing investigates whether compounds with existing development, approval or pharmacological information may have additional therapeutic applications. It can provide a complementary route to drug discovery by using existing chemical and biological information to prioritize compounds for further investigation.

Computational approaches can support this process by integrating bioactivity databases, molecular representations, machine-learning models and structure-based screening. Public bioactivity resources such as ChEMBL provide compound identifiers, activity measurements and assay annotations that can be curated for quantitative structure–activity modelling (Zdrazil et al., 2024).

Molecular fingerprints provide numerical representations of chemical structures and can be used as input features for machine-learning algorithms. Random Forest regression is one approach for modelling relationships between molecular features and measured biological activity. However, model performance depends on dataset quality, chemical diversity, assay consistency and the validation strategy. Random data splitting may produce more optimistic estimates than scaffold-based evaluation when related chemical structures occur in both training and test sets.

Model explainability and applicability-domain assessment can provide additional context for interpreting predictions. Explainability methods describe which model features contribute to predictions, while applicability-domain approaches assess whether a query compound resembles the chemical space represented by the training dataset. Neither approach independently establishes biological activity.

Structure-based molecular docking provides a complementary method for investigating the possible placement of compounds within a protein binding site. Docking scores and predicted interactions can help prioritize compounds, but they are approximate computational estimates and require experimental validation. Redocking a crystallographic ligand can assess recovery of a known binding pose under a selected protocol, although successful pose recovery does not by itself demonstrate accurate ranking of unrelated compounds.

A 2025 study by Zhao and colleagues combined machine learning and structure-based virtual screening to identify potential TNIK inhibitors from DrugBank. The authors reported experimental TNIK enzyme inhibition for several candidates and identified osimertinib as a compound with measurable TNIK inhibitory activity and effects in a cellular fibrogenesis model. This study demonstrates the relevance of integrating computational prioritization with biochemical and cellular validation.

## 5. Research rationale and knowledge gap

Despite progress in TNIK-focused antifibrotic research, the evaluation of additional approved drugs as potential TNIK-binding compounds remains relevant. Computational predictions from different datasets, molecular representations, screening libraries and structural protocols may produce different candidate priorities.

An independent assessment can help examine how curated TNIK bioactivity data, molecular fingerprint-based regression, chemical-space evaluation and structure-based analysis can be combined to prioritize approved drugs for further investigation.

Importantly, a predicted docking pose or a similarity-based prediction should not be considered equivalent to direct experimental evidence of TNIK inhibition. Existing bioactivity records must also be interpreted in their original assay context, particularly when measurements originate from cell-lysate competition assays or are reported as censored affinity values.

The present study therefore uses an integrated computational workflow to examine TNIK-associated chemical activity and prioritize selected approved drugs for further investigation. The study is intended to generate testable hypotheses rather than establish therapeutic efficacy in IPF.

## 6. Aim and objectives

### Aim

To develop and apply an integrated computational workflow for prioritizing approved drugs as potential TNIK-targeting compounds relevant to idiopathic pulmonary fibrosis.

### Specific objectives

1. Retrieve and curate TNIK bioactivity records from ChEMBL while preserving activity endpoints and assay information.
2. Generate RDKit Morgan molecular fingerprints (Rogers and Hahn, 2010) and develop a Random Forest regression model (Breiman, 2001) to predict TNIK pIC50 values.
3. Evaluate model performance using independent test-set, cross-validation and scaffold-based validation strategies.
4. Investigate model feature contributions using SHAP analysis and assess prediction reliability using an applicability-domain approach.
5. Evaluate the TNIK crystal structure-based docking workflow through redocking of the co-crystallized reference ligand.
6. Screen an approved-drug collection and prioritize selected compounds using docking and structural interaction analyses.
7. Compare computational findings with available experimental TNIK affinity records while preserving their assay provenance and measurement limitations.
8. Examine preliminary restrained molecular-dynamics equilibration outputs for selected TNIK–drug systems.
9. Identify candidates and experimental questions that warrant further biochemical and cellular validation.

The overall objective is to establish a transparent and reproducible computational prioritization framework that can support subsequent experimental investigation of TNIK-targeting compounds in the context of IPF.

# Materials and Methods

## 1. Study design and computational workflow

A computational drug-repurposing workflow was developed to investigate potential interactions between approved drugs and TRAF2- and NCK-interacting kinase (TNIK) in the context of idiopathic pulmonary fibrosis (IPF) (Ren et al., 2025; Moua et al., 2024).

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

A Random Forest regression model (Breiman, 2001) was developed using the Morgan fingerprint features as input variables and pIC50 as the continuous prediction target.

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

SHAP analysis (Lundberg and Lee, 2017) was used to investigate the contributions of molecular fingerprint features to Random Forest predictions.

Feature-importance and summary outputs were generated to describe the influence of fingerprint features on model behaviour.

Because Morgan fingerprint bits represent encoded molecular environments, their contributions were interpreted as model-level explanations rather than direct evidence of specific biological mechanisms.

## 7. Applicability-domain assessment

A nearest-neighbour molecular similarity approach was used to assess whether test compounds were represented within the chemical space of the training dataset.

Tanimoto similarity was used to quantify molecular fingerprint similarity. A threshold of 0.40 was used in the documented applicability-domain assessment.

Test compounds were classified as inside or outside the applicability domain according to the selected similarity criterion.

Prediction errors were compared between the two groups using MAE. The small number of outside-domain compounds was considered when interpreting this comparison.

## 8. Protein structure preparation and docking

The TNIK kinase-domain crystal structure PDB 6RA7 was used for the structure-based investigation (Worldwide Protein Data Bank, 2019).

The co-crystallized reference ligand JWK was used to evaluate the docking protocol through redocking. The redocked pose was compared with the crystallographic reference pose using heavy-atom root mean square deviation (RMSD).

AutoDock Vina was used for molecular docking and virtual screening (Trott and Olson, 2010; Eberhardt et al., 2021). Meeko was used in the receptor and ligand preparation workflow to generate docking-compatible molecular files.

An approved-drug library obtained from e-Drug3D contained 2,197 compounds and was used for virtual screening (Pihan et al., 2012).

Four compounds were selected for integrated structural investigation:
- Ponatinib.
- Tepotinib.
- Entrectinib.
- Nilotinib.

Docking scores were treated as approximate scoring-function outputs rather than experimentally measured binding affinities.

Repeated-seed docking was used to examine score consistency for selected compounds.

## 9. Protein–ligand interaction analysis

Predicted ligand poses were examined using protein–ligand interaction analysis and three-dimensional molecular visualisation.

PLIP (Salentin et al., 2015) was used to identify and summarise predicted protein–ligand interactions, including hydrogen-bond and hydrophobic interaction assignments.

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

AmberTools and tleap were used for system preparation and topology generation; however, the exact AmberTools version could not be verified from the available software records. GROMACS (version 2026.3; Abraham et al., 2015) was used for energy minimisation, molecular-dynamics equilibration and trajectory analysis. The computational environment used Python 3.10.15, RDKit 2026.03.1, AutoDock Vina 1.2.7, Meeko 0.8.0, SHAP 0.49.1, scikit-learn 1.7.2, NumPy 1.26.4, pandas 2.3.3 and Biopython 1.88.

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


---

# Results

## 1. TNIK bioactivity dataset and chemical diversity

A total of 2,177 TNIK-associated activity records were retrieved from ChEMBL. Among the reported activity types, 644 records were annotated as IC50 measurements. After restricting the dataset to numeric IC50 values reported in nM with an exact equality relation, 554 records were retained for pIC50 calculation. These measurements represented 453 unique compounds, with repeated measurements present for a subset of compounds.

The resulting unique-compound dataset contained 453 compounds and 2,048-bit Morgan fingerprints. The pIC50 values had a mean of 7.3656, a standard deviation of 1.0001, and a range of 4.60–9.92. Assay metadata were available for all 554 retained records. Of these, 324 records were associated with Homo sapiens and confidence score 9, while 230 records had confidence score 8 but no organism annotation.

## 2. Random Forest model development and predictive performance

A Random Forest regression model (Breiman, 2001) was developed using Morgan molecular fingerprints to predict TNIK pIC50 values. The 453 unique compounds were divided into training and independent test sets, containing 362 and 91 compounds, respectively. No compounds were shared between the two sets.

On the independent test set, the model achieved an R² of 0.4914, an RMSE of 0.7976 and an MAE of 0.5917. These results indicate moderate predictive performance, with prediction errors remaining substantial for some compounds.

Five-fold cross-validation using the full dataset produced a mean R² of 0.5439 ± 0.0931, RMSE of 0.6664 ± 0.0803 and MAE of 0.4877 ± 0.0585. The cross-validation results were broadly consistent with the independent test-set evaluation.

Validation using chemical scaffolds gave weaker performance than random splitting. The all-assay scaffold-validation R² was approximately 0.4024. In the high-confidence subset of 224 compounds, random cross-validation gave an R² of 0.5044 ± 0.1819, RMSE of 0.8054 ± 0.1217 and MAE of 0.6023 ± 0.0791. High-confidence scaffold validation gave an R² of 0.3726 ± 0.2745, RMSE of 0.8523 ± 0.1351 and MAE of 0.6693.

The lower and more variable scaffold-validation performance suggests that predictive generalisation to structurally distinct chemical series is more challenging than prediction under random splitting. The model should therefore be interpreted as a prioritisation tool rather than a reliable predictor for every new chemical scaffold.

## 3. Model interpretation and applicability domain

SHAP analysis was performed to examine the contribution of molecular fingerprint features to the Random Forest predictions. The resulting feature-importance and summary outputs provide a description of which fingerprint features influence model predictions. Because Morgan fingerprint bits do not automatically correspond to uniquely interpretable chemical mechanisms, these results should not be treated as direct evidence of specific molecular interactions or biological mechanisms.

Applicability-domain analysis classified 86 of the 91 test compounds as inside the domain and 5 as outside it. The MAE was 0.5466 for compounds inside the domain and 1.3664 for compounds outside it. This difference suggests poorer prediction reliability for compounds outside the model's represented chemical space. However, the outside-domain group contained only five compounds, so this comparison should be interpreted cautiously.

## 4. Docking validation and approved-drug screening

The TNIK crystal structure PDB 6RA7 was used for the docking workflow. Redocking of the reference ligand JWK produced a docking score of -9.770 kcal/mol and a heavy-atom RMSD of 0.782 Å. This supports recovery of the reference binding pose under the selected protocol. It does not, by itself, demonstrate enrichment of known active compounds or establish the accuracy of virtual-screening rankings.

An approved-drug library obtained from e-Drug3D was screened, and four compounds were selected for integrated follow-up analysis: Ponatinib, Tepotinib, Entrectinib and Nilotinib. Their integrated docking scores were -11.139, -11.149, -10.736 and -10.523 kcal/mol, respectively.

Tepotinib and Ponatinib therefore produced the most favourable docking scores among these four compounds, with a difference of only 0.010 kcal/mol. Repeated-seed docking gave mean scores of -11.166 ± 0.027 kcal/mol for Ponatinib and -11.128 ± 0.041 kcal/mol for Tepotinib. These small differences do not support a definitive potency ranking between the two compounds.

Structural interaction analyses were also performed to compare the predicted binding environments. The resulting interaction diagrams and contact analyses provide structural hypotheses for follow-up investigation. They do not establish that the predicted contacts occur experimentally.

## 5. Existing experimental TNIK activity evidence for shortlisted drugs

The ChEMBL records contain a reported TNIK Kd value of 385 nM for Ponatinib. This value is associated with a Kinobeads cell-lysate competition assay and should not be described as a purified-TNIK enzyme IC50 or as proof of functional TNIK inhibition.

Tepotinib and Entrectinib have reported Kd values greater than 30,000 nM in the cited 2017 ChEMBL assay record. Nilotinib has reported values greater than 10,000 nM in a 2011 record and greater than 30,000 nM in a 2017 record. These are censored measurements and should be described using their reported inequality signs rather than converted into exact values.

The experimental observations do not produce the same apparent ordering as the docking scores. In particular, the slightly more favourable docking score for Tepotinib does not indicate stronger experimental TNIK binding. Ponatinib is therefore a reasonable candidate for further experimental investigation because of its reported Kd observation, but the assay context does not establish purified-enzyme inhibition or therapeutic relevance to idiopathic pulmonary fibrosis.

Similarity comparisons to reference TNIK compounds were also performed. The maximum Tanimoto similarities reported for the four candidates were low, ranging from approximately 0.098 to 0.204. Activities of their nearest reference compounds must not be transferred to the shortlisted drugs; those measurements belong to the reference compounds, not to the candidates.

## 6. Restrained molecular-dynamics equilibration

Short restrained equilibration simulations were performed for the selected TNIK–drug systems. The workflow included 500 ps of NVT equilibration followed by 1 ns of NPT equilibration, using a 2 fs timestep and heavy-atom restraints.

The available temperature and structural analyses provide preliminary information about the equilibration behaviour of the prepared systems. However, the simulations were short and restrained. Low backbone RMSD values under these conditions cannot be interpreted as evidence of long-term stability in unrestrained production simulations. The MD results are therefore treated as preliminary supporting observations rather than independent validation of stable drug binding.

## 7. Overall interpretation

The computational workflow identified a set of approved drugs for further investigation of possible TNIK engagement. The Random Forest model showed moderate predictive performance, with weaker generalisation under scaffold-based validation. Docking and structural analyses prioritised Ponatinib and Tepotinib, while the available experimental Kd records provided a different and assay-dependent picture.

Taken together, the findings support further experimental evaluation rather than a claim that any shortlisted compound is a confirmed TNIK inhibitor or an effective treatment for idiopathic pulmonary fibrosis.


---

# Discussion

## 1. Overview of the principal findings

This study used an integrated computational workflow to investigate the potential of approved drugs to interact with TNIK, a target of interest in the context of idiopathic pulmonary fibrosis (IPF) (Ren et al., 2025; Moua et al., 2024). The workflow combined ChEMBL bioactivity data curation, molecular fingerprint-based machine learning, model interpretation, applicability-domain assessment, molecular docking, structural interaction analysis and preliminary restrained molecular-dynamics equilibration.

The Random Forest model demonstrated moderate predictive performance on the independent test set. Docking and structural analyses prioritised Ponatinib and Tepotinib among the four shortlisted drugs. However, the available experimental TNIK Kd records did not support a simple ranking based on docking scores. These findings highlight the importance of combining computational predictions with experimental evidence and careful consideration of assay context.

## 2. Interpretation of the TNIK bioactivity dataset

The curated dataset contained 453 unique compounds with exact numeric IC50 measurements suitable for pIC50 calculation. This provided a basis for developing a TNIK activity prediction model using Morgan molecular fingerprints.

The dataset also showed limitations relevant to model development. Compounds had repeated activity records, assay metadata varied, and organism annotations were unavailable for a substantial proportion of measurements. Although assay confidence scores were available for all retained records, differences in assay conditions and measurement contexts may contribute to variability in the observed activity values.

The model results should therefore be interpreted as predictions based on the available curated ChEMBL dataset rather than as universal estimates of TNIK activity.

## 3. Interpretation of Random Forest model performance

The Random Forest model achieved an independent test-set R² of 0.4914, RMSE of 0.7976 and MAE of 0.5917. Five-fold random cross-validation produced a mean R² of 0.5439 ± 0.0931.

These results indicate that the molecular fingerprints captured useful relationships between chemical structure and reported TNIK activity. Nevertheless, the prediction errors indicate that substantial uncertainty remains for individual compounds.

An important observation was the reduction in performance under scaffold-based validation. The all-assay scaffold-validation R² was approximately 0.4024, while high-confidence scaffold validation produced an R² of 0.3726 ± 0.2745. This weaker performance suggests that the model has greater difficulty predicting compounds belonging to chemical scaffolds that are not adequately represented in the training data.

Consequently, the model is more appropriately used to support prioritisation within its represented chemical space than to make definitive predictions for structurally novel compounds. Future model improvements could investigate additional chemical representations, larger curated datasets and validation strategies that explicitly assess generalisation to new scaffolds.

## 4. Model interpretation and applicability-domain considerations

SHAP analysis (Lundberg and Lee, 2017) was used to investigate the influence of molecular fingerprint features on model predictions. This provides information about the features associated with the model's predictions, but fingerprint-level contributions do not automatically establish specific molecular mechanisms.

Applicability-domain analysis classified 86 of 91 test compounds as inside the domain and five as outside it. The inside-domain MAE was 0.5466, compared with 1.3664 for the outside-domain compounds. This difference suggests that prediction reliability may decrease for compounds that are less represented in the model's chemical space.

However, the outside-domain group contained only five compounds. The difference should therefore be treated as an indication requiring further investigation rather than as a definitive statistical conclusion.

## 5. Interpretation of molecular docking and structural interactions

The redocking of the reference ligand JWK into TNIK structure 6RA7 produced a heavy-atom RMSD of 0.782 Å and a docking score of -9.770 kcal/mol. The RMSD supports recovery of the reference ligand pose under the selected docking protocol. However, successful redocking does not establish virtual-screening enrichment or guarantee that the docking scores accurately rank unrelated compounds.

Among the four shortlisted approved drugs, Tepotinib and Ponatinib produced the most favourable integrated docking scores, at -11.149 and -11.139 kcal/mol, respectively. Their difference was very small. Repeated-seed docking also showed closely similar mean scores for Ponatinib and Tepotinib. Therefore, the docking results do not justify a definitive ranking of their TNIK-binding potency.

Structural interaction analysis provided additional information about the predicted positioning and contacts of the compounds within the binding pocket. These observations can help formulate hypotheses for experimental investigation. Nevertheless, predicted contacts and interaction assignments should not be treated as experimentally established binding interactions.

## 6. Importance of experimental TNIK Kd evidence

The available ChEMBL records provide an important qualification to the docking-based prioritisation.

Ponatinib has a reported TNIK Kd value of 385 nM. However, this measurement originated from a Kinobeads cell-lysate competition assay. It is not equivalent to a purified-TNIK enzyme IC50 measurement and does not independently demonstrate functional inhibition of TNIK.

Tepotinib and Entrectinib have censored Kd records greater than 30,000 nM in the cited 2017 assay record. Nilotinib has two censored records, greater than 10,000 nM in 2011 and greater than 30,000 nM in 2017. These values must retain their reported inequality signs and should not be interpreted as exact measurements.

The difference between docking-based prioritisation and the available experimental records illustrates why docking scores alone are insufficient for selecting candidate inhibitors. In particular, Tepotinib's slightly more favourable docking score does not establish stronger experimental TNIK binding than Ponatinib.

Ponatinib is a reasonable candidate for further experimental investigation because of its reported Kd observation and favourable docking results. However, its assay provenance must remain explicit, and direct biochemical testing is required before describing it as a confirmed TNIK inhibitor.

The nearest-neighbour similarity analysis also requires careful interpretation. The maximum Tanimoto similarities between the shortlisted drugs and the reference compounds were low. The Kd values belonging to nearest-neighbour compounds cannot be assigned to the shortlisted drugs. They provide reference information only.

## 7. Interpretation of restrained molecular-dynamics equilibration

The molecular-dynamics workflow included 500 ps of NVT equilibration followed by 1 ns of NPT equilibration, with a 2 fs timestep and heavy-atom restraints.

These simulations provide preliminary information about the behaviour of the prepared systems during restrained equilibration. However, the use of restraints and the short simulation duration limit the conclusions that can be drawn. Low backbone RMSD under heavy-atom restraints does not demonstrate long-term stability of the ligand-bound complex under unrestrained conditions.

The available MD analyses should therefore be presented as preliminary supporting observations. Longer unrestrained production simulations, with appropriate replicate trajectories and convergence assessment, would be required for stronger conclusions about dynamic stability.

## 8. Relevance to the TNIK–IPF research hypothesis

The computational findings provide a basis for investigating whether selected approved drugs can engage TNIK. They do not establish that the shortlisted compounds modify TNIK-dependent biological pathways in lung cells, reduce fibrotic responses or provide therapeutic benefit in IPF.

The distinction between predicted target engagement and disease-level efficacy is essential. Experimental evaluation would need to establish direct TNIK binding or inhibition, followed by appropriate cellular studies examining relevant pathway responses and fibrotic phenotypes. Further pharmacological and safety evaluation would also be necessary before considering therapeutic translation.

## 9. Limitations

Several limitations should be considered when interpreting this study:

1. The curated TNIK activity dataset contained a relatively limited number of unique compounds, which constrains model generalisation.
2. Random-split validation produced stronger performance than scaffold-based validation, indicating reduced reliability for structurally novel chemical series.
3. The SHAP analysis describes model feature contributions and does not independently establish molecular mechanisms.
4. Applicability-domain assessment identified only five outside-domain test compounds, limiting the strength of comparisons between inside- and outside-domain errors.
5. Successful reference-ligand redocking demonstrates pose recovery under the selected protocol but does not establish virtual-screening enrichment.
6. Docking scores are computational estimates and should not be interpreted as measured binding affinities or definitive potency rankings.
7. The experimental Kd records have assay-specific limitations, particularly the cell-lysate Kinobeads context of the Ponatinib measurement.
8. The molecular-dynamics simulations were short restrained equilibration runs rather than long unrestrained production simulations.
9. Computational findings do not establish cellular TNIK inhibition, antifibrotic activity, clinical efficacy or safety in IPF.

## 10. Future directions

The most important next step is experimental validation of the shortlisted compounds. Ponatinib is a reasonable priority for direct TNIK biochemical testing, while Tepotinib may also be investigated because of its favourable docking results despite its censored weak-binding Kd record.

Future studies should use an appropriate purified-TNIK biochemical assay, followed by orthogonal binding or target-engagement experiments. Cellular studies should then investigate whether confirmed TNIK modulation produces relevant effects on fibrotic signalling and cellular phenotypes.

Additional computational work could include longer unrestrained molecular-dynamics simulations, more extensive validation against known active and inactive compounds, and evaluation of alternative molecular representations. These extensions should be guided by the results of experimental validation.

## 11. Overall conclusion

This study establishes an integrated computational framework for TNIK-focused approved-drug prioritisation in the context of IPF research. The Random Forest model showed moderate predictive performance, while scaffold-based validation revealed limitations in generalisation to structurally distinct compounds. Docking and structural analyses prioritised Ponatinib and Tepotinib, but experimental assay records highlighted the need for cautious interpretation.

Ponatinib merits further investigation because of its reported TNIK Kd observation and computational prioritisation. However, the available evidence does not establish purified-enzyme TNIK inhibition or therapeutic efficacy in IPF. Direct biochemical and cellular validation remains essential.


---

# TNIK–IPF Computational Study: Figure Captions

## Main manuscript figures

### Figure 1. Performance of the TNIK activity prediction model
Actual versus predicted pIC50 values for the Random Forest regression model (Breiman, 2001) trained using Morgan molecular fingerprints. The model was evaluated on an independent test set of 91 compounds. Test-set performance was R² = 0.4914, RMSE = 0.7976 and MAE = 0.5917. The results indicate moderate predictive performance.

File: 09_results/figures/TNIK_RF_actual_vs_predicted.png

### Figure 2. SHAP-based interpretation of the TNIK prediction model
SHAP summary plot showing the contribution of molecular fingerprint features to Random Forest predictions. These contributions describe model behaviour and should not be interpreted directly as experimentally confirmed molecular mechanisms.

File: 09_results/figures/TNIK_SHAP_summary.png

### Figure 3. Predicted binding pose of Ponatinib in the TNIK binding pocket
Molecular docking pose of Ponatinib within the TNIK crystal structure (PDB: 6RA7). The figure illustrates the predicted ligand orientation and its spatial relationship with the binding pocket. Docking provides a computational binding hypothesis and does not independently establish biochemical inhibition.

File: 09_results/figures/TNIK_ponatinib_final_binding_pose.png

### Figure 4. Comparative docking and interaction analysis of shortlisted drugs
Comparison of Ponatinib, Tepotinib, Entrectinib and Nilotinib using integrated docking and structural-interaction analysis. Their docking scores were -11.139, -11.149, -10.736 and -10.523 kcal/mol, respectively. Small differences in docking scores should not be interpreted as definitive differences in experimental TNIK potency.

File: 09_results/figures/TNIK_four_drug_comparative_interactions.png

### Figure 5. Predicted interaction profiles of Ponatinib and Tepotinib
Protein–ligand interaction analysis of Ponatinib and Tepotinib in the TNIK binding pocket. Panel A: Ponatinib. Panel B: Tepotinib. The diagrams show computationally identified contacts, not experimentally confirmed interactions.

Panel A: 06_docking/interaction_analysis/comparative_analysis/3D_figures/TNIK_Ponatinib_PLIP.png
Panel B: 06_docking/interaction_analysis/comparative_analysis/3D_figures/TNIK_Tepotinib_PLIP.png

## Supplementary figures

### Supplementary Figure S1. Restrained equilibration backbone RMSD
Backbone RMSD profiles for the TNIK–drug systems during restrained equilibration. The simulations included NVT and NPT equilibration stages. These short, restrained simulations do not establish long-term unrestrained complex stability.

File: 06_docking/molecular_dynamics/analysis/npt_stability/plots/TNIK_three_complexes_backbone_RMSD.png

### Supplementary Figure S2. PLIP interaction-count comparison
Comparison of interaction counts identified during the structural analysis of shortlisted TNIK–drug complexes.

File: 06_docking/interaction_analysis/comparative_analysis/figures/TNIK_PLIP_interaction_counts.png

### Supplementary Figure S3. Ligand–protein contact occupancy
Comparison of ligand–protein contact occupancy in the available restrained equilibration trajectories. These observations are preliminary and should not be interpreted as production-MD evidence.

File: 06_docking/interaction_analysis/comparative_analysis/figures/TNIK_PLIP_MD_contact_occupancy.png

The applicability-domain analysis is reported in the Results text and its source CSV. No corresponding plot was located during the figure-file check, so it is not assigned a figure number.

Source: 07_explainability/applicability_domain/TNIK_applicability_domain_results.csv

# References


 (2025). **A small-molecule TNIK inhibitor targets fibrosis in preclinical and clinical models.** Nature Biotechnology. DOI: https://doi.org/10.1038/s41587-024-02143-0

 (2025). **A generative AI-discovered TNIK inhibitor for idiopathic pulmonary fibrosis: a randomized phase 2a trial.** Nature Medicine. DOI: https://doi.org/10.1038/s41591-025-03743-2

 (2025). **Integrated Machine Learning and Structure-Based Virtual Screening Identify Osimertinib as a TNIK Inhibitor for Idiopathic Pulmonary Fibrosis.** Journal of Chemical Information and Modeling. DOI: https://doi.org/10.1021/acs.jcim.5c01521

 (2024). **The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods.** Nucleic Acids Research. DOI: https://doi.org/10.1093/nar/gkad1004

 (2012). **e-Drug3D: 3D structure collections dedicated to drug repurposing and fragment-based drug design.** Bioinformatics. DOI: https://doi.org/10.1093/bioinformatics/bts186

{Worldwide Protein Data Bank (2019). **Human TNIK in complex with compound 9: PDB entry 6RA7.** . DOI: https://doi.org/10.2210/pdb6ra7/pdb

Breiman, Leo (2001). **Random Forests.** Machine Learning. DOI: https://doi.org/10.1023/A:1010933404324

Pedregosa, Fabian and Varoquaux, Gael and Gramfort, Alexandre and others (2011). **Scikit-learn: Machine Learning in Python.** Journal of Machine Learning Research.

Lundberg, Scott M. and Lee, Su-In (2017). **A Unified Approach to Interpreting Model Predictions.** Advances in Neural Information Processing Systems.

Trott, Oleg and Olson, Arthur J. (2010). **AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading.** Journal of Computational Chemistry. DOI: https://doi.org/10.1002/jcc.21334

Eberhardt, Jan and Santos-Martins, Diogo and Tillack, Andreas F. and Forli, Stefano (2021). **AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings.** Journal of Chemical Information and Modeling. DOI: https://doi.org/10.1021/acs.jcim.1c00203

Salentin, Sebastian and Schreiber, Sven and Haupt, V. Joachim and Adasme, Melissa F. and Schroeder, Michael (2015). **PLIP: fully automated protein-ligand interaction profiler.** Nucleic Acids Research. DOI: https://doi.org/10.1093/nar/gkv315

Abraham, Mark James and Murtola, Teemu and Schulz, Roland and others (2015). **GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers.** SoftwareX. DOI: https://doi.org/10.1016/j.softx.2015.06.001

Moua, Tao and Baqir, Muneeb and Ryu, Jay H. (2024). **What Is on the Horizon for Treatments in Idiopathic Pulmonary Fibrosis?.** Journal of Clinical Medicine. DOI: https://doi.org/10.3390/jcm13216304

Rogers, David and Hahn, Mathew (2010). **Extended-Connectivity Fingerprints.** Journal of Chemical Information and Modeling. DOI: https://doi.org/10.1021/ci100050t
