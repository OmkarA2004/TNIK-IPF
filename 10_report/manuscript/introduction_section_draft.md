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

A 2024 study by the developers of a small-molecule TNIK inhibitor identified TNIK as a potential antifibrotic target and reported preclinical activity of INS018_055, subsequently named rentosertib. The study described antifibrotic effects across experimental models and reported safety, tolerability and pharmacokinetic findings from phase I evaluation.

Subsequently, a randomized, placebo-controlled phase 2a study of rentosertib in patients with IPF was published in 2025. The study reported a dose-dependent pattern in lung-function outcomes over 12 weeks, including a mean forced vital capacity change of +98.4 ml in the 60 mg once-daily group compared with −20.3 ml in the placebo group. The authors concluded that the findings warranted further investigation in larger and longer clinical studies.

These developments strengthen the rationale for studying TNIK in the context of IPF. However, the clinical development of a specific TNIK inhibitor does not establish that other drugs can engage TNIK or reproduce its biological effects. Independent evaluation of potential compounds remains important.

## 4. Drug repurposing and computational approaches

Drug repurposing investigates whether compounds with existing development, approval or pharmacological information may have additional therapeutic applications. It can provide a complementary route to drug discovery by using existing chemical and biological information to prioritize compounds for further investigation.

Computational approaches can support this process by integrating bioactivity databases, molecular representations, machine-learning models and structure-based screening. Public bioactivity resources such as ChEMBL provide compound identifiers, activity measurements and assay annotations that can be curated for quantitative structure–activity modelling.

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
2. Generate RDKit Morgan molecular fingerprints and develop a Random Forest regression model to predict TNIK pIC50 values.
3. Evaluate model performance using independent test-set, cross-validation and scaffold-based validation strategies.
4. Investigate model feature contributions using SHAP analysis and assess prediction reliability using an applicability-domain approach.
5. Evaluate the TNIK crystal structure-based docking workflow through redocking of the co-crystallized reference ligand.
6. Screen an approved-drug collection and prioritize selected compounds using docking and structural interaction analyses.
7. Compare computational findings with available experimental TNIK affinity records while preserving their assay provenance and measurement limitations.
8. Examine preliminary restrained molecular-dynamics equilibration outputs for selected TNIK–drug systems.
9. Identify candidates and experimental questions that warrant further biochemical and cellular validation.

The overall objective is to establish a transparent and reproducible computational prioritization framework that can support subsequent experimental investigation of TNIK-targeting compounds in the context of IPF.

# References cited in this Introduction

1. Xu, Z., Ren, F., Wang, P., et al. (2025). A generative AI-discovered TNIK inhibitor for idiopathic pulmonary fibrosis: a randomized phase 2a trial. Nature Medicine, 31, 2602–2610. https://doi.org/10.1038/s41591-025-03743-2

2. A small-molecule TNIK inhibitor targets fibrosis in preclinical and clinical models. (2025). Nature Biotechnology, 43, 63–75. https://doi.org/10.1038/s41587-024-02143-0

3. Moua, T., Baqir, M., and Ryu, J. H. (2024). What Is on the Horizon for Treatments in Idiopathic Pulmonary Fibrosis? Journal of Clinical Medicine, 13(21), 6304. https://doi.org/10.3390/jcm13216304

4. Zhao, L., Liu, H., Yao, X., Ma, X., Liu, B., Li, B., Tong, H. H. Y., and Zhang, Q. (2025). Integrated Machine Learning and Structure-Based Virtual Screening Identify Osimertinib as a TNIK Inhibitor for Idiopathic Pulmonary Fibrosis. Journal of Chemical Information and Modeling, 65(19), 10673–10687. https://doi.org/10.1021/acs.jcim.5c01521
