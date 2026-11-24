# Results

## 1. TNIK bioactivity dataset and chemical diversity

A total of 2,177 TNIK-associated activity records were retrieved from ChEMBL. Among the reported activity types, 644 records were annotated as IC50 measurements. After restricting the dataset to numeric IC50 values reported in nM with an exact equality relation, 554 records were retained for pIC50 calculation. These measurements represented 453 unique compounds, with repeated measurements present for a subset of compounds.

The resulting unique-compound dataset contained 453 compounds and 2,048-bit Morgan fingerprints. The pIC50 values had a mean of 7.3656, a standard deviation of 1.0001, and a range of 4.60–9.92. Assay metadata were available for all 554 retained records. Of these, 324 records were associated with Homo sapiens and confidence score 9, while 230 records had confidence score 8 but no organism annotation.

## 2. Random Forest model development and predictive performance

A Random Forest regression model was developed using Morgan molecular fingerprints to predict TNIK pIC50 values. The 453 unique compounds were divided into training and independent test sets, containing 362 and 91 compounds, respectively. No compounds were shared between the two sets.

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
