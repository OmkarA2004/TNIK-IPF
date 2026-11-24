# Discussion

## 1. Overview of the principal findings

This study used an integrated computational workflow to investigate the potential of approved drugs to interact with TNIK, a target of interest in the context of idiopathic pulmonary fibrosis (IPF). The workflow combined ChEMBL bioactivity data curation, molecular fingerprint-based machine learning, model interpretation, applicability-domain assessment, molecular docking, structural interaction analysis and preliminary restrained molecular-dynamics equilibration.

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

SHAP analysis was used to investigate the influence of molecular fingerprint features on model predictions. This provides information about the features associated with the model's predictions, but fingerprint-level contributions do not automatically establish specific molecular mechanisms.

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
