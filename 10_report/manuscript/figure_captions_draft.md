# TNIK–IPF Computational Study: Figure Captions

## Main manuscript figures

### Figure 1. Performance of the TNIK activity prediction model
Actual versus predicted pIC50 values for the Random Forest regression model trained using Morgan molecular fingerprints. The model was evaluated on an independent test set of 91 compounds. Test-set performance was R² = 0.4914, RMSE = 0.7976 and MAE = 0.5917. The results indicate moderate predictive performance.

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

## Applicability-domain result
The applicability-domain analysis is reported in the Results text and its source CSV. No corresponding plot was located during the figure-file check, so it is not assigned a figure number.

Source: 07_explainability/applicability_domain/TNIK_applicability_domain_results.csv
