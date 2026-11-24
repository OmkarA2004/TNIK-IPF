# Integrated TNIK–IPF Drug Repurposing Project
## Master Results Report

**Report preparation date:** 2026-11-07

**Purpose:** Consolidation of the existing computational and experimental-evidence records for the TNIK–IPF drug repurposing project.

**Scope:** This document brings together previously generated machine-learning, molecular-docking and molecular-dynamics results. It does not represent a new computational analysis.

---

## 1. Project overview

The project investigates TNIK as a potential target for drug repurposing in idiopathic pulmonary fibrosis (IPF). The computational workflow includes ChEMBL activity-data processing, machine-learning prediction, structure-based virtual screening, docking validation, ligand interaction analysis and molecular-dynamics equilibration.

The findings below must be interpreted as computational prioritisation and supporting evidence, not as proof of therapeutic efficacy.


## 2. Integrated project summary


# Integrated TNIK–IPF Project Results Summary

Generated: 2026-11-07

## 1. Project overview

This project investigates the potential of drug repurposing against
TNIK (TRAF2 and NCK-interacting kinase) in the context of idiopathic
pulmonary fibrosis (IPF), using cheminformatics, machine learning,
molecular docking and molecular dynamics simulations.

The computational workflow combines:
- ChEMBL bioactivity data processing
- Molecular fingerprint generation
- Random Forest activity prediction
- Molecular docking and pose analysis
- Ligand parameterisation
- Molecular dynamics simulations
- Protein–ligand interaction analysis

## 2. ChEMBL dataset and activity processing

Target: Human TNIK, CHEMBL4527.

- Total activity records: 2,177
- IC50 records: 644
- Numerical IC50 values in nM: 622
- Exact IC50 records retained: 554
- Unique compounds: 453

The activity endpoint was converted to pIC50 using:

pIC50 = 9 - log10(IC50 in nM)

Unique-compound pIC50 statistics:
- Mean: 7.365578
- Standard deviation: 1.000105
- Minimum: 4.600000
- Median: 7.366532
- Maximum: 9.920819

## 3. Machine learning results

### 3.1 Molecular representation

Morgan fingerprints:
- Radius: 2
- Fingerprint length: 2,048 bits
- Unique-compound dataset: 453 compounds

### 3.2 Train/test evaluation

- Training compounds: 362
- Test compounds: 91
- Shared compounds between sets: 0

Random Forest test-set performance:
- R²: 0.4914
- RMSE: 0.7976
- MAE: 0.5917

### 3.3 Cross-validation

Five-fold cross-validation:
- Mean R²: 0.5439 ± 0.0931
- Mean RMSE: 0.6664 ± 0.0803
- Mean MAE: 0.4877 ± 0.0585

Scaffold cross-validation:
- R²: 0.4024
- RMSE: 0.7502
- MAE: 0.5630

### 3.4 Additional model validation

Y-randomisation (100 permutations):
- Mean R²: -0.2142
- Maximum R²: -0.0774
- Empirical p-value: 0.0099

Applicability-domain assessment using a Tanimoto threshold of 0.40:
- Test compounds inside domain: 86/91
- Inside-domain MAE: 0.5466
- Outside-domain MAE: 1.3664 (5 compounds)

Interpretation:
The Random Forest model demonstrates moderate predictive performance.
Its predictions should be interpreted within the applicability domain
and alongside the limitations of the dataset and validation design.

## 4. Molecular docking results

### 4.1 Target structure and redocking validation

- TNIK crystal structure: PDB 6RA7
- Protein chain: A
- Resolution: 1.20 Å
- Co-crystallised ligand: JWK

JWK redocking:
- Docking score: -9.770 kcal/mol
- Heavy-atom RMSD: 0.782 Å

The redocking RMSD supports reproduction of the reference ligand pose
within the docking workflow.

### 4.2 Virtual screening

- e-Drug3D library: 2,197 compounds

Top docking scores:
- Zavegepant: -11.49 kcal/mol
- Ponatinib: -11.21 kcal/mol
- Tepotinib: -11.15 kcal/mol

Shortlisted candidates:
- Ponatinib
- Tepotinib
- Entrectinib
- Nilotinib

### 4.3 Repeated-seed comparison

Ponatinib:
- Mean score: approximately -11.166 kcal/mol
- Standard deviation: approximately 0.027

Tepotinib:
- Mean score: approximately -11.128 kcal/mol
- Standard deviation: approximately 0.041

The difference between these mean scores is small and should not
be interpreted as evidence of a meaningful difference in binding
affinity or biological potency.

## 5. Experimental bioactivity evidence

ChEMBL evidence identified an exact canonical compound match for
Ponatinib (CHEMBL1171837), with a reported Kd of 385 nM in a
cell-lysate competition assay.

This result is not equivalent to a purified-protein TNIK inhibition
measurement and does not independently establish functional inhibition.

Other reported values:
- Tepotinib: >30,000 nM
- Entrectinib: >30,000 nM
- Nilotinib: >10,000 nM / >30,000 nM

A separate patent source reported 16 Ponatinib Kd values between
3 and 95 nM. These results originate from a different source and
assay context and should not be directly combined with the ChEMBL
measurement.

## 6. Molecular dynamics results

### 6.1 Systems prepared

The MD workflow used the TNIK 6RA7 protein structure and
GAFF2-parameterised ligands.

Systems included:
- Apo TNIK
- TNIK–Ponatinib, N32-H protonation state
- TNIK–Ponatinib, N37-H protonation state
- TNIK–Tepotinib, +1 protonation state

### 6.2 Simulation preparation and equilibration

- Protein force field: ff14SB with phosphoresidue support
- Water model: TIP3P
- Ligand force field: GAFF2
- Energy minimisation completed for all three ligand complexes
- NVT equilibration: 500 ps
- NPT equilibration: 1 ns

The ligand parameterisation checks confirmed atom counts,
approximate net charges and successful parmchk2 validation.

### 6.3 Protein and ligand RMSD

Protein backbone RMSD mean:
- Ponatinib N32-H: 0.0220 nm
- Ponatinib N37-H: 0.0223 nm
- Tepotinib: 0.0220 nm

Ligand heavy-atom RMSD mean, aligned to protein backbone:
- Ponatinib N32-H: 0.0317 nm
- Ponatinib N37-H: 0.0295 nm
- Tepotinib: 0.0291 nm

These simulations used heavy-atom restraints. Consequently, the
low RMSD values do not independently establish unrestrained
protein–ligand stability.

### 6.4 Protein–ligand contacts

Mean contact counts at 0.35 nm:
- Ponatinib N32-H: 31.198
- Ponatinib N37-H: 32.010
- Tepotinib: 33.376

These are geometric contact measurements, not binding-affinity
measurements. Differences may also reflect ligand size and atom count.

### 6.5 Hydrogen bonds

- Ponatinib N32-H: no hydrogen bonds detected under the applied
  GROMACS criteria.
- Ponatinib N37-H: no hydrogen bonds detected under the applied
  GROMACS criteria.
- Tepotinib: mean 2.3168 hydrogen bonds per frame.

Tepotinib showed recurring interactions involving GLU29 and CYS108
in the analysed trajectory.

### 6.6 Residue contact occupancy

Shared residues at the 50% occupancy threshold:
VAL31, GLY32, TYR36, VAL39, ALA52, LYS54, PHE107, CYS108,
GLY111, SER112, LEU160, VAL170 and ASP171.

Ponatinib-specific residues:
LEU15, GLY109 and ASP115.

Tepotinib-specific residues:
SER16, GLU29, LEU30, ASN33, GLU69, MET105, GLU106 and ARG313.

These assignments describe the analysed trajectories and should
not be interpreted as proof of stronger binding.

## 7. Docking interaction analysis

PLIP analysis of representative docking poses identified:

Ponatinib:
- 10 interaction records
- 8 hydrophobic
- 1 hydrogen bond
- 1 halogen-bond assignment

Tepotinib:
- 20 interaction records
- 14 hydrophobic
- 5 hydrogen bonds
- 1 salt-bridge assignment

Some geometric and chemical assignments require cautious
interpretation. In particular, several reported distances and
angles are borderline, and the Ponatinib fluorine-associated
halogen-bond assignment should not be treated as definitive.

Docking-pose interactions and MD trajectory contacts are distinct
analyses and should be reported separately.

## 8. Integrated interpretation

The project combines ligand-based activity modelling, structure-based
virtual screening and preliminary restrained MD analysis.

The ML model provides moderate predictive performance on the
available TNIK IC50 dataset. Docking identified several candidate
compounds with favourable predicted scores, while the MD analysis
describes short-timescale behaviour of selected restrained complexes.

Ponatinib has additional reported TNIK-associated Kd evidence,
but assay context differs between sources. Tepotinib showed recurring
hydrogen-bond observations in the analysed MD trajectory.

These computational findings support prioritisation for further
investigation. They do not establish clinical efficacy, therapeutic
benefit in IPF, or definitive TNIK inhibition.

## 9. Important limitations

1. The activity dataset contains heterogeneous assay records and
   only 453 unique compounds after filtering.
2. Random Forest performance is moderate and should not be
   interpreted as highly accurate prediction.
3. Docking scores are approximate scoring-function outputs,
   not experimental binding free energies.
4. Some experimental affinity records originate from different
   assay contexts.
5. MD equilibration was limited to 500 ps NVT and 1 ns NPT.
6. The complex simulations used heavy-atom restraints.
7. The observed RMSD, contact and hydrogen-bond patterns do not
   demonstrate long-term unrestrained stability.
8. Further experimental validation is required to establish
   TNIK inhibition and relevance to IPF.

## 10. Supporting project records

The following records have been created and audited:

- `09_results/molecular_dynamics/TNIK_GAFF2_parameter_validation.md`
- `09_results/molecular_dynamics/TNIK_MD_reproducibility_record.md`
- `09_results/molecular_dynamics/TNIK_MD_results_index.md`

The MD results index audit found 93 indexed paths and zero missing files.

## 11. Remaining project finalisation

- Visually inspect the seven individual MD heatmaps.
- Select report- and presentation-quality figures.
- Finalise figure captions.
- Integrate selected figures into the manuscript and presentation.
- Review all source tables and outputs before final submission.


## 3. Machine-learning results


# TNIK Random Forest – Final Evaluation Metrics

Metrics calculated from saved evaluation CSVs and test predictions.

| Evaluation | Metric | Mean | SD | Folds/runs |
|---|---:|---:|---:|---:|
| High-confidence random 5-fold CV | R2 | 0.5044 | 0.1819 | 5 |
| High-confidence random 5-fold CV | RMSE | 0.8054 | 0.1217 | 5 |
| High-confidence random 5-fold CV | MAE | 0.6023 | 0.0791 | 5 |
| High-confidence scaffold 5-fold CV | R2 | 0.3725 | 0.2745 | 5 |
| High-confidence scaffold 5-fold CV | RMSE | 0.8523 | 0.1351 | 5 |
| High-confidence scaffold 5-fold CV | MAE | 0.6693 | 0.0873 | 5 |
| General scaffold 5-fold CV | R2 | 0.4024 | 0.1559 | 5 |
| General scaffold 5-fold CV | RMSE | 0.7502 | 0.0640 | 5 |
| General scaffold 5-fold CV | MAE | 0.5630 | 0.0618 | 5 |
| Held-out test set | R2 | 0.4914 | — | 91 |
| Held-out test set | RMSE | 0.7976 | — | 91 |
| Held-out test set | MAE | 0.5917 | — | 91 |
| Y-randomisation (100 runs) | Mean_R2 | -0.2142 | 0.0573 | 100 |

## Y-randomisation interpretation

- Randomisation runs: 100
- Mean randomised R²: -0.2142
- Standard deviation: 0.0573
- Minimum R²: -0.3714
- Maximum R²: -0.0774
- Observed held-out test R²: 0.4914
- Randomised runs with R² at least as high as observed: 0
- Empirical permutation p-value: 0.0099

## Test prediction error distribution

- Test compounds: 91
- Median absolute error: 0.3857
- Maximum absolute error: 2.3396
- Compounds with absolute error below 0.5: 51
- Compounds with absolute error below 1.0: 73

## Interpretation note

General scaffold CV and high-confidence scaffold CV are reported separately. The high-confidence analyses must not be presented as the evaluation of the full dataset. Y-randomisation is a model validation check, not a measure of predictive accuracy.


### 3.1 Machine-learning interpretation

The held-out test-set R² is approximately 0.49, indicating moderate predictive performance on the selected test compounds. General scaffold cross-validation gives a mean R² of approximately 0.40, while high-confidence scaffold cross-validation gives a mean R² of approximately 0.37 with substantial fold-to-fold variability.

The Y-randomisation results provide a model-validation check: the 100 randomised runs have a mean R² below zero, with an empirical p-value of 0.0099 as calculated in the saved metrics table.

These results support cautious use of the model for prioritisation. They do not establish reliable activity prediction for every new chemical structure.

## 4. Molecular-docking results


# TNIK Docking and Experimental Evidence Summary

| Evidence | Compound | Result | Additional_information | Interpretation |
| --- | --- | --- | --- | --- |
| Crystal ligand redocking | JWK | -9.770 kcal/mol | Heavy-atom RMSD = 0.782 Å | Redocking pose recovery; supports the docking protocol |
| Virtual screening | Ponatinib | -11.21 kcal/mol | Initial screening score | Predicted docking score; not experimental affinity |
| Virtual screening | Tepotinib | -11.15 kcal/mol | Initial screening score | Predicted docking score; not experimental affinity |
| Repeated-seed docking | Ponatinib | -11.166 kcal/mol | 5 seeds; SD = 0.027 kcal/mol; representative seed 11 | Consistent docking score across tested seeds |
| Repeated-seed docking | Tepotinib | -11.128 kcal/mol | 5 seeds; SD = 0.041 kcal/mol; representative seed 22 | Consistent docking score across tested seeds |
| Experimental bioactivity | Ponatinib | Kd = 385 nM | ChEMBL Kinobeads cell-lysate competition assay | Binding evidence; not proof of purified-enzyme inhibition |
| Experimental bioactivity | Tepotinib | Kd > 30,000 nM | ChEMBL-reported result | Weak or undetected binding under reported assay conditions |
| Experimental bioactivity | Entrectinib | Kd > 30,000 nM | ChEMBL-reported result | Weak or undetected binding under reported assay conditions |
| Experimental bioactivity | Nilotinib | Kd > 10,000 / > 30,000 nM | Reported values from different records | Assay-dependent evidence; not a confirmed TNIK inhibitor |

## Interpretation and limitations

- JWK redocking produced a score of -9.770 kcal/mol and heavy-atom RMSD of 0.782 Å.
- Ponatinib and Tepotinib showed similar repeated-seed docking scores. Their small score difference does not establish a potency ranking.
- Docking scores are computational predictions, not measured binding affinities.
- The Ponatinib Kinobeads Kd result comes from a cell-lysate competition assay and does not establish purified-TNIK inhibition.
- These results do not establish therapeutic efficacy against IPF.


### 4.1 Docking interpretation

The JWK redocking result (heavy-atom RMSD 0.782 Å) supports recovery of the crystallographic ligand pose under the applied docking setup.

Ponatinib and Tepotinib have similar repeated-seed docking scores. The difference between their mean scores is small and should not be interpreted as a meaningful potency ranking.

The reported Ponatinib Kd of 385 nM comes from a Kinobeads cell-lysate competition assay. It is evidence of TNIK-related binding in that assay context, not proof of inhibition of purified TNIK. Docking scores are predictions and are not experimental binding affinities.

## 5. Molecular-dynamics results


# TNIK Molecular-Dynamics Final Evidence Table

| Analysis | System_or_residue | Primary_result | Additional_result | Interpretation_and_notes |
| --- | --- | --- | --- | --- |
| Simulation protocol | All systems | NVT: 500 ps; NPT: 1 ns | 2 fs timestep; heavy-atom restraints | Restrained equilibration, not production MD |
| NPT equilibration | All systems | Approximately 300 K | Initial 0–200 ps volume/density transient; later values relatively stable | Pressure fluctuations remained large; equilibration is preliminary |
| Protein backbone RMSD | Ponatinib N32-H | Mean 0.0220 nm; SD 0.0022 nm | Final 0.0209 nm | Low RMSD under heavy-atom restraints |
| Protein backbone RMSD | Ponatinib N37-H | Mean 0.0223 nm; SD 0.0023 nm | Final 0.0219 nm | Low RMSD under heavy-atom restraints |
| Protein backbone RMSD | Tepotinib | Mean 0.0220 nm; SD 0.0022 nm | Final 0.0224 nm | Low RMSD under heavy-atom restraints |
| Ligand heavy-atom RMSD | Ponatinib N32-H | Mean 0.0317 nm; SD 0.0059 nm | Maximum 0.0458 nm; final 0.0395 nm | Protein-backbone-aligned; restraints affect interpretation |
| Ligand heavy-atom RMSD | Ponatinib N37-H | Mean 0.0295 nm; SD 0.0049 nm | Maximum 0.0401 nm; final 0.0276 nm | Protein-backbone-aligned; restraints affect interpretation |
| Ligand heavy-atom RMSD | Tepotinib | Mean 0.0291 nm; SD 0.0056 nm | Maximum 0.0413 nm; final 0.0291 nm | Protein-backbone-aligned; restraints affect interpretation |
| Protein-ligand contacts (0.35 nm) | Ponatinib N32-H | Mean 31.198; SD 1.265 | Range 28–34; final 30 | Geometric contacts, not binding affinity |
| Protein-ligand contacts (0.35 nm) | Ponatinib N37-H | Mean 32.010; SD 1.473 | Range 28–35; final 33 | Geometric contacts, not binding affinity |
| Protein-ligand contacts (0.35 nm) | Tepotinib | Mean 33.376; SD 0.926 | Range 30–35; final 34 | Geometric contacts, not binding affinity |
| Mean minimum contact distance | Ponatinib N32-H | 0.2438 nm | Final 0.2478 nm | Minimum atom-pair distance |
| Mean minimum contact distance | Ponatinib N37-H | 0.2388 nm | Final 0.2423 nm | Minimum atom-pair distance |
| Mean minimum contact distance | Tepotinib | 0.1926 nm | Final 0.2058 nm | Minimum atom-pair distance |
| Hydrogen bonds | Ponatinib N32-H | Zero detected | Under the applied GROMACS criteria | Does not exclude other interaction types |
| Hydrogen bonds | Ponatinib N37-H | Zero detected | Under the applied GROMACS criteria | Does not exclude other interaction types |
| Hydrogen bonds | Tepotinib | Mean 2.3168 per frame; SD 0.734 | Maximum 4; at least one in 100% of frames | Geometry-based hydrogen-bond analysis |
| Tepotinib hydrogen-bond occupancy | GLU29 | 79.21% | GLU29 OE2 to ligand N3-H14 | Residue numbering mapped to original PDB |
| Tepotinib hydrogen-bond occupancy | CYS108 | 83.17% | CYS108 N-H to ligand O2 | Residue numbering mapped to original PDB |
| Shared contact residues (>=50%) | All three systems | 13 residues | VAL31, GLY32, TYR36, VAL39, ALA52, LYS54, PHE107, CYS108, GLY111, SER112, LEU160, VAL170, ASP171 | Residue-level contact occupancy |
| Ponatinib-specific contact residues (>=50%) | Ponatinib | 3 residues | LEU15, GLY109, ASP115 | Specific to the analysed Ponatinib systems |
| Tepotinib-specific contact residues (>=50%) | Tepotinib | 8 residues | SER16, GLU29, LEU30, ASN33, GLU69, MET105, GLU106, ARG313 | Specific to the analysed Tepotinib system |
| Ponatinib protonation comparison | N32-H versus N37-H | Different residue-contact occupancies | Examples: GLN157 37.62% vs 6.93%; MET105 37.62% vs 10.89%; ARG313 40.59% vs 48.51% | Protonation state influences observed contact patterns |

## Overall interpretation

The three restrained MD systems maintained low protein-backbone and ligand RMSD values during the simulated equilibration period. Tepotinib showed recurring hydrogen bonds with GLU29 and CYS108. Ponatinib systems had no hydrogen bonds detected under the applied GROMACS criteria. Contact counts and residue occupancies describe observed geometry and should not be interpreted as binding affinity.

## Important limitations

- The simulations comprised 500 ps NVT and 1 ns NPT equilibration.
- Heavy-atom restraints were applied; low RMSD values therefore do not demonstrate stability in an unrestrained production simulation.
- NPT volume and density showed an initial transient.
- Pressure fluctuations were substantial.
- No binding free-energy calculation was performed.
- These MD observations do not establish TNIK inhibition or therapeutic efficacy against idiopathic pulmonary fibrosis.


### 5.1 MD interpretation and limitations

The analysed systems underwent 500 ps NVT and 1 ns NPT equilibration with heavy-atom restraints. The low protein and ligand RMSD values must therefore be interpreted in the context of those restraints.

The NPT analysis showed an initial volume and density transient, and pressure fluctuations remained substantial. These simulations are preliminary restrained equilibration runs; they are not long, unrestrained production simulations.

Tepotinib showed recurring hydrogen bonds with GLU29 and CYS108 under the applied analysis criteria. No hydrogen bonds were detected for the two Ponatinib protonation-state systems under those same criteria. Contact counts, minimum distances and residue occupancies describe geometry and do not directly quantify binding affinity.

The comparison of Ponatinib N32-H and N37-H shows that the selected protonation state can affect residue-contact patterns. These observations do not establish which protonation state is dominant in the protein-bound state.

## 6. Supporting reproducibility and parameterisation records

### 6.1 MD reproducibility record


# TNIK–Ligand Molecular Dynamics Reproducibility Record

**Record generated:** 2026-11-07

## 1. Systems investigated

Three protein–ligand systems were investigated:

1. TNIK–Ponatinib N32-H
2. TNIK–Ponatinib N37-H
3. TNIK–Tepotinib +1

An apo-TNIK system was also prepared and simulated separately.

## 2. Molecular structures

- Protein crystal structure: PDB ID 6RA7
- Protein chain: A
- Original protein sequence: 302 residues
- Protein sequence mapping: validated against the original PDB chain A
- Ligand protonation candidates: prepared for nominal pH 7.4

The validated residue mapping was used to relate GROMACS residue numbering to original PDB residue numbering.

## 3. Force fields and molecular system preparation

| Component | Method or force field |
|---|---|
| Protein | AMBER ff14SB |
| Phosphorylated residues | phosaa14SB |
| Ligands | GAFF2 |
| Water | TIP3P |
| Ligand parameterisation | AmberTools / ACPYPE |
| Long-range electrostatics | Particle Mesh Ewald (PME) |

The systems were solvated using a 12 Å water buffer and neutralised with counterions.

## 4. Ligand parameter validation

Three ligand parameter sets were examined:

| Ligand | Protonated nitrogen | Total MOL2 charge |
|---|---|---:|
| Tepotinib +1 | N3 | +0.9990 |
| Ponatinib N32-H | N5 | +1.0000 |
| Ponatinib N37-H | N6 | +0.9980 |

All three MOL2 files were processed successfully using parmchk2 with GAFF2.

The regenerated FRCMOD files were compared with the original FRCMOD files. All three comparisons produced exact content matches.

The validation record is available at:
`09_results/molecular_dynamics/TNIK_GAFF2_parameter_validation.md`

## 5. Energy minimisation

Energy minimisation was performed using the steepest-descent algorithm.

| System | Steps | Final potential energy (kJ/mol) | Maximum force (kJ/mol/nm) |
|---|---:|---:|---:|
| Ponatinib N32-H | 1263 | -707420.69 | 984.86 |
| Ponatinib N37-H | 1364 | -709737.69 | 885.73 |
| Tepotinib +1 | 734 | -691523.31 | 933.75 |

All three systems reached the specified maximum-force criterion of less than 1000 kJ/mol/nm.

## 6. NVT equilibration

| Parameter | Setting |
|---|---|
| Ensemble | NVT |
| Integrator | Molecular dynamics |
| Time step | 0.002 ps |
| Number of steps | 250000 |
| Simulation duration | 500 ps |
| Target temperature | 300 K |
| Thermostat | V-rescale |
| Thermostat time constant | 1.0 ps |
| Pressure coupling | Disabled |
| Bond constraints | Hydrogen-containing bonds |
| Electrostatics | PME |
| Position restraints | Protein and ligand heavy atoms |

All three complexes completed the 500 ps NVT simulations.

The trajectories contained 501 frames at 1 ps intervals.

## 7. NPT equilibration

| Parameter | Setting |
|---|---|
| Ensemble | NPT |
| Time step | 0.002 ps |
| Number of steps | 500000 |
| Simulation duration | 1000 ps (1 ns) |
| Target temperature | 300 K |
| Thermostat | V-rescale |
| Thermostat time constant | 1.0 ps |
| Barostat | C-rescale |
| Pressure coupling | Isotropic |
| Reference pressure | 1 bar |
| Pressure time constant | 5 ps |
| Compressibility | 4.5 × 10⁻⁵ bar⁻¹ |
| Bond constraints | Hydrogen-containing bonds |
| Electrostatics | PME |
| Position restraints | Continued from the restrained equilibration setup |

All three systems completed their 1 ns NPT simulations.

The compressed trajectories contained 101 frames at 10 ps intervals, spanning 0–1000 ps.

## 8. NPT stability observations

The average temperatures remained close to 300 K.

The initial 0–200 ps block showed box contraction and an increase in density. Later 200 ps blocks showed comparatively stable volumes and densities.

Approximate volume and density ranges across the later blocks were:

| System | Later-block volume range (nm³) | Later-block density range (kg/m³) |
|---|---:|---:|
| Ponatinib N32-H | 493.827–494.049 | 1007.237–1007.685 |
| Ponatinib N37-H | 493.499–493.851 | 1007.638–1008.356 |
| Tepotinib +1 | 493.389–494.186 | 1007.061–1008.688 |

These observations indicate that the simulation boxes became comparatively stable after the initial volume adjustment.

## 9. Structural and interaction analyses

The following analyses were performed on the NPT trajectories:

- Protein backbone RMSD
- Ligand heavy-atom RMSD after protein-backbone alignment
- Protein–ligand contacts using a 0.35 nm cutoff
- Hydrogen-bond analysis
- Residue-level contact occupancy
- PLIP interaction profiling of representative docking poses

The corresponding analysis outputs and figures are stored in:
`06_docking/molecular_dynamics/analysis/npt_stability/`

PLIP comparison outputs are stored in:
`06_docking/interaction_analysis/comparative_analysis/`

## 10. Software

| Software | Recorded version |
|---|---|
| GROMACS | 2026.3 |
| AmberTools | 26 |
| ACPYPE | 2026.10.15 |

## 11. Important methodological limitations

1. The NPT simulations were 1 ns in duration and should be treated as preliminary.
2. Heavy-atom position restraints were applied during the reported equilibration simulations.
3. Low RMSD values in restrained trajectories do not establish unrestrained binding stability.
4. Contact counts and docking scores are not direct measurements of binding affinity.
5. Ligand protonation states represent computational assumptions.
6. The simulations do not establish TNIK inhibition or therapeutic efficacy against IPF.
7. Total potential energies from different systems should not be used to rank ligand affinity.

## 12. Overall record

The three TNIK–ligand systems completed energy minimisation, 500 ps NVT equilibration and 1 ns NPT simulation.

The ligand parameter files passed the reported GAFF2 parameter-coverage and FRCMOD content-comparison checks.

The resulting trajectories were used for structural and interaction analyses. All findings should be interpreted within the limitations of the simulation duration and applied position restraints.


### 6.2 GAFF2 parameter validation


# TNIK Ligand GAFF2 Parameter Validation

**Validation date:** 2026-11-07

## 1. Objective

To assess the consistency of the GAFF2 ligand parameter files
used in the three TNIK–ligand molecular dynamics systems.

## 2. Ligand systems examined

| Ligand system | MOL2 file |
|---|---|
| Tepotinib +1 | `tepotinib_pH74_plus1/tepotinib_pH74_plus1.mol2` |
| Ponatinib N32-H | `ponatinib_N32_H/ponatinib_N32_H.mol2` |
| Ponatinib N37-H | `ponatinib_N37_H/ponatinib_N37_H.mol2` |

## 3. Protonation-site observations

| Ligand system | Protonated nitrogen | GAFF2 type | Total MOL2 charge |
|---|---|---|---:|
| Tepotinib +1 | N3 | nx | +0.9990 |
| Ponatinib N32-H | N5 | nx | +1.0000 |
| Ponatinib N37-H | N6 | nx | +0.9980 |

The additional hydrogen was observed bonded to the specified
nitrogen atom in each MOL2 structure.

## 4. parmchk2 validation

| Ligand system | Exit code | Original FRCMOD size | Check FRCMOD size |
|---|---:|---:|---:|
| Tepotinib +1 | 0 | 4484 bytes | 4484 bytes |
| Ponatinib N32-H | 0 | 4656 bytes | 4656 bytes |
| Ponatinib N37-H | 0 | 4304 bytes | 4304 bytes |

All three parmchk2 executions completed successfully.

## 5. FRCMOD content comparison

The original FRCMOD files were compared line by line with the
corresponding files regenerated using parmchk2.

| Ligand system | Comparison result |
|---|---|
| Tepotinib +1 | Exact content match |
| Ponatinib N32-H | Exact content match |
| Ponatinib N37-H | Exact content match |

## 6. Conclusion

All three ligand parameter files passed the reported
parmchk2 processing and FRCMOD content-comparison checks.
The regenerated parameter files were identical in content
to the original files.

These checks support the internal consistency of the
existing GAFF2 parameter files. They do not independently
prove the accuracy of every force-field parameter or
validate the biological protonation state.

No original ligand parameter files or simulation outputs
were modified during these validation checks.


## 7. Integrated interpretation

The evidence should be considered in distinct layers:

1. **Machine learning:** The baseline model has moderate predictive performance. Scaffold-based validation is more challenging and shows variability, particularly in the high-confidence subset.
2. **Docking:** JWK redocking supports the docking setup. Ponatinib and Tepotinib have similar docking scores, so docking alone does not distinguish them reliably.
3. **Experimental bioactivity:** Ponatinib has a reported TNIK-related Kd result in a cell-lysate competition assay. The assay context must be retained when discussing this result.
4. **Molecular dynamics:** The short restrained equilibration runs provide preliminary structural observations. They do not establish long-term stability or binding free energy.
5. **Interaction analysis:** PLIP and trajectory-based contacts provide complementary geometric descriptions. Automated interaction assignments require chemical and geometric interpretation.

Taken together, the current results can support prioritisation of candidates for further experimental investigation. They do not demonstrate that Ponatinib or Tepotinib is an effective TNIK inhibitor or an IPF treatment.

## 8. Important project-wide limitations

- The ChEMBL activity dataset contains a limited number of unique compounds and heterogeneous assay records.
- Machine-learning performance is moderate and may be affected by dataset size, chemical-series distribution and assay heterogeneity.
- Docking scores are approximate scoring-function outputs and are not equivalent to measured binding constants.
- The Ponatinib Kinobeads result is from a cell-lysate competition assay rather than a purified-TNIK enzyme assay.
- MD simulations are short and restrained; longer unrestrained production simulations and additional validation would be needed.
- No binding free-energy calculation has been included in this consolidated evidence.
- Computational binding predictions do not establish target inhibition, cellular efficacy, safety or therapeutic benefit in IPF.

## 9. Suggested future validation

1. Confirm TNIK inhibition using a purified-protein biochemical assay.
2. Evaluate concentration-response behaviour and determine an experimental IC50 or Ki.
3. Confirm cellular target engagement and downstream pathway effects.
4. Assess selectivity against related kinases.
5. Conduct longer, appropriately equilibrated unrestrained MD simulations if justified.
6. Consider binding free-energy calculations only after the structural models and simulation protocols are sufficiently validated.
7. Evaluate relevant ADMET and toxicity properties before considering further development.

## 10. Evidence file locations

The accompanying CSV index lists the consolidated evidence files and selected supporting records. Original analysis outputs remain in their existing project folders.

---

**End of master results report.**
