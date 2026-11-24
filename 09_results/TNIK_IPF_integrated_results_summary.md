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
