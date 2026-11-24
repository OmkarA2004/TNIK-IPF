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
