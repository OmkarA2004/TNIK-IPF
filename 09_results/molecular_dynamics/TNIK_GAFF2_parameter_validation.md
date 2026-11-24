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
