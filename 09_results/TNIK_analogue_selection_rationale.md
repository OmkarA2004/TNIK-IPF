# TNIK Analogue Selection Rationale

## Summary

Five BRICS-generated TNIK analogue candidates were evaluated using
random-forest activity predictions, applicability-domain assessment,
AutoDock Vina docking, distance-based interaction analysis and
preliminary ADMET predictions.

## Candidate interpretation

### TNIK_BRICS_0189 — structural/docking priority
- Docking score: -9.846 kcal/mol.
- Applicability domain: inside.
- Predicted pIC50: 7.6229.
- One plausible geometrically screened hydrogen bond with GLY32.
- Predicted AMES, hERG II and hepatotoxicity flags are positive.

### TNIK_BRICS_0112 — complementary ADMET priority
- Docking score: -9.759 kcal/mol.
- Applicability domain: inside.
- Predicted pIC50: 7.8244.
- Predicted hepatotoxicity: negative.
- Predicted AMES and hERG II flags remain positive.
- Geometric hydrogen-bond analysis has not yet been performed.

### TNIK_BRICS_0179 — close structural comparator
- Docking score: -9.750 kcal/mol.
- Applicability domain: inside.
- Predicted pIC50: 7.6365.
- One plausible geometrically screened hydrogen bond with GLY32.
- Predicted AMES, hERG II and hepatotoxicity flags are positive.

### TNIK_BRICS_0110 — activity-prediction comparator
- Highest predicted pIC50 among the four candidates inside the
  applicability domain: 7.8318.
- Docking score: -9.471 kcal/mol.
- Predicted AMES, hERG II and hepatotoxicity flags are positive.

### TNIK_BRICS_0001 — lower-confidence prediction
- Highest predicted pIC50 overall: 8.4840.
- Outside the applicability domain.
- Docking score: -9.029 kcal/mol.
- pkCSM toxicity predictions were not assessed for this candidate.
- Its activity prediction should be interpreted cautiously.

## Provisional prioritisation

1. BRICS_0189: structural/docking priority for further investigation.
2. BRICS_0112: complementary candidate because of its negative
   predicted hepatotoxicity result.
3. BRICS_0179: close structural comparator to BRICS_0189.
4. BRICS_0110: activity-prediction comparator.
5. BRICS_0001: lower-confidence, outside-domain candidate.

This is a prioritisation for further research, not a validated potency
or safety ranking.

## Important limitations

- Vina docking scores are computational estimates, not measured
  binding affinities. Small score differences should not be treated
  as meaningful potency differences.
- The GLY32 contacts for BRICS_0179 and BRICS_0189 passed the selected
  geometric screening thresholds. Each represents one unique
  donor–acceptor pair, not two independent hydrogen bonds.
- Other ligand–protein N/O distance contacts are not automatically
  hydrogen bonds.
- ADMET results are in-silico predictions and do not establish
  toxicity or safety.
- Positive AMES and hERG II predictions require follow-up assessment.
- Experimental TNIK inhibition, selectivity, cellular activity,
  pharmacokinetics and toxicity remain unverified.

## Recommended next research stage

Prioritise experimental or orthogonal computational validation of
BRICS_0189 and BRICS_0112. Retain BRICS_0179 as a close structural
comparator. Do not describe any candidate as a confirmed TNIK
inhibitor or safe drug without experimental evidence.
