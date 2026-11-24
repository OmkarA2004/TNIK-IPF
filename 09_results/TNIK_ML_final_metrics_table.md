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
