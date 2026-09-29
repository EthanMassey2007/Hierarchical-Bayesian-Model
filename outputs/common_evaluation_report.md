# Corrected common-sample evaluation

Run: `20260929T093220_37466`

All 17 models were trained on available 2017–2022 observations and scored on the same **324 municipality-week observations** (81 municipalities), dated **2023-01-02 to 2023-01-23**.

The previous cross-model predictive metrics used unequal observation sets and must be replaced. Model-specific training samples still differ; this comparison evaluates fitted pipelines, rather than isolating the effect of a feature on an identical training cohort.

## Replacement predictive metrics

| Model | MAE | RMSE | WAPE | R² | N |
|---|---:|---:|---:|---:|---:|
| M0 | 6.5100 | 16.8722 | 0.6989 | 0.7020 | 324 |
| M1 | 6.5922 | 17.8408 | 0.7077 | 0.6667 | 324 |
| M2 | 6.7220 | 20.0390 | 0.7216 | 0.5796 | 324 |
| M3 | 6.5980 | 17.7953 | 0.7083 | 0.6684 | 324 |
| M4 | 4.5590 | 12.4284 | 0.4894 | 0.8383 | 324 |
| M5 | 4.3084 | 11.5998 | 0.4625 | 0.8591 | 324 |
| S1 | 4.3262 | 11.6686 | 0.4644 | 0.8574 | 324 |
| S2 | 4.1431 | 10.8601 | 0.4448 | 0.8765 | 324 |
| S3 | 4.2967 | 11.6404 | 0.4613 | 0.8581 | 324 |
| S4 | 4.1521 | 10.9162 | 0.4457 | 0.8752 | 324 |
| S5 | 4.1375 | 10.8601 | 0.4442 | 0.8765 | 324 |
| S6 | 4.0752 | 10.5736 | 0.4375 | 0.8829 | 324 |
| S7 | 4.1923 | 11.2267 | 0.4501 | 0.8680 | 324 |
| S8 | 4.1095 | 10.7301 | 0.4412 | 0.8795 | 324 |
| S9 | 4.3308 | 11.6376 | 0.4649 | 0.8582 | 324 |
| S10 | 4.7120 | 12.9211 | 0.5059 | 0.8252 | 324 |
| S11 | 4.4701 | 11.8982 | 0.4799 | 0.8518 | 324 |

Lowest common-sample WAPE: **S6**. Lowest common-sample RMSE: **S6**.

## Information-criterion comparison groups

DIC and WAIC came from separate full-data descriptive fits. Compare these criteria only within groups having identical fitted response observations.

- fit_1: R_M0; 33488 fitted observations.
- fit_2: R_M1; 31588 fitted observations.
- fit_3: R_M2; 30393 fitted observations.
- fit_4: R_M3; 32111 fitted observations.
- fit_5: R_M4; 31143 fitted observations.
- fit_6: R_M5, S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11; 30302 fitted observations.

## Sample exclusions

| Model | Eligible January rows | Shared rows | Excluded |
|---|---:|---:|---:|
| R_M0 | 368 | 324 | 44 |
| R_M1 | 324 | 324 | 0 |
| R_M2 | 365 | 324 | 41 |
| R_M3 | 324 | 324 | 0 |
| R_M4 | 324 | 324 | 0 |
| R_M5 | 365 | 324 | 41 |
| S1 | 365 | 324 | 41 |
| S2 | 365 | 324 | 41 |
| S3 | 365 | 324 | 41 |
| S4 | 365 | 324 | 41 |
| S5 | 365 | 324 | 41 |
| S6 | 365 | 324 | 41 |
| S7 | 365 | 324 | 41 |
| S8 | 365 | 324 | 41 |
| S9 | 365 | 324 | 41 |
| S10 | 365 | 324 | 41 |
| S11 | 365 | 324 | 41 |

## Manuscript implications

- Replace all predictive metrics and rankings using the table above.
- Describe a four-week January 2023 evaluation, not a full-year or two-year holdout.
- Do not compare DIC/WAIC across different fitted-response groups.
- Keep full-data rainfall-effect maps and fitted-value figures distinct from held-out prediction.
- Models using contemporaneous weather condition on observed weather; this is not an operational weather forecast evaluation.
- Broader seasonal and epidemic generalization requires a separate rolling-origin evaluation.
