# Quantile Trend Filtering Experiments

R experiments comparing quantile trend filtering estimates on simulated data.

## Contents

- `Fast_Gibbs_RMSE_Evaluation.R`: Gibbs sampling with P-MALA updates and an estimated lambda.
- `Individual_PMALA_RMSE_Evaluation.R`: individual P-MALA sampling with fixed lambda.
- `Quantile_TF_with_proximal_MCMC*.pdf`: accompanying reference documents.

Both scripts use seed `100`, simulate 1,000 observations, evaluate quantiles from 0.1 to 0.9, and save RMSE summaries as CSV files in the working directory.

## Requirements and execution

Requires R and the packages `glmgen`, `ggplot2`, `tidyr`, `dplyr`, and `zoo` (Gibbs script only).
From the repository directory on `owls.stat.rice.edu`, run:

```sh
Rscript Fast_Gibbs_RMSE_Evaluation.R
Rscript Individual_PMALA_RMSE_Evaluation.R
```

Outputs: `fast_gibbs_quantile_rmse_evaluation.csv` and `individual_pmala_quantile_rmse_evaluation.csv`.

## Status

Repository contents and script configuration have been inspected; runtime behavior and numerical accuracy have not been validated. The scripts retain their historical settings: Gibbs uses Gamma shape 9 and rate 3/n, while individual P-MALA fixes lambda at 1,000. Future comparisons should explicitly configure the current lambda-prior default, Gamma shape 1 and rate sqrt(n), while preserving historical snapshots.
