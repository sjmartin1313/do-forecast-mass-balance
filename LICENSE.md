# License

This repository is licensed under the **Creative Commons Attribution 3.0 Unported License (CC BY 3.0)**.

You may share and adapt the materials in this repository, including the code and documentation, provided that you give appropriate credit, provide a link to the license, and indicate whether changes were made.

License text: https://creativecommons.org/licenses/by/3.0/

SPDX-License-Identifier: CC-BY-3.0

## Attribution

This repository includes a modified version of code from the original **BASEmetab** R package:

- Original project: **BASEmetab**
- Original repository: https://github.com/dgiling/BASEmetab
- Original authors: Darren Giling and Ralph Mac Nally
- Original license: Creative Commons Attribution 3.0 Unported License (CC BY 3.0)
- Associated publication: Grace, M. R., Giling, D. P., Hladyz, S., Caron, V., Thompson, R. M., and Mac Nally, R. (2015). *Fast processing of diel oxygen curves: estimating stream metabolism with BASE (BAyesian Single-station Estimation).* Limnology and Oceanography: Methods, 13, 103–114. https://doi.org/10.1002/lom3.10011

The modified BASEmetab function included in this repository is adapted from the original `bayesmetab` function in BASEmetab and is distributed under the same CC BY 3.0 license.

## Modifications to BASEmetab

The modified script, `bayesmetab_mod.R`, was adapted for the project:

**Short-Term Dissolved Oxygen Forecasting in Aquaculture Systems Using a Process-Based Mass-Balance Model**

Major modifications include:

- Added `R.mean`, `R.sd`, and `R.median` to the output table.
- Added `DO.meas` to the monitored JAGS model parameters.
- Added safer extraction of model outputs and Rhat values.
- Added validation checks before extracting model outputs.
- Added fallback output rows when model output is invalid or incomplete.
- Wrapped result generation and model-fit metric calculations in `tryCatch` blocks.
- Ensured missing output columns are added as `NA` rather than causing the full run to fail.
- Replaced direct indexing of model outputs with controlled validation to reduce dimension mismatch errors.
- Replaced the original instantaneous output logic with mean-based aggregation.

## Citation Request

If you use this repository, please cite both the original BASEmetab work and this modified repository/project, where appropriate.

Recommended citation for the original method:

Grace, M. R., Giling, D. P., Hladyz, S., Caron, V., Thompson, R. M., and Mac Nally, R. (2015). Fast processing of diel oxygen curves: estimating stream metabolism with BASE (BAyesian Single-station Estimation). *Limnology and Oceanography: Methods*, 13, 103–114. https://doi.org/10.1002/lom3.10011

## No Warranty

The materials in this repository are provided as-is, without warranties or conditions of any kind. See the Creative Commons Attribution 3.0 Unported License for the full license terms.
