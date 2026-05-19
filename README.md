# Short-Term Dissolved Oxygen Forecasting with Modified BASEmetab

This repository contains code used to estimate daily dissolved oxygen (DO) metabolism parameters with a modified version of `BASEmetab` and then use those parameters to generate 24-hour DO predictions with a process-based Monte Carlo model.

The workflow was developed for the project:

> **Short-Term Dissolved Oxygen Forecasting in Aquaculture Systems Using a Process-Based Mass-Balance Model**  
> Sonny, Martin, Joseph Dvorak, Ken Semmens and Bill Ford

## Project background

Dissolved oxygen is one of the most important water-quality variables in aquaculture systems. Low DO events can stress fish, reduce growth, or cause mortality, especially during the early morning when photosynthesis has not yet resumed and respiration continues to remove oxygen from the water. A useful management tool should therefore be able to warn producers during the day when DO is likely to fall to critical levels overnight.

This project uses a process-based modeling approach instead of a site-specific machine learning model. The main idea is to estimate metabolism and reaeration coefficients from the previous 24 hours of measurements, then use those coefficients to predict the next 24 hours of DO at 10-minute resolution.

The model uses:

- measured DO,
- water temperature,
- photosynthetically active radiation (PAR),
- atmospheric pressure,
- salinity,
- daily metabolism coefficients,
- and a daily reaeration coefficient.

In the manuscript, future DO is computed from a mass-balance equation using coefficients estimated by `BASEmetab`. The prediction period is treated as a 24-hour management window, such as 18:00 to 18:00, rather than necessarily a midnight-to-midnight calendar day. The goal is to support early decisions about whether additional aeration may be needed before an overnight DO drop.

For model testing, observed water temperature, atmospheric pressure, PAR, and salinity were used in place of forecasted inputs. This was done to isolate errors caused by the model structure and coefficient selection before adding additional uncertainty from weather and water-temperature forecasts. Using true forecasted inputs is the next intended step for field deployment.

## Repository contents

```text
.
├── bayesmetab_mod.R
├── DO Prediction.py
├── Measurement Datasets/
│   └── Final_Format_P7T.csv                  # example expected location/name; user-provided
├── Metabolism Datasets/
│   └── Final_Format_P7T_Metabolism.csv       # example expected location/name; generated/processed from BASEmetab output
└── Results/
    ├── DO_predictions_montecarlo.csv
    ├── daily_r2_montecarlo.csv
    └── DO_comparison_<period_id>.png
```

Only the scripts are included here. Measurement and metabolism datasets are expected to be supplied by the user. Data from the paper may be available upon request.

## Script descriptions

### `bayesmetab_mod.R`

`bayesmetab_mod.R` is a modified version of the original `bayesmetab` function from the `BASEmetab` R package. It estimates single-station metabolism and reaeration coefficients from diel DO curves using Bayesian MCMC through JAGS.

The modified function keeps the same general purpose as the original `bayesmetab` function, but it adds output fields and error handling needed for the DO prediction workflow.

The input data files must be CSV files containing the following case-sensitive columns:

```text
Date
Time
I
tempC
DO.meas
atmo.pressure
salinity
```

For 10-minute data, use:

```r
interval = 600
```

Each modeled 24-hour period should contain 144 rows when using a 600-second interval.

### `DO Prediction.py`

`DO Prediction.py` loads measurement data and metabolism-parameter data, aligns each measurement with the correct 24-hour metabolism period, samples daily coefficients from normal distributions, and recursively computes future DO using the mass-balance model.

For each prediction period, the script:

1. uses the previous period's metabolism-parameter means and standard deviations,
2. samples `K`, `theta`, `R`, `A`, and `p`,
3. converts `K` from a daily scale to a 10-minute timestep scale,
4. uses the final measured DO from the previous period as the initial condition,
5. simulates 1000 Monte Carlo DO trajectories,
6. saves the median predicted DO trajectory,
7. calculates daily performance metrics,
8. and creates measured-versus-predicted DO plots.

## Installation

### R requirements

Install R, JAGS, and the required R packages.

```r
install.packages(c("remotes", "R2jags", "zoo", "coda", "lattice"))

remotes::install_github("dgiling/BASEmetab")
```

You also need JAGS installed on your computer because `R2jags` calls JAGS during parameter estimation.

### Python requirements

The Python script was written for Python 3.11. Install the required packages with:

```bash
pip install numpy numba pandas matplotlib scikit-learn
```

## How to run the modified BASEmetab script

1. Install the original `BASEmetab` package and R dependencies.

2. Place your formatted input CSV files in a folder such as:

```text
Measurement Datasets/
```

Each CSV should contain complete 24-hour periods and the required columns:

```text
Date, Time, I, tempC, DO.meas, atmo.pressure, salinity
```

3. In R, source the modified script and run `bayesmetab_mod()`:

```r
library(BASEmetab)
library(R2jags)

source("bayesmetab_mod.R")

data.dir <- file.path(getwd(), "Measurement Datasets")
results.dir <- file.path(getwd(), "Metabolism Datasets")

if (!dir.exists(results.dir)) {
  dir.create(results.dir)
}

results <- bayesmetab_mod(
  data.dir = data.dir,
  results.dir = results.dir,
  interval = 600,
  n.iter = 20000,
  n.burnin = 10000,
  K.est = TRUE,
  p.est = FALSE,
  theta.est = FALSE,
  instant = FALSE
)
```

4. The modified BASEmetab output will be written to `results.dir`. The main output file is named similar to:

```text
BASE_results_<timestamp>.csv
```

5. Before using the output in the Python prediction script, make sure the metabolism dataset includes the period-start fields expected by `DO Prediction.py`:

```text
Start Date
Start Time
K.mean
K.sd
theta.mean
theta.sd
R.mean
R.sd
A.mean
A.sd
p.mean
p.sd
```

If your raw BASEmetab output uses a different date format, create `Start Date` and `Start Time` columns that define the beginning of each 24-hour period. For example, if the model periods are 18:00 to 18:00, `Start Time` should be `18:00:00`.

## How to run the DO prediction script

1. Put the measurement dataset in:

```text
Measurement Datasets/
```

2. Put the metabolism-parameter dataset in:

```text
Metabolism Datasets/
```

3. Open `DO Prediction.py` and update the file paths at the bottom of the script if needed:

```python
measurement_data = 'Measurement Datasets/Final_Format_P7T.csv'
metabolism_data = 'Metabolism Datasets/Final_Format_P7T_Metabolism.csv'
```

4. Run the script from the repository root:

```bash
python "DO Prediction.py"
```

The filename contains a space, so keep the quotation marks. Alternatively, rename the file to `do_prediction.py` and run:

```bash
python do_prediction.py
```

5. The script creates a `Results/` folder and writes:

```text
Results/DO_predictions_montecarlo.csv
Results/daily_r2_montecarlo.csv
Results/DO_comparison_<period_id>.png
```

## Expected input formats

### Measurement dataset

The measurement dataset should contain at least:

```text
Date
Time
I
tempC
DO.meas
atmo.pressure
salinity
```

The Python script combines `Date` and `Time` into a `datetime` column.

### Metabolism dataset

The metabolism dataset should contain at least:

```text
Start Date
Start Time
K.mean
K.sd
theta.mean
theta.sd
R.mean
R.sd
A.mean
A.sd
p.mean
p.sd
```

The Python script combines `Start Date` and `Start Time` into `period_start`, then assigns each measurement row to a 24-hour metabolism period.

## What changed from the original BASEmetab version?

This repository does not replace the original `BASEmetab` package. Instead, `bayesmetab_mod.R` is an adapted version of the original `bayesmetab` function designed to provide the additional outputs and robustness needed for this DO prediction workflow.

Main changes include:

- Added `R.mean`, `R.sd`, and `R.median` to the output table.
- Added `DO.meas` to the monitored JAGS parameters.
- Added safer extraction of model outputs and Rhat values.
- Added checks for valid model output structure before extracting results.
- Added fallback empty output rows when model output is invalid or incomplete.
- Ensured required output columns are present, with `NA` added where needed.
- Replaced direct indexing of model outputs with controlled validation.
- Wrapped R², RMSE, and related metric calculations in `tryCatch()` blocks.
- Modified instantaneous output handling to use mean-based aggregation.
- Preserved the original BASEmetab model call and JAGS model structure from the installed `BASEmetab` package.

These changes were made so that long batch runs would not fail completely when one day or one file produced incomplete or invalid model output.

## Notes on measured vs. forecasted inputs

The broader model is intended for short-term forecasting, where future water temperature, atmospheric pressure, PAR, and salinity would come from forecastable sources or coupled models. In the current testing workflow, measured future inputs are used in place of forecasted inputs.

This is intentional. It allows the model to be tested under idealized input conditions so that errors from the coefficient-estimation and mass-balance framework can be evaluated separately from errors caused by weather, PAR, or water-temperature forecasts. A future version of this workflow should replace the observed future inputs with forecasted inputs.

## Citation and attribution

This repository includes a modified version of code adapted from the original `BASEmetab` project:

> Giling, D. and Mac Nally, R. `BASEmetab`: R package for single-station stream metabolism estimation. GitHub repository: https://github.com/dgiling/BASEmetab

The original `BASEmetab` code accompanies:

> Grace, M. R., Giling, D. P., Hladyz, S., Caron, V., Thompson, R. M., and Mac Nally, R. (2015). Fast processing of diel oxygen curves: estimating stream metabolism with BASE (BAyesian Single-station Estimation). *Limnology and Oceanography: Methods*, 13, 103–114. https://doi.org/10.1002/lom3.10011

The modified `bayesmetab_mod.R` script retains attribution to the original authors and notes that the original code was released under the Creative Commons Attribution 3.0 license.

If you use this repository, cite both the original `BASEmetab`/BASE paper and the associated dissolved oxygen forecasting study.

## Suggested repository citation

```bibtex
@software{martin_do_prediction_basemetab,
  author = {Martin, Sonny},
  title = {Short-Term Dissolved Oxygen Forecasting with Modified BASEmetab},
  year = {2026},
  url = {https://github.com/sjmartin1313/do-forecast-mass-balance
}
```

Replace the URL with the final GitHub repository URL.

## License

The modified `bayesmetab_mod.R` script is adapted from `BASEmetab`, whose original code is identified in the script header as being licensed under Creative Commons Attribution 3.0 (CC BY 3.0). Any redistribution of this modified version should preserve the original attribution and license notice.
