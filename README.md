# The evolving epidemiology of scrub typhus in Thailand (2003–2024)

Code and data for *"The evolving epidemiology of scrub typhus in Thailand
(2003–2024): insights from latent process modelling of national surveillance
data"* (Wongnak et al., PLOS Neglected Tropical Diseases, under revision).

A Bayesian latent process model (Stan) jointly fits monthly reported cases and
deaths, their age distributions, and province-level seroprevalence to estimate
the force of infection, reporting proportion, age-specific susceptibility,
case-fatality and seasonality for all 77 provinces of Thailand.

## Repository structure

| Path | Contents |
|---|---|
| `Scripts/main_analysis.qmd` | Full workflow: data preparation, model fitting, sensitivity analysis, all figures and the figure source data |
| `Scripts/data_prep_functions.R` | Builds the Stan data from the input files |
| `Scripts/plot_functions.R` | Summaries, figures and `export_figure()` (PNG + PLOS TIFF) |
| `Scripts/run_models_local.R`, `Scripts/init_function.R` | Fits a prepared scenario: `Rscript --vanilla Scripts/run_models_local.R <scenario> <seed>` |
| `Stan_models/Model10.stan` | Final model |
| `Analysis_Data/` | Processed input data (surveillance, population, seroprevalence, map) |
| `Plots/` | Manuscript figures (`Fig1`–`Fig5`), supplementary figures (`FigureS1`–`FigureS5`), `tiff/` submission files, diagnostics |
| `Source_data/` | Numerical data underlying every figure (CSV per figure, and `S1_Data.xlsx` with a README sheet) |

## Reproducing the analysis

1. Install R (≥ 4.3) with `rstan`, `tidyverse`, `here`, `sf`, `ggspatial`,
   `ggh4x`, `cowplot`, `patchwork`, `biscale`, `ggrepel`, `ggalluvial`,
   `cluster`, `ragg`, `writexl` and `knitr`, and Quarto.
2. Open `Scripts/main_analysis.qmd`. Model fits are stored in `Rout/`, which is
   not tracked because of its size (about 1 GB per fit). To refit:
   - baseline model: set `RUN_MODEL <- TRUE` in the `model-fitting` chunk;
   - sensitivity analysis (clinical fraction and waning rate halved and
     doubled): set `PREP_SENSITIVITY <- TRUE` and `RUN_SENSITIVITY <- TRUE`.

   Each fit takes several hours (4 chains, 2,000 iterations).
3. Render the qmd (with a UTF-8 locale) to regenerate all figures in `Plots/`
   and the figure source data in `Source_data/`.

Manuscript figures correspond to: Fig 1 study overview; Fig 2 reporting
proportion; Fig 3 force of infection; Fig 4 age-specific susceptibility and
case-fatality; Fig 5 seasonal patterns. Figure S1 is the posterior correlation
between force of infection and reporting proportion; Figure S2 the choice of
the number of seasonal clusters; Figures S3–S5 the sensitivity analysis.

## Data sources

- Reported scrub typhus cases and deaths: National Disease Surveillance
  (Report 506), Division of Epidemiology, Department of Disease Control,
  Ministry of Public Health, Thailand (publicly available).
- Population by province and age: Bureau of Registration Administration,
  Ministry of Interior, Thailand (publicly available).
- Seroprevalence: Gonwong et al. (2022), *Am J Trop Med Hyg* 106(5):1363–1369.

Posterior predictive checks and MCMC trace plots are archived at
https://doi.org/10.5281/zenodo.23941245.
