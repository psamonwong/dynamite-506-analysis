# Archive

Code from earlier model development (Models 1–9 and exploratory scripts). None
of it is needed to reproduce the manuscript results; it is kept for reference.

The current workflow is `Scripts/main_analysis.qmd`, which uses
`Stan_models/Model10.stan` via `Scripts/run_models_local.R`.

- `Scripts/` – run scripts, trace plots and sandboxes for Models 1–9, and
  superseded plotting helpers (`plot_posterior_functions.R` loads `fit_09`).
- `Stan_models/` – Stan code and compiled caches for Models 1–9 and the
  additive false-positive model.
- `OPTIMIZATION_NOTES.md`, `lambda_trend_*.png` – notes and figures from the
  Model 7 development stage.
