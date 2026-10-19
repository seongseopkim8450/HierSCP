# HierSCP

R and C++ code for

> Kim, S. and Choi, T. Shape-Defined Robust Changepoint Detection for Hierarchical Time Series. Submitted to *Computational Statistics & Data Analysis*.

HierSCP is a Bayesian changepoint model for multiple time series. A regime is defined by the signs and intensities of the first two derivatives of the mean function, so a changepoint is a change of shape, with a jump in level determined by the two shape atoms that meet at the boundary. Series in a prespecified cluster share the sequence of shape regimes, and each series places its own changepoints. An observation-level variance-inflation mixture keeps isolated outliers from creating boundaries.

## Requirements

- R with a C++ toolchain for `Rcpp` (Rtools on Windows, Xcode command line tools on macOS).
- Model: `Rcpp`, `RcppArmadillo`, `mvtnorm`.
- Competing methods: `ppmSuite`, `bcp`.
- Convergence diagnostics: `posterior`.
- Data handling and figures: `dplyr`, `lubridate`, `ggplot2`.

```r
install.packages(c("Rcpp", "RcppArmadillo", "mvtnorm", "ppmSuite", "bcp",
                   "posterior", "dplyr", "lubridate", "ggplot2"))
```

Each script sets its working directory with `setwd()` at the top; `realdata_CLI_source.R` also sets `MODEL_DIR`, the folder that contains `OECD_CLI.csv`. Change these paths before running.

## Model code (`corefunctions_code/`)

The files are sourced in numerical order. `15_rcpp_bridge.R` compiles `hiercpd_core.cpp` and replaces the R functions of the same names with compiled versions; if compilation fails, the R versions are used.

| File | Contents | Paper |
|---|---|---|
| `00_utils.R` | Numerical helpers, samplers and the base measure $F_0$ of the shape atoms | Section 2.2 |
| `01_precompute.R` | Cosine basis and $D_{mm'}(x)$, so that $H_j(x) = \theta' D(x)\theta / \theta'\theta$ | Section 2.2 |
| `02_shape_function.R` | Shape function $f(x_t) = \{\delta_1 g + b(\delta_1+\delta_2)/2\}x_t - \delta_2 b H_j(x_t)$ and $\mu_{j,t} = \alpha + \beta_c + \gamma^*_j + f_{j,k}(x_t)$ | Sections 2.2, 2.5 |
| `03_state_management.R` | Regime paths, changepoints and free sets under (C1) and (C2) | Section 2.3 |
| `04_sampling.R` | Slice and elliptical slice samplers; truncated univariate and bivariate normal samplers | |
| `05_likelihood.R` | Log-likelihoods and the unnormalized log-posterior | Section 2.5 |
| `06_local_adjusting.R` | Local Adjusting: each boundary drawn from its position distribution over the free set | Section 3.2, S6 |
| `07_interval_adjusting.R` | Interval Adjusting: split and merge passes on the cluster-level regimes | Section 3.1, S3.4, S5 |
| `08_param_update_continuous.R` | Parameter Updating: level parameters, Gaussian process coefficients, shape atoms, series effects, robust error block, stick-breaking weights | Section 3.3, S4 |
| `10_initialization.R` | Initial values, hyperparameters and `create_model_spec()` | |
| `11_mcmc_main.R` | `run_mcmc()`: the schedule of Algorithm 1 and storage of the draws | Section 3 |
| `15_rcpp_bridge.R`, `hiercpd_core.cpp` | Compiled kernels | |
| `compute_atom_loglik_fallback.R` | R version of the atom log-likelihood, used when the compiled one is not loaded | |

## Analysis scripts

| File | Contents | Paper |
|---|---|---|
| `simulation_source.R` | Data generation for S1–S3, posterior summaries, HierSCP and competitor fits (ICP-PPM, CCP-PPM, BCP), evaluation metrics, replication drivers and boxplots | Sections 3.4, 4; S7.1–S7.7 |
| `sensitivity_S2.R` | One-at-a-time sensitivity analysis on S2 (seeds 2026–2045): $m_{\min} = 8, 12$ (default 10), $\lambda_0 = 0.5, 2$ (default 1), $m_0 = 0.01, 0.10$ (default 0.03); driver `run_sensitivity_S2.R` | S7.8 |
| `realdata_CLI_source.R` | OECD CLI application: data, model fit, figures, changepoint table and convergence diagnostics | Section 5, S8 |

### Data References
Scripts and datasets for reproducing the empirical application on macroeconomic indicators.

* **`OECD_CLI.csv`**: The OECD Composite Leading Indicators (CLI) dataset for nine economies (North America, Europe, East Asia). 
  * *Data Source:* Retrieved from the official [OECD Data Portal - Composite Leading Indicator (CLI)](https://data.oecd.org/leadind/composite-leading-indicator-cli.htm).

### Simulation Study (`simulation_code_run/`)
Scripts to reproduce the simulation scenarios described in the paper.

* **`simulation_source.R`**: Generates synthetic data and runs the MCMC sampler across three primary scenarios: benign shape detection (S1), detection under heavy-tailed contamination (S2), and detection under model misspecification (S3).

