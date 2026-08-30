# HierSCP
R and C++ implementation of HierSCP (Hierarchical Shape-Defined Robust Changepoint Detection for Multivariate Time Series).

# HierSCP: Hierarchical Shape-Defined Robust Changepoint Detection (2026)

**Author:** Seongseop Kim 

## Overview
This repository contains the R and C++ implementation of the **Hierarchical Shape-defined Change Point (HierSCP)** detection model. Unlike traditional moment-based detectors, HierSCP characterizes a regime entirely by the geometric shape of its trajectory (the signs and magnitudes of the first two derivatives) and identifies changepoints as transitions between these shapes. 

The model employs a two-level hierarchical state process for multiple time series, a continuous shape function based on a Karhunen-Loève expansion, and a robust normal-inverse-gamma scale mixture to prevent heavy-tailed outliers from distorting shape inferences.

## Repository Structure

The codebase is highly modularized into three main directories: core MCMC functions, simulation scripts, and the real-world data application.

### 1. Core Functions (`corefunctions_code/`)
This directory contains the primary MCMC sampler engine, separated by algorithm steps.

* **`00_utils.R`**: General utility and helper functions for data manipulation and basic mathematical operations.
* **`01_precompute.R`**: Precomputes the cosine basis double integrals and normalizers to avoid redundant numerical integration inside the MCMC loop.
* **`02_shape_function.R`**: Evaluates the core shape function curves using the spectral construction and the four-parameter shape-determining atoms.
* **`03_state_management.R`**: Manages the two-level hierarchical latent states, strictly enforcing the structural constraints (adjacent forward transition and minimum segment length) and defining the free observation sets.
* **`04_sampling.R`**: Core statistical sampling routines, including categorical draws and stick-breaking probability formulations.
* **`05_likelihood.R`**: Computes the segment likelihoods, integrating the robust scale-mixture error structure to discount outliers.
* **`06_local_adjusting.R`**: Implements **Direct Independence Sampling (DIS)**. It relocates existing changepoint boundaries within their free observation sets exactly in a single pass without relying on sequential restricted scans.
* **`07_interval_adjusting.R`**: Implements the dimension-changing **Split and Merge** moves. It dynamically infers the number of shape regimes using an **sequential categorical draw with no accept-reject step**, completely bypassing standard Reversible Jump MCMC inefficiencies.
* **`08_param_update_continuous.R`**: Updates all continuous parameters conditional on the current segmentation, utilizing Gibbs sampling, elliptical slice sampling (for GP coefficients), and standard slice sampling.
* **`10_initialization.R`**: Sets the initial parameter values and hyperparameter starting points to launch the MCMC chains.
* **`11_mcmc_main.R`**: The main wrapper script that orchestrates the entire MCMC loop (Interval Adjusting $\rightarrow$ Local Adjusting $\rightarrow$ Parameter Updating).
* **`15_rcpp_bridge.R` & `hiercpd_core.cpp`**: C++ integration via Rcpp to accelerate severe computational bottlenecks (e.g., intensive likelihood evaluations and sequential MCMC scans).
* **`compute_atom_loglik_fallback.R`**: A pure R fallback routine for atom log-likelihood computations in case the C++ bridge encounters issues.

### 2. Real Data Application (`realdata_data_code/`)
Scripts and datasets for reproducing the empirical application on macroeconomic indicators.

* **`OECD_CLI.csv`**: The OECD Composite Leading Indicators (CLI) dataset for nine economies (North America, Europe, East Asia). 
  * *Data Source:* Retrieved from the official [OECD Data Portal - Composite Leading Indicator (CLI)](https://data.oecd.org/leadind/composite-leading-indicator-cli.htm).
* **`realdata_CLI_source.R`**: Execution script that applies the HierSCP model to the OECD dataset to extract cluster-level phase transitions and analyze the 2020 pandemic's impact.

### 3. Simulation Study (`simulation_code_run/`)
Scripts to reproduce the simulation scenarios described in the paper.

* **`simulation_source.R`**: Generates synthetic data and runs the MCMC sampler across three primary scenarios: benign shape detection (S1), detection under heavy-tailed contamination (S2), and detection under model misspecification (S3).

## Usage
To run the model, source the core functions from the `corefunctions_code` directory and execute either the simulation or real-data scripts. Ensure that the `Rcpp` package is installed and properly configured to compile `hiercpd_core.cpp`.
