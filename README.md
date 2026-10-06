# QEPF: Quantile-Based Effectiveness Persistence Function

R code for

> Sankaran, P.G., Prasanth, V.P. and Midhu, N.N. *Quantile-Based Effectiveness Persistence Function: A Tail-Focused Metric with Theory, Estimation, and Application to Biosimilar Evaluation.* (Under review.)

The QEPF of a positive outcome $X$ with quantile function $Q(u)$ is

$$P(u)=\frac{V(u)}{Q(u)}=\frac{1}{1-u}\int_u^1\frac{Q(p)}{Q(u)}\,dp,\qquad 0<u<1,$$

the average multiple of its own entry threshold $Q(u)$ reached by the top $(1-u)$ fraction of responders. The repository contains the nonparametric estimator, a permutation test of equal persistence, an equivalence test for tail shape, the tail-mean ratio, and the code for every simulation and the real-data analysis in the paper.

## Using the methods on your own data

```r
source("R/qepf_functions.R")

qepf_hat(x, u = c(0.6, 0.7, 0.8, 0.9))            # estimate P(u)
qepf_perm_test(x, y, U = c(0.60, 0.90), B = 2000)   # test of equal persistence
qepf_equiv_test(x, y, Delta = 0.15)                 # equivalence of tail shape
qepf_tail_mean_ratio(x, y, u1 = 0.60)               # level component V_y(u1)/V_x(u1)
```

`qepf_functions.R` also contains closed-form $P(u)$ for the uniform, exponential, Pareto I, power, Weibull, gamma and LMRQD distributions.

## Reproducing the paper

Open `QEPF.Rproj` in RStudio, or set the working directory to this folder, then run the scripts in `R/`. Each script is self-contained and writes to `results/` or `figures/`.

| Script | Reproduces | Output | Approximate run time |
|---|---|---|---|
| `R/01_sim_estimator.R` | Table 2 (bias and MSE of the estimator) | `results/sim_estimator_bias_mse.csv` | minutes |
| `R/02_sim_equality_test.R` | Tables 3–4 (size and power of the test of equal persistence) | `results/sim_equality_size.csv`, `results/sim_equality_power.csv` | several hours (parallel) |
| `R/03_sim_equivalence_test.R` | Table 5 (equivalence test) | `results/sim_equivalence.csv` | about 15 minutes (parallel) |
| `R/04_fig3_estimator_curves.R` + `figures/fig3.tex` | Figure 3 | `figures/qepf_*.csv`, `figures/fig3.pdf` | about a minute |
| `R/05_smarter_analysis.R` | Figure 4, Section 8 test, Table 6 | `figures/fig4_smarter_curves.pdf`, `results/smarter_*.csv` | several hours |
| `R/06_smarter_level_component.R` | Section 8: tail-mean ratio and $\Delta_{\min}$ | `results/smarter_level_component.csv` | minutes |
| `R/07_run_time.R` | Section 6.5 (timing) | printed | seconds |
| `figures/fig1.tex` | Figure 1 (closed-form $P(u)$) | `figures/fig1.pdf` | seconds |

The simulation scripts have a `QUICK` switch at the end: `QUICK <- TRUE` runs a one-minute check, `QUICK <- FALSE` (the default) runs the full design used in the paper. The `results/` folder already contains the outputs used in the paper. Figures 3 and 1 are built with pgfplots: compile `fig3.tex` and `fig1.tex` from inside `figures/` with pdflatex or lualatex.

## Requirements

R (≥ 4.1) with the base package `parallel`. `R/05_smarter_analysis.R` also needs `readr`, `dplyr`, `purrr`, `ggplot2` and `latex2exp`; the Figure 3 preview uses `ggplot2` if installed.

## Data

`data/smarter_anonymised_data.csv` is the anonymised individual-level dataset of the SMARTER cluster-randomised trial, obtained from the Dryad repository (Li, 2025; doi:[10.5061/dryad.tmpg4f58w](https://doi.org/10.5061/dryad.tmpg4f58w)). Please cite the original trial (Zhang et al., 2025, *BMJ* 389, e082765) and the Dryad dataset when using it. See `data/README.md`.

## Licence

Code: MIT licence (see `LICENSE`). The SMARTER data are subject to the terms of the Dryad repository.

## Contact

Prasanth V.P. (prasanth.stat@gmail.com), Department of Statistics, Cochin University of Science and Technology, Kochi, India.
