# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

R package for maximum likelihood estimation of Weibull component parameters in series systems from masked failure data. This package implements the methodology from the master's paper:

> **Towell, A.** (2023). *Reliability Estimation in Series Systems: Maximum Likelihood Techniques for Right-Censored and Masked Failure Data*. Master's Project, Southern Illinois University Edwardsville.

- **Paper Repository**: https://github.com/queelius/reliability-estimation-in-series-systems
- **This Package**: https://github.com/queelius/wei.series.md.c1.c2.c3

The package handles:
- Right-censored system lifetimes
- Masked component failure causes (candidate sets satisfying conditions C1, C2, C3)

## Build & Development Commands

```bash
# Install package locally (from R)
devtools::install()

# Build and check package
R CMD build .
R CMD check wei.series.md.c1.c2.c3_*.tar.gz

# Generate documentation (roxygen2)
devtools::document()

# Run tests
devtools::test()

# Check test coverage
covr::package_coverage()

# Build pkgdown site
pkgdown::build_site()
```

## Dependencies

Install required GitHub packages before building:
```r
devtools::install_github("queelius/md.tools")
devtools::install_github("queelius/algebraic.mle")
devtools::install_github("queelius/algebraic.dist")
```

## Architecture

### Core Likelihood Functions (R/md_wei_series_c1_c2_c3.R)
- `loglik_wei_series_md_c1_c2_c3()` - Log-likelihood for masked Weibull series data
- `score_wei_series_md_c1_c2_c3()` - Analytical score function (gradient)
- `hessian_wei_series_md_c1_c2_c3()` - Numerical Hessian via `numDeriv::jacobian`

### MLE Solvers (R/md_mle_wei_series_c1_c2_c3.R)
Multiple optimization approaches wrapping `optim`:
- `mle_lbfgsb_*` - L-BFGS-B with analytical gradient (recommended)
- `mle_newton_*` - Newton-Raphson with line search
- `mle_nelder_*` - Nelder-Mead (gradient-free)
- `mle_sann_*` - Simulated annealing
- `mle_grad_*` - Gradient ascent

### Distribution Functions (R/wei_series.R)
Standard R distribution API for Weibull series:
- `dwei_series`, `pwei_series`, `qwei_series`, `rwei_series`
- `surv_wei_series`, `hazard_wei_series`
- `wei_series_mttf`, `wei_series_cause`

### Masked Data Generation (R/md_series_bernoulli_masked_component_cause.R)
- `md_bernoulli_cand_c1_c2_c3()` - Assigns masking probabilities to components
- `md_cand_sampler()` - Generates candidate sets from probabilities

## Data Format

Parameter vector `theta` is ordered as: `(shape1, scale1, shape2, scale2, ..., shapem, scalem)`

Data frames require:
- `t` - System lifetime column
- `x1, x2, ..., xm` - Boolean candidate set indicators (which components could have failed)
- `delta` - Right-censoring indicator (TRUE = failure observed, FALSE = censored)

## Key Theoretical Conditions

The likelihood model assumes candidate sets satisfy (see paper Section 3):
- **C1**: True failure cause is always in candidate set
- **C2**: Masking probability is independent of which candidate caused failure
- **C3**: Masking probabilities are independent of the parameter vector theta

## Related Paper Sections

Key correspondences between this package and the paper:
- **Likelihood derivation**: Paper Section 3 → `loglik_wei_series_md_c1_c2_c3()`
- **Score function**: Paper Section 3.2 → `score_wei_series_md_c1_c2_c3()`
- **Simulation studies**: Paper Section 5 uses this package for MLE computation
- **Bernoulli masking model**: Paper Section 4.1 → `md_bernoulli_cand_c1_c2_c3()`
- **Guo et al. (2013) replication**: Paper Section 5.1 → `guo_weibull_series_md` dataset
