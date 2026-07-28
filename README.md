# dgtf

> Variational and MCMC inference for Bayesian state-space models of count data.

`dgtf` provides a unified R interface to variational Bayes and MCMC for a general class of dynamic generalised transfer function models for count data.

## Installation

1. Download the git repo via terminal:
```
git clone git@github.com:minimeini/dgtf.git
```
2. Change directory to the repository: `cd dgtf`
3. Install the R package
    - In terminal: `R CMD build .`
    - Or in R: `devtools::install(".")`


### System requirements

- C++17 compiler
- LAPACK / BLAS (provided by R)

OpenMP is intentionally disabled in the current release while the parallel SMC
sections are audited for thread safety.

### R-package dependencies

- `LinkingTo`: `Rcpp`, `RcppArmadillo`, `BH`, `RcppProgress`, `pg`.
- `Imports`: `Rcpp`, `ggplot2` (plus base `stats`, `utils`, `graphics`, and
  `grDevices`).

## License

MIT (see [`LICENSE`](LICENSE)).
