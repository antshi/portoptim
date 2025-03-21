# portoptim

<!-- badges: start -->
<!-- badges: end -->

**portoptim** is an R package designed for portfolio optimization and evaluation.
It includes functions for calculating expected returns, estimation of expected returns, optimizing portfolios, and calculating performance and performing statistical inference.

## Installation

To install the package from GitHub, use the following commands:

```r
install.packages("devtools")
devtools::install_github("antshi/portoptim")
```

## Datasets

The package includes the following datasets:

- `prices_d`: Daily stock returns.
- `prices_spx_d`: Daily index returns.
- `prices_m`: Monthly stock returns.
- `prices_spx_m`: Monthly index returns.
- `ff_5factors`: 5 Fama-French Research Factors.

## Features

### Data Preparation

- `find_reps`: Find repeated values in a data frame.
- `calc_rets`: Calculate stock returns.
- `mu_estim_wrapper`: Estimate the expected returns (various ways).

### Portfolio Optimization

- `port_optim_naive`: Naive portfolio diversification.
- `port_optim_gmv`: Global Minimum-Variance portfolio.
- `port_optim_tang`: Tangency portfolio.
- `port_optim_markowitz`: Markowitz portfolio.
- `port_optim_eff`: Efficient frontier portfolios.
- `port_optim_cml`: Capital Market Line portfolios.
- and many more.

### Portfolio Performance

### Statistical Inference

- `boot_nonparam`:  Nonparametric Bootstrap.
- `boot_param`: Parametric Botstrap.
- `f_test_rets`: F-test for variance comparison.
- `hac_infer`: HAC-based inference.
- `boot_infer`: Bootstrap-based inference.

## Example Usage

```r
library(portoptim)
data(prices_m)
rets_m <- calc_rets(prices_m)
sigma_mat <- cov(rets_m)
port_naive <- port_optim_naive(sigma_mat)
port_gmv <- port_optim_gmv(sigma_mat)

matplot(cbind(port_naive, port_gmv), type = "l", col = c("black", "green"), lty = 1, lwd = 1, ylab = "Weights")
legend("topleft", legend = c("Naive", "GMV"), col = c("black", "green"), lty = 1, lwd = 1)
```
