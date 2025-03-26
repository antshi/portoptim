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

There are multiple portfolio optimization functions included. Some deal with the traditional, theoretical methodologies, such as

- `port_optim_naive`: Naive portfolio diversification.
- `port_optim_gmv`: Global Minimum-Variance portfolio.
- `port_optim_tang`: Tangency portfolio.
- `port_optim_markowitz`: Markowitz portfolio.
- `port_optim_eff`: Efficient frontier portfolios.
- `port_optim_cml`: Capital Market Line portfolios.

Moreover, there is a numerical portfolio optimization function that can calculate the GMV or the tangency portfolios with constraints: `port_optim_solver`.

### Portfolio Performance

We include functions for calculating and estimating various portfolio performance parameters, such as

- `calc_port_rets`: Portfolio returns.
- `calc_port_mu_is` and `calc_port_var_is`: Expected return and variance of portfolio returns in-sample.
- `calc_mean`, `calc_var`, and `calc_sr`: Expected return, variance, and Sharpe ratio of portfolio returns out-of-sample.
- `calc_port_turnover`, `calc_port_grosslev`, and `calc_port_proplev`: Portfolio turnover, gross leverage, and proportional leverage.

### Statistical Inference

- `boot_nonparam`:  Nonparametric Bootstrap.
- `boot_param`: Parametric Botstrap.
- `f_test`: F-test for variance comparison of two time series.
- `hac_infer`: HAC-based inference.
- `boot_infer`: Bootstrap-based inference.

## Example Usage

```r
library(portoptim)
library(tidyverse)
library(ggplot2)

data(prices_m)

rets_m <- calc_rets(prices_m)
sigma_mat <- cov(rets_m)
mu_vec <- mu_estim_wrapper(rets_m, mu_estim_ml)

port_naive <- port_optim_naive(sigma_mat)
port_gmv <- port_optim_wrapper(sigma_mat, port_optim_gmv)
port_tang <- port_optim_wrapper(sigma_mat, port_optim_tang, mu_vec)
port_gmv_noshort <- port_optim_wrapper(sigma_mat, port_optim_solver, porttype = "GMV", A = diag(1, ncol(sigma_mat)), b = rep(0, ncol(sigma_mat)))

port_weights <- data.frame(
  Asset = colnames(sigma_mat),
  Naive = port_naive,
  GMV = port_gmv,
  Tangency = port_tang,
  GMV_NoShort = port_gmv_noshort
) %>%
  pivot_longer(cols = -Asset, names_to = "Portfolio", values_to = "Weight")

ggplot(port_weights, aes(x = Asset, y = Weight, fill = Portfolio)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Portfolio Weights", x = "Assets", y = "Weights") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
```
