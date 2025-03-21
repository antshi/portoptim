% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/port_optim_trad.R
\name{port_optim_cml}
\alias{port_optim_cml}
\title{CML Portfolio Optimization}
\usage{
port_optim_cml(sigma_mat, mu_vec, mu_grid, rf = 0, res_all = FALSE)
}
\arguments{
\item{sigma_mat}{a pxp covariance matrix of asset returns.}

\item{mu_vec}{a numeric vector of length p, the expected returns.}

\item{mu_grid}{a vector with length m as the grid of expected returns,
along which the capital market line is to be estimated.}

\item{rf}{a double, the assumed risk-free return. Default value is 0.}

\item{res_all}{a logical.
If TRUE, the result includes the calculated weights and the standard deviations for the CML.
If FALSE, only the weights. Default value is FALSE.}
}
\value{
a pxm matrix with the weights of the CML portfolio along mu_grid.

a vector of length m with the corresponding standard deviations.
}
\description{
Calculates the weights and the standard deviations along the Capital Market Line (CML).
}
\examples{
data(prices_m)
rets_m <- calc_rets(prices_m)
sigma_mat <- cov(rets_m)
mu_vec <- mu_estim_wrapper(rets_m, mu_estim_ml)
mu_grid <- seq(0.8 * min(mu_vec), 1.2 * max(mu_vec), length.out = 50)
port_cml <- port_optim_cml(sigma_mat, mu_vec, mu_grid)
port_cml_all <- port_optim_cml(sigma_mat, mu_vec, mu_grid, res_all = TRUE)

}
