% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/port_optim_trad.R
\name{port_optim_gmv}
\alias{port_optim_gmv}
\title{Global Minimum-Variance (GMV) Portfolio Optimization}
\usage{
port_optim_gmv(sigma_mat)
}
\arguments{
\item{sigma_mat}{a pxp covariance matrix of asset returns.}
}
\value{
a numeric vector with the GMV portfolio weights.
}
\description{
Calculates the weights of a GMV portfolio strategy.
}
\examples{
data(prices_m)
rets_m <- calc_rets(prices_m)
sigma_mat <- cov(rets_m)
port_gmv <- port_optim_gmv(sigma_mat)

}
