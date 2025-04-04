% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/port_perf_rets.R
\name{calc_port_rets}
\alias{calc_port_rets}
\title{Portfolio Returns}
\usage{
calc_port_rets(port_weights, rets)
}
\arguments{
\item{port_weights}{a numeric vector with the portfolio weights.}

\item{rets}{an nxp matrix of stock returns.}
}
\value{
a numeric vector of the portfolio returns.
}
\description{
Calculates the portfolio returns.
}
\examples{
data(prices_m)
rets_m <- calc_rets(prices_m)
sigma_mat <- cov(rets_m)
port_naive <- port_optim_naive(sigma_mat)
port_naive_rets <- calc_port_rets(port_naive, rets_m)
}
