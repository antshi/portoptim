% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/port_perf_turn.R
\name{calc_port_grosslev}
\alias{calc_port_grosslev}
\title{Portfolio Gross Leverage}
\usage{
calc_port_grosslev(port_weights)
}
\arguments{
\item{port_weights}{a numerical nxp data matrix with portfolio weights.}
}
\value{
a numeric vector with the gross leverage over the specified observation period.
}
\description{
Calculates the gross leverage rate of a portfolio strategy.
}
\examples{
set.seed(123)
# Simulated weights matrix (each row is a time period, columns are assets)
port_sim <- matrix(runif(100, min = 0, max = 1), nrow = 10, ncol = 10)
port_sim <- port_sim / rowSums(port_sim)
port_grosslev <- calc_port_grosslev(port_sim)

}
