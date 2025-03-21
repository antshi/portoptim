% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/port_optim_trad.R
\name{port_optim_naive_slim}
\alias{port_optim_naive_slim}
\title{Naive Portfolio Diversification II}
\usage{
port_optim_naive_slim(p)
}
\arguments{
\item{p}{an integer, specifying the number of stocks.}
}
\value{
a numeric vector with the naive portfolio weights.
}
\description{
Calculates the weights of a naive portfolio strategy.
}
\examples{
data(prices_m)
rets_m <- calc_rets(prices_m)
port_naive <- port_optim_naive_slim(dim(rets_m)[2])

}
