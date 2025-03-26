#' Naive Portfolio Diversification I
#'
#' Calculates the weights of a naive portfolio strategy.
#'
#' @param sigma_mat a pxp covariance matrix of asset returns.
#'
#' @return a numeric vector with the naive portfolio weights.
#'
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' sigma_mat <- cov(rets_m)
#' port_naive <- port_optim_naive(sigma_mat)
#'
#' @export port_optim_naive
#'
port_optim_naive <- function(sigma_mat) {
  p <- dim(sigma_mat)[1]
  weights <- as.numeric(rep.int(1 / p, p))
  weights
}

#' Naive Portfolio Diversification II
#'
#' Calculates the weights of a naive portfolio strategy.
#'
#' @param p an integer, specifying the number of stocks.
#'
#' @return a numeric vector with the naive portfolio weights.
#'
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' port_naive <- port_optim_naive_slim(dim(rets_m)[2])
#'
#' @export port_optim_naive_slim
#'
port_optim_naive_slim <- function(p) {
  weights <- as.numeric(rep.int(1 / p, p))
  weights
}

#' Global Minimum-Variance (GMV) Portfolio Optimization
#'
#' Calculates the weights of a GMV portfolio strategy.
#'
#' @param sigma_mat a pxp covariance matrix of asset returns.
#'
#' @return a numeric vector with the GMV portfolio weights.
#'
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' sigma_mat <- cov(rets_m)
#' port_gmv <- port_optim_gmv(sigma_mat)
#'
#' @export port_optim_gmv
#'
port_optim_gmv <- function(sigma_mat) {
  if (any(diag(sigma_mat) <= 0)) stop("Covariance matrix must be positive definite.")

  p <- dim(sigma_mat)[1]
  ones <- rep.int(1, p)
  sigma_mat_inv <- solve(sigma_mat)

  weights <- as.numeric(sigma_mat_inv %*% ones / as.numeric(ones %*% sigma_mat_inv %*% ones))
  weights
}

#' Tangency Portfolio Optimization
#'
#' Calculates the weights of a tangency portfolio strategy.
#'
#' @param sigma_mat a pxp matrix, the covariance matrix of asset returns.
#' @param mu_vec a numeric vector of length p, the expected returns.
#' @param rf a double, the assumed risk-free return. Default value is 0.
#'
#' @return a numeric vector with the Tangency portfolio weights.
#'
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' sigma_mat <- cov(rets_m)
#' mu_vec <- mu_estim_wrapper(rets_m, mu_estim_ml)
#' port_tang <- port_optim_tang(sigma_mat, mu_vec)
#'
#' @export port_optim_tang
#'
port_optim_tang <- function(sigma_mat, mu_vec, rf = 0) {
  if (any(diag(sigma_mat) <= 0)) stop("Covariance matrix must be positive definite.")

  p <- dim(sigma_mat)[1]
  if (length(mu_vec) != p) stop("Expected returns vector must match covariance matrix dimensions.")

  mu_vec_excess <- mu_vec - rf
  ones <- rep.int(1, p)
  sigma_mat_inv <- solve(sigma_mat)

  weights <- as.numeric(sigma_mat_inv %*% mu_vec_excess / (as.numeric(ones %*% sigma_mat_inv %*% mu_vec_excess)))
  weights
}

#' Markowitz Portfolio Optimization
#'
#' Calculates the weights of an optimal Markowitz portfolio strategy.
#'
#' @param sigma_mat a pxp covariance matrix of asset returns.
#' @param mu_vec a numeric vector of length p, the expected returns.
#' @param rf a double, the assumed risk-free return. Default value is 0.
#' @param gamma an integer, the risk aversion parameter. Default value is 2.
#'
#' @return a numeric vector with the Markowitz portfolio weights.
#'
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' sigma_mat <- cov(rets_m)
#' mu_vec <- mu_estim_wrapper(rets_m, mu_estim_ml)
#' port_markowitz <- port_optim_markowitz(sigma_mat, mu_vec)
#'
#' @export port_optim_markowitz
#'
port_optim_markowitz <- function(sigma_mat, mu_vec, rf = 0, gamma = 2) {
  sigma_mat <- as.matrix(sigma_mat)
  if (any(diag(sigma_mat) <= 0)) stop("Covariance matrix must be positive definite.")

  sigma_mat_inv <- solve(sigma_mat)
  mu_vec_excess <- mu_vec - rf
  if (length(mu_vec) != dim(sigma_mat)[1]) stop("Expected returns vector must match covariance matrix dimensions.")

  weights <- as.numeric((1 / gamma) * sigma_mat_inv %*% mu_vec_excess)
  weights
}

#' Efficient Portfolio Optimization
#'
#' Calculates the portfolios along the Efficient Frontier.
#'
#' @param sigma_mat a pxp covariance matrix of asset returns.
#' @param mu_vec a numeric vector of length p, the expected returns.
#' @param mu_grid a vector of length m as the grid of expected returns,
#' along which the capital market line is to be estimated.
#' @param res_all a logical. If TRUE,
#' the result includes the calculated weights, the expected returns,
#' and the standard deviations for the efficient frontier.
#' If FALSE, only the weights. Default value is FALSE.
#'
#' @return a pxm matrix with the weights of the efficient portfolios along mu_grid.
#' @return a numeric vector of length m with the corresponding expected returns, actually mu_grid.
#' @return a numeric vector of length m with the corresponding standard deviations.
#'
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' sigma_mat <- cov(rets_m)
#' mu_vec <- mu_estim_wrapper(rets_m, mu_estim_ml)
#' mu_grid <- seq(0.8 * min(mu_vec), 1.2 * max(mu_vec), length.out = 50)
#' port_eff <- port_optim_eff(sigma_mat, mu_vec, mu_grid)
#' port_eff_all <- port_optim_eff(sigma_mat, mu_vec, mu_grid, res_all = TRUE)
#'
#' @export port_optim_eff
port_optim_eff <- function(sigma_mat, mu_vec, mu_grid, res_all = FALSE) {
  if (any(diag(sigma_mat) <= 0)) stop("Covariance matrix must be positive definite.")

  p <- dim(sigma_mat)[1]
  if (length(mu_vec) != p) stop("Expected returns vector must match covariance matrix dimensions.")

  ones <- rep.int(1, p)
  sigma_mat_inv <- solve(sigma_mat)
  A <- as.numeric(ones %*% sigma_mat_inv %*% mu_vec)
  B <- as.numeric(mu_vec %*% sigma_mat_inv %*% mu_vec)
  C <- as.numeric(ones %*% sigma_mat_inv %*% ones)
  g <- as.numeric(1 / (C * B - A^2) * (B * sigma_mat_inv %*% ones - A * sigma_mat_inv %*% mu_vec))
  h <- 1 / (C * B - A^2) * (C * sigma_mat_inv %*% mu_vec - A * sigma_mat_inv %*% ones)

  weights <- h %*% mu_grid + g
  mus <- mu_grid
  sds <- as.numeric(sqrt(C / (C * B - A^2) * (mu_grid - A / C)^2 + 1 / C))

  if (res_all) {
    return(list(weights, mus, sds))
  } else {
    return(weights)
  }
}
