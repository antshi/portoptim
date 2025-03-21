#' Portfolio Returns
#'
#' Calculates the portfolio returns.
#'
#' @param port_weights a numeric vector with the portfolio weights.
#' @param rets an nxp matrix of stock returns.
#'
#' @return a numeric vector of the portfolio returns.
#'
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' sigma_mat <- cov(rets_m)
#' port_naive <- port_optim_naive(sigma_mat)
#' port_naive_rets <- calc_port_rets(port_naive, rets_m)
#' @export calc_port_rets
#'
calc_port_rets <- function(port_weights, rets) {
  rets <- as.matrix(rets)
  port_rets <- rets %*% port_weights
  as.numeric(port_rets)
}

#' Portfolio Expected Returns (in-sample)
#'
#' Calculates the in-sample portfolio expected returns.
#'
#' @param port_weights a numeric vector, the portfolio weights.
#' @param mu_vec a numeric vector of expected returns.
#'
#' @return a double, the portfolio expected return.
#'
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' sigma_mat <- cov(rets_m)
#' mu_vec <- mu_estim_ml(rets_m)
#' port_naive <- port_optim_naive(sigma_mat)
#' port_mu_is <- calc_port_mu_is(port_naive, mu_vec)
#'
#' @export calc_port_mu_is
calc_port_mu_is <- function(port_weights, mu_vec) {
  port_mu_is <- t(port_weights) %*% mu_vec
  as.numeric(port_mu_is)
}

#' Portfolio Variance (in-sample)
#'
#' Calculates the in-sample portfolio variance.
#' @param port_weights a numeric vector, the portfolio weights.
#' @param sigma_mat a pxp covariance matrix of returns.
#'
#' @return a double, the portfolio variance.
#'
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' sigma_mat <- cov(rets_m)
#' port_naive <- port_optim_naive(sigma_mat)
#' port_var_is <- calc_port_var_is(port_naive, sigma_mat)
#'
#' @export calc_port_var_is
#'
calc_port_var_is <- function(port_weights, sigma_mat) {
  port_var_is <- t(port_weights) %*% sigma_mat %*% port_weights
  as.numeric(port_var_is)
}

#' Expected Return of Returns
#'
#' Calculates the expected return of stock returns.
#'
#' @param rets a numerical vector of stock or portfolio returns.
#' @param ann_factor a double, the annualization factor. If ann_factor=1 (default), no annualization is performed.
#' For monthly returns, set ann_factor=12. For daily returns, set ann_factor=252, etc.
#'
#' @return a double the expected return.
#'
#' @examples
#' data(prices_spx_m)
#' rets_spx_m <- calc_rets(prices_spx_m)
#' spx_mean <- calc_mean(rets_spx_m)
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' stocks_mean <- apply(rets_m, 2, calc_mean)
#'
#' @export calc_mean
#'
calc_mean <- function(rets,
                      ann_factor = 1) {
  ann_factor * mean(rets, na.rm = TRUE)
}

#' Variance of Returns
#'
#' Calculates the variance of stock returns.
#'
#' @param rets a numerical vector of stock or portfolio returns.
#' @param ann_factor a double, the annualization factor. If ann_factor=1 (default), no annualization is performed.
#' For monthly returns, set ann_factor=12. For daily returns, set ann_factor=252, etc.
#'
#' @return a double the variance.
#'
#' @examples
#' data(prices_spx_m)
#' rets_spx_m <- calc_rets(prices_spx_m)
#' spx_var <- calc_var(rets_spx_m)
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' stocks_var <- apply(rets_m, 2, calc_var)
#'
#' @export calc_var
#'
calc_var <- function(rets,
                     ann_factor = 1) {
  ann_factor * stats::var(rets, na.rm = TRUE)
}

#' Sharpe Ratio of Returns
#'
#' Calculates the Sharpe ratio of returns.
#'
#' @param rets a numerical vector of stock or portfolio returns.
#' @param rf a double, the assumed risk-free return. Default=0.
#' @param ann_factor a double, the annualization factor. If ann_factor=1 (default), no annualization is performed.
#' For monthly returns, set ann_factor=12. For daily returns, set ann_factor=252, etc.
#'
#' @return a double the Sharpe ratio.
#'
#' @examples
#' data(prices_spx_m)
#' rets_spx_m <- calc_rets(prices_spx_m)
#' spx_sharpe <- calc_sr(rets_spx_m)
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' stocks_sharpe <- apply(rets_m, 2, calc_sr)
#'
#' @export calc_sr
#'
calc_sr <- function(rets,
                    rf = 0,
                    ann_factor = 1) {
  (ann_factor * (mean(rets, na.rm = TRUE) - rf)) / (sqrt(ann_factor) * stats::sd(rets, na.rm = TRUE))
}
