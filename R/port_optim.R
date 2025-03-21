#' Wrapper Function for Portfolio Optimization I
#'
#' Allows the execution of all included optimization functions.
#'
#' @param sigma_mat a pxp covariance matrix of asset returns.
#' @param estim_func a function for portfolio optimization.
#' @param ... additional arguments to be passed to estim_func
#'
#' @return the value of the executed function estim_func.
#'
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' sigma_mat <- cov(rets_m)
#' mu_vec <- mu_estim_wrapper(rets_m, mu_estim_ml)
#' port_gmv <- port_optim_wrapper(sigma_mat, port_optim_gmv)
#' port_tang <- port_optim_wrapper(sigma_mat, port_optim_tang, mu_vec)
#'
#' @export port_optim_wrapper
#'
port_optim_wrapper <- function(sigma_mat, estim_func, ...) {
  estim_func(sigma_mat, ...)
}
