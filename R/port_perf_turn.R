#' Portfolio Turnover
#'
#' Calculates the average turnover of a portfolio strategy.
#'
#' @param port_weights a numerical nxp data matrix with portfolio weights over the whole observation period.
#'
#' @return a numeric vector with the turnover over the specified observation period.
#'
#' @examples
#' set.seed(123)
#' # Simulated weights matrix (each row is a time period, columns are assets)
#' port_sim <- matrix(runif(100, min = 0, max = 1), nrow = 10, ncol = 10)
#' port_sim <- port_sim / rowSums(port_sim)
#' port_turn <- calc_port_turnover(port_sim)
#'
#' @export calc_port_turnover
#'
calc_port_turnover <- function(port_weights) {
  rowSums(abs(diff(port_weights)))
}

#' Portfolio Gross Leverage
#'
#' Calculates the gross leverage rate of a portfolio strategy.
#'
#' @param port_weights a numerical nxp data matrix with portfolio weights.
#'
#' @return a numeric vector with the gross leverage over the specified observation period.
#'
#' @examples
#' set.seed(123)
#' # Simulated weights matrix (each row is a time period, columns are assets)
#' port_sim <- matrix(runif(100, min = 0, max = 1), nrow = 10, ncol = 10)
#' port_sim <- port_sim / rowSums(port_sim)
#' port_grosslev <- calc_port_grosslev(port_sim)
#'
#' @export calc_port_grosslev
#'
calc_port_grosslev <- function(port_weights) {
  apply(abs(port_weights), 1, sum)
}

#' Portfolio Proportional Leverage
#'
#' Calculates the proportional leverage (% short sales) of a portfolio strategy.
#'
#' @param port_weights a numerical nxp data matrix with portfolio weights.
#'
#' @return a numeric vector with the proportional leverage over the specified observation period.
#'
#' @examples
#' set.seed(123)
#' # Simulated weights matrix (each row is a time period, columns are assets)
#' port_sim <- matrix(runif(100, min = 0, max = 1), nrow = 10, ncol = 10)
#' port_sim <- port_sim / rowSums(port_sim)
#' port_proplev <- calc_port_proplev(port_sim)
#'
#' @export calc_port_proplev
#'
calc_port_proplev <- function(port_weights) {
  apply(port_weights < 0, 1, sum)
}
