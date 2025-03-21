#' Markowitz Portfolio Optimization (Unbiased)
#'
#' Calculates the weights of an optimal unbiased Markowitz portfolio strategy.
#'
#' @param sigma_mat a pxp covariance matrix of asset returns.
#' @param mu_vec a numeric vector of length p, the expected returns.
#' @param rf a double, the assumed risk-free return. Default value is 0.
#' @param gamma an integer, the risk aversion parameter. Default value is 2.
#' @param is_size an integer, the number of observations of the returns data.
#'
#' @return a numeric vector with the Markowitz unbiased portfolio weights.
#'
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' n <- nrow(rets_m)
#' sigma_mat <- cov(rets_m)
#' mu_vec <- mu_estim_wrapper(rets_m, mu_estim_ml)
#' port_markowitz_unb <- port_optim_markowitz_unb(sigma_mat, mu_vec, is_size = n)
#'
#' @export port_optim_markowitz_unb
#'
port_optim_markowitz_unb <-
  function(sigma_mat,
           mu_vec,
           rf = 0,
           gamma = 2,
           is_size) {
    n <- is_size
    p <- dim(sigma_mat)[1]
    unb <- (n - p - 2) / n
    port_markowitz <- port_optim_markowitz(sigma_mat, mu_vec, rf, gamma)

    weights <-
      as.numeric(unb * port_markowitz)

    weights
  }

#' KanZhou-2Fund Portfolio Optimization
#'
#' Implements the KanZhou-2Fund portfolio optimization strategy.
#'
#' @param sigma_mat a pxp covariance matrix of asset returns.
#' @param mu_vec a numeric vector of length p, the expected returns.
#' @param rf a double, the assumed risk-free return. Default value is 0.
#' @param gamma an integer, the risk aversion parameter. Default value is 2.
#' @param is_size an integer, the sample size of the returns data.
#'
#' @return a numeric vector with the KanZhou-2Fund portfolio weights.
#'
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' n <- nrow(rets_m)
#' sigma_mat <- cov(rets_m)
#' mu_vec <- mu_estim_wrapper(rets_m, mu_estim_ml)
#' port_kz2f <- port_optim_kz2f(sigma_mat, mu_vec, is_size = n)
#'
#' @export
port_optim_kz2f <- function(sigma_mat,
                            mu_vec,
                            rf = 0,
                            gamma = 2,
                            is_size) {
  n <- is_size
  p <- dim(sigma_mat)[1]
  mu_vec_excess <- mu_vec - rf
  sigma_mat_inv <- solve(sigma_mat)
  sharpesquared <- as.numeric(mu_vec_excess %*% sigma_mat_inv %*% mu_vec_excess)
  c3 <- (n - p - 1) * (n - p - 4) / (n * (n - 2))
  cstar <- c3 * sharpesquared / (sharpesquared + (p / n))

  weights <- as.numeric((1 / gamma) * cstar * sigma_mat_inv %*% mu_vec_excess)

  return(weights)
}

#' KanZhou-3Fund Portfolio Optimization
#'
#' Implements the KanZhou-3Fund portfolio optimization strategy.
#'
#' @param sigma_mat a pxp covariance matrix of asset returns.
#' @param mu_vec a numeric vector of length p, the expected returns.
#' @param rf a double, the assumed risk-free return. Default value is 0.
#' @param gamma an integer, the risk aversion parameter. Default value is 2.
#' @param is_size an integer, the sample size of the returns data.
#'
#' @return a numeric vector with the KanZhou-3Fund portfolio weights.
#'
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' n <- nrow(rets_m)
#' sigma_mat <- cov(rets_m)
#' mu_vec <- mu_estim_wrapper(rets_m, mu_estim_ml)
#' port_kz3f <- port_optim_kz3f(sigma_mat, mu_vec, is_size = n)
#'
#' @export
port_optim_kz3f <- function(sigma_mat,
                            mu_vec,
                            rf = 0,
                            gamma = 2,
                            is_size) {
  n <- is_size
  p <- dim(sigma_mat)[1]
  ones <- rep.int(1, p)
  mu_vec_excess <- mu_vec - rf
  sigma_mat_inv <- solve(sigma_mat)
  c3 <- (n - p - 1) * (n - p - 4) / (n * (n - 2))
  mug <-
    as.numeric(mu_vec_excess %*% sigma_mat_inv %*% ones) / as.numeric(ones %*% sigma_mat_inv %*%
      ones)
  psisquared <-
    as.numeric(mu_vec_excess %*% sigma_mat_inv %*% mu_vec_excess) - (as.numeric((mu_vec_excess %*% sigma_mat_inv %*%
      ones)^2) / as.numeric(ones %*% sigma_mat_inv %*% ones))

  weights <-
    as.numeric((c3 / gamma) * ((psisquared / (psisquared + p / n)) * sigma_mat_inv %*%
      mu_vec_excess + ((p / n) / (psisquared + p / n)) * mug * sigma_mat_inv %*% ones
    ))

  weights
}

#' TuZhou-Markowitz Portfolio Optimization
#'
#' Implements the TuZhou-Markowitz portfolio optimization strategy.
#'
#' @param sigma_mat a pxp covariance matrix of asset returns.
#' @param mu_vec a numeric vector of length p, the expected returns.
#' @param rf a double, the assumed risk-free return. Default value is 0.
#' @param gamma an integer, the risk aversion parameter. Default value is 2.
#' @param is_size an integer, the sample size of the returns data.
#'
#' @return a numeric vector with the TuZhou-Markowitz portfolio weights.
#'
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' n <- nrow(rets_m)
#' sigma_mat <- cov(rets_m)
#' mu_vec <- mu_estim_wrapper(rets_m, mu_estim_ml)
#' port_tz_markowitz <- port_optim_tz_markowitz(sigma_mat, mu_vec, is_size = n)
#'
#' @export
port_optim_tz_markowitz <-
  function(sigma_mat,
           mu_vec,
           rf = 0,
           gamma = 2,
           is_size) {
    n <- is_size
    p <- dim(sigma_mat)[1]
    mu_vec_excess <- mu_vec - rf
    sigma_mat_inv <- solve(sigma_mat)
    port_naive <- port_optim_naive(sigma_mat)
    port_markowitz_unb <-
      port_optim_markowitz_unb(sigma_mat, mu_vec, rf, gamma, n)
    sharpesquared <- as.numeric(mu_vec_excess %*% sigma_mat_inv %*% mu_vec_excess)
    pi1 <-
      as.numeric(
        as.numeric(port_naive %*% sigma_mat %*% port_naive) - ((2 / gamma) * as.numeric(port_naive %*% mu_vec_excess)) + ((1 /
          (gamma^2)) * sharpesquared)
      )
    c1 <- ((n - 2) * (n - p - 2)) / ((n - p - 1) * (n - p - 4))
    pi2 <-
      as.numeric(((1 / (gamma^2)) * (c1 - 1) * sharpesquared) + (c1 / (gamma^
        2)) * (p / n))
    delta_tz_markowitz <- pi1 / (pi1 + pi2)

    weights <-
      as.numeric((1 - delta_tz_markowitz) * port_naive + delta_tz_markowitz * port_markowitz_unb)

    weights
  }


#' TuZhou-KanZhou-3Fund Portfolio Optimization
#'
#' Combines the TuZhou and KanZhou-3Fund portfolio optimization strategies.
#'
#' @param sigma_mat a pxp covariance matrix of asset returns.
#' @param mu_vec a numeric vector of length p, the expected returns.
#' @param rf a double, the assumed risk-free return. Default value is 0.
#' @param gamma an integer, the risk aversion parameter. Default value is 2.
#' @param is_size an integer, the sample size of the returns data.
#'
#' @return a numeric vector with the TuZhou-KanZhou-3Fund portfolio weights.
#'
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' n <- nrow(rets_m)
#' sigma_mat <- cov(rets_m)
#' mu_vec <- mu_estim_wrapper(rets_m, mu_estim_ml)
#' port_tz_kz3f <- port_optim_tz_kz3f(sigma_mat, mu_vec, is_size = n)
#'
#' @export
port_optim_tz_kz3f <- function(sigma_mat,
                               mu_vec,
                               rf = 0,
                               gamma = 2,
                               is_size) {
  n <- is_size
  p <- dim(sigma_mat)[1]
  ones <- rep.int(1, p)
  mu_vec_excess <- mu_vec - rf
  sigma_mat_inv <- solve(sigma_mat)
  port_naive <- port_optim_naive_slim(p)
  port_KZ3f <- port_optim_kz3f(sigma_mat, mu_vec, rf, gamma, n)
  mug <-
    as.numeric(mu_vec_excess %*% sigma_mat_inv %*% ones) / as.numeric(ones %*% sigma_mat_inv %*%
      ones)
  psisquared <-
    as.numeric(mu_vec_excess %*% sigma_mat_inv %*% mu_vec_excess) - (as.numeric((mu_vec_excess %*% sigma_mat_inv %*%
      ones)^2) / as.numeric(ones %*% sigma_mat_inv %*% ones))
  sharpesquared <- as.numeric(mu_vec_excess %*% sigma_mat_inv %*% mu_vec_excess)
  c1 <- ((n - 2) * (n - p - 2)) / ((n - p - 1) * (n - p - 4))
  pi1 <-
    as.numeric(
      as.numeric(port_naive %*% sigma_mat %*% port_naive) -
        ((2 / gamma) * as.numeric(port_naive %*% mu_vec_excess)) + ((1 / (gamma^2)) * sharpesquared)
    )
  pi3 <-
    as.numeric(((1 / (gamma^2)) * sharpesquared) - (1 / (c1 * (gamma^2))) *
      (sharpesquared - ((p / n) * psisquared)))
  pi13 <- as.numeric(((1 / (gamma^2)) * sharpesquared) -
    ((1 / gamma) * as.numeric(port_naive %*% mu_vec_excess)) +
    (1 / (gamma * c1)) * ((
      psisquared * as.numeric(port_naive %*% mu_vec_excess) + (1 - psisquared) * mug * as.numeric(port_naive %*% ones)
    ) - (1 / gamma) * (
      psisquared * as.numeric(mu_vec_excess %*% sigma_mat_inv %*% mu_vec_excess) +
        (1 - psisquared) * mug * as.numeric(mu_vec_excess %*% sigma_mat_inv %*% ones)
    )
    ))
  delta_TZKZ3f <- (pi1 - pi13) / (pi1 + pi3 - 2 * pi13)

  weights <-
    as.numeric((1 - delta_TZKZ3f) * port_naive + delta_TZKZ3f * port_KZ3f)

  return(weights)
}
