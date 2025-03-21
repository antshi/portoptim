#' Simulated Expected Returns
#'
#' Simulates a vector of expected returns with specific properties.
#'
#' @param p an integer, indicating the number of assets.
#' @param min_exp_val a double, indicating the minimum possible expected value.
#' Default value is 0.00.
#' @param max_exp_val a double, indicating the maximum possible expected value.
#' Default value is 0.50.
#'
#' @return a numeric vector with the simulated expected returns.
#'
#' @details The expected values are drown from an uniform distribution.
#' @examples
#' mu_vec <- mu_sim(p = 100, min_exp_val = 0.00, max_exp_val = 0.50)
#'
#' @export mu_sim
mu_sim <- function(
    p,
    min_exp_val = 0.00,
    max_exp_val = 0.50) {
  mu_vec <- stats::runif(p, min = min_exp_val, max = max_exp_val)

  mu_vec
}

#' Maximum-Likelihood Expected Returns
#'
#' Estimates the expected returns with the Maximum-Likelihood (ML) estimator.
#'
#' @param rets an nxp matrix of the stock returns.
#'
#' @return a numeric vector with the estimated expected returns.
#'
#' @details The ML estimator for the expected returns of an \eqn{n\times p} matrix with stock returns \eqn{R}
#' is calculated with the following formula:
#' \deqn{E[R]=\frac{1}{n}\sum_{i=1}^{n}R_{i}}.
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' mu_ml <- mu_estim_ml(rets_m)
#' head(mu_ml)
#' @export mu_estim_ml
mu_estim_ml <- function(rets) {
  rets <- as.matrix(rets)
  mu_vec_ml <- colMeans(rets, na.rm = TRUE)
  names(mu_vec_ml) <- colnames(rets)

  mu_vec_ml
}

#'
#' Single-Factor Expected Returns
#'
#' Estimates the expected returns with a single-factor model, specifically the Capital Asset Pricing Model (CAPM).
#'
#' @param rets an nxp matrix, the stock returns. n stands for the observation length, p for the number of stocks.
#' @param market_rets_excess an nx1 vector, the excess market return.
#' @param rf a risk-free rate of return. Default value is rf=0.
#'
#' @return a numeric vector with the estimated expected returns.
#'
#' @details The Single-Factor estimator for each expected return \eqn{i} of an \eqn{n\times p} matrix
#' with stock returns \eqn{R} is
#' calculated with the following formula:
#' \deqn{E[R_{i}]=r_{f} + \widehat{\beta}_{i}\left(E[R_{M}]-r_{f}\right),}
#' where \eqn{r_{f}} is the risk-free rate of return, \eqn{E[R_{M}]} is the expected market return and
#' \eqn{\widehat{{\beta}}_{i}} are the estimated beta coefficients from the linear regression
#' \deqn{R_{ij} - r_{f} = \alpha_{i} + \beta_{i}\left(R_{Mj}-r_{f}\right) + \epsilon_{ij}}
#' for \eqn{i=1\cdots p} and \eqn{j=1\cdots n}.
#' @examples
#' data(prices_m)
#' data(ff_5factors)
#' rets_m <- calc_rets(prices_m)
#' mu_capm <- mu_estim_capm(rets_m, ff_5factors$Mkt.RF, ff_5factors$RF)
#' head(mu_capm)
#' @export mu_estim_capm
mu_estim_capm <- function(rets, market_rets_excess, rf = 0) {
  if (length(market_rets_excess) != nrow(rets)) {
    stop("The length of market_rets_excess must be equal to the number of observations in rets.")
  }
  if (length(rf) != 1) {
    if (length(rf) != nrow(rets)) {
      stop("The length of rf must be equal to the number of observations in rets.")
    }
  }
  rets <- as.matrix(rets)
  rets_excess <- rets - rf

  capm_models <-
    stats::lm(market_rets_excess ~ rets_excess)
  capm_betas <- as.numeric(capm_models$coefficients)[-1]

  mu_vec_capm <-
    as.numeric(mean(rf)) + capm_betas * mean(rets_excess)

  names(mu_vec_capm) <- colnames(rets)

  mu_vec_capm
}

#' Factor-Model Expected Returns
#'
#' Estimates the expected returns with a factor model.
#' The factors can be, for example, the Fama-French 5 Research Factors.
#'
#' @param rets an nxp data matrix.
#' @param rf a numeric vector or double for the risk-free return.
#' @param factors an nxf data frame of factors, e.g. the Fama-French 5 Research Factors.
#'
#' @return a numeric vector with the estimated expected returns.
#'
#' @details The expected returns are calculated according to a factor model, estimated with an OLS regression.
#' For a set of factors F, the expected returns for an \eqn{n\times p} data matrix X are defined as
#' \deqn{\hat{\mu_vec}= \hat{B}'E[F],}
#' where \eqn{E[F]=\frac{1}{n}\sum_{i=1}^{n}F_{i}} and \eqn{\hat{B}} are the estimated
#' beta coefficients from the linear regression.
#'
#' @examples
#' data(prices_m)
#' data(ff_5factors)
#' rets_m <- calc_rets(prices_m)
#' mu_fm <- mu_estim_fm(rets_m, ff_5factors$RF, ff_5factors[, -ncol(ff_5factors)])
#' head(mu_fm)
#' @export mu_estim_fm
mu_estim_fm <- function(rets, rf, factors) {
  if (nrow(factors) != nrow(rets)) {
    stop("The observation length of factors must be equal to the number of observations in rets.")
  }
  if (length(rf) != 1) {
    if (length(rf) != nrow(rets)) {
      stop("The length of rf must be equal to the number of observations in rets.")
    }
  }
  rets <- as.matrix(rets)
  rets_excess <- rets - rf
  factors <- as.data.frame(factors)
  mu_vec_fm <- c()

  for (i in seq_len(ncol(rets_excess))) {
    data_lm <- as.matrix(cbind(rets_excess[, i], factors))
    model_lm <-
      stats::lm(data_lm[, 1] ~ data_lm[, -1])
    betas <-
      model_lm$coefficients[-1]
    mu_vec_fm[i] <- mean(rf) + as.numeric(betas %*% colMeans(factors))
  }

  mu_vec_fm
}

#' Linear-Shrinkage Expected Returns
#'
#' Estimates the expected returns with a linear shrinkage estimator towards the
#' expected return of a global minimum-variance (GMV) portfolio.
#'
#' @param rets an nxp matrix, the stock returns. n stands for the observation length, p for the number of stocks.
#'
#' @return a numeric vector with the estimated expected returns.
#'
#' @details The Linear-Shrinkage estimator for the expected return vector of an \eqn{n\times p} matrix
#' with stock returns \eqn{R} is
#' calculated with the following formula:
#' \deqn{E[R]= \left(1-\delta\right)\widehat{\mu_vec} + \delta\widehat{\mu_vec}_{g}{1}_{p},}
#' where \eqn{\widehat{\mu_vec}} is the sample mean vector,
#' \eqn{\widehat{\mu_vec}_{g}} is the shrinkage target (here, the grand mean), \eqn{{1}_{p}} is a 1xp vector of ones,
#' and \eqn{\delta} is the optimal shrinkage intensity as in \insertCite{Stein_1956;textual}{portoptim}.
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' mu_shrink <- mu_estim_shrink(rets_m)
#' head(mu_shrink)
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited
#'
#' @export mu_estim_shrink
mu_estim_shrink <- function(rets) {
  rets <- as.matrix(rets)
  n <- dim(rets)[1]
  p <- dim(rets)[2]

  mu_vec_ml <- mu_estim_ml(rets)
  grand_mean <- mean(mu_vec_ml)
  sigma2 <- grand_mean * (1 - grand_mean) / n
  mu_vec_shrink <- grand_mean +
    (1 - ((p - 3) * sigma2 / sum((mu_vec_ml - grand_mean)^2))) * (mu_vec_ml - grand_mean)

  mu_vec_shrink
}

#' Wrapper Function for Expected Returns Estimation I
#'
#' Estimates the expected returns of a dataset
#' according to the user-defined function.
#'
#' @param rets an nxp stock returns matrix.
#' @param est_func an estimation function.
#' @param ... additional parameters, parsed to est_func.
#'
#' @return a numeric vector with the estimated expected returns.
#'
#' @examples
#' data(prices_m)
#' data(ff_5factors)
#' rets_m <- calc_rets(prices_m)
#' mu_shrink <- mu_estim_wrapper(rets_m, mu_estim_shrink)
#' mu_fm <- mu_estim_wrapper(rets_m, mu_estim_fm, ff_5factors$RF, ff_5factors[, -ncol(ff_5factors)])
#'
#' @export mu_estim_wrapper
#'
mu_estim_wrapper <- function(rets, est_func, ...) {
  mu_vec <- est_func(rets, ...)
  mu_vec
}
