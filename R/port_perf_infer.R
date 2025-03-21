#' Nonparametric Bootstrap
#'
#' Performs a nonparametric bootstrap on returns time series and a predefined performance measure.
#'
#' @param rets an nxp stock returns matrix.
#' @param boot_reps an integer, the number of bootstrap repetitions.
#' @param sample_size an integer, the sample size for the nonparametric bootstrap.
#' @param type a function, the performance measure to be calculated. Default is calc_sr.
#'
#' @return a Bxp data matrix with the bootstrapped values of the performance measure from type.
#'
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' srs_boot <- boot_nonparam(rets_m, boot_reps = 10, sample_size = 50, type = calc_sr)
#'
#' @export boot_nonparam
#'
boot_nonparam <- function(rets, boot_reps, sample_size, type = calc_sr) {
  boot_mat <- matrix(NA, boot_reps, ncol(rets))
  for (i in 1:boot_reps) {
    ret_ind <- sample(seq_len(nrow(rets)), sample_size, replace = TRUE)
    boot_mat[i, ] <- apply(rets[ret_ind, ], 2, type)
  }

  return(boot_mat)
}

#' Parametric Bootstrap
#'
#' Performs a parametric bootstrap on returns time series (with an assumed normal distribution)
#' and predefined performance measure.
#'
#' @param rets an nxp stock returns matrix.
#' @param boot_reps an integer, the number of bootstrap repetitions.
#' @param sample_size an integer, the sample size for the parametric bootstrap.
#' @param type a function, the performance measure to be calculated. Default is sr_calc.
#'
#' @return a Bxp data matrix with the bootstrapped values of the performance measure from type.
#'
#'
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' srs_boot <- boot_param(rets_m, boot_reps = 10, sample_size = 50, type = calc_sr)
#'
#' @export boot_param
#'
boot_param <- function(rets, boot_reps, sample_size, type = calc_sr) {
  boot_mat <- matrix(NA, boot_reps, ncol(rets))

  for (b in 1:boot_reps) {
    boot_rets <- apply(rets, 2, function(x) {
      mu_vec <- mean(x)
      sigma <- stats::sd(x)
      stats::rnorm(sample_size, mean = mu_vec, sd = sigma)
    })
    boot_mat[b, ] <- apply(boot_rets, 2, type)
  }
  return(boot_mat)
}

#' F-Test
#' Performs a two-sided F-test to compare the variances of two samples from normal populations.
#'
#' @param rets an nx2 matrix with two returns time series.
#' n stands for the observation length.
#' Only two return time series can be compared simultaneously.
#'
#' @return a list
#' \itemize{
#' \item vars, the estimated variances.
#' \item vars-ratio, the ratio of the estimated variances.
#' \item p-value, the estimated p-value.
#' }
#'
#' @details The null hypothesis is that the ratio of the variances of the populations
#' from which the returns were drawn is equal to 1,
#' that is, the variances are equal.
#' The F-test requires that the data come from a bivariate normal distribution
#' without correlation and dependence over time.
#' If the data are correlated, have tails heavier than the normal distribution, or are dependent over time,
#' the test is not valid \insertCite{LedoitWolf_2011}{portoptim}.
#'
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' sigma_mat <- cov(rets_m)
#' mu_vec <- mu_estim_ml(rets_m)
#' port_naive <- port_optim_naive(sigma_mat)
#' port_tang <- port_optim_tang(sigma_mat, mu_vec)
#' port_rets <- cbind(calc_port_rets(port_naive, rets_m), calc_port_rets(port_tang, rets_m))
#' ftest_res <- f_test_rets(port_rets)
#'
#' @export f_test_rets
f_test_rets <- function(rets) {
  rets <- as.matrix(rets)
  result <- stats::var.test(rets[, 1], rets[, 2])
  results_list <-
    list(
      as.numeric(c(stats::var(rets[, 1]), stats::var(rets[, 2]))),
      as.numeric(result$estimate),
      as.numeric(result$p.value)
    )
  names(results_list) <- c("vars", "vars-ratio", "p-value")
  results_list
}

#' HAC Inference
#'
#' Performs a HAC (Heteroskedasticity and Autocorrelation Consistent) inference
#'
#' @param rets an nx2 matrix with two returns time series.
#' n stands for the observation length.
#' Only two return time series can be compared simultaneously.
#' @param type a character, indicating the performance measure (parameter) of interest.
#' Possible values are type=c("var", "sharpe"). Default value is type="var".
#' @param digits an integer, indicating how many digits the results are to be rounded to. Default value is 5.
#'
#' @return a list
#' \itemize{
#' \item param-diff, the difference between the estimated parameters in type
#' (first minus second time series).
#' For type="var", the function returns the differences in log-variances.
#' \item se, the estimated standard error of the parameter in type.
#' \item p-value, the estimated p-value.
#' }
#'
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' sigma_mat <- cov(rets_m)
#' mu_vec <- mu_estim_ml(rets_m)
#' port_naive <- port_optim_naive(sigma_mat)
#' port_tang <- port_optim_tang(sigma_mat, mu_vec)
#' port_rets <- cbind(calc_port_rets(port_naive, rets_m), calc_port_rets(port_tang, rets_m))
#' hac_infer_res <- hac_infer(port_rets)
#'
#' @export hac_infer
#'
hac_infer <- function(rets, type = "var", digits = 5) {
  if (type == "var") {
    result <- hac.inference.log.var(rets, digits)
  } else if (type == "sharpe") {
    result <- hac.inference(rets, digits)
  } else {
    stop("Type unknown.")
  }
  results_list <-
    list(
      as.numeric(result$Difference),
      as.numeric(result$Standard.Errors[2]),
      as.numeric(result$p.Values[2])
    )
  names(results_list) <-
    c("param-diff", "se", "p-value")
  results_list
}

#' Bootstrap Inference
#'
#' Performs a studentized Bootstrap inference
#'
#' @param rets an nx2 matrix with two returns time series.
#' n stands for the observation length.
#' Only two return time series can be compared simultaneously.
#' @param boot_reps an integer, the number of bootstrap repetitions.
#' @param blocks_n an integer, the number of blocks for the circular block bootstrap.
#' Default value is 1.
#' @param type a character, indicating the performance measure (parameter) of interest.
#' Possible values are type=c("var", "sharpe"). Default value is type="var".
#' @param digits an integer, indicating how many digits the results are to be rounded to. Default value is 5.
#'
#' @return a list
#' \itemize{
#' \item param-diff, the difference between the estimated parameters in type
#' (first minus second time series).
#' For type="var", the function returns the differences in log-variances.
#' \item se, the estimated standard error of the parameter in type.
#' \item p-value, the estimated p-value.
#' }
#'
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' sigma_mat <- cov(rets_m)
#' mu_vec <- mu_estim_ml(rets_m)
#' port_naive <- port_optim_naive(sigma_mat)
#' port_tang <- port_optim_tang(sigma_mat, mu_vec)
#' port_rets <- cbind(calc_port_rets(port_naive, rets_m), calc_port_rets(port_tang, rets_m))
#' boot_infer_res <- boot_infer(port_rets)
#'
#' @export boot_infer
#'
boot_infer <- function(
    rets,
    boot_reps = 100,
    blocks_n = 1,
    type = "var",
    digits = 5) {
  rets <- as.matrix(rets)
  if (type == "var") {
    result <- boot.time.inference.log.var(rets,
      b = blocks_n,
      M = boot_reps,
      digits = digits
    )
  } else if (type == "sharpe") {
    result <- boot.time.inference(rets,
      b = blocks_n,
      M = boot_reps,
      digits = digits
    )
  } else {
    stop("Type unknown.")
  }
  results_list <-
    list(
      as.numeric(result$Difference),
      as.numeric(result$Standard.Error),
      as.numeric(result$p.Value)
    )
  names(results_list) <-
    c("param-diff", "se", "p-value")
  results_list
}

#' @noRd
hac.inference <- function(ret, digits = 3) {
  ret1 <- ret[, 1]
  ret2 <- ret[, 2]
  mu1.hat <- mean(ret1)
  mu2.hat <- mean(ret2)
  sig1.hat <- stats::sd(ret1)
  sig2.hat <- stats::sd(ret2)
  SR1.hat <- mu1.hat / sig1.hat
  SR2.hat <- mu2.hat / sig2.hat
  SRs <- round(c(SR1.hat, SR2.hat), digits)
  diff <- SR1.hat - SR2.hat
  names(SRs) <- c("SR1.hat", "SR2.hat")
  se <- compute.se.Parzen(ret)
  se.pw <- compute.se.Parzen.pw(ret)
  SEs <- round(c(se, se.pw), digits)
  names(SEs) <- c("HAC", "HAC.pw")
  PV <- 2 * stats::pnorm(-abs(diff) / se)
  PV.pw <- 2 * stats::pnorm(-abs(diff) / se.pw)
  PVs <- round(c(PV, PV.pw), digits)
  names(PVs) <- c("HAC", "HAC.pw")
  list(
    Sharpe.Ratios = SRs,
    Difference = round(diff, digits),
    Standard.Errors = SEs,
    p.Values = PVs
  )
}

#' @noRd
hac.inference.log.var <- function(ret, digits = 3) {
  ret1 <- ret[, 1]
  ret2 <- ret[, 2]
  var1.hat <- stats::var(ret1)
  var2.hat <- stats::var(ret2)
  log.var1.hat <- log(var1.hat)
  log.var2.hat <- log(var2.hat)
  VARs <- round(c(var1.hat, var2.hat), digits)
  LogVARs <- round(c(log.var1.hat, log.var2.hat), digits)
  diff <- log.var1.hat - log.var2.hat
  names(VARs) <- c("Var1.hat", "Var2.hat")
  names(LogVARs) <- c("LogVar1.hat", "LogVar2.hat")
  Delta.hat <- diff
  se <- compute.se.Parzen.log.var(ret)
  se.pw <- compute.se.Parzen.pw.log.var(ret)
  SEs <- round(c(se, se.pw), digits)
  names(SEs) <- c("HAC", "HAC.pw")
  PV <- 2 * stats::pnorm(-abs(diff) / se)
  PV.pw <- 2 * stats::pnorm(-abs(diff) / se.pw)
  PVs <- round(c(PV, PV.pw), digits)
  names(PVs) <- c("HAC", "HAC.pw")
  list(
    Variances = VARs,
    Log.Variances = LogVARs,
    Difference = round(
      diff,
      digits
    ),
    Standard.Errors = SEs,
    p.Values = PVs
  )
}

#' @noRd
boot.time.inference <-
  function(ret,
           b,
           M,
           Delta.null = 0,
           digits = 4) {
    T <- length(ret[, 1])
    l <- floor(T / b)
    Delta.hat <- sharpe.ratio.diff(ret)
    d <- abs(Delta.hat - Delta.null) / compute.se.Parzen.pw(ret)
    p.value <- 1
    for (m in (1:M)) {
      ret.star <- ret[cbb.sequence(T, b), ]
      Delta.hat.star <- sharpe.ratio.diff(ret.star)
      ret1.star <- ret.star[, 1]
      ret2.star <- ret.star[, 2]
      mu1.hat.star <- mean(ret1.star)
      mu2.hat.star <- mean(ret2.star)
      gamma1.hat.star <- mean(ret1.star^2)
      gamma2.hat.star <- mean(ret2.star^2)
      gradient <- rep(0, 4)
      gradient[1] <- gamma1.hat.star / (gamma1.hat.star - mu1.hat.star^
        2)^1.5
      gradient[2] <- -gamma2.hat.star / (gamma2.hat.star - mu2.hat.star^
        2)^1.5
      gradient[3] <- -0.5 * mu1.hat.star / (gamma1.hat.star -
        mu1.hat.star^2)^1.5
      gradient[4] <- 0.5 * mu2.hat.star / (gamma2.hat.star - mu2.hat.star^
        2)^1.5
      y.star <- data.frame(
        ret1.star - mu1.hat.star,
        ret2.star -
          mu2.hat.star,
        ret1.star^2 - gamma1.hat.star,
        ret2.star^2 -
          gamma2.hat.star
      )
      Psi.hat.star <- matrix(0, 4, 4)
      for (j in (1:l)) {
        zeta.star <- b^0.5 * colMeans(y.star[((j - 1) * b + 1):(j *
          b), ])
        Psi.hat.star <- Psi.hat.star + zeta.star %*% t(zeta.star)
      }
      Psi.hat.star <- Psi.hat.star / l
      se.star <- as.numeric(sqrt(t(gradient) %*% Psi.hat.star %*%
        gradient / T))
      d.star <- abs(Delta.hat.star - Delta.hat) / se.star
      if (d.star >= d) {
        p.value <- p.value + 1
      }
    }
    p.value <- p.value / (M + 1)
    list(
      Difference = round(Delta.hat, digits),
      p.Value = round(
        p.value,
        digits
      )
    )
  }

#' @noRd
boot.time.inference.log.var <- function(ret, b, M, Delta.null = 0, digits = 3) {
  T <- length(ret[, 1])
  l <- floor(T / b)
  Delta.hat <- log.var.diff(ret)
  d <- abs(Delta.hat - Delta.null) / compute.se.Parzen.pw.log.var(ret)
  p.value <- 1
  for (m in (1:M)) {
    ret.star <- ret[cbb.sequence(T, b), ]
    Delta.hat.star <- log.var.diff(ret.star)
    ret1.star <- ret.star[, 1]
    ret2.star <- ret.star[, 2]
    mu1.hat.star <- mean(ret1.star)
    mu2.hat.star <- mean(ret2.star)
    gamma1.hat.star <- mean(ret1.star^2)
    gamma2.hat.star <- mean(ret2.star^2)
    gradient <- rep(0, 4)
    gradient[1] <- -2 * mu1.hat.star / (gamma1.hat.star - mu1.hat.star^2)
    gradient[2] <- 2 * mu2.hat.star / (gamma2.hat.star - mu2.hat.star^2)
    gradient[3] <- 1 / (gamma1.hat.star - mu1.hat.star^2)
    gradient[4] <- -1 / (gamma2.hat.star - mu2.hat.star^2)
    y.star <- data.frame(ret1.star - mu1.hat.star, ret2.star -
      mu2.hat.star, ret1.star^2 - gamma1.hat.star, ret2.star^2 -
      gamma2.hat.star)
    Psi.hat.star <- matrix(0, 4, 4)
    for (j in (1:l)) {
      zeta.star <- b^0.5 * colMeans(y.star[((j - 1) * b + 1):(j *
        b), ])
      Psi.hat.star <- Psi.hat.star + zeta.star %*% t(zeta.star)
    }
    Psi.hat.star <- Psi.hat.star / l
    Psi.hat.star <- (T / (T - 4)) * Psi.hat.star
    se.star <- as.numeric(sqrt(t(gradient) %*% Psi.hat.star %*%
      gradient / T))
    d.star <- abs(Delta.hat.star - Delta.hat) / se.star
    if (d.star >= d) {
      p.value <- p.value + 1
    }
  }
  p.value <- p.value / (M + 1)
  list(Difference = round(Delta.hat, digits), p.Value = round(
    p.value,
    digits
  ))
}

#' @noRd
block.size.calibrate <-
  function(ret,
           b.vec = c(1, 3, 6, 10),
           alpha = 0.05,
           M = 199,
           K = 1000,
           b.av = 5,
           T.start = 50) {
    b.len <- length(b.vec)
    emp.reject.probs <- rep(0, b.len)
    Delta.hat <- sharpe.ratio.diff(ret)
    ret1 <- ret[, 1]
    ret2 <- ret[, 2]
    T <- length(ret1)
    Var.data <- matrix(0, T.start + T, 2)
    Var.data[1, 1] <- ret[1, 1]
    Var.data[1, 2] <- ret[1, 2]
    Delta.hat <- sharpe.ratio.diff(ret)
    fit1 <- stats::lm(ret1[2:T] ~ ret1[1:(T - 1)] + ret2[1:(T - 1)])
    fit2 <- stats::lm(ret2[2:T] ~ ret1[1:(T - 1)] + ret2[1:(T - 1)])
    coef1 <- as.numeric(fit1$coef)
    coef2 <- as.numeric(fit2$coef)
    resid.mat <- cbind(as.numeric(fit1$resid), as.numeric(fit2$resid))
    for (k in (1:K)) {
      resid.mat.star <- rbind(c(0, 0), resid.mat[sb.sequence(T -
        1, b.av, T.start + T - 1), ])
      for (t in (2:(T.start + T))) {
        Var.data[t, 1] <- coef1[1] + coef1[2] * Var.data[t -
          1, 1] + coef1[3] * Var.data[t - 1, 2] + resid.mat.star[
          t,
          1
        ]
        Var.data[t, 2] <- coef2[1] + coef2[2] * Var.data[t -
          1, 1] + coef2[3] * Var.data[t - 1, 2] + resid.mat.star[
          t,
          2
        ]
      }
      Var.data.trunc <- Var.data[(T.start + 1):(T.start + T), ]
      for (j in (1:b.len)) {
        p.Value <- boot.time.inference(
          Var.data.trunc, b.vec[j],
          M, Delta.hat
        )$p.Value
        if (p.Value <= alpha) {
          emp.reject.probs[j] <- emp.reject.probs[j] + 1
        }
      }
    }
    emp.reject.probs <- emp.reject.probs / K
    b.order <- order(abs(emp.reject.probs - alpha))
    b.opt <- b.vec[b.order[1]]
    b.vec.with.probs <- rbind(b.vec, emp.reject.probs)
    colnames(b.vec.with.probs) <- rep("", length(b.vec))
    list(
      Empirical.Rejection.Probs = b.vec.with.probs,
      b.optimal = b.opt
    )
  }

#' @noRd
block.size.calibrate.log.var <-
  function(ret,
           b.vec = c(1, 3, 6, 10),
           alpha = 0.05,
           M = 199,
           K = 1000,
           b.av = 5,
           T.start = 50) {
    b.len <- length(b.vec)
    emp.reject.probs <- rep(0, b.len)
    Delta.hat <- log.var.diff(ret)
    ret1 <- ret[, 1]
    ret2 <- ret[, 2]
    T <- length(ret1)
    Var.data <- matrix(0, T.start + T, 2)
    Var.data[1, ] <- ret[1, ]
    Delta.hat <- log.var.diff(ret)
    fit1 <- stats::lm(ret1[2:T] ~ ret1[1:(T - 1)] + ret2[1:(T - 1)])
    fit2 <- stats::lm(ret2[2:T] ~ ret1[1:(T - 1)] + ret2[1:(T - 1)])
    coef1 <- as.numeric(fit1$coef)
    coef2 <- as.numeric(fit2$coef)
    resid.mat <- cbind(as.numeric(fit1$resid), as.numeric(fit2$resid))
    for (k in (1:K)) {
      resid.mat.star <- rbind(c(0, 0), resid.mat[sb.sequence(T -
        1, b.av, T.start + T - 1), ])
      for (t in (2:(T.start + T))) {
        Var.data[t, 1] <- coef1[1] + coef1[2] * Var.data[t -
          1, 1] + coef1[3] * Var.data[t - 1, 2] + resid.mat.star[
          t,
          1
        ]
        Var.data[t, 2] <- coef2[1] + coef2[2] * Var.data[t -
          1, 1] + coef2[3] * Var.data[t - 1, 2] + resid.mat.star[
          t,
          2
        ]
      }
      Var.data.trunc <- Var.data[(T.start + 1):(T.start + T), ]
      for (j in (1:b.len)) {
        p.Value <- boot.time.inference.log.var(
          Var.data.trunc,
          b.vec[j], M, Delta.hat
        )$p.Value
        if (p.Value <= alpha) {
          emp.reject.probs[j] <- emp.reject.probs[j] + 1
        }
      }
    }
    emp.reject.probs <- emp.reject.probs / K
    b.order <- order(abs(emp.reject.probs - alpha))
    b.opt <- b.vec[b.order[1]]
    b.vec.with.probs <- rbind(b.vec, emp.reject.probs)
    colnames(b.vec.with.probs) <- rep("", length(b.vec))
    list(
      Empirical.Rejection.Probs = b.vec.with.probs,
      b.optimal = b.opt
    )
  }

#' @noRd
sharpe.ratio.diff <- function(ret) {
  ret1 <- ret[, 1]
  ret2 <- ret[, 2]
  mu1.hat <- mean(ret1)
  mu2.hat <- mean(ret2)
  sig1.hat <- stats::sd(ret1)
  sig2.hat <- stats::sd(ret2)
  diff <- mu1.hat / sig1.hat - mu2.hat / sig2.hat
  diff
}

#' @noRd
log.var.diff <- function(ret) {
  log(stats::var(ret[, 1])) - log(stats::var(ret[, 2]))
}

#' @noRd
compute.Psi.hat <- function(V.hat) {
  T <- length(V.hat[, 1])
  alpha.hat <- compute.alpha.hat(V.hat)
  S.star <- 2.6614 * (alpha.hat * T)^0.2
  Psi.hat <- compute.Gamma.hat(V.hat, 0)
  j <- 1
  while (j < S.star) {
    Gamma.hat <- compute.Gamma.hat(V.hat, j)
    Psi.hat <- Psi.hat + kernel.Parzen(j / S.star) * (Gamma.hat + t(Gamma.hat))
    j <- j + 1
  }
  Psi.hat <- (T / (T - 4)) * Psi.hat
  Psi.hat
}

#' @noRd
compute.alpha.hat <- function(V.hat) {
  dimensions <- dim(V.hat)
  T <- dimensions[1]
  p <- dimensions[2]
  numerator <- 0
  denominator <- 0
  for (i in (1:p)) {
    fit <- stats::ar(V.hat[, i], 0, 1, method = "ols")
    rho.hat <- as.numeric(fit[2])
    sig.hat <- sqrt(as.numeric(fit[3]))
    numerator <- numerator + 4 * rho.hat^2 * sig.hat^4 / (1 - rho.hat)^
      8
    denominator <- denominator + sig.hat^4 / (1 - rho.hat)^4
  }
  numerator / denominator
}

#' @noRd
compute.Gamma.hat <- function(V.hat, j) {
  dimensions <- dim(V.hat)
  T <- dimensions[1]
  p <- dimensions[2]
  Gamma.hat <- matrix(0, p, p)
  if (j >= T) {
    stop("j must be smaller than the row dimension!")
  }
  for (i in ((j + 1):T)) {
    Gamma.hat <- Gamma.hat + V.hat[i, ] %*%
      t(V.hat[i - j, ])
  }
  Gamma.hat <- Gamma.hat / T
  Gamma.hat
}

#' @noRd
compute.V.hat <- function(ret) {
  ret1 <- ret[, 1]
  ret2 <- ret[, 2]
  V.hat <- cbind(
    ret1 - mean(ret1),
    ret2 - mean(ret2),
    ret1^2 -
      mean(ret1^2),
    ret2^2 - mean(ret2^2)
  )
  V.hat
}

#' @noRd
cbb.sequence <- function(T, b) {
  l <- floor(T / b)
  index.sequence <- c(1:T, 1:b)
  sequence <- rep(0, T)
  start.points <- sample(1:T, l, replace = T)
  for (j in (1:l)) {
    start <- start.points[j]
    sequence[((j - 1) * b + 1):(j * b)] <- index.sequence[start:(start +
      b - 1)]
  }
  sequence
}

#' @noRd
sb.sequence <- function(T, b.av, length = T) {
  index.sequence <- c(1:T, 1:T)
  sequence <- rep(0, length + T)
  current <- 0
  while (current < length) {
    start <- sample(1:T, 1)
    b <- stats::rgeom(1, 1 / b.av) + 1
    sequence[(current + 1):(current + b)] <- index.sequence[start:(start + b - 1)]
    current <- current + b
  }
  sequence[1:length]
}

#' @noRd
kernel.Parzen <- function(x) {
  if (abs(x) <= 0.5) {
    result <- 1 - 6 * x^2 + 6 * abs(x)^3
  } else if (abs(x) <= 1) {
    result <- 2 * (1 - abs(x))^3
  } else {
    result <- 0
  }
  result
}

#' @noRd
compute.se.Parzen <- function(ret) {
  ret1 <- ret[, 1]
  ret2 <- ret[, 2]
  T <- length(ret1)
  mu1.hat <- mean(ret1)
  mu2.hat <- mean(ret2)
  ret1.2 <- ret1^2
  ret2.2 <- ret2^2
  gamma1.hat <- mean(ret1.2)
  gamma2.hat <- mean(ret2.2)
  gradient <- rep(0, 4)
  gradient[1] <- gamma1.hat / (gamma1.hat - mu1.hat^2)^1.5
  gradient[2] <- -gamma2.hat / (gamma2.hat - mu2.hat^2)^1.5
  gradient[3] <- -0.5 * mu1.hat / (gamma1.hat - mu1.hat^2)^1.5
  gradient[4] <- 0.5 * mu2.hat / (gamma2.hat - mu2.hat^2)^1.5
  V.hat <- cbind(
    ret1 - mu1.hat,
    ret2 - mu2.hat,
    ret1.2 - gamma1.hat,
    ret2.2 - gamma2.hat
  )
  Psi.hat <- compute.Psi.hat(V.hat)
  se <- as.numeric(sqrt(t(gradient) %*% Psi.hat %*% gradient / T))
  se
}

#' @noRd
compute.se.Parzen.pw <- function(ret) {
  ret1 <- ret[, 1]
  ret2 <- ret[, 2]
  mu1.hat <- mean(ret1)
  mu2.hat <- mean(ret2)
  ret1.2 <- ret1^2
  ret2.2 <- ret2^2
  gamma1.hat <- mean(ret1.2)
  gamma2.hat <- mean(ret2.2)
  gradient <- rep(0, 4)
  gradient[1] <- gamma1.hat / (gamma1.hat - mu1.hat^2)^1.5
  gradient[2] <- -gamma2.hat / (gamma2.hat - mu2.hat^2)^1.5
  gradient[3] <- -0.5 * mu1.hat / (gamma1.hat - mu1.hat^2)^1.5
  gradient[4] <- 0.5 * mu2.hat / (gamma2.hat - mu2.hat^2)^1.5
  T <- length(ret1)
  V.hat <- cbind(
    ret1 - mu1.hat,
    ret2 - mu2.hat,
    ret1.2 - gamma1.hat,
    ret2.2 - gamma2.hat
  )
  A.ls <- matrix(0, 4, 4)
  V.star <- matrix(0, T - 1, 4)
  reg1 <- V.hat[1:T - 1, 1]
  reg2 <- V.hat[1:T - 1, 2]
  reg3 <- V.hat[1:T - 1, 3]
  reg4 <- V.hat[1:T - 1, 4]
  for (j in (1:4)) {
    fit <- stats::lm(V.hat[2:T, j] ~ -1 + reg1 + reg2 + reg3 + reg4)
    A.ls[j, ] <- as.numeric(fit$coef)
    V.star[, j] <- as.numeric(fit$resid)
  }
  svd.A <- svd(A.ls)
  d <- svd.A$d
  d.adj <- d
  for (i in (1:4)) {
    if (d[i] > 0.97) {
      d.adj[i] <- 0.97
    } else if (d[i] < -0.97) {
      d.adj[i] <- -0.97
    }
  }
  A.hat <- svd.A$u %*% diag(d.adj) %*% t(svd.A$v)
  D <- solve(diag(4) - A.hat)
  reg.mat <- rbind(reg1, reg2, reg3, reg4)
  for (j in (1:4)) {
    V.star[, j] <- V.hat[2:T, j] - A.hat[j, ] %*% reg.mat
  }
  Psi.hat <- compute.Psi.hat(V.star)
  Psi.hat <- D %*% Psi.hat %*% t(D)
  se <- as.numeric(sqrt(gradient %*% Psi.hat %*% gradient / T))
  se
}

#' @noRd
compute.se.Parzen.log.var <- function(ret) {
  ret1 <- ret[, 1]
  ret2 <- ret[, 2]
  T <- length(ret1)
  mu1.hat <- mean(ret1)
  mu2.hat <- mean(ret2)
  gamma1.hat <- mean(ret1^2)
  gamma2.hat <- mean(ret2^2)
  gradient <- rep(0, 4)
  gradient[1] <- -2 * mu1.hat / (gamma1.hat - mu1.hat^2)
  gradient[2] <- 2 * mu2.hat / (gamma2.hat - mu2.hat^2)
  gradient[3] <- 1 / (gamma1.hat - mu1.hat^2)
  gradient[4] <- -1 / (gamma2.hat - mu2.hat^2)
  V.hat <- compute.V.hat(ret)
  Psi.hat <- compute.Psi.hat(V.hat)
  se <- as.numeric(sqrt(t(gradient) %*% Psi.hat %*% gradient / T))
  se
}

#' @noRd
compute.se.Parzen.pw.log.var <- function(ret) {
  ret1 <- ret[, 1]
  ret2 <- ret[, 2]
  mu1.hat <- mean(ret1)
  mu2.hat <- mean(ret2)
  gamma1.hat <- mean(ret1^2)
  gamma2.hat <- mean(ret2^2)
  gradient <- rep(0, 4)
  gradient[1] <- -2 * mu1.hat / (gamma1.hat - mu1.hat^2)
  gradient[2] <- 2 * mu2.hat / (gamma2.hat - mu2.hat^2)
  gradient[3] <- 1 / (gamma1.hat - mu1.hat^2)
  gradient[4] <- -1 / (gamma2.hat - mu2.hat^2)
  T <- length(ret1)
  V.hat <- compute.V.hat(ret)
  A.ls <- matrix(0, 4, 4)
  V.star <- matrix(0, T - 1, 4)
  reg1 <- V.hat[1:T - 1, 1]
  reg2 <- V.hat[1:T - 1, 2]
  reg3 <- V.hat[1:T - 1, 3]
  reg4 <- V.hat[1:T - 1, 4]
  for (j in (1:4)) {
    fit <- stats::lm(V.hat[2:T, j] ~ -1 + reg1 + reg2 + reg3 + reg4)
    A.ls[j, ] <- as.numeric(fit$coef)
    V.star[, j] <- as.numeric(fit$resid)
  }
  svd.A <- svd(A.ls)
  d <- svd.A$d
  d.adj <- d
  for (i in (1:4)) {
    if (d[i] > 0.97) {
      d.adj[i] <- 0.97
    } else if (d[i] < -0.97) {
      d.adj[i] <- -0.97
    }
  }
  A.hat <- svd.A$u %*% diag(d.adj) %*% t(svd.A$v)
  D <- solve(diag(4) - A.hat)
  reg.mat <- rbind(reg1, reg2, reg3, reg4)
  for (j in (1:4)) {
    V.star[, j] <- V.hat[2:T, j] - A.hat[j, ] %*% reg.mat
  }
  Psi.hat <- compute.Psi.hat(V.star)
  Psi.hat <- D %*% Psi.hat %*% t(D)
  se <- as.numeric(sqrt(t(gradient) %*% Psi.hat %*% gradient / T))
  se
}
