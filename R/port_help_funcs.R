#' Calculate Stock Returns
#'
#' This function computes the stock returns of a given vector or time series of prices.
#' Log returns are calculated as the natural logarithm of the ratio of consecutive prices.
#'
#' @param prices a time series object representing the prices.
#' @param type a character, identifying the type of returns. Possible values are "log" for log returns
#' and "simple" for simple returns. Default value is set to "log".
#' @return a time series object with log returns.
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' head(rets_m[, 1:5])
#' rets_spx_m <- calc_rets(prices_spx_m)
#' head(rets_spx_m)
#' @export
calc_rets <- function(prices, type = "log") {
  prices <- as.matrix(prices)
  if (type == "log") {
    rets <- diff(log(prices), lag = 1)
  } else if (type == "simple") {
    rets <- diff(prices, lag = 1) / prices[-nrow(prices), ]
  } else {
    print("Invalid type. Please choose 'log' or 'simple'.")
  }
  return(rets)
}

#' Find Repeated Values in a Data Frame
#'
#' This function identifies repeated values in a data frame and returns a summary
#' of the repetitions, up to a specified maximum number of repetitions.
#'
#' @param df a data frame in which to search for repeated values.
#' @param max_rep an integer, specifying the maximum number of repetitions to consider.
#'   Defaults to 10.
#'
#' @return a summary of the positions of repeated values in the data frame, up to the specified maximum.
#'
#' @examples
#' df <- data.frame(a = c(1, 2, 2, 3, 3, 3), b = c(4, 5, 5, 6, 6, 6))
#' find_rep_vals(df, max_rep = 2)
#'
#' @export
find_rep_vals <- function(df, max_rep = 10) {
  sapply(df, function(col) {
    if (!is.numeric(col)) {
      return(NA)
    }
    rle_out <- rle(col)
    max(rle_out$lengths, na.rm = TRUE) >= max_rep
  }) |> which()
}
