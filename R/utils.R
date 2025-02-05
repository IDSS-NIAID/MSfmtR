
#' tukey_biweight_mean
#' Tukey Biweight Mean
#'
#' @param x numeric vector
#' @param c tuning constant
#'
#' @source Original function definition provided by Gemini 2.0 Flash
#'
#' @return numeric
#' @export
#' @importFrom stats mad
tukey_biweight_mean <- function(x, c = 4.685) { # c is the tuning constant
  if(length(x) == 1)
    return(x)

  x_med <- median(x)
  u <- (x - x_med) / (c * mad(x))  # mad is Median Absolute Deviation
  w <- ifelse(abs(u) < 1, (1 - u^2)^2, 0)
  sum(w * x) / sum(w)
}
