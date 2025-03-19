# filters.R
# functions for filtering data


#' check_missingness
#'
#' Check for missingness in data for a group of columns and flag rows with too many missing or low values.
#'
#' @param dat data.frame or matrix to check
#' @param cols Vector of columns to check for missingness
#' @param prop_good Numeric proportion of good values to keep.
#' @param llod Numeric lower limit of detection. Any value below this will be considered missing.
#'
#' @details
#' By default, this function ignores the `llod` argument (since `llod = Inf`) and requires all values to be non-missing.
#' If the lower limit of detection for your assay is 1e6, for example, you can set `llod = 1e6` to treat any value below this as missing.
#' To allow for missing values, set `prop_good` to something less than 1.
#' For example, 0.5 indicates at least 50 percent of values on a row must be non-missing and above the lower limit of detection.
#'
#' @return vector of logicals indicating which rows to filter / drop
#'
#' @examples
#' set.seed(294873)
#' dat <- data.frame(pep = c(1, 2, 3),
#'                   A = runif(3) * 1e7,
#'
#'                   # middle value is missing and last value is too small
#'                   B = runif(3) * c(1e7, NA, 1e4),
#'
#'                   # middle value is missing
#'                   C = runif(3) * c(1e7, NA, 1e7),
#'
#'                   # last value is too small
#'                   D = runif(3) * c(1e7, 1e7, 1e4))
#'
#' # all rows have at least two good values, so no rows are flagged to drop
#' check_missingness(dat, 2:5, prop_good = 0.5)
#'
#' # the second and third rows have less than 60% good values, so they are flagged to drop
#' check_missingness(dat, 2:5, prop_good = 0.6)
#'
#' @export
check_missingness <- function(dat, cols, prop_good = 1, llod = 0)
{
  # count number of rows passing the check
  ngood <- (dat[,cols] > llod) |>
    rowSums(na.rm = TRUE)

  # flag rows to filter / drop
  retval <- ngood < prop_good * length(cols)

  return(retval)
}
