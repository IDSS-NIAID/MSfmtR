# filters.R
# functions for filtering data


#' check_missingness
#'
#' Check for missingness in data for a group of columns and flag rows with too many missing or low values.
#'
#' @rdname check_missingness
#'
#' @param dat data.frame or matrix to check
#' @param cols Vector of columns to check for missingness
#' @param prop_good Numeric proportion of good values to keep.
#' @param llod Numeric lower limit of detection. Any value below this will be considered missing.
#' @param valid_chr Character vector of valid values for character columns. If `NULL` (default), all non-missing values are considered valid.
#' @param ... Additional arguments passed to `check_missingness_num()` or `check_missingness_chr()`.
#'
#' @details
#' By default, this function ignores the `llod` argument (since `llod = 0`) and requires all values to be non-missing.
#' If the lower limit of detection for your assay is 1e6, for example, you can set `llod = 1e6` to treat any value below this as missing.
#' To allow for missing values, set `prop_good` to something less than 1.
#' For example, 0.5 indicates at least 50 percent of values on a row must be non-missing and above the lower limit of detection.
#'
#' @return vector of logicals indicating which rows to filter / drop. In other words, rows with `TRUE` should be dropped.
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
#' check_missingness(dat, 2:5, prop_good = 0.5, llod = 1e6)
#'
#' # the second and third rows have less than 60% good values, so they are flagged to drop
#' check_missingness(dat, 2:5, prop_good = 0.6, llod = 1e6)
#'
#' @export
check_missingness <- function(dat, cols, ...)
{
  # check if the columns are numeric or character
  is_num <- sapply(dat[,cols], is.numeric)
  is_chr <- sapply(dat[,cols], is.character)

  # make sure we have all desired columns covered - if not, throw an error
  if(any(!is_num & !is_chr)){
    stop("Columns identified by `cols` must be either numeric or character")
  }

  # check for missingness in numeric columns
  if (any(is_num)) {
    retval_num <- check_missingness_num(dat, cols[is_num], ...)
  }else{
    retval_num <- rep(FALSE, nrow(dat))
  }

  # check for missingness in character columns
  if (any(is_chr)) {
    retval_chr <- check_missingness_chr(dat, cols[is_chr], ...)
  }else{
    retval_chr <- rep(FALSE, nrow(dat))
  }

  return(retval_num | retval_chr)
}

#' check_missingness_chr
#' @rdname check_missingness
#' @export
check_missingness_chr <- function(dat, cols, prop_good = 1, valid_chr = NULL, ...)
{
  if(is.null(valid_chr)){
    # if no valid values are provided, use all non-missing values
    valid_chr <- unique(na.omit(unlist(dat[,cols])))
  }

  # count number of rows passing the check
  ngood <- apply(dat[,cols], 1, function(.x) sum(.x %in% valid_chr))

  # flag rows to filter / drop
  retval <- ngood < prop_good * length(cols)

  return(retval)
}

#' check_missingness_num
#' @rdname check_missingness
#' @export
check_missingness_num <- function(dat, cols, prop_good = 1, llod = 0, ...)
{
  # count number of rows passing the check
  ngood <- (dat[,cols] > llod) |>
    rowSums(na.rm = TRUE)

  # flag rows to filter / drop
  retval <- ngood < prop_good * length(cols)

  return(retval)
}
