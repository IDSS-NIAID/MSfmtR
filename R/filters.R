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
#' @importFrom stats na.omit
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


#' filter_completeness
#'
#' Filter specified columns based on completeness for each condition.
#'
#' @param dat data.frame to filter
#' @param cols Integer, columns to filter
#' @param filter_cols Integer, columns of `dat` with missingness information (see Details).
#' @param conditions Character vector specifying conditions. By default we assume names of `filter_cols` are the conditions prepended with "drop_".
#' @param ratio Logical value indicating whether the columns to be filtered are ratio statistics, as ratio statistics are calculated from two conditions and non-ratio statistics only pertain to one.
#' @param max_ratio Numeric vector indicating the maximum ratio to return for each column. Anything above this (or below `1/max_ratio`) will be truncated at the limit. If a single value is provided, it is used for all columns. If a vector is provided, it must be the same length as `cols`. (ignored when `ratio` is FALSE)
#' @param name_with Function for naming new, filtered columns.
#'
#' @details
#' This function filters the specified columns based on the completeness of the data, replacing filtered values with `NA`.
#' The `filter_cols` argument specifies the columns of `dat` that contain missingness information (e.g. from `check_missingness()`).
#' These columns should be logical vectors with `TRUE` indicating too much missingness and that row should be filtered for the specified condition.
#' Also, `filter_cols` should correspond directly to the values in `conditions`, with one condition per column.
#'
#' @return data frame with filtered data.
#' @export
#' @importFrom dplyr filter mutate select starts_with
#' @importFrom stringr str_replace str_split
#' @importFrom tidyr pivot_longer pivot_wider
filter_completeness <- function(dat, cols, filter_cols,
                                conditions = gsub('drop_', '', names(dat)[filter_cols]),
                                ratio = FALSE, max_ratio = Inf,
                                name_with = function(x) str_replace(x, ":", " (Filtered):"))
{
  # check length of cols
  if(length(cols) == 0)
  {
    warning("Skipping filtering, no columns specified")
    return(dat)
  }

  # check length of filter_cols
  if(length(max_ratio) == 1)
    max_ratio <- rep(max_ratio, length(cols))

  if(length(max_ratio) != length(cols))
  {
    stop("Length of max_ratio must be 1 or the same length as cols")
  }

  # filter cols
  for(i in 1:length(cols))
  {
    # figure out which conditions we are filtering for
    j <- map_lgl(conditions, ~ grepl(.x, names(dat)[cols[i]], fixed = TRUE)) |>
      which()

    if(ratio)
    {
      # check that we have exactly two matches
      if(length(j) != 2)
      {
        warning("Incorrect number of matches for ", names(dat)[cols[i]], " in conditions - expecting 2 - skipping")
        next
      }

      # perform filtering
      # for ratios we need to check both conditions
      # if either is good, we keep the value
      # if only one is good, this usually means the observed ratio is either 0 or Inf and we will truncate it later
      dat$new <- ifelse(dat[[ filter_cols[j[1]] ]] & dat[[ filter_cols[j[2]] ]],
                        NA, dat[[ cols[i] ]])

      # perform ratio truncation
      # if the ratio is too large, set it to the max_ratio
      dat$new <- ifelse(dat$new > max_ratio[i], max_ratio[i], dat$new)

      # if the ratio is too small, set it to 1/max_ratio
      dat$new <- ifelse(dat$new < 1/max_ratio[i], 1/max_ratio[i], dat$new)

    }else{
      # check that we have exactly one match
      if(length(j) == 0)
      {
        warning("No match for ", names(dat)[cols[i]], " in conditions - skipping")
        next
      }

      if(length(j) > 1)
      {
        warning("Multiple matches for ", names(dat)[cols[i]], " in conditions - skipping")
        next
      }

      # perform filtering
      dat$new <- ifelse(dat[[ filter_cols[j] ]], NA, dat[[ cols[i] ]])
    }

    names(dat)[which(names(dat) == 'new')] <- name_with(names(dat)[cols[i]])
  }

  return(dat)
}
