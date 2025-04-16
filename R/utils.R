
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

#' load_MSfmtR
#'
#' Load a specific version of MSfmtR
#'
#' @param lib.loc Character, specifying the location of a local library to use for installation (i.e. to avoid overwriting the system-wide installation)
#' @param commit Character, specifying the GitHub commit to use
#'
#' @return Logical value indicating the correct version of the package was successfully loaded
#' @export
#' @importFrom devtools install_github
load_MSfmtR <- function(lib.loc, commit)
{
  # this is where we expect the installed DESCRIPTION file to be located
  description <- file.path(lib.loc, 'MSfmtR', 'DESCRIPTION')

  # check if we have the correct version installed
  if(!file.exists(description)) # if there is no DESCRIPTION file, it isn't installed where we want it
  {
    update_install <- TRUE
  }else if(read.dcf(description)[,'GithubSHA1'] != commit) # if the commit is different, install
  {
    update_install <- TRUE
  }else{
    update_install <- FALSE
  }

  if(update_install)
  {
    devtools::install_github('IDSS-NIAID/MSfmtR', ref = commit, lib = lib.loc)

    # double check that we didn't install over the top of an already-loaded version
    if('package:MSfmtR' %in% search())
      warning('MSfmtR is already loaded. Restarting R is advised after a fresh install')
  }

  require(MSfmtR, lib.loc = lib.loc)
}
