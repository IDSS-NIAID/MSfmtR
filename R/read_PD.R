#' read_PD
#'
#' Read excel files exported from Proteome Discoverer
#'
#' @param config List containing configuration parameters
#' @param file Character, path to input file. If NULL (default), parameters in `config` are used.
#' @param ... additional agruments to pass to `updt_config`
#'
#' @details The input file should be an excel spreadsheet <...describe formatting...>
#' Protein and peptide rows are written to a temporary directory and read in separately to maintain correct column types.
#'
#'  Configuration parameters relevant to this function include: `input_dir`, `in_file`, `in_sheet`.
#'  See `?updt_config` for more configuration details.
#'
#' @return A list with two data.frames, `proteins` and `peptides`
#' @export
#' @importFrom dplyr filter select
#' @importFrom readr read_csv write_csv
#' @importFrom readxl read_excel
#' @importFrom tidyselect starts_with
read_PD <- function(config, file = NULL, ...)
{
  # take care of annoying no visible binding notes
  if(FALSE)
    Master <- NULL
  
  # update config
  config <- updt_config(config, ...)

  # validate parameters
  if(is.null(file))
    file <- file.path(config$input_dir, config$in_file)

  # temporary directory to work from
  tmp_dir <- tempdir()

  # read in and separate full sheet
  dat_full <- read_excel(file, na = c('', 'n/a'), sheet = config$in_sheet)

  # read in proteins (write only protein data and read again to get auto-typing of columns correct)
  dat_full |>
    filter(!is.na(Master)) |>
    select(-starts_with('...')) |>
    write_csv(file = file.path(tmp_dir, 'proteins.csv'))

  proteins <- suppressMessages(read_csv(file.path(tmp_dir, 'proteins.csv')))

  # filter out Master protein column
  peptides <- dat_full |>
    filter(is.na(Master)) |>
    select(-Master)
  
  # remove any extra columns from proteins (i.e. if there are more columns of proteins)
  peptides <- peptides[,!(is.na(peptides[1,]) |> unlist() |> as.vector())]

  # read in peptides (write only peptide data and read again to get auto-typing of columns correct)
  peptides[c(TRUE, peptides[[1]][-1] != peptides[[1]][1]),] |>              # remove column headers (except first)
    write_csv(file = file.path(tmp_dir, 'peptides.csv'), col_names = FALSE) # don't want protein headers

  peptides <- suppressMessages(read_csv(file.path(tmp_dir, 'peptides.csv')))

  return(list(proteins = proteins,
              peptides = peptides))
}
