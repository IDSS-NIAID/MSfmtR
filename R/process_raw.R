#' process_raw
#' Process raw data into MSstats format
#'
#' @param config list of configuration parameters
#' @param stage character path to checkpoint file
#' @param save_intermediate logical save intermediate data
#'
#' @details This function processes raw data into MSstats format. If save_intermediate is TRUE, the processed data are also saved to the checkpoint file. This also has the side effect of changing/updating `config` in the calling environment.
#' @return data frame of processed raw data
#' @export
#' @importFrom Biostrings readAAStringSet
#' @importFrom dplyr filter select distinct pull
#' @importFrom MSstats dataProcess
#' @importFrom MSstatsConvert SpectronauttoMSstatsFormat
#' @importFrom purrr map_chr
#' @importFrom readr read_delim
#' @importFrom stringr fixed str_split
#' @importFrom utils combn
process_raw <- function(config, stage = file.path(config$output_dir, config$processed_checkpoint),
                        save_intermediate = TRUE)
{
  # for those pesky no visible binding warnings
  if(FALSE)
    PG.ProteinAccessions <- R.Condition <- NULL

  # contaminants
  contam <- readAAStringSet(file.path(config$fasta_dir, config$cont_fasta))@ranges@NAMES |>
    str_split(fixed('|')) |>
    map_chr(~ .x[2])

  # Load and process the data
  raw <- with(config, file.path(input_dir, in_file)) |>
    read_delim(delim = "\t", col_names = TRUE) |>
    dplyr::filter(!PG.ProteinAccessions %in% contam &
                    !grepl('Cont_', PG.ProteinAccessions, fixed = TRUE))

  ########## do protein group analysis here ##########

  # format ratios
  if(is.null(config$ratios))
  {
    config$ratios <- raw |>
      dplyr::select(R.Condition) |>
      distinct() |>
      pull() |>
      combn(2, simplify = TRUE) |>
      t()
  }else{
    tmp <- config$ratios |>
      str_split('/')

    config$ratios <- cbind(map_chr(tmp, ~ .x[1]), map_chr(tmp, ~ .x[2]))

    # make sure these are all in raw$R.Condition
    if(any(!config$ratios %in% raw$R.Condition))
      stop("Not all conditions in ratios are in the raw data")
  }


  # process raw data
  data <- raw |>

    # convert raw data into MSstats format
    SpectronauttoMSstatsFormat(intensity = 'NormalizedPeakArea',
                               use_log_file = FALSE)

  # remove out-of-spec peptides (MSstats assumes one row for every run for each feature, even if the intensity value is missing
  #                              therefore, use mutate and convert to NA instead of filter...actually mutate doesn't work)
  data$Intensity[data$Intensity < config$lloq |
                 data$Intensity > config$uloq] <- NA

  # peptide-level data (see documentation of `MSstats::dataProcess` at https://www.bioconductor.org/packages/devel/bioc/vignettes/MSstats/inst/doc/MSstats.html)
  data <- dataProcess(data,
                      normalization = FALSE,
                      use_log_file = FALSE)

  # save data
  if(FALSE)
    save(raw, file = file.path(config$output_dir, 'raw.RData'))

  # checkpoint
  if(save_intermediate)
  {
    config_bak <- config
    save(data, config_bak, file = stage)
  }

  config <<- config
  return(data)
}
