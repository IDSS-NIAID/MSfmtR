#' process_raw
#' Process raw data into MSstats format
#'
#' @param config list of configuration parameters
#' @param stage character path to checkpoint file
#' @param save_intermediate logical save intermediate data
#' @param preprocess logical preprocess data using MSstats
#'
#' @details This function processes raw data into MSstats format. If save_intermediate is TRUE, the processed data are also saved to the checkpoint file. This also has the side effect of changing/updating `config` in the calling environment.
#'
#' There are currently two tested functions for `convert_function`: `MSstatsConvert::SpectronauttoMSstatsFormat` and `MSfmtR::convert_spectronaut`. The `MSfmtR` version uses the abundance score from Spectronaut directly and does not perform any preprocessing.
#'
#' @return data frame of processed raw data
#' @export
#' @importFrom Biostrings readAAStringSet
#' @importFrom dplyr filter select distinct pull n
#' @importFrom MSstats dataProcess
#' @importFrom MSstatsConvert SpectronauttoMSstatsFormat
#' @importFrom purrr map_chr
#' @importFrom readr read_delim
#' @importFrom stringr fixed str_split
#' @importFrom utils combn
process_raw <- function(config, stage = file.path(config$output_dir, config$processed_checkpoint),
                        save_intermediate = TRUE, preprocess = TRUE)
{
  # for those pesky no visible binding warnings
  if(FALSE)
    ABUNDANCE <- censored <- EG.PrecursorId <- FEATURE <- FRACTION <- FrgIon.uid <- F.Charge <-
      F.FrgIon <- F.FrgLossType <- F.NormalizedPeakArea <- F.NormalizedPeakHeight <- GROUP <-
      INTENSITY <- originalRUN <- LABEL <- LogIntensities <- MissingPercentage <- more50missing <-
      newABUNDANCE <- NumImputedFeature <- NumMeasuredFeature <- PEPTIDE <- PG.ProteinAccessions <-
      PG.Quantity <- PROTEIN <- RUN <- R.Condition <- R.FileName <- R.Replicate <- SUBJECT <-
      TotalGroupMeasurements <- TRANSITION <- NULL

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
  if(preprocess)
  {
    # convert raw data into MSstats format
    data <- SpectronauttoMSstatsFormat(raw,
                                       intensity = 'NormalizedPeakArea',
                                       use_log_file = FALSE)

    # remove out-of-spec peptides (MSstats assumes one row for every run for each feature, even if the intensity value is missing
    #                              therefore, use mutate and convert to NA instead of filter...actually mutate doesn't work)
    data$Intensity[data$Intensity < config$lloq |
                   data$Intensity > config$uloq] <- NA

    # peptide-level data (see documentation of `MSstats::dataProcess` at https://www.bioconductor.org/packages/devel/bioc/vignettes/MSstats/inst/doc/MSstats.html)
    data <- dataProcess(data,
                        normalization = FALSE,
                        use_log_file = FALSE)
  }else{
    # peptide-level data (see documentation of `MSstats::dataProcess` at https://www.bioconductor.org/packages/devel/bioc/vignettes/MSstats/inst/doc/MSstats.html)
    raw <- raw |>

      # Figure out unique ID for each FrgIon (not unique)
      group_by(PG.ProteinAccessions, EG.PrecursorId, F.FrgIon) |>
      mutate(FrgIon.uid = paste(EG.PrecursorId, F.FrgIon, F.Charge, F.FrgLossType, sep = "_") |>
               factor() |>
               as.numeric()) |>
      ungroup() |>

      # format information for FeatureLevelData
      mutate(PROTEIN = factor(PG.ProteinAccessions),
             PEPTIDE = str_replace(EG.PrecursorId, "_", "") |>
               str_replace("_.", "_") |>
               factor(),
             TRANSITION = paste(F.FrgIon, FrgIon.uid, sep = '_') |>
               factor(),
             FEATURE = paste(PEPTIDE, TRANSITION, sep = '_') |>
               factor(),
             LABEL = factor('L'),                                                               ##### don't know what this is
             GROUP = factor(R.Condition),
             SUBJECT = factor(R.Replicate),
             FRACTION = as.integer(1),                                                          ##### don't know what this is
             originalRUN = factor(R.FileName, levels = unique(R.FileName[order(R.Condition)])), # this is the ordering MSstats uses
             RUN = as.numeric(originalRUN) |> as.factor(),
             censored = FALSE,
             INTENSITY = ifelse(rep(config$normMeasure == 'NormalizedPeakArea', nrow(raw)),
                                F.NormalizedPeakArea,
                                F.NormalizedPeakHeight),
             ABUNDANCE = log10(INTENSITY),
             newABUNDANCE = ABUNDANCE,
             predicted = as.numeric(NA))

    # create FeatureLevelData and filter
    FeatureLevelData <- raw |>

      dplyr::select(PROTEIN, PEPTIDE, TRANSITION, FEATURE, LABEL, GROUP, RUN, SUBJECT,
                    FRACTION, originalRUN, censored, INTENSITY, ABUNDANCE, newABUNDANCE, PG.Quantity) |> # only using PG.Quantity for protein-level data. Drop it later.

      dplyr::filter(INTENSITY > config$lloq, # remove out-of-spec peptides
                    INTENSITY < config$uloq)


    # protein-level data
    ProteinLevelData = FeatureLevelData |>

      group_by(RUN, PROTEIN) |>
      mutate(NumMeasuredFeature = n(),
             MissingPercentage = sum(INTENSITY < 1) / n(),
             more50missing = MissingPercentage > 0.5,
             NumImputedFeature = sum(is.na(INTENSITY))) |>
      ungroup() |>

      group_by(GROUP, PROTEIN) |>
      mutate(TotalGroupMeasurements = n()) |>
      ungroup() |>

      dplyr::select(RUN, PROTEIN, originalRUN, GROUP, SUBJECT, PG.Quantity,                          # only keep protein-specific info
                    TotalGroupMeasurements, NumMeasuredFeature, MissingPercentage, more50missing) |> # plus some summary stats
      unique() |>

      mutate(LogIntensities = log10(PG.Quantity),
             NumImputedFeature = 0) |>

      dplyr::select(RUN, PROTEIN, LogIntensities, originalRUN, GROUP, SUBJECT, TotalGroupMeasurements,
                    NumMeasuredFeature, MissingPercentage, more50missing, NumImputedFeature) |>

      dplyr::rename(Protein = PROTEIN)


    # final data object
    data <- list(FeatureLevelData = dplyr::select(FeatureLevelData, -PG.Quantity),
                 ProteinLevelData = ProteinLevelData,
                 SummaryMethod = 'linear')
  }

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
