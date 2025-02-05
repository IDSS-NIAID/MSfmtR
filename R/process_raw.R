#' process_raw
#' Process raw data into MSstats format
#'
#' @param config list of configuration parameters
#' @param stage character path to checkpoint file
#' @param save_intermediate logical save intermediate data
#' @param preprocess logical preprocess data using MSstats
#' @param format character format of data (see `?raw_to_fld` for options)
#' @param normMeasure character normalization measure (ignored when `preprocess` is TRUE or when `peptide_summary != 'none'`)
#' @param peptide_summary character flagging level of peptide summarization (ignored when `preprocess` is TRUE). Options: 'PEP' for peptide level (default), 'FG' for fragment group level, or 'none' to leave at the fragment level.
#' @param input_dir character path to directory containing raw data
#' @param in_file character name of raw data file
#' @param output_dir character path to directory for output
#' @param processed_checkpoint character name of checkpoint file for processed data
#' @param fasta_dir character path to directory containing fasta files
#' @param cont_fasta character path to fasta file containing contaminants
#' @param ratios character vector of the contrasts to be used in the MSstats analysis.
#'  Default is all combinations of all levels in of `R.Condition` in `in_file`.
#'  Each ratio in the list should be of the form "<group1>/<group2>", so for group1=Case and group2=Control, the ratio would be "Case/Control".
#' @param groups Need to review how this is used
#' @param lloq numeric lower limit of quantification
#' @param uloq numeric upper limit of quantification
#'
#'
#' @details This function processes raw data into an MSstats formatted object in R.
#' If save_intermediate is TRUE, the processed data are also saved to the checkpoint file.
#' This also has the side effect of changing/updating `config` in the calling environment.
#'
#' Any parameters that are `NULL` will be taken from the package defaults.
#'
#'
#' @return data frame of processed raw data
#' @export
#' @importFrom Biostrings readAAStringSet
#' @importFrom dplyr filter select distinct pull n rename group_by mutate ungroup arrange
#' @importFrom MSstats dataProcess
#' @importFrom MSstatsConvert SpectronauttoMSstatsFormat
#' @importFrom purrr map_chr
#' @importFrom readr read_delim
#' @importFrom stringr fixed str_split
#' @importFrom utils combn
process_raw <- function(config = configure_formatR(),
                        stage = file.path(config$output_dir, config$processed_checkpoint),
                        save_intermediate = TRUE,
                        preprocess = NULL, format = NULL, normMeasure = NULL, peptide_summary = NULL,
                        input_dir = NULL, in_file = NULL, output_dir = NULL, processed_checkpoint = NULL,
                        fasta_dir = NULL, cont_fasta = NULL,
                        ratios = NULL, groups = NULL,
                        lloq = NULL, uloq = NULL)
{
  # for those pesky no visible binding warnings
  if(FALSE)
    ABUNDANCE <- censored <- EG.PrecursorId <- FEATURE <- FRACTION <- FrgIon.uid <- F.Charge <-
      F.FrgIon <- F.FrgLossType <- F.NormalizedPeakArea <- F.NormalizedPeakHeight <- GROUP <-
      INTENSITY <- originalRUN <- LABEL <- LogIntensities <- MissingPercentage <- more50missing <-
      newABUNDANCE <- NumImputedFeature <- NumMeasuredFeature <- PEPTIDE <- PG.ProteinAccessions <-
      PG.Quantity <- PROTEIN <- RUN <- R.Condition <- R.FileName <- R.Replicate <- SUBJECT <-
      TotalGroupMeasurements <- TRANSITION <- NULL

  # update config with parameters
  config <- updt_config(config,
                        preprocess = preprocess, format = format,
                        normMeasure = normMeasure, peptide_summary = peptide_summary,
                        input_dir = input_dir, in_file = in_file,
                        output_dir = output_dir, processed_checkpoint = processed_checkpoint,
                        fasta_dir = fasta_dir, cont_fasta = cont_fasta,
                        ratios = ratios, groups = groups,
                        lloq = lloq, uloq = uloq)


  # figure out where contaminants file is located - if it is already a valid path we can move on
  if(!file.exists(config$cont_fasta))
  {
    # if it isn't a valid path, look for it in the fasta_dir
    if(file.exists(file.path(config$fasta_dir, config$cont_fasta)))
    {
      config$cont_fasta <- file.path(config$fasta_dir, config$cont_fasta)
    }else{
      # if it isn't in the fasta_dir, give up
      paste("Can't find contaminants fasta file:", config$cont_fasta) |>
        stop()
    }
  }

  # contaminants
  contam <- readAAStringSet(config$cont_fasta)@ranges@NAMES |>
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
    if(!is.null(config$groups))
    {
      config$ratios <- config$groups |>
        combn(2, simplify = TRUE) |>
        t()
    }else{
      config$ratios <- raw |>
        dplyr::select(R.Condition) |>
        distinct() |>
        pull() |>
        combn(2, simplify = TRUE) |>
        t()
    }
  }else{
    tmp <- config$ratios |>
      str_split('/')

    config$ratios <- cbind(map_chr(tmp, ~ .x[1]), map_chr(tmp, ~ .x[2]))

    # make sure these are all in raw$R.Condition
    if(any(!config$ratios %in% raw$R.Condition))
      stop("Not all conditions in ratios are in the raw data")
  }


  # process raw data
  if(config$preprocess & config$format == 'MSstats')
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

    # convert raw data to FeatureLevelData (this keeps PG.Quantity for use in formatting protein-level data)
    FeatureLevelData <- raw_to_fld(raw, config$format, config)

    # final data object
    data <- list(FeatureLevelData = dplyr::select(FeatureLevelData, -PG.Quantity),
                 ProteinLevelData = fld_to_pld(FeatureLevelData),
                 SummaryMethod = 'linear')

  }

  # save data
  if(FALSE)
    save(raw, file = file.path(config$output_dir, 'raw.RData'))

  # checkpoint
  if(save_intermediate)
  {
    config_bak <- config
    save(data, config_bak, file = file.path(config$output_dir, config$processed_checkpoint))
  }

  config <<- config
  return(data)
}


#' raw_to_fld
#' Convert raw data to FeatureLevelData
#'
#' @param raw data.frame
#' @param format character describing the format of the raw data
#' @param config list of configuration parameters
#'
#' @details The only config values used are `lloq` and `uloq`. These are used to filter out-of-spec data.
#'
#' Currently supported formats are
#' 'MSstats', the MSstats report exported by Spectronaut, and
#' 'other', the ... report exported by Spectronaut.
#'
#' @references See documentation of `MSstats::dataProcess` at https://www.bioconductor.org/packages/devel/bioc/vignettes/MSstats/inst/doc/MSstats.html)
#'
#' @return FeatureLevelData data.frame
#' @export
raw_to_fld <- function(raw, format = 'MSstats',
                       config = list(lloq = 0, uloq = Inf))
{
  if(format == 'MSstats')
  {
    # Figure out unique ID for each FrgIon (not unique)
    raw <- raw |>
      group_by(PG.ProteinAccessions, EG.PrecursorId, F.FrgIon) |>

      mutate(FrgIon.uid = paste(EG.PrecursorId, F.FrgIon, F.Charge, F.FrgLossType, sep = "_") |>
               factor() |>
               as.numeric()) |>

      ungroup() |>

      # format information for FeatureLevelData at the fragment level
      mutate(TRANSITION = paste(F.FrgIon, FrgIon.uid, sep = "_") |>
               factor(),
             FEATURE = paste(EG.PrecursorId, TRANSITION, sep = "_") |>
               factor())


    # roll up peptides - we are assuming MSstats format from Spectronaut for now
    if(config$peptide_summary == 'none')
    {
      if(config$normMeasure == 'NormalizedPeakArea')
      {
        raw <- rename(raw, Intensity_measure = F.NormalizedPeakArea)
      }else{
        raw <- rename(raw, Intensity_measure = F.NormalizedPeakHeight)
      }

    }else{

      # first group appropriately according to config$peptide_summary
      if(config$peptide_summary == 'FG'){
        # roll up to Fragment Group level
        raw <- rename(raw, Intensity_measure = FG.Quantity) |>
          group_by(PG.ProteinAccessions, R.FileName, R.Replicate, PEP.GroupingKey, EG.ModifiedSequence)

      }else if(config$peptide_summary == 'PEP'){
        # roll up to Peptide level
        raw <- rename(raw, Intensity_measure = PEP.Quantity) |>
          group_by(PG.ProteinAccessions, R.FileName, R.Replicate, PEP.GroupingKey)
      }else{
        stop("Unknown value for `peptide_summary`: ", config$peptide_summary)
      }

      raw <- raw |>

        # concatenate TRANSITION and FEATURE by group
        mutate(TRANSITION = as.character(TRANSITION) |> unique() |> paste(collapse = ','),
               FEATURE    = as.character(FEATURE   ) |> unique() |> paste(collapse = ',')) |>
        ungroup() |>

        # remove extra columns and duplicate rows that we just summarized
        dplyr::select(R.FileName, R.Condition, R.Replicate,
                      PG.ProteinAccessions, PG.Quantity,
                      EG.ModifiedSequence,
                      TRANSITION, FEATURE,
                      Intensity_measure) |>
        unique()
    }


    raw <- raw |>

      # format information for FeatureLevelData
      mutate(PROTEIN = factor(PG.ProteinAccessions),
             PEPTIDE = str_replace_all(EG.ModifiedSequence, "_", "") |> factor(),
             # TRANSITION and FEATURE are defined above
             LABEL = factor('L'),                                                               ##### for label-free data, this is always 'L' - need to update this when we have a good example TMT data set
             GROUP = factor(R.Condition),
             SUBJECT = factor(R.Replicate),
             FRACTION = as.integer(1),                                                          ##### For fractionation - need to update this when we get a good example data set
             originalRUN = factor(R.FileName, levels = unique(R.FileName[order(R.Condition)])), # this is the ordering MSstats uses
             RUN = as.numeric(originalRUN) |> as.factor(),
             censored = FALSE,
             INTENSITY = Intensity_measure,
             ABUNDANCE = log10(INTENSITY),
             newABUNDANCE = ABUNDANCE,
             predicted = as.numeric(NA))
  }else{
    stop("Unknown format")
  }

  # create FeatureLevelData and filter
  FeatureLevelData <- raw |>

    dplyr::select(PROTEIN, PEPTIDE, TRANSITION, FEATURE, LABEL, GROUP, RUN, SUBJECT,
                  FRACTION, originalRUN, censored, INTENSITY, ABUNDANCE, newABUNDANCE, PG.Quantity) |> # only using PG.Quantity for protein-level data. Drop it later.

    dplyr::filter(INTENSITY > config$lloq, # remove out-of-spec peptides
                  INTENSITY < config$uloq)

  if(dim(FeatureLevelData)[1] == 0)
    stop("No features left after filtering")

  return(FeatureLevelData)
}


#' fld_to_pld
#' Convert FeatureLevelData to ProteinLevelData
#'
#' @param fld FeatureLevelData
#'
#' @return ProteinLevelData data.frame
#' @export
#' @importFrom dplyr group_by mutate ungroup select rename
fld_to_pld <- function(fld)
{
  fld |>

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
}

