
#' configuration options
#' @param config_file The name of the configuration file to use. This must be specified on
#'  the command line. Default is `config.yml`. Edit the `config.yml` file that comes with
#'  this repository to change the configuration options, or include parameters via the command line.
#'  Command line parameters will override the `config.yml` file. See `config.yml` for defaults.

### Directory parameters
#' @param input_dir The directory containing the data files (default is the current working directory).
#'  This script will only run one file at a time, so either the directory should contain only one file,
#'  or the `in_file` variable should be set to the file to be processed. If `in_file` is not provided it
#'  will look for a `.tsv` file to analyze inside of `input_dir`.
#' @param fasta_dir The path to the directory containing the fasta files to use for annotation. If not
#'  provided, the script will look for fasta files in the current working directory.
#' @param output_dir The directory to save the output files to (default is the current working directory).
#'  If it does not exist, it will be created.

### File parameters
#' @param in_file The name of the input file to be processed. This should be an MSstats-formatted file from Spectronaut.
#'  If not provided, the script will look for a `.tsv` file in `input_dir`, and if there is more than one it will throw
#'  an error.
#' @param out_xlsx The name of the output file. Default is `str_replace(config$in_file, fixed('.tsv'), '.xlsx')`.
#' @param out_sqlite The name of the output SQLite database. Default is `str_replace(config$out_file, fixed('.xlsx'), '.sqlite')`.
#' @param sheet The name of the sheet to write the output to. Default is `str_replace(config$out_file, fixed('.xlsx'), '')`.

### Filtering parameters
#' @param uloq Upper limit of quantification. Proteins with abundances greater than this will be removed from the analysis. Default is `Inf`.
#' @param lloq Lower limit of quantification. Proteins with abundances less than this will be removed from the analysis. Default is `0`.
#' @param cont_fasta The name of the contaminant fasta file to use (assumed to be in `fasta_dir`). Default
#'  is the Universal Contaminant file that comes with this package.

### Metadata parameters
#' @param fasta_meta The name of the fasta file(s) to use for annotation. If no annotation files are provided,
#'  `fasta_dir` will be checked for fasta files different from `cont_fasta`. If none are found, UniProt
#'   will be used for annotation.
#' @param taxId The taxonomy ID of the organism being analyzed. Default is `9606` (human). This is used to
#'  search for anything outside of the files in `fasta`.

### MSStats parameters
#' @param ratios The contrasts to be used in the MSStats analysis. Default is all combinations of all levels
#'  in `in_file`. Each ratio in the list should be of the form "<group1>/<group2>", so for group1=Case and
#'  group2=Control, the ratio would be "Case/Control". These labels should match values in the `R.Condition`
#'  column of `in_file`.
#'  @param normMeasure The normalization measure to use - pick from `NormalizedPeakArea` and `NormalizedPeakHeight`. Default is `NormalizedPeakArea`.

### Style parameters
#' @param protein_header_fill The fill color for the protein header in the output Excel file.
#' @param protein_rows_fill The fill color for the protein rows in the output Excel file.
#' @param peptide_header_fill The fill color for the peptide header in the output Excel file.
#' @param peptide_rows_fill The fill color for the peptide rows in the output Excel file.

### Checkpoint parameters
#' @param checkpoints A list of files to generate (if they already exist, this will be ignored). Options
#'  include: `xlsx,sql,processed,protein,peptide,wb`. Specifying `all` will result in all checkpoints being
#'  generated. Default is `xlsx,sql`.
#' @param processed_checkpoint The name of the checkpoint file to save the processed data to. Default is `paste0(config$sheet, '_processed.RData')`.
#' @param protein_checkpoint The name of the checkpoint file to save the protein statistics to. Default is `paste0(config$sheet, '_protein.RData')`.
#' @param peptide_checkpoint The name of the checkpoint file to save the peptide statistics to. Default is `paste0(config$sheet, '_peptide.RData')`.
#' @param wb_checkpoint The name of the checkpoint file to save the R object containing the formatted Excel workbook. Default is `paste0(config$sheet, '_wb.RData')`.

#' @details This script is designed to be run from the command line. It will read the `config.yml` file in the
#' current directory and use the parameters to process the data.
#'
#' @examples
#' Rscript driver.R --args input_dir=data output_dir=output


##################
# Load libraries #
##################

library(MSstats)
library(MSstatsConvert)

library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)

library(openxlsx)
library(RSQLite)

library(Biostrings)
library(muscle)
library(UniProt.ws)
library(Peptides)


###############
# Load config #
###############

# read command line arguments
args <- commandArgs(trailingOnly = TRUE) |>
  str_split('=') |>
  map_chr(~ .x[2]) |>
  as.list()

names(args) <- commandArgs(trailingOnly = TRUE) |>
  str_split('=') |>
  map_chr(~ .x[1])

if(is.null(args$config_file))
  args$config_file <- 'config.yml'

# read yaml file (requires config package)
config <- config::get(file = args$config_file)

# replace any defaults with command line arguments
if(length(args) > 0)
  for(i in 1:length(args))
  {
    config[names(args)[i]] <- args[i]
  }

## Directory parameters
if(is.null(config$input_dir))
  config$input_dir <- '.'

if(is.null(config$fasta_dir))
  config$fasta_dir <- '.'

if(is.null(config$output_dir))
  config$output_dir <- '.'


## File parameters
if(is.null(config$in_file))
{
  config$in_file <- list.files(config$input_dir, pattern = 'tsv')

  if(length(config$in_file) != 1)
    stop('Expecting to find exactly 1 .tsv file in input_dir, found ', length(config$in_file))
}

if(is.null(config$out_xlsx))
  config$out_xlsx <- str_replace(config$in_file, fixed('.tsv'), '.xlsx')

if(is.null(config$out_sqlite))
  config$out_sqlite <- str_replace(config$out_xlsx, fixed('.xlsx'), '.sqlite')

if(is.null(config$sheet))
  config$sheet <- str_replace(config$out_xlsx, fixed('.xlsx'), '')


## Filtering parameters
if(is.null(config$uloq))
  config$uloq <- Inf
if(is.null(config$lloq))
  config$lloq <- 1
if(is.null(config$cont_fasta))
{
  config$cont_fasta <- 'inst/extdata/Universal Contaminants.fasta'
}else if(!file.exists(file.path(config$fasta_dir, config$cont_fasta)))
{
  warning('cont_fasta does not exist, using default: ', config$cont_fasta)
  config$cont_fasta <- 'inst/extdata/Universal Contaminants.fasta'
}


## Metadata parameters
if(is.null(config$fasta_meta))
{
  config$fasta_meta <- list.files(config$fasta_dir, pattern = 'fasta') |>
    grep(pattern = config$cont_fasta, invert = TRUE, value = TRUE)
}
if(is.null(config$taxId))
  config$taxId <- 9606


## MSStats parameters
# if ratios is null, we'll fill it in after reading in the raw data

if(is.null(config$normMeasure))
  config$normMeasure <- 'NormalizedPeakArea'


## Style parameters
if(is.null(config$protein_header_fill))
  config$protein_header_fill <- "#A7CDF0"
if(is.null(config$protein_rows_fill))
  config$protein_rows_fill <- "#DDEBF7"
if(is.null(config$peptide_header_fill))
  config$peptide_header_fill <- "#F0CBA8"
if(is.null(config$peptide_rows_fill))
  config$peptide_rows_fill <- "#FCE4D6"


## Checkpoint parameters
if(!is.null(config$checkpoints))
{
  if(config$checkpoints == 'all')
  {
    config$checkpoints <- c('xlsx', 'sql', 'processed', 'protein', 'peptide', 'wb')
  }else{
    config$checkpoints <- str_split(config$checkpoints, ',') %>%
      unlist()
  }
}else{
  config$checkpoints <- c('xlsx', 'sql')
}

if(is.null(config$processed_checkpoint))
  config$processed_checkpoint <- paste0(config$sheet, '_processed.RData')

if(is.null(config$peptide_checkpoint))
  config$peptide_checkpoint <- paste0(config$sheet, '_peptide.RData')

if(is.null(config$protein_checkpoint))
  config$protein_checkpoint <- paste0(config$sheet, '_protein.RData')

if(is.null(config$wb_checkpoint))
  config$wb_checkpoint <- paste0(config$sheet, '_wb.RData')


##################
# check progress #
##################

# check that output_dir exists
if(!dir.exists(config$output_dir))
  dir.create(config$output_dir)


# start at the end and move backward (i.e. higher files on the list depend on files later in the list)
all_files <- c(xlsx = 'out_xlsx',
               wb = 'wb_checkpoint',
               sql = 'out_sqlite',
               protein = 'protein_checkpoint',
               peptide = 'peptide_checkpoint',
               processed = 'processed_checkpoint',
               input = 'in_file')
all_files <- factor(all_files, levels = all_files, ordered = TRUE)

# file dependencies
depends <- list(xlsx = c('wb_checkpoint'),
                wb = c('protein_checkpoint', 'peptide_checkpoint'),
                sql = c('protein_checkpoint', 'peptide_checkpoint'),
                protein = c('processed_checkpoint', 'peptide_checkpoint'),
                peptide = c('processed_checkpoint'),
                processed = c('in_file'),
                input = character(0))

# files to be generated (ordered by dependency)
files <- all_files[config$checkpoints]

# check what files exist, which need to be loaded, and which need to be run
progress <- as.character(all_files) |>

  map_df(~ ifelse(.x == 'in_file',
                  file.path(config$input_dir,  config[[.x]]),
                  file.path(config$output_dir, config[[.x]])) |>
           file.info()) |>

  mutate(conf_name = all_files,
         generate = conf_name %in% files,
         run = generate & is.na(size),
         load = FALSE)

for(i in 1:nrow(progress))
{
  # does anything depend on this?
  need <- which(map_lgl(depends, ~ progress$conf_name[i] %in% .x))

  if(any(progress[need,'run']))
  {
    # if so, can we load a checkpoint?
    if(!is.na(progress$size[i]))
    {
      # if so, we don't need to run this
      progress$load[i] <- TRUE
    }else{
      # if not, we need to run this
      progress$run[i] <- TRUE
    }
  }
}


###########
# do work #
###########

##### Load and process raw data #####

stage <- with(config, file.path(output_dir, processed_checkpoint))

if(progress[stage,'load'])
{
  load(stage)
  config$ratios <- config_bak$ratios
}else if(progress[stage,'run']){

  # contaminants
  contam <- Biostrings::readAAStringSet(file.path(config$fasta_dir, config$cont_fasta))@ranges@NAMES |>
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


  # save data
  data <- list(FeatureLevelData = dplyr::select(FeatureLevelData, -PG.Quantity),
               ProteinLevelData = ProteinLevelData,
               SummaryMethod = 'linear')

  if(FALSE)
    save(raw, file = file.path(config$output_dir, 'raw.RData'))
  # checkpoint
  if(progress[stage, 'generate'])
  {
    config_bak <- config
    save(data, config_bak, file = stage)
  }
  
  # clean up / memory management
  rm(raw, FeatureLevelData, ProteinLevelData)
}


##### Process peptide stats #####

stage <- with(config, file.path(output_dir, peptide_checkpoint))

if(progress[stage,'load'])
{
  load(stage)
}else if(progress[stage,'run']){

  # calculate peptide stats
  peptides_long <- data$FeatureLevelData %>%

    group_by(PROTEIN, FEATURE, GROUP) %>%

    summarize(INTENSITY = median(INTENSITY, na.rm = TRUE),
              cv = sd(ABUNDANCE, na.rm = TRUE) / mean(ABUNDANCE, na.rm = TRUE) * 100) %>%

    ungroup()


  peptides <- peptides_long |>
    dplyr::rename(id = GROUP) |>

    pivot_wider(names_from = id, values_from = c(INTENSITY, cv))


  # extract peptide data
  peptides <- data$FeatureLevelData |>

    mutate(id = paste0(GROUP, ', ', originalRUN, ' (', SUBJECT, ')'),
           
           PEPTIDE = map_chr(PEPTIDE, ~ strsplit(as.character(.x),
                                                  split = '_', fixed = TRUE)[[1]][1] |>
                               gsub(pattern = '\\[.*?\\]', replacement = '')),

           # split modifications out into a different column
           Modification = map_chr(FEATURE, ~
                                    {
                                      openMod <- str_locate_all(.x, fixed('[')) %>%
                                        unlist() %>%
                                        as.vector() %>%
                                        unique()
                                      closeMod <- str_locate_all(.x, fixed(']')) %>%
                                        unlist() %>%
                                        as.vector() %>%
                                        unique()

                                      if(length(openMod) > 0)
                                      {
                                        mod <- substr(.x, openMod[1], closeMod[1]) %>%
                                          paste(collapse = ';')
                                      }else{
                                        mod <- ''
                                      }

                                      return(mod)
                                    })) %>%


    dplyr::select(PROTEIN, PEPTIDE, Modification, TRANSITION, FEATURE, id, INTENSITY) %>%

    arrange(id) |> # sort by id - this keeps column names in the same order as in `proteins`

    unique() |> # remove duplicates - these do appear rarely in the data

    pivot_wider(names_from = id, values_from = INTENSITY) %>%

    # merge stats into peptides
    left_join(peptides, by = c('PROTEIN', 'FEATURE'))

  if(FALSE)
  {
    # some diagnostic plots for looking at peptides_long
    library(ggplot2)
    library(cowplot)
    theme_set(theme_cowplot())

    peptides_long %>%
      ggplot(aes(x = cv, y = INTENSITY)) +
      geom_point() +
      facet_wrap(~GROUP)
  }else{
    rm(peptides_long)
  }

  # checkpoint
  if(progress[stage, 'generate'])
    save(peptides, file = stage)
}


##### Process protein stats #####

stage <- with(config, file.path(output_dir, protein_checkpoint))

if(progress[stage,'load'])
{
  load(stage)
}else if(progress[stage,'run']){

  # contrasts
  contrasts <- matrix(0, nrow = nrow(config$ratios), ncol = length(levels(data$FeatureLevelData$GROUP)),
                      dimnames = list(paste(config$ratios[,1], '/', config$ratios[,2]),
                                      levels(data$FeatureLevelData$GROUP)))

  for(i in 1:nrow(config$ratios))
  {
    contrasts[i,c(config$ratios[i,1],
                  config$ratios[i,2])] <- c(1,-1)
  }

  ### Metadata to add ###

  # if we have fasta files provided, use them to get protein metadata
  if(length(config$fasta_meta) > 0)
  {
    # read fasta
    fasta <- map(config$fasta_meta, ~ Biostrings::readAAStringSet(file.path(config$fasta_dir, .x))) |>
      AAStringSetList() |>
      unlist()

    # get named meta data on the end of each fasta entry
    fasta_meta <- names(fasta) |>
      str_locate_all(fixed('=')) |>
      map2(names(fasta), ~
             {
               names_start <-   .x[  ,1] - 2
               names_end   <-   .x[  ,1] - 1

               dat_start   <-   .x[  ,1] + 1
               dat_end     <- c(.x[-1,1] - 4, nchar(.y))

               dat <- str_sub_all(.y, dat_start, dat_end)[[1]]
               names(dat) <- str_sub_all(.y, names_start, names_end)[[1]]

               dat
             })

    # compile metadata
    meta <- tibble(Protein = names(fasta) |>
                     str_split(fixed('|')) |>
                     map_chr(~ .x[2]),
                   Description = names(fasta) |>
                     str_split(fixed('=')) |>
                     map_chr(~ .x[1]) |>
                     str_split(fixed('|')) |>
                     map_chr(~ .x[3]),
                   Organism = map_chr(fasta_meta, ~ .x['OS']),
                   Sequence = as.character(fasta))

    # otherwise get UniProt metadata - manually search individual IDs here: https://www.uniprot.org/uniprotkb
  }else{
    # all protein IDs
    proteinIDs <- data$ProteinLevelData$Protein |>
      str_split(fixed(';')) |>
      unlist()

    meta <- UniProt.ws::select(UniProt.ws(taxId=config$taxId),
                               proteinIDs,
                               c('organism_name', 'protein_name', 'sequence')) %>%
      dplyr::select(Entry, Organism, Protein.names, Sequence) %>%
      dplyr::rename(Protein = Entry,
                    Description = Protein.names)
  }


  ### calculate protein stats ###
  proteins <- groupComparison(contrast.matrix = contrasts, 
                              data = data,
                              use_log_file = FALSE)$ComparisonResult |>
    
    dplyr::filter(!is.na(Protein)) |>

    dplyr::select(Protein, Label, log2FC, pvalue, adj.pvalue) |>

    pivot_wider(names_from = Label, values_from = c(log2FC, pvalue, adj.pvalue))

  ### Extract protein data ###
  tmp <- data$ProteinLevelData |>

    mutate(id = paste0(GROUP, ', ', originalRUN, ' (', SUBJECT, ')')) |>   # create a unique ID for each measurement

    # add summary (by median) of log intensities for each protein
    bind_rows({data$ProteinLevelData |>
        group_by(Protein, GROUP) |>
        summarize(LogIntensities = median(LogIntensities, na.rm = TRUE)) |>
        ungroup() |>
        dplyr::rename(id = GROUP)}) |>

    mutate(Intensity = 10^LogIntensities,                                  # report linear scale
           primary_id = str_split(Protein, fixed(';')) |>                  # pick a primary protein ID for now - still need to deal with protein groups
             map_chr(~ .x[1])) |>

    arrange(is.na(GROUP), id) |>                                           # sort by ID - this keeps column names in the same order as in `peptides`
    
    dplyr::select(-RUN, -LogIntensities, -originalRUN, -GROUP, -SUBJECT, -NumMeasuredFeature,
                  -MissingPercentage, -more50missing, -NumImputedFeature, -TotalGroupMeasurements) |>

    pivot_wider(names_from = id, values_from = Intensity) |>

    left_join(meta, by = c('primary_id' = 'Protein')) |>                   # merge metadata into proteins
    
    mutate(nAA = str_length(Sequence),                                     # calculate protein length
           `coverage%` = NA,                                               # placeholder for coverage
           `mass (kDa)` = NA)                                              # placeholder for mass

  # calculate masses
  skip <- grepl('X', tmp$Sequence) & !is.na(tmp$Sequence) # skip these that have unknown amino acids
  tmp$`mass (kDa)`[!skip] <- Peptides::mw(as.character(tmp$Sequence[!skip]))

  # look up any masses with unknown amino acids
  if(any(skip))
  {
    meta_mass <- UniProt.ws::select(UniProt.ws(taxId=config$taxId),
                       tmp$primary_id[skip],
                       c('mass'))

    tmp$`mass (kDa)`[skip] <- as.numeric(meta_mass$Mass) / 1000
  }


  # merge protein data
  proteins <- tmp |>
    dplyr::select(Protein, Description, Organism,
                names(tmp)[!names(tmp) %in% c('Protein', 'Description', 'Organism', 
                                              'primary_id')]) |>
    
    mutate(Description = str_split(Description, fixed('_')) |>
             map_chr(~ .x[[2]]) |>
             str_replace(' OS$', '')) |>

    left_join(proteins, by = c('Protein' = 'Protein'))


  # Calculate coverage
  for(i in 1:nrow(proteins))
  {
    if(is.na(proteins$Sequence[i]))
       next

    seqs <- c(Biostrings::AAStringSet(proteins$Sequence[i]),
              Biostrings::AAStringSet(filter(peptides, PROTEIN == as.character(proteins$Protein[i]))$PEPTIDE %>% unique()))

    # align sequences and convert to matrix
    aligned <- muscle::muscle(seqs) %>%
      as.matrix()

    # calculate coverage (proportion of non-dashes in the alignment)
    coverage <- {aligned != '-'} %>%
      colSums()

    proteins$`coverage%`[i] <- sum(coverage > 1) / length(coverage) * 100
  }

  # checkpoint
  if(progress[stage, 'generate'])
    save(proteins, file = stage)
}



##### Write sqlite file #####

stage <- with(config, file.path(output_dir, out_sqlite))

if(progress[stage,'run']){
  # create a new sqlite database
  con <- dbConnect(SQLite(), dbname = stage)

  # write the data to the database
  dbWriteTable(con, 'peptides', peptides, overwrite = TRUE)
  dbWriteTable(con, 'proteins', proteins, overwrite = TRUE)

  # close the connection
  dbDisconnect(con)
}

##### Format Excel wb #####

stage <- with(config, file.path(output_dir, wb_checkpoint))

if(progress[stage,'load'])
{
  load(stage)
}else if(progress[stage,'run']){
  # infer number of columns to indent peptides
  peptide_indent <- 1
  
  # connect to the database
  con <- dbConnect(SQLite(), dbname = with(config, file.path(output_dir, out_sqlite)))


  # create a new workbook
  wb <- createWorkbook()

  addWorksheet(wb, config$sheet)

  # keep track of what row we're on
  nextRow <- 1

  # keep track of where we put things
  indices <- list(protein_rows    = integer(),
                  peptide_headers = integer(),
                  peptide_rows    = integer())

  for(i in 1:3)#dim(proteins)[1])
  {
    proteins[i,] |>
      
      mutate_all(~ ifelse(is.nan(.), NA, .)) |> # convert a few NaN's to NA
      
      writeData(wb = wb,
                sheet = config$sheet,
                startRow = nextRow,
                startCol = 1,
                rowNames = FALSE,
                colNames = i == 1)

    # account for header row on the first protein
    if(i == 1)
    {
      nextRow <- nextRow + 1
    }
    indices$protein_rows <- c(indices$protein_rows, nextRow)

    # add peptides
    tmp <- filter(peptides, PROTEIN == as.character(proteins$Protein[i]))

    tmp %>%
      dplyr::select(-PROTEIN) |>
      
      mutate_all(~ ifelse(is.nan(.), NA, .)) |> # convert a few NaN's to NA
      
      writeData(wb = wb, 
                sheet = config$sheet,
                startRow = nextRow + 1, 
                startCol = peptide_indent, 
                rowNames = FALSE)
    
    indices$peptide_headers <- c(indices$peptide_headers, nextRow + 1)

    # group and hide peptides
    groupRows(wb, config$sheet, nextRow + 1:(dim(tmp)[1] + 1), hidden = TRUE)

    # update next row
    nextRow <- nextRow + dim(tmp)[1] + 2
  }

  indices$peptide_rows <- with(indices,
                               c((1:max(peptide_headers)+1)[-c(protein_rows, peptide_headers) + 1],
                                 max(peptide_headers) + 1:dim(peptides)[1]))

  # format protein rows
  addStyle(wb = wb,
           sheet = config$sheet,
           style = createStyle(fgFill = config$protein_header_fill,
                               borderColour = openxlsx_getOp("borderColour", "black"),
                               border = "LeftRight"),
           rows = 1,
           cols = 1:dim(proteins)[2],
           gridExpand = TRUE)
  addStyle(wb = wb,
           sheet = config$sheet,
           style = createStyle(fgFill = config$protein_rows_fill,
                               borderColour = openxlsx_getOp("borderColour", "grey70"),
                               border = "TopBottomLeftRight"),
           rows = indices$protein_rows,
           cols = 1:dim(proteins)[2],
           gridExpand = TRUE)

  # format peptide rows
  addStyle(wb = wb,
           sheet = config$sheet,
           style = createStyle(fgFill = config$peptide_header_fill,
                               borderColour = openxlsx_getOp("borderColour", "black"),
                               border = "LeftRight"),
           rows = indices$peptide_headers,
           cols = 3+2:max(which(!is.na(proteins[2,]))),
           gridExpand = TRUE)
  addStyle(wb = wb,
           sheet = config$sheet,
           style = createStyle(fgFill = config$peptide_rows_fill,
                               borderColour = openxlsx_getOp("borderColour", "grey70"),
                               border = "TopBottomLeftRight"),
           rows = indices$peptide_rows,
           cols = 3+2:max(which(!is.na(proteins[2,]))),
           gridExpand = TRUE)

  #openXL(wb)
  if(progress[stage, 'generate'])
    save(wb, indices, file = stage)
}


##### Write Excel file #####

stage <- with(config, file.path(output_dir, out_xlsx))

if(progress[stage,'run']){
  saveWorkbook(wb, file = stage)
}
