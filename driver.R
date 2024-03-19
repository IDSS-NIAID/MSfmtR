
#' configuration options
#' Edit the `config.yml` file in this directory to change the configuration options, or include parameters via the command line.
#' Command line parameters will override the `config.yml` file. See `config.yml` for defaults.

## required parameters

#' @param input_dir The directory containing the data files. This script will only run one file at a time,
#'  so either the directory should contain only one file, or the `in_file` variable should be set to the file to be
#'  processed. If `in_file` is not provided it will look for a `.tsv` file to analyze inside of `input_dir`.
#' @param output_dir The directory to save the output files to. If it does not exist, it will be created.

## optional parameters

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
#' @param cont_fasta The name of the contaminant fasta file to use (assumed to be in `input_dir`). Default is the Universal Contaminant file that comes with this package.

### Metadata parameters
#' @param fasta The name of the fasta file(s) to use for annotation. If no annotation files are provided, `input_dir` will be checked for fasta files different from `cont_fasta`. If none are found, UniProt will be used for annotation.
#' @param taxId The taxonomy ID of the organism being analyzed. Default is `9606` (human). This is used to search for anything outside of the files in `fasta`.

### MSStats parameters
#' @param ratios The contrasts to be used in the MSStats analysis. Default is all combinations of all levels in `in_file`. Each ratio in the list should be of the form "<group1>/<group2>", so for group1=Case and group2=Control, the ratio would be "Case/Control". These labels should match values in the `R.Condition` column of `in_file`.

### Style parameters
#' @param protein_header_fill The fill color for the protein header in the output Excel file.
#' @param protein_rows_fill The fill color for the protein rows in the output Excel file.
#' @param peptide_header_fill The fill color for the peptide header in the output Excel file.
#' @param peptide_rows_fill The fill color for the peptide rows in the output Excel file.

### Checkpoint parameters
#' @param checkpoints A list of files to generate (if they already exist, this will be ignored). Options include: `xlsx,sql,processed,protein,peptide,wb`.
#' @param processed_checkpoint The name of the checkpoint file to save the processed data to. Default is `paste0(config$sheet, '_processed.RData')`.
#' @param protein_checkpoint The name of the checkpoint file to save the protein statistics to. Default is `paste0(config$sheet, '_protein.RData')`.
#' @param peptide_checkpoint The name of the checkpoint file to save the peptide statistics to. Default is `paste0(config$sheet, '_peptide.RData')`.
#' @param wb_checkpoint The name of the checkpoint file to save the R object containing the formatted Excel workbook. Default is `paste0(config$sheet, '_wb.RData')`.

#' @details This script is designed to be run from the command line. It will read the `config.yml` file in the
#' current directory and use the parameters to process the data.
#'
#' `ratios` should be a comma-separated list of contrasts to be used in the MSStats analysis, with the numerator and
#' denominator separated by "/" and each contrast separated by a comma. For example, `EME 0H/UN 0H,EME 2H/EME 0H`. The
#' values in each ratio (e.g. `EME 0H` and `UN 0H`) should match the values in the `R.Condition` column of `in_file`.
#'
#' @examples
#' Rscript driver.R --args input_dir=data output_dir=output


##################
# Load libraries #
##################

library(MSstats)

library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)

library(openxlsx)
library(config)
library(RSQLite)

library(Biostrings)
library(muscle)
library(UniProt.ws)


###############
# Load config #
###############

# read yaml file (requires config package)
config <- config::get()

# read command line arguments
args <- commandArgs(trailingOnly = TRUE) |>
  str_split('=') |>
  map_chr(~ .x[2])

names(args) <- commandArgs(trailingOnly = TRUE) |>
  str_split('=') |>
  map_chr(~ .x[1])

# replace any defaults with command line arguments
if(length(args) > 0)
  for(i in 1:length(args))
  {
    config[names(args)[i]] <- args[i]
  }

# check for required parameters
if(is.null(config$input_dir) | is.null(config$output_dir))
  stop('input_dir and output_dir are required parameters')

# calculate unspecified defaults

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
  config$lloq <- 0
if(is.null(config$cont_fasta))
{
  config$cont_fasta <- 'inst/extdata/Universal Contaminants.fasta'
}else if(!file.exists(file.path(config$input_dir, config$cont_fasta)))
{
  warning('cont_fasta does not exist, using default: ', config$cont_fasta)
  config$cont_fasta <- 'inst/extdata/Universal Contaminants.fasta'
}


## Metadata parameters
if(is.null(config$fasta))
{
  config$fasta <- list.files(config$input_dir, pattern = 'fasta') |>
    grep(pattern = config$cont_fasta, invert = TRUE, value = TRUE)
} 
if(is.null(config$taxId))
  config$taxId <- 9606


## MSStats parameters
# if ratios is null, we'll fill it in after reading in the raw data


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
                protein = c('processed', 'peptide_checkpoint'),
                peptide = c('processed'),
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
}else if(progress[stage,'run']){
  # Load and process the data
  raw <- with(config, file.path(input_dir, in_file)) |>
    read_delim(delim = "\t", col_names = TRUE)

  # format ratios
  if(is.null(config$ratios))
  {
    config$ratios <- raw |>
      select(R.Condition) |>
      distinct() |>
      pull() |>
      combn(2, simplify = TRUE) |>
      t()
  }else{
    tmp <- config$ratios |>
      str_split(',') |>
      unlist() |>
      str_split('/')

    config$ratios <- cbind(map_chr(tmp, ~ .x[1]), map_chr(tmp, ~ .x[2]))
  }


  data <- SpectronauttoMSstatsFormat(raw, use_log_file = FALSE) %>%
    dataProcess()

  # checkpoint
  if(progress[stage, 'generate'])
    save(data, config, file = stage)
}


##### Process peptide stats #####

stage <- with(config, file.path(output_dir, peptide_checkpoint))

if(progress[stage,'load'])
{
  load(stage)
}else if(progress[stage,'run']){
  # calculate peptide stats
  peptides <- data$FeatureLevelData %>%

    rename(abund = newABUNDANCE) %>% # newABUNDANCE is the normalized abundances

    group_by(PROTEIN, FEATURE, GROUP) %>%

    summarize(abund = median(abund, na.rm = TRUE),
              cv = sd(abund, na.rm = TRUE) / mean(abund, na.rm = TRUE) * 100) %>%

    ungroup() %>%

    rename(id = GROUP) %>%

    pivot_wider(names_from = id, values_from = c(abund, cv))


  # extract peptide data
  peptides <- data$FeatureLevelData %>%

    rename(abund = newABUNDANCE) %>%

    mutate(id = paste0(GROUP, ', ', originalRUN, ' (', SUBJECT, ')'),

           # split modifications out into a different column
           Modification = map_chr(PEPTIDE, ~
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
                                    }),

           PEPTIDE = map_chr(PEPTIDE, ~ strsplit(as.character(.x),
                                                 split = '_', fixed = TRUE)[[1]][1] %>%
                               gsub(pattern = '\\[.*?\\]', replacement = ''))) %>%


    select(PROTEIN, PEPTIDE, Modification, TRANSITION, FEATURE, id, abund) %>%

    pivot_wider(names_from = id, values_from = abund) %>%

    # merge stats into peptides
    left_join(peptides, by = c('PROTEIN', 'FEATURE'))

  # checkpoint
  if(progress[stage, 'generate'])
    save(peptides, config, file = stage)
}


##### Process protein stats #####

stage <- with(config, file.path(output_dir, protein_checkpoint))

if(progress[stage,'load'])
{
  load(stage)
}else if(progress[stage,'run']){
  contrasts <- matrix(0, nrow = nrow(config$ratios), ncol = length(levels(data$FeatureLevelData$GROUP)),
                      dimnames = list(paste(config$ratios[,1], '/', config$ratios[,2]),
                                      levels(data$FeatureLevelData$GROUP)))

  for(i in 1:nrow(config$ratios))
  {
    contrasts[i,c(config$ratios[i,1],
                  config$ratios[i,2])] <- c(1,-1)
  }

  ### calculate protein stats ###
  proteins <- groupComparison(contrast.matrix = contrasts, data = data)$ComparisonResult %>%

    dplyr::select(Protein, Label, log2FC, pvalue, adj.pvalue) %>%

    pivot_wider(names_from = Label, values_from = c(log2FC, pvalue, adj.pvalue))


  ### Extract protein data ###
  proteins <- data$ProteinLevelData %>%

    mutate(id = paste0(GROUP, ', ', originalRUN, ' (', SUBJECT, ')'),
           `coverage%` = NA) %>%

    select(Protein, `coverage%`, TotalGroupMeasurements, LogIntensities, id) %>%

    pivot_wider(names_from = id, values_from = LogIntensities) %>%

    left_join(proteins, by = c('Protein' = 'Protein'))


  ### Add metadata ###

  # all protein IDs
  proteinIDs <- proteins$Protein |>
    str_split(fixed(';')) |>
    unlist() |>
    str_replace(fixed('Cont_'), '') |>
    unique()

  # get UniProt metadata - manually search individual IDs here: https://www.uniprot.org/uniprotkb
  UPmeta <- UniProt.ws::select(UniProt.ws(taxId=config$taxId),
                             proteinIDs,
                             c('organism_name', 'protein_name', 'mass', 'sequence')) %>%
    dplyr::select(Entry, Organism, Protein.names, Mass, Sequence) %>%
    dplyr::rename(Protein = Entry,
                  Description = Protein.names)

  # merge metadata into proteins
  proteins <- proteins |>
    mutate(primary_id = str_split(Protein, fixed(';')) |>
             map_chr(~ .x[1]) |>
             str_replace(fixed('Cont_', ''))) |>
    left_join(UPmeta, proteins, by = c('primary_id' = 'Protein'))

  # Calculate coverage
  for(i in 1:nrow(proteins))
  {
    if(is.na(proteins$Sequence[i]))
       next

    # write FASTA file
    cat('>protein', proteins$Sequence[i], sep = '\n', file = "sequence.txt")

    pep <- filter(peptides, PROTEIN == proteins$Protein[i])$PEPTIDE %>% unique()
    cat(paste0('>pep', 1:length(pep), '\n', pep), sep = '\n', file = "sequence.txt", append = TRUE)

    # read sequences into Biostrings::StringSet object
    seqs <- Biostrings::readAAStringSet("sequence.txt")

    # align sequences and convert to matrix
    aligned <- muscle::muscle(seqs) %>%
      as.matrix()

    # calculate coverage
    coverage <- {aligned != '-'} %>%
      colSums()

    proteins$`coverage%`[i] <- sum(coverage > 1) / length(coverage) * 100
  }

  # clean up
  file.remove("sequence.txt")

  # checkpoint
  if(progress[stage, 'generate'])
    save(proteins, config, proteinIDs, UPmeta, file = stage)
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

  for(i in 1:dim(proteins)[1])
  {
    proteins[i,] %>%
      writeData(wb = wb, sheet = config$sheet, startRow = nextRow, startCol = 1, rowNames = FALSE,
                colNames = i == 1)

    # account for header row on the first protein
    if(i == 1)
    {
      nextRow <- nextRow + 1
    }
    indices$protein_rows <- c(indices$protein_rows, nextRow)

    # add peptides
    tmp <- filter(peptides, PROTEIN == proteins$Protein[i])

    tmp %>%
      select(-PROTEIN) %>%
      writeData(wb = wb, sheet = config$sheet, startRow = nextRow + 1, startCol = 2, rowNames = FALSE)
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
  addStyle(wb = wb, sheet = config$sheet, 
           style = createStyle(fgFill = config$protein_header_fill),
           rows = 1, cols = 1:dim(proteins)[2],
           gridExpand = TRUE)
  addStyle(wb = wb, sheet = config$sheet,
           style = createStyle(fgFill = config$protein_rows_fill),
           rows = indices$protein_rows,
           cols = 1:dim(proteins)[2],
           gridExpand = TRUE)

  # format peptide rows
  addStyle(wb = wb, sheet = config$sheet,
           style = createStyle(fgFill = config$peptide_header_fill),
           rows = indices$peptide_headers,
           cols = 2:max(which(!is.na(proteins[2,]))),
           gridExpand = TRUE)
  addStyle(wb = wb, sheet = config$sheet,
           style = createStyle(fgFill = config$peptide_rows_fill),
           rows = indices$peptide_rows,
           cols = 2:max(which(!is.na(proteins[2,]))),
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
