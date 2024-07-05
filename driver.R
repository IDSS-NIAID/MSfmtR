
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

source('R/configure_formatR.R')
source('R/check_progress.R')
source('R/process_raw.R')
source('R/process_peptides.R')
source('R/process_proteins.R')
source('R/process_wb.R')

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

# load configurations
config <- configure_formatR(config_file = 'config_aspinwall2.yml',
                            args = args)

##################
# check progress #
##################

progress <- check_progress(config)


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
  data <- process_raw(config, save_intermediate = progress[stage, 'generate'])
}


##### Process peptide stats #####

stage <- with(config, file.path(output_dir, peptide_checkpoint))

if(progress[stage,'load'])
{
  load(stage)
}else if(progress[stage,'run']){
  peptides <- process_peptides(data, config, save_intermediate = progress[stage, 'generate'])
}


##### Process protein stats #####

stage <- with(config, file.path(output_dir, protein_checkpoint))

if(progress[stage,'load'])
{
  load(stage)
}else if(progress[stage,'run']){
  proteins <- process_proteins(data, peptides, config, save_intermediate = progress[stage, 'generate'])
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
  wb <- process_wb(proteins, peptides, config,
                   save_intermediate = progress[stage, 'generate'])
}


##### Write Excel file #####

stage <- with(config, file.path(output_dir, out_xlsx))

if(progress[stage,'run']){
  saveWorkbook(wb, file = stage)
}
