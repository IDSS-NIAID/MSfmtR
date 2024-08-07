
configure_formatR <- function(config_file = NULL, args = NULL)
{
  if(!is.null(args$config_file))
  {
    config_file <- args$config_file
  }

  # read yaml file (requires config package)
  if(!is.null(config_file))
  {
    config <- config::get(file = config_file)
  }else{
    config <- list()
  }

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
  if(is.null(config$max_ratio))
    config$max_ratio <- 100


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

  return(config)
}