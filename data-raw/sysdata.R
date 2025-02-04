
#' defaults
#' Default configuration for MSfmtR
#' See README for details

defaults <- list(
  # Directory parameters
  input_dir = '.',
  fasta_dir = '.',
  output_dir = '.',

  # File parameters
  in_file = parse(text = "list.files(input_dir, pattern = 'tsv')"),
  out_xlsx = parse(text = "gsub('.tsv', '.xlsx', in_file, fixed = TRUE)"),
  out_sqlite = parse(text = "gsub('.tsv', '.sqlite', in_file, fixed = TRUE)"),
  sheet = parse(text = "gsub('.tsv', '', in_file, fixed = TRUE)"),

  # Sample parameters
  smp_grp_rep = 'R.FileName',
  merge_method = 'median',

  # Filtering parameters
  uloq = Inf,
  lloq = 0,
  cont_fasta = parse(text = "system.file('extdata/Universal_Contaminants.fasta', package = 'MSfmtR')"),
  max_ratio = 100,

  # Metadata parameters
  fasta_meta = parse(text = "list.files(fasta_dir, pattern = 'fasta') |> grep(pattern = cont_fasta, invert = TRUE, value = TRUE)"),
  taxId = 9606,

  # MSStats parameters
  ratios = NULL, # if ratios is null, we'll fill it in after reading in the raw data
  normMeasure = 'NormalizedPeakArea',
  preprocess = FALSE,
  format = 'MSstats',

  # Style parameters
  protein_header_fill = "#A7CDF0",
  protein_rows_fill = "#DDEBF7",
  peptide_header_fill = "#F0CBA8",
  peptide_rows_fill = "#FCE4D6",

  # Checkpoint parameters
  checkpoints = c('xlsx', 'sql', 'processed', 'protein', 'peptide', 'wb'),
  processed_checkpoint = parse(text = "paste0(config$sheet, '_processed.RData')"),
  peptide_checkpoint = parse(text = "paste0(config$sheet, '_peptide.RData')"),
  protein_checkpoint = parse(text = "paste0(config$sheet, '_protein.RData')"),
  wb_checkpoint = parse(text = "paste0(config$sheet, '_wb.RData')")
)

usethis::use_data(defaults, overwrite = TRUE, internal = TRUE)
