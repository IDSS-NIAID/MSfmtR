#' check_progress
#' Check progress of previous runs and determine what needs to be done
#'
#' @param config list of configuration options
#'
#' @return data.frame with progress information
#' @export
#' @importFrom dplyr mutate
#' @importFrom purrr map_df map_lgl
check_progress <- function(config)
{
  # for those pesky no visible binding warnings
  if(FALSE)
    conf_name <- generate <- size <- NULL

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

  return(progress)
}