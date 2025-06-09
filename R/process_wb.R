#' process_wb
#' Format peptides and proteins for export to Excel
#'
#' @param proteins data frame of proteins data
#' @param peptides data frame of peptides data
#' @param config list of configuration parameters
#' @param save_intermediate logical save intermediate data
#' @param n_proteins integer number of proteins to process (use a small number to test output prior to running a large data set)
#' @param wb an existing workbook object (optional). If provided a new sheet will be added to the existing workbook with the name specified by the `sheet` configuration parameter.
#' @param sort_cols logical sort columns, assuming `proteins` and `peptides` have many overlapping columns, potentially with some blank columns needed to maintain registration of like-columns. If FALSE the columns will be in the order they were provided.
#' @param protein_alignment character, `proteins` column name to use for alignment. Default is to align `proteins$Modifications` with `peptides$Modification` in the excel worksheet (i.e. default: "Modifications").
#' @param peptide_alignment character, `peptides` column name to use for alignment. Default is to align `proteins$Modifications` with `peptides$Modification` in the excel worksheet (i.e. default: "Modification").
#' @param protein_pid character, column of `proteins` containing protein IDs. Default is "Protein".
#' @param peptide_pid character, column of `peptides` containing protein IDs. Default is "PROTEIN".
#' @param peptide_sep character, separator for protein accession numbers in `peptides`, for peptides with multiple possible proteins of origin. Default is ";".
#' @param overwrite logical, force `process_wb` to run again, even when saved intermediate data already exist
#' @param start_col_matching character, start matching column names at the first column that starts with `start_col_matching`. Default is "Abundance".
#' @param end_col_matching character, start matching column names at the first column that starts with `end_col_matching`. Default is "Found in Sample".
#' @param ... additional arguments passed to `updt_config`. Configuration parameters provided in this function will override those in `config`.
#'
#' @details This function processes formatted protein and peptide data and returns an excel workbook object ready for export to Excel. If save_intermediate is TRUE, the processed data are also saved to the checkpoint file.
#'
#' @return An excel workbook object
#' @export
#' @importFrom dplyr tibble mutate mutate_all starts_with
#' @importFrom openxlsx createWorkbook addWorksheet writeData groupRows addStyle createStyle
#' @importFrom stringr str_replace_all str_trim
process_wb <- function(proteins, peptides, config, save_intermediate = TRUE,
                       n_proteins = dim(proteins)[1], wb = NULL, sort_cols = TRUE,
                       protein_alignment = 'Modifications', peptide_alignment = 'Modification',
                       protein_pid = 'Protein', peptide_pid = 'PROTEIN', peptide_sep = ';',
                       overwrite = FALSE, 
                       start_col_matching = 'Abundance', end_col_matching = 'Found in Sample', ...)
{
  # for those pesky no visible binding warnings
  if(FALSE)
    Protein <- Description <- Organism <- nAA <- `coverage%` <- `mass (kDa)` <- Modifications <-
      PROTEIN <- PEPTIDE <- TRANSITION <- FEATURE <- Modification <- peptide_headers <- protein_rows <-
      peptide_rows <- NULL

  config <- updt_config(config, ...)

  # if we are using checkpoints (i.e. when save_intermediate is TRUE) load and return saved data
  checkpoint <- file.path(config$output_dir, config$wb_checkpoint)
  if(file.exists(checkpoint) & save_intermediate == TRUE & overwrite == FALSE)
  {
    load(checkpoint)

    paste("Loading saved data from", checkpoint) |>
      warning()

    return(wb)
  }


  # check sheet name for length requirements
  if(nchar(config$sheet) > 31)
  {
    warning('Sheet name is too long. Truncating to 31 characters.')
    config$sheet <- substr(config$sheet, 1, 31)
  }

  # check that alignment parameters are characters
  # if they provided a column number, convert to character
  if(is.numeric(protein_alignment))
    protein_alignment <- names(proteins)[protein_alignment]

  if(is.numeric(peptide_alignment))
    peptide_alignment <- names(peptides)[peptide_alignment]


  # sort columns?
  if(sort_cols)
  {
    peptide_cols_picked <- rep(FALSE, dim(peptides)[2])
    
    # find the first and last columns in the region where all column names should match
    # i.e. the columns that start with start_col_matching and end with end_col_matching
    start_matching <- starts_with(start_col_matching, vars = names(proteins)) |> min()
    end_matching   <- starts_with(  end_col_matching, vars = names(proteins)) |> max()
    
    # calculate number of columns to indent peptide columns
    buffer <- which(names(proteins) == protein_alignment) -
              which(names(peptides) == peptide_alignment)
    
    # may need to address this at some point, but don't expect it will be an issue most of the time
    if(buffer < 0)
    {
      stop('Peptide alignment column is to the right of protein alignment column')
    }
    if(buffer + 1 > start_matching)
    {
      stop('Start matching column is to the left of protein alignment column')
    }
    
    # we will keep the output order of protein columns unchanged
    
    # this is the order we'll pick the columns of peptides for the first few columns
    # i.e. until we reach the matching region (match what we can here, but allow mismatches)
    pick <- numeric(0)
    
    for(i in (buffer + 1):(start_matching - 1))
    {
      if(names(proteins)[i] %in% names(peptides))
      {
        # if the column is in peptides, pick it
        tmp <- which(names(peptides) == names(proteins)[i])
        peptide_cols_picked[tmp] <- TRUE
        
        pick <- c(pick, tmp)
      }else{
        # if not, pick the next column unpicked column in peptides
        # that isn't in proteins
        tmp <- which(!peptide_cols_picked & 
                     !names(peptides) %in% names(proteins))
        
        pick <- c(pick, tmp[1])
        peptide_cols_picked[tmp[1]] <- TRUE
      }
    }
    
    peptides_sorted <- peptides[, pick]
    blank_col <- tibble(` ` = rep(' ', dim(peptides)[1]))

    # continue picking columns in the matching region
    # i.e. all columns should have matching names in this region
    for(i in start_matching:end_matching)
    {
      if(names(proteins)[i] %in% names(peptides))
      {
        # if the column is in peptides, pick it
        tmp <- which(names(peptides) == names(proteins)[i])
        peptide_cols_picked[tmp] <- TRUE
        
        peptides_sorted[[ names(peptides)[tmp] ]] <- peptides[[tmp]]
      }else{
        # if not, add a blank column
        peptides_sorted <- cbind(peptides_sorted, blank_col)
        names(blank_col) <- paste0(names(blank_col), ' ') # keep names unique and invisible when in the excel sheet
      }
    }
    
    # continue picking columns after the matching region
    # these can be in any order and don't need to match
    if(end_matching < dim(proteins)[2])
      for(i in (end_matching + 1):dim(proteins)[2])
      {
        # if we run off the end of the peptides, just skip to the end
        if(all(peptide_cols_picked))
          break
      
        if(names(proteins)[i] %in% names(peptides))
        {
          # if the column is in peptides, pick it
          tmp <- which(names(peptides) == names(proteins)[i])
          peptide_cols_picked[tmp] <- TRUE
        
          peptides_sorted[[ names(peptides)[tmp] ]] <- peptides[[tmp]]
        }else{
          # if not, add whatever is left
          tmp <- which(!peptide_cols_picked & 
                       !names(peptides) %in% names(proteins))
        
          if(length(tmp) == 0)
          {
            # add a blank column if there is nothing left
            # (presumably there is another matching column later on, or we would have broken the loop above)
            peptides_sorted <- cbind(peptides_sorted, blank_col)
            names(blank_col) <- paste0(names(blank_col), ' ') # keep names unique and invisible when in the excel sheet
          }else{
            peptide_cols_picked[ tmp[1] ] <- TRUE
            peptides_sorted[[ names(peptides)[ tmp[1] ] ]] <- peptides[[tmp[1]]]
          }
        }
      }
    
    # add the remaining, unpicked column on the end if there are any
    if(any(!peptide_cols_picked))
    {
      peptides_sorted <- cbind(peptides_sorted, peptides[, !peptide_cols_picked])
    }
    
    peptides <- peptides_sorted
  }


  # figure out which columns contain modifications
  prot_mod_col <- which(names(proteins) == protein_alignment)
  pep_mod_col <- which(names(peptides) == peptide_alignment)

  # infer number of columns to indent peptides       (account for dropping PROTEIN column from peptides)
  peptide_indent <- prot_mod_col - pep_mod_col + 1 + ('PROTEIN' %in% names(peptides))

  # # connect to the database
  # con <- dbConnect(SQLite(), dbname = with(config, file.path(output_dir, out_sqlite)))

  # create a new workbook
  if(is.null(wb))
    wb <- createWorkbook()

  addWorksheet(wb, config$sheet)

  # keep track of what row we're on
  nextRow <- 1

  # how many peptide rows are we expecting?
  n_peptide_rows <- proteins[[protein_pid]] |>
    map_int(~ grepl(.x, peptides[[peptide_pid]]) |>
              sum()) |>
    sum()

  # keep track of where we put things
  # (first row is protein_headers,
  #  one row for each protein plus a peptide header row for each protein,
  #  one row for each peptide)
  indices <- tibble(row = 1:(1 + nrow(proteins)*2 + n_peptide_rows),
                    protein_rows    = FALSE,
                    peptide_headers = FALSE,
                    peptide_rows    = FALSE)

  for(i in 1:min(n_proteins, dim(proteins)[1]))
  {
    # start with the protein we are currently working with
    prot_curr <- proteins[i,]

    # make sure the protein ID is a character (not a factor)
    prot_curr[[protein_pid]] <- as.character(prot_curr[[protein_pid]])

    # write this row
    prot_curr |>

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
    indices$protein_rows[nextRow] <- TRUE

    # add peptides that match the current protein
    tmp <- peptides[grepl(prot_curr[[protein_pid]], peptides[[peptide_pid]]),]

    if(!is.null(tmp))
    {
      # this is the default for Spectronaut output, but not for PD
      # a little hacky, so we should revisit this at a later point
      if('PROTEIN' %in% names(tmp))
        tmp$PROTEIN <- NULL

      tmp |>

        mutate_all(~ ifelse(is.nan(.), NA, .)) |> # convert a few NaN's to NA

        writeData(wb = wb,
                  sheet = config$sheet,
                  startRow = nextRow + 1,
                  startCol = peptide_indent,
                  rowNames = FALSE)

      indices$peptide_headers[nextRow + 1] <- TRUE

      # group and hide peptides
      groupRows(wb, config$sheet, nextRow + 1:(dim(tmp)[1] + 1), hidden = TRUE)

      # update next row
      nextRow <- nextRow + dim(tmp)[1] + 2

      indices$peptide_rows[max(filter(indices, peptide_headers)$row + 1):(nextRow - 1)] <- TRUE
    }else{
      nextRow <- nextRow + 1
    }
  }

  # format protein header row
  addStyle(wb = wb,
           sheet = config$sheet,
           style = createStyle(fgFill = config$protein_header_fill,
                               borderColour = rep('black', 4),
                               border = "LeftRight"),
           rows = 1,
           cols = 1:dim(proteins)[2],
           gridExpand = TRUE)
  # format protein rows (non-abundance rows)
  addStyle(wb = wb,
           sheet = config$sheet,
           style = createStyle(fgFill = config$protein_rows_fill,
                               borderColour = rep('grey70', 4),
                               border = "TopBottomLeftRight"),
           rows = filter(indices, protein_rows)$row,
           cols = (1:dim(proteins)[2])[!grepl('Abundance', names(proteins))],
           gridExpand = TRUE)
  # format protein abundances
  addStyle(wb = wb,
           sheet = config$sheet,
           style = createStyle(fgFill = config$protein_rows_fill,
                               borderColour = rep('grey70', 4),
                               border = "TopBottomLeftRight",
                               numFmt = 'SCIENTIFIC'),
           rows = filter(indices, protein_rows)$row,
           cols = (1:dim(proteins)[2])[grepl('Abundance', names(proteins))],
           gridExpand = TRUE)

  if(!is.null(peptides))
  {
    # format peptide header rows
    addStyle(wb = wb,
             sheet = config$sheet,
             style = createStyle(fgFill = config$peptide_header_fill,
                                 borderColour = rep('black', 4),
                                 border = "LeftRight"),
             rows = filter(indices, peptide_headers)$row,
             cols = (peptide_indent - 1) + 1:(dim(peptides)[2] - 1),
             gridExpand = TRUE)
    # format peptide rows (non-abundance rows)
    addStyle(wb = wb,
             sheet = config$sheet,
             style = createStyle(fgFill = config$peptide_rows_fill,
                                 borderColour = rep('grey70', 4),
                                 border = "TopBottomLeftRight"),
             rows = filter(indices, peptide_rows)$row,
             cols = (peptide_indent - 1) + (1:(dim(peptides)[2] - 1))[!grepl('Abundance', names(peptides)[-1])],
             gridExpand = TRUE)
    # format peptide abundances
    addStyle(wb = wb,
             sheet = config$sheet,
             style = createStyle(fgFill = config$peptide_rows_fill,
                                 borderColour = rep('grey70', 4),
                                 border = "TopBottomLeftRight",
                                 numFmt = 'SCIENTIFIC'),
             rows = filter(indices, peptide_rows)$row,
             cols = (peptide_indent - 1) + (1:(dim(peptides)[2] - 1))[grepl('Abundance', names(peptides)[-1])],
             gridExpand = TRUE)
  }

  #openXL(wb)
  if(save_intermediate)
    save(wb, indices, file = file.path(config$output_dir, config$wb_checkpoint))

  return(wb)
}