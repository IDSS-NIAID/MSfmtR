
process_wb <- function(proteins, peptides, config,
                       stage = file.path(config$output_dir, config$wb_checkpoint),
                       save_intermediate = TRUE)
{
  # infer number of columns to indent peptides
  peptide_indent <- 1

  # # connect to the database
  # con <- dbConnect(SQLite(), dbname = with(config, file.path(output_dir, out_sqlite)))

  # create a new workbook
  wb <- createWorkbook()

  addWorksheet(wb, config$sheet)

  # keep track of what row we're on
  nextRow <- 1

  # keep track of where we put things
  indices <- tibble(row = 1:(nrow(proteins)*2 + nrow(peptides) + 1),
                    protein_rows    = FALSE,
                    peptide_headers = FALSE,
                    peptide_rows    = FALSE)

  for(i in 1:dim(proteins)[1])
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
    indices$protein_rows[nextRow] <- TRUE

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

    indices$peptide_headers[nextRow + 1] <- TRUE

    # group and hide peptides
    groupRows(wb, config$sheet, nextRow + 1:(dim(tmp)[1] + 1), hidden = TRUE)

    # update next row
    nextRow <- nextRow + dim(tmp)[1] + 2

    indices$peptide_rows[max(filter(indices, peptide_headers)$row + 1):(nextRow - 1)] <- TRUE
  }

  # format protein rows
  addStyle(wb = wb,
           sheet = config$sheet,
           style = createStyle(fgFill = config$protein_header_fill,
                               borderColour = rep('black', 4),
                               border = "LeftRight"),
           rows = 1,
           cols = 1:dim(proteins)[2],
           gridExpand = TRUE)
  addStyle(wb = wb,
           sheet = config$sheet,
           style = createStyle(fgFill = config$protein_rows_fill,
                               borderColour = rep('grey70', 4),
                               border = "TopBottomLeftRight"),
           rows = filter(indices, protein_rows)$row,
           cols = 1:dim(proteins)[2],
           gridExpand = TRUE)

  # format peptide rows
  addStyle(wb = wb,
           sheet = config$sheet,
           style = createStyle(fgFill = config$peptide_header_fill,
                               borderColour = rep('black', 4),
                               border = "LeftRight"),
           rows = filter(indices, peptide_headers)$row,
           cols = (peptide_indent - 1) + 1:(dim(peptides)[2] - 1),
           gridExpand = TRUE)
  addStyle(wb = wb,
           sheet = config$sheet,
           style = createStyle(fgFill = config$peptide_rows_fill,
                               borderColour = rep('grey70', 4),
                               border = "TopBottomLeftRight"),
           rows = filter(indices, peptide_rows)$row,
           cols = (peptide_indent - 1) + 1:(dim(peptides)[2] - 1),
           gridExpand = TRUE)

  #openXL(wb)
  if(save_intermediate)
    save(wb, indices, file = stage)

  return(wb)
}