
process_wb <- function(proteins, peptides, config,
                       stage = file.path(config$output_dir, config$wb_checkpoint),
                       save_intermediate = TRUE,
                       n_proteins = NULL)
{
  # output order of protein columns
  proteins <- dplyr::select(proteins,
                            Protein, Description, Organism, nAA, `coverage%`, `mass (kDa)`,
                            Modifications,
                            starts_with('Abundance'),
                            starts_with('Group Abundance'),
                            starts_with('log2FC'),
                            starts_with('pvalue'),
                            starts_with('adj.pvalue'))

  # output order of peptide columns
  peptides <- dplyr::select(peptides,
                            PROTEIN, PEPTIDE, TRANSITION, FEATURE,
                            Modification,
                            starts_with('Abundance'),
                            starts_with('Group Abundance'),
                            starts_with('cv'))

  # fix a bug in 63456e without rerunning the whole thing
  names(proteins) <- str_replace_all(names(proteins), '  ', ' ')
  
  # make sure columns line up
  if(any(names(proteins)[starts_with('Abundance', vars = names(proteins))] !=
         names(peptides)[starts_with('Abundance', vars = names(peptides))]))
  {
    stop('Abundance columns do not match between proteins and peptides')
  }

  if(any(names(proteins)[starts_with('Group Abundance', vars = names(proteins))] !=
         names(peptides)[starts_with('Group Abundance', vars = names(peptides))]))
  {
    stop('Group Abundance columns do not match between proteins and peptides')
  }


  # figure out which columns contain modifications
  prot_mod_col <- which(names(proteins) == 'Modifications')
  pep_mod_col <- which(names(peptides) == 'Modification')

  # infer number of columns to indent peptides (account for dropping PROTEIN column from peptides)
  peptide_indent <- prot_mod_col - pep_mod_col + 2

  # # connect to the database
  # con <- dbConnect(SQLite(), dbname = with(config, file.path(output_dir, out_sqlite)))

  # create a new workbook
  wb <- createWorkbook()

  addWorksheet(wb, config$sheet)

  # keep track of what row we're on
  nextRow <- 1

  # keep track of where we put things
  # (first row is protein_headers,
  #  one row for each protein plus a peptide header row for each protein,
  #  one row for each peptide)
  indices <- tibble(row = 1:(1 + nrow(proteins)*2 + nrow(peptides)),
                    protein_rows    = FALSE,
                    peptide_headers = FALSE,
                    peptide_rows    = FALSE)

  for(i in 1:min(n_proteins, dim(proteins)[1]))
  {
    proteins[i,] |>
      
      mutate(Protein = as.character(Protein)) |>

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

  #openXL(wb)
  if(save_intermediate)
    save(wb, indices, file = stage)

  return(wb)
}