
process_proteins <- function(data, config,
                             stage = file.path(config$output_dir, config$peptide_checkpoint))
{
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

    left_join({group_by(peptides, PROTEIN) |>
               summarize(Modifications = paste0(unique(Modification), collapse = '; ')) |>
               ungroup()},
              by = c('Protein' = 'PROTEIN')) |>               # merge in modification summary

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
  if(save_intermediate)
    save(proteins, file = stage)

  return(proteins)
}