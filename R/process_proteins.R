#' process_proteins
#' Process proteins from raw data
#'
#' @param data MSstats formatted data
#' @param peptides processed and formatted peptides data
#' @param config list of configuration parameters
#' @param save_intermediate logical save intermediate data
#' @param ... additional arguments to pass to `updt_config`
#'
#' @details This function processes MSstats formatted data and returns proteins ready for ProtResDash. If save_intermediate is TRUE, the processed data are also saved to the checkpoint file.
#' @return data frame of processed proteins data
#' @export
#' @importFrom Biostrings AAStringSet readAAStringSet AAStringSetList
#' @importFrom dplyr group_by summarize ungroup arrange left_join mutate tibble bind_rows starts_with
#' @importFrom MSstats groupComparison
#' @importFrom lme4 fixef lmer lmerControl .makeCC
#' @importFrom msa msa
#' @importFrom purrr map map2 map_chr map_dbl map_int map_lgl
#' @importFrom stats median
#' @importFrom stringr fixed str_locate_all str_split str_sub_all str_length str_replace
#' @importFrom tidyr pivot_wider
#' @importFrom UniProt.ws UniProt.ws
process_proteins <- function(data, peptides, config, save_intermediate = TRUE, ...)
{
  # for those pesky no visible binding warnings
  if(FALSE)
    Entry <- Organism <- Protein.names <- Sequence <- Protein <- Label <- log2FC <- pvalue <- adj.pvalue <-
      GROUP <- originalRUN <- SUBJECT <- LogIntensities <- id <- run <- NumMeasuredFeature <- MissingPercentage <-
      more50missing <- NumImputedFeature <- TotalGroupMeasurements <- Intensity <- Description <-
      PROTEIN <- Modification <- RUN <- group <- models <- primary_id <- NULL

  # update config and pull package defaults if needed
  config <- updt_config(config, ...)

  # contrasts
  contrasts <- matrix(0, nrow = nrow(data$ratios), ncol = length(levels(data$FeatureLevelData$GROUP)),
                      dimnames = list(paste(data$ratios[,1], '/', data$ratios[,2]),
                                      levels(data$FeatureLevelData$GROUP)))

  # if config$groups is defined, we will need this mapping between Group and orginalRUN
  if(!is.null(config$groups))
  {
    group_mapping <- data$FeatureLevelData |>
      dplyr::select(GROUP, originalRUN) |>
      distinct()
  }


  for(i in 1:nrow(data$ratios))
  {
    if(!data$ratios[i,1] %in% levels(data$FeatureLevelData$GROUP))
    {
      num <- filter(group_mapping, grepl(data$ratios[i,1], originalRUN)) |>
        dplyr::select(GROUP) |>
        distinct() |>
        unlist() |>
        as.character()
    }else{
      num <- data$ratios[i,1]
    }

    if(!data$ratios[i,2] %in% levels(data$FeatureLevelData$GROUP))
    {
      den <- filter(group_mapping, grepl(data$ratios[i,2], originalRUN)) |>
        dplyr::select(GROUP) |>
        distinct() |>
        unlist() |>
        as.character()
    }else{
      den <- data$ratios[i,2]
    }

    contrasts[i, num] <-  1
    contrasts[i, den] <- -1
  }

  ### Metadata to add ###

  # if we have fasta files provided, use them to get protein metadata
  if(length(config$fasta_meta) > 0)
  {
    # read fasta
    fasta <- map(config$fasta_meta, ~ readAAStringSet(file.path(config$fasta_dir, .x))) |>
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
  }else{
    meta <- tibble(Protein = NA,
                   Description = NA,
                   Organism = NA,
                   Sequence = NA) |>
      dplyr::filter(!is.na(Protein))
  }

  # get UniProt metadata for anything not in `meta`
  # manually search individual IDs here: https://www.uniprot.org/uniprotkb

  # missing protein IDs
  proteinIDs <- data$ProteinLevelData$Protein |>
    str_split(fixed(';')) |>
    unlist()
  proteinIDs <- proteinIDs[!proteinIDs %in% meta$Protein]

  if(length(proteinIDs) > 0)
  {
    meta <- UniProt.ws::select(UniProt.ws(taxId=config$taxId),
                               proteinIDs,
                               c('organism_name', 'protein_name', 'sequence')) |>
      dplyr::select(Entry, Protein.names, Organism, Sequence) |>
      dplyr::rename(Protein = Entry,
                    Description = Protein.names) |>

      bind_rows(meta)
  }


  ### calculate protein stats ###
  proteins <- groupComparison(contrast.matrix = contrasts,
                              data = data,
                              use_log_file = FALSE)$ComparisonResult |>

    dplyr::filter(!is.na(Protein)) |>

    dplyr::select(Protein, Label, log2FC, pvalue, adj.pvalue) |>

    pivot_wider(names_from = Label, values_from = c(log2FC, pvalue, adj.pvalue))

  ### Extract protein data ###

  # add sample/group information if provided
  data$ProteinLevelData <- map_samples(data$ProteinLevelData, config) |>

    map_groups(config)

  # calculate group abundance statistics for peptides
  tmp <- data$ProteinLevelData |>

    mutate(id = paste0('Abundance: ', GROUP, ', ', originalRUN, ' (', SUBJECT, ')'))   # create a unique ID for each measurement

  if(config$merge_method == 'median')
  {
    # add summary (by median) of log intensities for each protein
    tmp_summ <- data$ProteinLevelData |>
      group_by(Protein, group) |>
      summarize(LogIntensities = median(LogIntensities, na.rm = TRUE)) |>
      ungroup() |>
      dplyr::rename(id = group) |>
      mutate(id = paste0('Group Abundance: ', id))
  }else if(config$merge_method == 'mean'){
    # add summary (by mean) of log intensities for each protein
    tmp_summ <- data$ProteinLevelData |>
      group_by(Protein, group) |>
      summarize(LogIntensities = mean(LogIntensities, na.rm = TRUE)) |>
      ungroup() |>
      dplyr::rename(id = group) |>
      mutate(id = paste0('Group Abundance: ', id))
  }else if(config$merge_method == 'lmer'){
    # add summary (by lmer) of log intensities for each protein
    tmp_summ <- data$ProteinLevelData |>
      filter(!is.na(LogIntensities)) |>

      group_by(Protein, group) |>

      summarize(models = list(try(lmer(LogIntensities ~ 1 + (1 | sample),
                                       control = lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))))), # we use lmerControl to silence singular fit warnings. We expect to see a few of these when there are few technical replicates (e.g. if there are 3 replicates). Since we are only concerned with the fixed effects estimates, we should be fine even with a singular fit.
                LogIntensities = map_dbl(models, ~
                                           {
                                             if('try-error' %in% class(.x))
                                               return(as.numeric(NA))
                                             fixef(.x)}
                                         )) |>
      dplyr::select(-models) |>
      ungroup() |>
      dplyr::rename(id = group) |>
      mutate(id = paste0('Group Abundance: ', id))
  }else{
    stop('Invalid merge_method')
  }


  tmp <- bind_rows(tmp, tmp_summ) |>

    mutate(Intensity = 10^LogIntensities,                                  # report linear scale
           primary_id = str_split(Protein, fixed(';')) |>                  # pick a primary protein ID for now - still need to deal with protein groups
             map_chr(~ .x[1])) |>

    arrange(is.na(group), id) |>                                           # sort by ID - this keeps column names in the same order as in `peptides`

    dplyr::select(Protein, primary_id, id, Intensity) |>

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
             map_chr(~
                    {
                      if(any(is.na(.x)))
                        return(NA)

                      if(length(.x) > 1)
                        return(.x[[2]])

                      .x
                    }) |>
             str_replace(' OS$', '')) |>

    left_join({group_by(peptides, PROTEIN) |>
               summarize(Modifications = unique(Modification) |>
                           grep(pattern = '^\\s*$', invert = TRUE, value = TRUE) |> # get rid of empty strings
                           paste0(collapse = '; ')) |>
               ungroup()},
              by = c('Protein' = 'PROTEIN')) |>               # merge in modification summary

    left_join(proteins, by = c('Protein' = 'Protein'))


  # Calculate coverage
  for(i in 1:nrow(proteins))
  {
    if(is.na(proteins$Sequence[i]))
      next

    seqs <- c(AAStringSet(proteins$Sequence[i]),
              AAStringSet(filter(peptides, PROTEIN == as.character(proteins$Protein[i]))$PEPTIDE |> unique()))

    # align sequences and convert to matrix
    aligned <- try({msa(seqs)@unmasked |>
                    as.matrix()})

    # calculate coverage (proportion of non-dashes in the alignment)
    if(! 'try-error' %in% class(aligned))
    {
      coverage <- {aligned != '-'} |>
        colSums()

      proteins$`coverage%`[i] <- sum(coverage > 1) / length(coverage) * 100
    }
  }

  # Cut off ratios at config$max_ratio
  for(i in names(proteins)[starts_with('log2FC', vars = names(proteins))])
  {
    proteins[[i]][proteins[[i]] >  log2(config$max_ratio)] <-  log2(config$max_ratio)
    proteins[[i]][proteins[[i]] < -log2(config$max_ratio)] <- -log2(config$max_ratio)
  }

  # checkpoint
  if(save_intermediate)
    save(proteins, file = file.path(config$output_dir, config$protein_checkpoint))

  return(proteins)
}