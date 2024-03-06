
library(MSstats)

library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)

library(openxlsx)


dirs <- c('aspinwallja-20230825-PI1', 'aspinwallja-20230825-PI2', 'aspinwallja-20230825-PI3')
files <- map_chr(dirs, ~ file.path('data', .x) |>
                   list.files() |>
                   grep(pattern = 'tsv', value = TRUE))
paths <- file.path('data', dirs, files)

# format styles
protein_header_style <- createStyle(fgFill = "#A7CDF0")
protein_rows_style   <- createStyle(fgFill = "#DDEBF7")
peptide_header_style <- createStyle(fgFill = "#F0CBA8")
peptide_rows_style   <- createStyle(fgFill = "#FCE4D6")


# do work
if(file.exists(paste0('output/', dirs, '_processed.RData')))
{
  load(paste0('output/', dirs, '_processed.RData'))
}else{
  if(file.exists(paste0('output/', dirs, '_msstats.RData')))
  {
    load(paste0('output/', dirs, '_msstats.RData'))
  }else{
    ########## MSStats step ##########
    # Load and process the data
    data <- read_delim(paths,
                       delim = "\t", col_names = TRUE) %>%
      SpectronauttoMSstatsFormat() %>%
      dataProcess()

    save(data, file = paste0('output/', dirs, '_msstats.RData'))
  }


  ########## Remainder of processing ##########
  # contrasts for protein statistics
  ratios <- rbind(c('EME; 0H',  'UN; 0H' ),
                  c('EME; 2H', 'EME; 0H' ),
                  c( 'UN; 2H',  'UN; 0H' ),
                  c('EME; 2H',  'UN; 2H' ),
                  c('EME; 24H','EME; 2H' ),
                  c( 'UN; 24H', 'UN; 2H' ),
                  c('EME; 72H','EME; 2H' ),
                  c( 'UN; 72H', 'UN; 2H' ),
                  c('EME; 24H', 'UN; 24H'),
                  c('EME; 48H','EME; 24H'),
                  c( 'UN; 48H', 'UN; 24H'),
                  c('EME; 48H', 'UN; 48H'),
                  c('EME; 72H','EME; 48H'),
                  c( 'UN; 72H', 'UN; 48H'),
                  c('EME; 72H', 'UN; 72H'))

  contrasts <- matrix(0, nrow = nrow(ratios), ncol = length(levels(data$FeatureLevelData$GROUP)),
                      dimnames = list(paste(ratios[,1], '/', ratios[,2]),
                                      levels(data$FeatureLevelData$GROUP)))

  for(i in 1:nrow(ratios))
  {
    contrasts[i,c(ratios[i,1],
                  ratios[i,2])] <- c(1,-1)
  }


  ###################################################################
  # Merge peptide and protein data into a single, nested data frame #
  ###################################################################

  # calculate peptide stats
  peptides <- data$FeatureLevelData %>%

    rename(abund = newABUNDANCE) %>% # newABUNDANCE is the normalized abundances

    group_by(PROTEIN, FEATURE, GROUP) %>%

    summarize(abund = median(abund, na.rm = TRUE),
              cv = sd(abund, na.rm = TRUE) / mean(abund, na.rm = TRUE) * 100) %>%

    ungroup() %>%

    rename(id = GROUP) %>%

    pivot_wider(names_from = id, values_from = c(abund, cv))


  # Extract peptide data
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


  # calculate protein stats
  proteins <- groupComparison(contrast.matrix = contrasts, data = data)$ComparisonResult %>%

    dplyr::select(Protein, Label, log2FC, pvalue, adj.pvalue) %>%

    pivot_wider(names_from = Label, values_from = c(log2FC, pvalue, adj.pvalue))


  # Extract protein data
  proteins <- data$ProteinLevelData %>%

    mutate(id = paste0(GROUP, ', ', originalRUN, ' (', SUBJECT, ')'),
           blank1 = '',
           blank2 = '',
           blank3 = '') %>%

    select(Protein, blank1, blank2, blank3, TotalGroupMeasurements, LogIntensities, id) %>%

    pivot_wider(names_from = id, values_from = LogIntensities) %>%

    mutate(peptides = map(Protein, ~ peptides %>% filter(PROTEIN == .x))) %>%

    left_join(proteins, by = c('Protein' = 'Protein'))

  # remove "blank*" column names
  names(proteins)[grep('blank', names(proteins))] <- ''
  save(proteins, file = paste0('output/', dirs, '_processed.RData'))

  # release some memory
  rm(data, peptides)
}


####################
# Write Excel file #
####################

if(file.exists(paste0('output/', dirs, '_wb.RData')))
{
  load(paste0('output/', dirs, '_wb.RData'))
}else{
  wb <- createWorkbook()

  addWorksheet(wb, dirs)

  # keep track of what row we're on
  nextRow <- 1

  # keep track of where we put things
  indices <- list(protein_rows    = integer(),
                  peptide_headers = integer(),
                  peptide_rows    = integer())

  for(i in 1:dim(proteins)[1])
  {
    proteins[i,-which('peptides' == names(proteins))] %>%
      writeData(wb = wb, sheet = dirs, startRow = nextRow, startCol = 1, rowNames = FALSE,
                colNames = i == 1)

    # account for header row on the first protein
    if(i == 1)
    {
      nextRow <- nextRow + 1
    }
    indices$protein_rows <- c(indices$protein_rows, nextRow)

    # add peptides
    proteins$peptides[[i]] %>%
      select(-PROTEIN) %>%
      writeData(wb = wb, sheet = dirs, startRow = nextRow + 1, startCol = 2, rowNames = FALSE)
    indices$peptide_headers <- c(indices$peptide_headers, nextRow + 1)

    # group and hide peptides
    groupRows(wb, dirs, nextRow + 1:(dim(proteins$peptides[[i]])[1] + 1), hidden = TRUE)

    # update next row
    nextRow <- nextRow + dim(proteins$peptides[[i]])[1] + 2
  }

  indices$peptide_rows <- with(indices,
                               c((1:max(peptide_headers)+1)[-c(protein_rows, peptide_headers) + 1],
                                 max(peptide_headers) + 1:dim(proteins$peptides[[i]])[1]))

  # Write data to excel
  #openXL(wb)
  save(wb, indices, file = paste0('output/', dirs, '_wb.RData'))
}

###############
# Format Rows #
###############

# format proteins
addStyle(wb = wb, sheet = dirs, style = protein_header_style,
         rows = 1, cols = 1:dim(proteins)[2],
         gridExpand = TRUE)
addStyle(wb = wb, sheet = dirs, style = protein_rows_style,
         rows = indices$protein_rows,
         cols = 1:dim(proteins)[2],
         gridExpand = TRUE)

# format peptides
addStyle(wb = wb, sheet = dirs, style = peptide_header_style,
         rows = indices$peptide_headers,
         cols = 2:max(which(!is.na(proteins[2,]))),
         gridExpand = TRUE)
addStyle(wb = wb, sheet = dirs, style = peptide_rows_style,
         rows = indices$peptide_rows,
         cols = 2:max(which(!is.na(proteins[2,]))),
         gridExpand = TRUE)


# save this monster
saveWorkbook(wb, file = paste0('output/', dirs, '.xlsx'))
