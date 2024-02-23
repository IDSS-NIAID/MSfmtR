
library(MSstats)

library(dplyr)
library(tidyr)
library(purrr)
library(readr)

library(openxlsx)

dirs <- c('aspinwallja-20230825-PI1', 'aspinwallja-20230825-PI2', 'aspinwallja-20230825-PI3')
files <- file.path('data', dirs) |>
  list.files() |>
  grep(pattern = 'MSStats.tsv', value = TRUE)
paths <- file.path('data', dirs, files)

# Write data to excel
file_out <- "output/20240117_101231_aspinwallja-20230825.xlsx"

wb <- createWorkbook()

for(j in 1:length(paths))
{
  # Load and process the data
  data <- read_delim(paths[j],
                     delim = "\t", col_names = TRUE) %>%
    SpectronauttoMSstatsFormat() %>%
    dataProcess()
  
  
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
    
    group_by(PROTEIN, PEPTIDE, TRANSITION, GROUP) %>%
    
    summarize(abund = median(abund, na.rm = TRUE),
              cv = sd(abund, na.rm = TRUE) / mean(abund, na.rm = TRUE) * 100) %>%
    
    ungroup() %>%
    
    rename(id = GROUP) %>%
    
    pivot_wider(names_from = id, values_from = c(abund, cv))
  
  
  # Extract peptide data
  peptides <- data$FeatureLevelData %>%
    
    rename(abund = newABUNDANCE) %>%
    
    mutate(id = paste0(GROUP, ', ', originalRUN, ' (', SUBJECT, ')')) %>%
    
    select(PROTEIN, PEPTIDE, TRANSITION, FEATURE, id, abund) %>%
    
    pivot_wider(names_from = id, values_from = abund) %>%
    
    # merge stats into peptides
    right_join(peptides, by = c('PROTEIN', 'PEPTIDE', 'TRANSITION'))
  
  
  # calculate protein stats
  proteins <- groupComparison(contrast.matrix = contrasts, data = data)$ComparisonResult %>%
    
    dplyr::select(Protein, Label, log2FC, pvalue, adj.pvalue) %>%
    
    pivot_wider(names_from = Label, values_from = c(log2FC, pvalue, adj.pvalue))
  
  
  # Extract protein data
  proteins <- data$ProteinLevelData %>%
    
    mutate(id = paste0(GROUP, ', ', originalRUN, ' (', SUBJECT, ')'),
           blank1 = '',
           blank2 = '') %>%
    
    select(Protein, blank1, blank2, TotalGroupMeasurements, LogIntensities, id) %>%
    
    pivot_wider(names_from = id, values_from = LogIntensities) %>%
    
    mutate(peptides = map(Protein, ~ peptides %>% filter(PROTEIN == .x))) %>%
    
    left_join(proteins, by = c('Protein' = 'Protein'))
  
  names(proteins)[2:3] <- ''
  
  addWorksheet(wb, dirs[j])
  
  nextRow <- 1 # keep track of what row we're on
  
  for(i in 1:dim(proteins)[1])
  {
    proteins[i,-which('peptides' == names(proteins))] %>%
      writeData(wb = wb, sheet = dirs[j], startRow = nextRow, startCol = 1, rowNames = FALSE,
                colNames = i == 1)
    
    # adjust for column names
    nextRow <- nextRow + as.numeric(i == 1)
    
    proteins$peptides[[i]] %>%
      select(-PROTEIN) %>%
      writeData(wb = wb, sheet = dirs[j], startRow = nextRow + 1, startCol = 2, rowNames = FALSE)
    groupRows(wb, dirs[j], nextRow + 0:dim(proteins$peptides[[i]])[1], hidden = TRUE)
    
    nextRow <- nextRow + dim(proteins$peptides[[i]])[1] + 2
  }
}

#openXL(wb)
saveWorkbook(wb, file_out)
