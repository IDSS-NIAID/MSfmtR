
library(MSstats)

library(dplyr)
library(tidyr)
library(purrr)
library(readr)

library(openxlsx)

library(foreach)

dirs <- c('aspinwallja-20230825-PI1', 'aspinwallja-20230825-PI2', 'aspinwallja-20230825-PI3')
files <- file.path('data', dirs) |>
  list.files() |>
  grep(pattern = 'MSStats.tsv', value = TRUE)
paths <- file.path('data', dirs, files)

# format sytles
protein_header <- createStyle(fgFill = "#A7CDF0")
protein_rows   <- createStyle(fgFill = "#DDEBF7")
peptide_header <- createStyle(fgFill = "#F0CBA8")
peptide_rows   <- createStyle(fgFill = "#FCE4D6")

foreach(j=1:length(paths)) %do%
{
  if(!file.exists(paste0('output/', dirs[j], '.RData')))
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
    save(data, proteins, peptides, file = paste0('output/', dirs[j], '.RData'))
  }else{
    load(paste0('output/', dirs[j], '.RData'))
  }
  
  if(!file.exists(paste0('output/', dirs[j], '.xlsx')))
  {
  
    wb <- createWorkbook()
  
    addWorksheet(wb, dirs[j])
  
    nextRow <- 1 # keep track of what row we're on
  
    for(i in 1:dim(proteins)[1])
    {
      proteins[i,-which('peptides' == names(proteins))] %>%
        writeData(wb = wb, sheet = dirs[j], startRow = nextRow, startCol = 1, rowNames = FALSE,
                  colNames = i == 1)
    
      # adjustment and formatting for column names on the first protein
      if(i == 1)
      {
        nextRow <- nextRow + 1
        
        addStyle(wb = wb, sheet = dirs[j], style = protein_header,
                 rows = 1, cols = 1:(dim(proteins)[2] - 1),
                 gridExpand = TRUE)
      }
      
      # format new protein row
      addStyle(wb = wb, sheet = dirs[j], style = protein_rows,
               rows = nextRow, cols = 1:(dim(proteins)[2] - 1),
               gridExpand = TRUE)
      
      # add peptides
      proteins$peptides[[i]] %>%
        select(-PROTEIN) %>%
        writeData(wb = wb, sheet = dirs[j], startRow = nextRow + 1, startCol = 2, rowNames = FALSE)
      
      # format peptides
      addStyle(wb = wb, sheet = dirs[j], style = peptide_header,
               rows = nextRow + 1,
               cols = 2:dim(peptides)[2],
               gridExpand = TRUE)
      addStyle(wb = wb, sheet = dirs[j], style = peptide_rows,
               rows = nextRow + 2:(dim(proteins$peptides[[i]])[1] + 1),
               cols = 2:dim(proteins$peptides[[i]])[2],
               gridExpand = TRUE)
      
      # group and hide peptides
      groupRows(wb, dirs[j], nextRow + 1:(dim(proteins$peptides[[i]])[1] + 1), hidden = TRUE)
    
      # update next row
      nextRow <- nextRow + dim(proteins$peptides[[i]])[1] + 2
    }
  
    # Write data to excel
    #openXL(wb)
    saveWorkbook(wb, paste0('output/', dirs[j], '.xlsx'))
  }
}
