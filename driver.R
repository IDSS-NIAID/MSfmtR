
library(MSstats)
library(dplyr)
library(readr)


# Load data
data <- read_delim("data/aspinwallja-20230825-PI1/20240117_101231_aspinwallja-20230825-PI1_3reps_Report-MSStats.tsv", 
                   delim = "\t", col_names = TRUE) %>%
  SpectronauttoMSstatsFormat() %>%
  dataProcess()

#processed <- dataProcess(data)

