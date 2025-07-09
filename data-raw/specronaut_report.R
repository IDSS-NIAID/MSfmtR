# generation of sample data from a simulated MS experiment as output by Spectronaut

library(readr)

root <- here::here()

# This is the sample data from the MSstatsConvert package with a minor change to
# use the actual accession number, rather than the expanded information in their file.
# We also save the file using tab as the delimiter, as this is what gets exported from Spectronaut.
spectronaut_report <- read_csv('https://raw.githubusercontent.com/Vitek-Lab/MSstatsConvert/refs/heads/devel/inst/tinytest/raw_data/Spectronaut/spectronaut_input.csv') |>
  mutate(PG.ProteinAccessions = 'O75475',
         PG.ProteinGroups = 'O75475')

write_tsv(spectronaut_report, file.path(root, "inst", "extdata", "O75475.tsv"))
