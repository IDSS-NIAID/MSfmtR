# test full Spectronaut workflow

library(readr)

##### check running with raw input #####
test_that("Specronaut workflow with raw input", {

  config <- configure_formatR(args = list(fasta_dir = system.file('extdata', package = 'MSfmtR')))

  raw <- read_delim(system.file('extdata', 'O75475.tsv', package = 'MSfmtR'),
                    delim = "\t", col_names = TRUE, show_col_types = FALSE)

  data <- process_raw(config, raw, save_intermediate = FALSE)

  expect_length(data, 4)
})


##### with no peptide roll up #####
test_that("Specronaut workflow with no peptide roll up", {

  config <- configure_formatR(args = list(input_dir = system.file('extdata', package = 'MSfmtR'),
                                          fasta_dir = system.file('extdata', package = 'MSfmtR'),
                                          peptide_summary = 'none'))

  data <- process_raw(config, save_intermediate = FALSE)

  expect_length(data, 4)
  expect_equal(dim(data$FeatureLevelData), c(372,16))
  expect_equal(dim(data$ProteinLevelData), c(6, 11))


  peptides <- process_peptides(data, config, save_intermediate = FALSE)

  expect_equal(dim(peptides), c(62, 21))
  expect_equal(peptides$FEATURE[1], "_ASNEDVTK_.2_y4_1")
  expect_equal(as.character(peptides$TRANSITION[1]), "y4_1")


  proteins <- process_proteins(data, peptides, config, save_intermediate = FALSE)

  expect_equal(dim(proteins), c(1, 19))
})


##### with peptide roll up #####
test_that("Specronaut workflow with peptide roll up", {

  config <- configure_formatR(args = list(input_dir = system.file('extdata', package = 'MSfmtR'),
                                          fasta_dir = system.file('extdata', package = 'MSfmtR'),
                                          peptide_summary = 'PEP'))

  data <- process_raw(config, save_intermediate = FALSE)

  expect_length(data, 4)
  expect_equal(dim(data$FeatureLevelData), c(108,16))
  expect_equal(dim(data$ProteinLevelData), c(6, 11))


  peptides <- process_peptides(data, config, save_intermediate = FALSE)

  expect_equal(dim(peptides), c(17, 21))
  expect_equal(peptides$FEATURE[1], "_ASNEDVTK_")
  expect_equal(peptides$TRANSITION[1], "y6_1,y7_1,y4_1,y7_2")


  proteins <- process_proteins(data, peptides, config, save_intermediate = FALSE)

  expect_equal(dim(proteins), c(1, 19))


  wb <- process_wb(proteins, peptides, config, save_intermediate = FALSE)

  expect_s4_class(wb, "Workbook")
  expect_equal(wb$sheet_names, "O75475")
  expect_length(wb$styleObjects, 6)
})
