# test raw data processing

test_that("Raw data processing", {

  raw <- read_delim(system.file('extdata/spectronaut_report.tsv', package = 'MSfmtR'),
                    delim = '\t')

  ##### with no peptide roll up #####
  data <- configure_formatR(args = list(input_dir = system.file('extdata', package = 'MSfmtR'),
                                        in_file = 'spectronaut_report.tsv',
                                        peptide_summary = 'none')) |>
    process_raw()

  # should have one row per fragment from raw in FeatureLevelData
  expect_equal(nrow(data$FeatureLevelData), nrow(raw))

  # and one row per protein and replicate in ProteinLevelData
  expect_equal(nrow(data$ProteinLevelData),
               select(raw, R.FileName, PG.ProteinAccessions) |>
                 distinct() |>
                 nrow())


  ##### with peptide roll up to fragment group level #####

  data <- configure_formatR(args = list(input_dir = system.file('extdata', package = 'MSfmtR'),
                                        in_file = 'spectronaut_report.tsv',
                                        peptide_summary = 'FG')) |>
    process_raw()

  # should have one row per fragment group from raw in FeatureLevelData
  expect_equal(nrow(data$FeatureLevelData),
               select(raw, R.FileName, PG.ProteinAccessions, EG.PrecursorId) |>
                 distinct() |>
                 nrow())

  # and one row per protein and replicate in ProteinLevelData
  expect_equal(nrow(data$ProteinLevelData),
               select(raw, R.FileName, PG.ProteinAccessions) |>
                 distinct() |>
                 nrow())


  ##### with peptide roll up to peptide level #####

  data <- configure_formatR(args = list(input_dir = system.file('extdata', package = 'MSfmtR'),
                                        in_file = 'spectronaut_report.tsv',
                                        peptide_summary = 'PEP')) |>
    process_raw()

  # should have one row per peptide from raw in FeatureLevelData
  expect_equal(nrow(data$FeatureLevelData),
               select(raw, R.FileName, PG.ProteinAccessions, EG.ModifiedSequence) |>
                 distinct() |>
                 nrow())
})
