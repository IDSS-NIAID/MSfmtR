# test raw data processing

test_that("Raw data processing", {

  raw <- read_delim(system.file('extdata/O75475.tsv', package = 'MSfmtR'),
                    delim = '\t')

  ##### with no peptide roll up #####
  data <- process_raw(input_dir = system.file('extdata', package = 'MSfmtR'),
                      in_file = 'O75475.tsv',
                      peptide_summary = 'none',
                      save_intermediate = FALSE)

  # should have one row per fragment from raw in FeatureLevelData (that aren't excluded from quantification)
  expect_equal(nrow(data$FeatureLevelData), sum(!raw$F.ExcludedFromQuantification))

  # and one row per protein and replicate in ProteinLevelData
  expect_equal(nrow(data$ProteinLevelData),
               raw |>
                 select(R.FileName, PG.ProteinAccessions) |>
                 distinct() |>
                 nrow())


  ##### with peptide roll up to fragment group level #####

  data <- process_raw(input_dir = system.file('extdata', package = 'MSfmtR'),
                      in_file = 'O75475.tsv',
                      peptide_summary = 'FG',
                      save_intermediate = FALSE)

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

  data <- process_raw(input_dir = system.file('extdata', package = 'MSfmtR'),
                      in_file = 'O75475.tsv',
                      peptide_summary = 'PEP',
                      save_intermediate = FALSE)

  # should have one row per peptide from raw in FeatureLevelData
  # expect_equal(nrow(data$FeatureLevelData),
  #              raw |>
  #                select(R.FileName, PG.ProteinAccessions, EG.ModifiedSequence) |>
  #                distinct() |>
  #                nrow())
})
