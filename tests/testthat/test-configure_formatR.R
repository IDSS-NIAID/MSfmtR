# test configuration and package defaults

test_that("Configuration management", {

  # test with minimal configuration - should have some blanks, but that is OK
  config <- configure_formatR()

  expect_equal(names(config), names(defaults))

  # test with package data
  config <- configure_formatR(args = list(input_dir = system.file('extdata', package = 'MSfmtR'),
                                          in_file = 'spectronaut_report.tsv'))

  expect_equal(config$in_file, 'spectronaut_report.tsv')
  expect_equal(config$sheet, 'spectronaut_report')
  expect_equal(basename(config$cont_fasta), 'Universal_Contaminants.fasta')

  # test an updated configuration
  config <- updt_config(config, uloq = 1e9, lloq = 1e6)

  expect_equal(config$uloq, 1e9)
  expect_equal(config$lloq, 1e6)
})
