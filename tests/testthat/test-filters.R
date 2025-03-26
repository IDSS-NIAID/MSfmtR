test_that("check_missingness works", {

  # test check_missingness function on a data set with three proteins
  # there are 4 bio-replicates with abundance values,
  # and 4 other columns with character values indicating if the protein was found
  set.seed(294873)
  dat <- tibble(pep = c(1, 2, 3),
                `Abundance: BR1` = runif(3) * 1e7,

                # middle protein abundance is missing and last value is too small
                `Abundance: BR2` = runif(3) * c(1e7, NA, 1e4),

                # middle protein abundance is missing
                `Abundance: BR3` = runif(3) * c(1e7, NA, 1e7),

                # last protein abundance is too small
                `Abundance: BR4` = runif(3) * c(1e7, 1e7, 1e4),

                `Found in Sample: BR1` = rep('High', 3),

                # middle value is missing and last value is too small
                `Found in Sample: BR2` = c('High', NA, 'Low'),

                # middle value is missing
                `Found in Sample: BR3` = c('High', NA, 'High'),

                # last value is too small
                `Found in Sample: BR4` = c('High', 'Low', 'Low'))


  ### numeric columns ###
  # all rows have at least two good values, so no rows are flagged to drop
  expect_equal(check_missingness(dat, 2:5,
                                 prop_good = 0.5,
                                 llod = 1e6),
               c(FALSE, FALSE, FALSE))

  # the second and third rows have less than 60% good values, so they are flagged to drop
  expect_equal(check_missingness(dat, 2:5,
                                 prop_good = 0.6,
                                 llod = 1e6),
               c(FALSE, TRUE, TRUE))


  ### character columns ###
  # middle row has three bad values
  expect_equal(check_missingness(dat, 6:9,
                                 prop_good = 0.5,
                                 valid_chr = 'High'),
               c(FALSE, TRUE, FALSE))

  # the second and third rows have less than 60% good values, so they are flagged to drop
  expect_equal(check_missingness(dat, 6:9,
                                 prop_good = 0.6,
                                 valid_chr = 'High'),
               c(FALSE, TRUE, TRUE))

  # test with no provided valid values - should use all non-missing values
  expect_equal(check_missingness(dat, 6:9,
                                 prop_good = 0.5),
               c(FALSE, FALSE, FALSE))

  # only the second rows has missing data, so it is flagged to drop
  expect_equal(check_missingness(dat, 6:9,
                                 prop_good = 0.6),
               c(FALSE, TRUE, FALSE))


  ### mixed columns ###
  # all rows pass numeric column check, but not the character column check, so it is flagged
  expect_equal(check_missingness(dat, 2:9,
                                 prop_good = 0.5,
                                 llod = 1e6,
                                 valid_chr = 'High'),
               c(FALSE, TRUE, FALSE))

  # the second and third rows have less than 60% good values, so they are flagged to drop
  expect_equal(check_missingness(dat, 2:9,
                                 prop_good = 0.6,
                                 llod = 1e6,
                                 valid_chr = 'High'),
               c(FALSE, TRUE, TRUE))
})

test_that("filter_completeness works", {

  # test filter_completeness function on a data set with three proteins
  # there are two conditions with grouped abundance,
  # and three other columns with ratios, log ratios, and p-values
  set.seed(294873)
  dat <- tibble(pep = c(1, 2, 3),
                `Group Abundance: trt A` = runif(3) * c(1e7, 1e4, 1e7), # one low abundance in A
                `Group Abundance: trt B` = runif(3) * c(1e7,   0, 1e4), # one not observed in treatment B and one too small

                `Abundance Ratio: (trt A) / (trt B)` = `Group Abundance: trt A` / `Group Abundance: trt B`,
                `Abundance Ratio (log2): (trt A) / (trt B)` = log2(`Abundance Ratio: (trt A) / (trt B)`),
                `Abundance Ratio Adj p-value: (trt A) / (trt B)` = runif(3),

                `drop_trt A` = c(FALSE, TRUE, FALSE),
                `drop_trt B` = c(FALSE, TRUE, TRUE))

  # filter rows
  dat_filtered <- dat |>
    filter_completeness(cols = tidyselect::starts_with("Group Abundance", vars = names(dat)),
                        filter_cols = tidyselect::starts_with("drop_", vars = names(dat)),
                        ratio = FALSE) |>

    filter_completeness(cols = tidyselect::starts_with("Abundance Ratio:", vars = names(dat)),
                        filter_cols = tidyselect::starts_with("drop_", vars = names(dat)),
                        ratio = TRUE, max_ratio = 100) |>

    filter_completeness(cols = tidyselect::starts_with("Abundance Ratio (log2):", vars = names(dat)),
                        filter_cols = tidyselect::starts_with("drop_", vars = names(dat)),
                        ratio = TRUE, max_ratio = log2(100)) |>

    filter_completeness(cols = tidyselect::starts_with("Abundance Ratio Adj p-value:", vars = names(dat)),
                        filter_cols = tidyselect::starts_with("drop_", vars = names(dat)),
                        ratio = TRUE)

  # check output
  expect_equal(dat_filtered$`Group Abundance (Filtered): trt A` |> round(),
               c(5789401, NA, 5885076))
  expect_equal(dat_filtered$`Group Abundance (Filtered): trt B` |> round(),
               c(4309969, NA, NA))
  expect_equal(dat_filtered$`Abundance Ratio (Filtered): (trt A) / (trt B)` |> round(2),
               c(1.34, NA, 100))
  expect_equal(dat_filtered$`Abundance Ratio (log2) (Filtered): (trt A) / (trt B)` |> round(3),
               c(0.426, NA, 6.644))
  expect_equal(dat_filtered$`Abundance Ratio Adj p-value (Filtered): (trt A) / (trt B)` |> round(2),
               c(0.93, NA, 0.28))

  # make sure we get a warning if we don't provide the correct columns
  expect_warning(filter_completeness(dat, cols = tidyselect::starts_with("Not here", vars = names(dat)),
                                     filter_cols = tidyselect::starts_with("drop_", vars = names(dat))))
  expect_warning(filter_completeness(dat, cols = tidyselect::starts_with("Abundance Ratio:", vars = names(dat)),
                                     filter_cols = tidyselect::starts_with("drop_", vars = names(dat)),
                                     ratio = FALSE))
})
