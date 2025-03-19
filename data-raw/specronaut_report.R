# generation of sample data from a simulated MS experiment as output by Spectronaut

library(dplyr)
library(stringr)
library(readr)
devtools::load_all()

# A tibble with sample data from a simulated MS experiment as output by Spectronaut
# There will be (total of 2 * 5 * 1 * 3 * 5 * 2 * 3 = 900) rows:
#    2 conditions/grouping factors,
#    5 replicates for each condition,
#    1 MS file per replicate,
#    3 proteins,
#    5 peptides per protein
#    2 fragment groups
#    3 fragments per fragment group

# Columns:
#    # replicate columns
#    R.Condition = col_character(), grouping factor, levels: "treatment", "control"
#    R.FileName = col_character(), MS file name, one for each replicate
#    R.Replicate = col_double(), replicate number, 1-5
#
#    # protein group columns
#    PG.ProteinAccessions = col_character(), protein accession number
#    PG.ProteinGroups = col_character(), protein group (leave this the same as ProteinAccessions for now)
#    PG.Qvalue = col_double(), q-value for protein identification: `0.01 * runif()`
#    PG.Quantity = col_double(), protein group quantity (`rnorm` with mean 6 and sd 2, then exponentiated)
#
#    # peptide columns
#    PEP.GroupingKey = col_character(), peptide grouping key (same as peptide sequence for now)
#    PEP.StrippedSequence = col_character(), peptide sequence
#    PEP.Quantity = col_double(), peptide quantity (generate normal distribution around `log10(protein)` quantity)
#
#    # experiment group columns
#    EG.iRTPredicted = col_double(), Indexed retion time predicted by Spectronaut (`rgamma` with shape 4 and scale 10)
#    EG.Library = col_character(), library name ("directDIA" for all of them)
#    EG.ModifiedSequence = col_character(), modified peptide sequence (needs to be the same for all conditions/replicates for a specific peptide)
#    EG.PrecursorId = col_character(), precursor ID (pick a random '.2' or '.3' to add on the end of each sequence modification for each peptide - should be the same for all conditions/replicates for a specific peptide - this is the charge)
#    EG.Qvalue = col_double(), q-value for peptide identification: `0.01 * runif()`
#
#    # fragment group columns
#    FG.Charge = col_double(), fragment charge (this is the number added to the end of the precursor ID)
#    FG.Id = col_character(), this is equal to the precursor ID
#    FG.PrecMz = col_double(), precursor m/z (`rnorm` with mean = 600 and sd = 100)
#    FG.Quantity = col_double(), fragment quantity (just use the peptide quantity for now)
#
#    # fragment columns
#    F.Charge = col_double(), fragment charge (mostly 1, but some 2 and a few 3's)
#    F.FrgIon = col_character(), fragment ion (`paste0(F.FrgType, F.FrgNum)`)
#    F.FrgLossType = col_character(), fragment loss type (call them all "noloss" for now)
#    F.FrgMz = col_double(), fragment m/z (`exp(rnorm())` with mean = 6.5 and sd = 0.35)
#    F.FrgNum = col_double(), fragment number (1-5)
#    F.FrgType = col_character(), fragment type (b, y, or p - mostly y, but some b and a few p's)
#    F.ExcludedFromQuantification = col_logical(), excluded from quantification (all FALSE)
#    F.NormalizedPeakArea = col_double(), normalized peak area (`rnorm` with mean = 8 and sd = 1.5)
#    F.NormalizedPeakHeight = col_double(), normalized peak height (`rnorm` with mean = 11 and sd = 2)
#    F.PeakArea = col_double(), peak area (F.NormalizedPeakArea + rnorm with mean of 0 and sd of 0.3)
#    F.PeakHeight = col_double(), peak height (F.NormalizedPeakHeight + rnorm with mean of 0 and sd of 0.3)

n <- 900

set.seed(123)
spectronaut_report <- tibble(
  R.Condition = c("treatment", "control") |>                 # 2 conditions
    rep(each = 5 * 2 * 3) |>                                 # 5 replicates per treatment, 2 * 3 fragments
    rep(5)        |>                                         # 5 peptides per protein
    rep(3),                                                  # 3 proteins
  R.FileName = paste0('file', 1:10) |>                       # 5 replicates, 2 conditions
    rep(each = 2 * 3) |>                                     # 2 * 3 fragments
    rep(3) |>                                                # 3 proteins
    rep(5),                                                  # 5 peptides per protein
  R.Replicate = (as.numeric(str_replace(R.FileName, 'file', '')) + 4) %% 5 + 1,

  PG.ProteinAccessions = c("P01308", "P68871", "P51681") |>  # 3 proteins
    rep(each = 5) |>                                         # 5 peptides per protein
    rep(each = 2) |>                                         # 2 conditions
    rep(each = 5 * 2 * 3),                                   # 5 replicates, 2 * 3 fragments
  PG.ProteinGroups = PG.ProteinAccessions,
  PG.Qvalue = NA,
  PG.Quantity = NA,

  PEP.GroupingKey = c('RLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTP', 'RLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTP', 'KTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQ', 'KRGIVEQCCTSICSLYQLENYCN', 'MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGE',
                      'MVHLTPEE', 'KSAVTALWGKVNVDEVGGEALG', 'KVLGAFSDGLAHLDNLKGTFATLSELHCD', 'KLHVDPENFRLLGNVLVCVLAHHFG', 'KEFTPPVQAAYQKVVAGVANALAH',
                      'RLLPPLYSLVFIFGFVGNMLVILILINC', 'KSMTDIYLLNLAISDLFFLLTVPFWAHYAAAQWDFGNTMCQLLTGLYFIGFFSGIFFIILLTID', 'RLDQAMQVTETLGMTHCCINPIIYAFVGE', 'KINVKQIAA', 'RLDQAMQVTETLGMTHCCINPIIYAFVGE') |>
    rep(each = 5 * 2 * 3) |>                                 # 5 replicates, 2 * 3 fragments
    rep(each = 2),                                           # 2 conditions
  PEP.StrippedSequence = PEP.GroupingKey,
  PEP.Quantity = NA,

  EG.iRTPredicted = NA,
  EG.Library = "directDIA™",
  EG.ModifiedSequence = c('_RLLPLLALLALWGPDPAAAFVNQHLC[Carbamidomethyl (C)]GSHLVEALYLVCGERGFFYTP_', '_RLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTP_', '_KTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQ_', '_KRGIVEQCCTSICSLYQLENYCN_', '_MALWM[Oxidation (M)]RLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGE_',
                          '_MVHLTPEE_', '_KSAVTALWGKVNVDEVGGEALG_', '_KVLGAFSDGLAHLDNLKGTFATLSELHCD_', '_KLHVDPENFRLLGNVLVCVLAHHFG_', '_KEFTPPVQAAYQKVVAGVANALAH_',
                          '_RLLPPLYSLVFIFGFVGNMLVILILINC_', '_KSMTDIYLLNLAISDLFFLLTVPFWAHYAAAQWDFGNTMCQLLTGLYFIGFFSGIFFIILLTID_', '_RLDQAMQVTETLGM[Oxidation (M)]THCCINPIIYAFVGE_', '_KINVKQIAA_', '_RLDQAM[Oxidation (M)]QVTETLGMTHCCINPIIYAFVGE_') |>
    rep(each = 5 * 2* 3) |>                                  # 5 replicates, 2 * 3 fragments
    rep(each = 2),                                           # 2 conditions
  EG.PrecursorId = paste0(EG.ModifiedSequence, sample(c(".2", ".3"), n, replace = TRUE, prob = c(1,2))),

  FG.Charge = sample(1:3, n, replace = TRUE, prob = 3:1),
  FG.Id = EG.PrecursorId,
  FG.PrecMz = NA,
  FG.Quantity = NA,

  F.Charge = sample(1:3, n, replace = TRUE, prob = 3:1),
  F.FrgIon = paste0(sample(c('y', 'b', 'p'), n, replace = TRUE, prob = 3:1),
                    sample(1:5, n, replace = TRUE, c(3,3,2,2,1))),
  F.FrgLossType = rep("noloss", n)) |>

  unique() |>


  # Fragment quantities

  mutate(F.FrgMz = exp(rnorm(n(), mean = 6.5, sd = 0.35)),
         F.FrgNum = substr(F.FrgIon, 2, 2) |> as.integer(),
         F.FrgType = substr(F.FrgIon, 1, 1),
         F.ExcludedFromQuantification = rep(FALSE, n()),
         F.NormalizedPeakArea = 10^(rnorm(n(), mean = 8, sd = 1.5)),
         F.NormalizedPeakHeight = 10^(log10(F.NormalizedPeakArea) + rnorm(n(), mean = 3, sd = 2)),
         F.PeakArea = exp(log(F.NormalizedPeakArea) + rnorm(n(), mean = 0, sd = 0.3)),
         F.PeakHeight = exp(log(F.NormalizedPeakHeight) + rnorm(n(), mean = 0, sd = 0.3))) |>


  # Fragment Group quantities

  group_by(R.FileName, PG.ProteinAccessions, PEP.GroupingKey, EG.ModifiedSequence, EG.PrecursorId) |>
  mutate(FG.Quantity = tukey_biweight_mean(F.NormalizedPeakHeight),
         FG.PrecMz = rnorm(1, mean = 600, sd = 100),
         EG.iRTPredicted = rgamma(1, shape = 4, scale = 10)) |>
  ungroup() |>


  # Peptide quantities

  group_by(R.FileName, PG.ProteinAccessions, PEP.GroupingKey, EG.ModifiedSequence) |>
  mutate(PEP.Quantity = median(FG.Quantity), # this isn't how Spectronaut does it, but it's close enough for made up data
         EG.Qvalue = 0.01 * runif(1)) |>
  ungroup() |>


  # Protein Group quantities

  group_by(R.FileName, PG.ProteinAccessions) |>
  mutate(PG.Quantity = median(PEP.Quantity), # this isn't how Spectronaut does it, but it's close enough for made up data
         PG.Qvalue = 0.01 * runif(n())) |>
  ungroup()


# repository root
root <- here::here()

# save to ext data
write_tsv(spectronaut_report, file.path(root, "inst", "extdata", "spectronaut_report.tsv"))
