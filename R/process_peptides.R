#' process_peptides
#' Process peptides from raw data
#'
#' @param data MSstats formatted data
#' @param config list of configuration parameters
#' @param merge_method character value defining the method to use for merging peptide data. Valid options are 'median', 'mean', 'lmer'.
#' @param stage character path to checkpoint file
#' @param save_intermediate logical save intermediate data
#'
#' @details This function processes MSstats formatted data and returns peptides ready for ProtResDash. If save_intermediate is TRUE, the processed data are also saved to the checkpoint file.
#'
#' If `merge_method` is 'median' (default), the function will use the median within `data$FeatureLevelData$GROUP` (or `config$groups` when defined) for group-level peptide abundance statistics.
#'
#' If `merge_method` is 'mean', the function will use the mean within `data$FratureLevelData$GROUP` (or `config$groups` when defined) for group-level peptide abundance statistics.
#'
#' If `merge_method` is 'lmer', the function will use a linear mixed effects model to estimate group-level peptide abundance statistics within `data$FratureLevelData$GROUP` (or `config$groups` when defined) and using random effects for `config$samples` (e.g. for technical replicates).
#'
#' @return data frame of processed peptide data
#' @export
#' @importFrom dplyr group_by summarize ungroup arrange left_join mutate
#' @importFrom purrr map_dfr map_lgl set_names
#' @importFrom lme4 lmer fixef
#' @importFrom stats median sd vcov
#' @importFrom stringr fixed str_replace_all
#' @importFrom tidyr pivot_wider
process_peptides <- function(data, config, merge_method = 'median',
                             stage = file.path(config$output_dir, config$peptide_checkpoint),
                             save_intermediate = TRUE)
{
  # for those pesky no visible binding warnings
  if(FALSE)
    PROTEIN <- FEATURE <- GROUP <- INTENSITY <- ABUNDANCE <- id <- cv <- originalRUN <- SUBJECT <-
      PEPTIDE <- TRANSITION <- Modification <- group <- model <- NULL

  # add sample/group information if provided
  data$FeatureLevelData <- map_samples(data$FeatureLevelData, config) |>
    
    map_groups(config)
  

  # calculate group abundance statistics for peptides
  if(merge_method == 'median'){

    peptides_long <- data$FeatureLevelData |>

      group_by(PROTEIN, FEATURE, group) |>

      summarize(INTENSITY = median(INTENSITY, na.rm = TRUE),
                cv = sd(ABUNDANCE, na.rm = TRUE) / mean(ABUNDANCE, na.rm = TRUE) * 100) |>

      ungroup()

  }else if(merge_method == 'mean'){

    peptides_long <- data$FeatureLevelData |>

      group_by(PROTEIN, FEATURE, group) |>

      summarize(INTENSITY = log(INTENSITY) |> mean(na.rm = TRUE) |> exp(),
                cv = sd(ABUNDANCE, na.rm = TRUE) / mean(ABUNDANCE, na.rm = TRUE) * 100) |>

      ungroup()

  }else if(merge_method == 'lmer'){
    peptides_long <- data$FeatureLevelData |>

      group_by(PROTEIN, FEATURE, group) |>

      summarize(models = list(lmer(log(INTENSITY) ~ 1 + (1 | sample),
                                   control = lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))), # we use lmerControl to silence singular fit warnings. We expect to see a few of these when there are few technical replicates (e.g. if there are 3 replicates). Since we are only concerned with the fixed effects estimates, we should be fine even with a singular fit.
                INTENSITY = map_dbl(models, ~ fixef(.x) |> exp()),
                cv = map_dbl(models, ~ 
                             {
                               retval <- try((sqrt(vcov(.x)) / fixef(.x) * 100)[1,1])
                               if(class(retval) == 'try-error')
                                 return(as.double(NA))
                               
                               return(as.double(retval))
                             })
                ) |>
      select(-models) |>

      ungroup()
  }else{
    stop('Invalid merge_method')
  }


  # summary data
  peptides <- peptides_long |>
    dplyr::rename(id = group) |>

    arrange(id) |> # sort by id - this keeps column names in the same order as in `proteins`

    pivot_wider(names_from = id, values_from = c(INTENSITY, cv))

  # update names to match protein data
  names(peptides) <- names(peptides) |>
    str_replace_all(pattern = fixed('INTENSITY_'), replacement = 'Group Abundance: ') |>
    str_replace_all(pattern = fixed('cv_'), replacement = 'CV: ')


  # extract peptide data
  peptides <- data$FeatureLevelData |>

    mutate(id = paste0('Abundance: ', GROUP, ', ', originalRUN, ' (', SUBJECT, ')'),

           PEPTIDE = map_chr(PEPTIDE, ~ strsplit(as.character(.x),
                                                 split = '_', fixed = TRUE)[[1]][1] |>
                               gsub(pattern = '\\[.*?\\]', replacement = '')),

           # split modifications out into a different column
           Modification = map_chr(FEATURE, ~
                                    {
                                      openMod <- str_locate_all(.x, fixed('[')) |>
                                        unlist() |>
                                        as.vector() |>
                                        unique()
                                      closeMod <- str_locate_all(.x, fixed(']')) |>
                                        unlist() |>
                                        as.vector() |>
                                        unique()

                                      if(length(openMod) > 0)
                                      {
                                        mod <- substr(.x, openMod[1], closeMod[1]) |>
                                          paste(collapse = ';')
                                      }else{
                                        mod <- ''
                                      }

                                      return(mod)
                                    })) |>


    dplyr::select(PROTEIN, PEPTIDE, TRANSITION, FEATURE, Modification, id, INTENSITY) |>

    arrange(id) |> # sort by id - this keeps column names in the same order as in `proteins`

    unique() |> # remove duplicates - these do appear rarely in the data

    pivot_wider(names_from = id, values_from = INTENSITY) |>

    # merge stats into peptides
    left_join(peptides, by = c('PROTEIN', 'FEATURE'))


  # # some diagnostic plots for looking at peptides_long
  # library(ggplot2)
  # library(cowplot)
  # theme_set(theme_cowplot())
  #
  # peptides_long |>
  #   ggplot(aes(x = cv, y = INTENSITY)) +
  #   geom_point() +
  #   facet_wrap(~GROUP)

  # checkpoint
  if(save_intermediate)
    save(peptides, file = stage)

  return(peptides)
}
