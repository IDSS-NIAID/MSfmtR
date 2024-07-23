
process_peptides <- function(data, config,
                             stage = file.path(config$output_dir, config$peptide_checkpoint),
                             save_intermediate = TRUE)
{
  # calculate peptide stats
  peptides_long <- data$FeatureLevelData %>%

    group_by(PROTEIN, FEATURE, GROUP) %>%

    summarize(INTENSITY = median(INTENSITY, na.rm = TRUE),
              cv = sd(ABUNDANCE, na.rm = TRUE) / mean(ABUNDANCE, na.rm = TRUE) * 100) %>%

    ungroup()


  # summary data
  peptides <- peptides_long |>
    dplyr::rename(id = GROUP) |>

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
                                      openMod <- str_locate_all(.x, fixed('[')) %>%
                                        unlist() %>%
                                        as.vector() %>%
                                        unique()
                                      closeMod <- str_locate_all(.x, fixed(']')) %>%
                                        unlist() %>%
                                        as.vector() %>%
                                        unique()

                                      if(length(openMod) > 0)
                                      {
                                        mod <- substr(.x, openMod[1], closeMod[1]) %>%
                                          paste(collapse = ';')
                                      }else{
                                        mod <- ''
                                      }

                                      return(mod)
                                    })) %>%


    dplyr::select(PROTEIN, PEPTIDE, TRANSITION, FEATURE, Modification, id, INTENSITY) %>%

    arrange(id) |> # sort by id - this keeps column names in the same order as in `proteins`

    unique() |> # remove duplicates - these do appear rarely in the data

    pivot_wider(names_from = id, values_from = INTENSITY) %>%

    # merge stats into peptides
    left_join(peptides, by = c('PROTEIN', 'FEATURE'))

  if(FALSE)
  {
    # some diagnostic plots for looking at peptides_long
    library(ggplot2)
    library(cowplot)
    theme_set(theme_cowplot())

    peptides_long %>%
      ggplot(aes(x = cv, y = INTENSITY)) +
      geom_point() +
      facet_wrap(~GROUP)
  }

  # checkpoint
  if(save_intermediate)
    save(peptides, file = stage)

  return(peptides)
}