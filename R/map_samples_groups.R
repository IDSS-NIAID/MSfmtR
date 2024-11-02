# map_samples_groups.R
# functions to map data$FeatureLevelData$originalRUN to samples and groups

#' map_samples
#' Map originalRUN IDs to sample IDs
#' 
#' @param data A data.frame from an MSstats object - either `MSstats_obj$FeatureLevelData` or `MSstats_obj$ProteinLevelData`
#' @param config A list containing the configuration parameters
#' 
#' @return An update to `data` with modified sample IDs. For sample mapping, a new variable, `sample` is added, and for group mapping a new variable `group` is added (not overwritting the variable `GROUP` that should be included in `data`).
#' @export
map_samples <- function(data, config)
{
  if(!is.null(config$samples))
  {
    # map originalRUN to config$samples
    sample_from <- unique(data$originalRUN) |> as.character()
    
    # create a mapping from originalRUN to sampleID
    sample_map_dfr <- map_dfr(config$samples, ~ data.frame(sample_from = grep(.x, sample_from, value = TRUE),
                                                           sample_to   = .x))
    
    # convert to named vector
    sample_map <- sample_map_dfr$sample_to |>
      set_names(sample_map_dfr$sample_from)
    
    # add sample information to data
    return(mutate(data, sample = sample_map[as.character(originalRUN)]))
  }
  
  if(config$merge_method == 'lmer')
      stop('merge_method is "lmer" but config$samples is not defined')
  
    
  return(mutate(data, sample = originalRUN))
}

#' map_groups
#' Map originalRUN IDs to group IDs
#' 
#' @rdname map_samples
#' @export
map_groups <- function(data, config)
{
  if(!is.null(config$groups))
  {
    # map originalRUN to config$groups
    group_from <- unique(data$originalRUN) |> as.character()
    
    # create a mapping from originalRUN to groupID
    group_map_dfr <- map_dfr(config$groups, ~ data.frame(group_from = grep(.x, group_from, value = TRUE),
                                                         group_to   = .x))
    
    # convert to named vector
    group_map <- group_map_dfr$group_to |>
      set_names(group_map_dfr$group_from)
    
    return(mutate(data, group = group_map[as.character(originalRUN)]))
  }
  
  return(mutate(data, group = GROUP))
}