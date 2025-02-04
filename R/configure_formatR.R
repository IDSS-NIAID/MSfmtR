#' configure_formatR
#' Default configuration for MSfmtR
#'
#' @param config_file path to a yaml file with configuration settings
#' @param args list of arguments to override defaults
#' @param config_profile character, name of profile to use from `config_file`
#' @param fill_default logical, fill in default values if using a profile other than 'default'
#'
#' @return list of configuration settings
#' @export
#' @importFrom stringr fixed str_replace
configure_formatR <- function(config_file = NULL, args = NULL,
                              config_profile = 'default', fill_default = TRUE)
{
  if(!is.null(args$config_file))
  {
    config_file <- args$config_file
  }

  # read yml file (requires config package)
  if(!is.null(config_file))
  {
    config <- config::get(file = config_file,
                          config = config_profile)

    # if we are using a different profile, fill in missing values with defaults from `config_file`
    if(fill_default & config_profile != 'default')
    {
      default <- config::get(file = config_file,
                             config = 'default')
      for(key in names(default))
      {
        if(is.null(config[[key]]))
        {
          config[[key]] <- default[[key]]
        }
      }
    }
  }else{
    # use package defaults if no config file is provided
    # additional defaults will be picked up when needed if they are not in config
    config <- defaults
  }

  # replace any defaults with command line arguments
  if(length(args) > 0)
  {
    args$config <- config
    config <- do.call('updt_config', args = args)
  }

  # go through and execute any remaining expressions from defaults
  for(key in names(config))
  {
    if(is.expression(config[[key]]))
    {
      config[[key]] <- with(config, eval(config[[key]])) # execute in the context of config
    }
  }

  ## Checkpoint parameters - this is due to be deprecated in a future version
  if(!is.null(config$checkpoints))
  {
    if(any(config$checkpoints == 'all'))
    {
      config$checkpoints <- c('xlsx', 'sql', 'processed', 'protein', 'peptide', 'wb')
    }

    if(any(grepl(',', config$checkpoints)))
    {
      config$checkpoints <- str_split(config$checkpoints, fixed(','))[[1]]
    }
  }else{
    config$checkpoints <- defaults$checkpoints
  }

  return(config)
}


#' updt_config
#' Update configuration settings (primarily used for local changes to defaults)
#'
#' @param config list of configuration settings
#' @param ... named arguments to update
updt_config <- function(config, ...)
{
  updates <- list(...)
  for(key in names(updates))
  {
    # update config with the new value
    if(!is.null(updates[[key]]))
      config[[key]] <- updates[[key]]

    # if the value is NULL, set it to the default
    if(is.null(config[[key]]) & !is.null(defaults[[key]]))
    {
      if(class(defaults[[key]]) == 'expression')
      {
        config[[key]] <- with(config, eval(defaults[[key]])) # execute in the context of config
      }else{
        config[[key]] <- defaults[[key]]
      }
    }
  }

  return(config)
}
