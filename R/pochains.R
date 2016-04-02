#' Extract posterior chains from MCMC output
#'
#' Extract posterior chains from MCMC output for specific parameters of interest.
#'
#' @param object Object containing MCMC output. See DETAILS below.
#' @param params Character string (or vector of character strings) denoting parameters of interest.
#' Partial names may be used to return all parameters containing that set of characters.
#'
#' Default \code{all} returns chains for all parameters.

#' @section Details:
#' \code{object} argument can be an \code{mcmc.list} object, an \code{R2jags} model object (output from the \code{R2jags}
#' package), or a matrix containing MCMC chains (each column representing MCMC output for a single parameter, rows
#' representing iterations in the chain).
#'
#' @return \code{pochains(params='all')} returns chains for all parameters.
#'
#' \code{pochains(params=c('beta[1]', 'beta[2]'))} returns chains for just parameters
#' \code{beta[1]} and \code{beta[2]}.
#'
#' \code{pochains(params=c('beta'))} returns chains for all parameters containing \code{beta}
#'  in their name.
#'
#' @export


pochains <- function(object,
                     params = 'all')
{
    if(coda::is.mcmc.list(object) == TRUE)
    {
      temp_in <- object
      names <- colnames(temp_in[[1]])

      temp <- do.call('rbind', temp_in)
    }

    if(typeof(object) == 'double')
    {
      temp <- object
      names <- colnames(temp)
    }

    if(typeof(object) == 'list' & coda::is.mcmc.list(object) == FALSE)
    {
      temp <- object$BUGSoutput$sims.matrix
      names <- colnames(temp)
    }

  if (length(params) == 1)
  {
    if (params == 'all')
    {
      OUT <- temp
    }else
    {
      get.cols <- grep(paste(params), names, fixed=TRUE)
      if (length(get.cols) < 1)
      {
        stop(paste0('"', params, '"', ' not found in MCMC ouput.'))
      }
      OUT <- temp[,get.cols]
    }
  }else
  {
    grouped <- c()
    for (i in 1:length(params))
    {
      #i <- 1
      get.cols <- grep(paste(params[i]), names, fixed=TRUE)
      if (length(get.cols) < 1)
      {
        stop(paste0('"', params[i], '"', ' not found in MCMC ouput.'))
      }
      grouped <- c(grouped, get.cols)
    }

    to.rm <- which(duplicated(grouped))
    if(length(to.rm) >0)
    {
      g_filt <- grouped[-to.rm]
    }else
    {
      g_filt <- grouped
    }

    OUT <- temp[,g_filt]
  }
}

