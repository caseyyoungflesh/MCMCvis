#' Extract posterior chains from MCMC output
#'
#' Extract posterior chains from MCMC output for specific parameters of interest.
#'
#' @param object Object containing MCMC output. See DETAILS below.
#' @param par Character string (or vector of character strings) denoting parameters of interest.
#' Partial names may be used to return all parameters containing that set of characters.
#'
#' Default \code{'all'} returns chains for all parameters.

#' @section Details:
#' Function returns matrix with one chain per column for specified parameters. Multiple input chains for each
#' parameter are combined to one posterior chain.
#'
#' \code{object} argument can be a \code{stanfit} object (\code{rstan} package), an \code{mcmc.list} object
#' (\code{coda} package), an \code{R2jags} model object (\code{R2jags} package), or a matrix containing MCMC
#' chains (each column representing MCMC output for a single parameter, rows representing iterations in the chain).
#'
#'
#' @examples
#' #Load data
#' data(MCMC_data)
#'
#' #Extract MCMC chains
#' ex <- pochains(MCMC_data)
#' apply(ex, 2, mean)
#'
#' #Extract MCMC chains for just 'beta' parameters
#' ex2 <- pochains(MCMC_data, par='beta')
#' apply(ex2, 2, mean)
#'
#' @export


pochains <- function(object,
                     par = 'all')
{

    if(typeof(object) == 'S4')
    {
      temp_in <- rstan::As.mcmc.list(object)
      names <- colnames(temp_in[[1]])

      temp <- do.call('rbind', temp_in)
    }

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


  if(coda::is.mcmc.list(object) != TRUE &
     typeof(object) != 'double' &
     typeof(object) != 'list' &
     typeof(object) != 'S4')
  {
    stop('Invalid object type. Input must be stanfit object (rstan), mcmc.list object (coda),
         rjags object (R2jags), or matrix with MCMC chains.')
  }

  if (length(par) == 1)
  {
    if (par == 'all')
    {
      OUT <- temp
    }else
    {
      get.cols <- grep(paste(par), names, fixed=TRUE)
      if (length(get.cols) < 1)
      {
        stop(paste0('"', par, '"', ' not found in MCMC ouput.'))
      }
      OUT <- temp[,get.cols]
    }
  }else
  {
    grouped <- c()
    for (i in 1:length(par))
    {
      #i <- 1
      get.cols <- grep(paste(par[i]), names, fixed=TRUE)
      if (length(get.cols) < 1)
      {
        stop(paste0('"', par[i], '"', ' not found in MCMC ouput.'))
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

  return(OUT)
}

