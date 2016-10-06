#' Extract posterior chains from MCMC output
#'
#' Extract posterior chains from MCMC output for specific parameters of interest.
#'
#' @param object Object containing MCMC output. See DETAILS below.
#' @param params Character string (or vector of character strings) denoting parameters of interest.
#' Partial names may be used to return all parameters containing that set of characters.
#'
#' Default \code{'all'} returns chains for all parameters.
#'
#' @param excl Character string (or vector of character strings) denoting parameters to exclude.
#' Partical names may be used to exclude all parameters contaiing that set of characters. Used in
#' conjunction with \code{params} argument to select parameters of interest.
#'
#' @section Details:
#' Function returns matrix with one chain per column for specified parameters. Multiple input chains for each
#' parameter are combined to one posterior chain.
#'
#' \code{object} argument can be a \code{stanfit} object (\code{rstan} package), an \code{mcmc.list} object
#' (\code{coda} package), an \code{R2jags} model object (\code{R2jags} package), or a matrix containing MCMC
#' chains (each column representing MCMC output for a single parameter, rows representing iterations in the chain).
#' The function automatically detects the object type and proceeds accordingly.
#'
#'
#' @examples
#' #Load data
#' data(MCMC_data)
#'
#' #Extract MCMC chains
#' ex <- MCMCchains(MCMC_data)
#' apply(ex, 2, mean)
#'
#' #Extract MCMC chains for just 'beta' parameters
#' ex2 <- MCMCchains(MCMC_data, params='beta')
#' apply(ex2, 2, mean)
#'
#' @export

MCMCchains <- function(object,
                     params = 'all',
                     excl = NULL)
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

  if(!is.null(excl))
  {
    to.rm1 <- c()
    for (i in 1:length(excl))
    {
      to.rm1 <- c(to.rm1, grep(excl[i], names, fixed = TRUE))
    }
    dups <- -which(duplicated(to.rm1))
    if(length(dups) > 0)
    {
      to.rm2 <- to.rm1[-dups]
    }else{
      to.rm2 <- to.rm1
    }
  }

  if (length(params) == 1)
  {
    if (params == 'all')
    {
      if(is.null(excl))
      {
        OUT <- temp
      }else{
        OUT <- temp[,-to.rm2]
      }
    }else
    {
      get.cols <- grep(paste(params), names, fixed=TRUE)
      if (length(get.cols) < 1)
      {
        stop(paste0('"', params, '"', ' not found in MCMC ouput.'))
      }

      if(!is.null(excl))
      {
        if(identical(get.cols, to.rm2))
        {
          stop('No parameters selected.')
        }

        matched <- which(get.cols == to.rm2)
        if (length(matched) > 0)
        {
          cols <- get.cols[-matched]
        }else {
          cols <- get.cols
        }

      }else{
        cols <- get.cols
      }

      OUT <- temp[,cols]
    }

  }else
  {
    grouped <- c()
    for (i in 1:length(params))
    {
      get.cols <- grep(paste(params[i]), names, fixed=TRUE)
      if (length(get.cols) < 1)
      {
        stop(paste0('"', params[i], '"', ' not found in MCMC ouput.'))
      }
      grouped <- c(grouped, get.cols)
    }

    if(!is.null(excl))
    {
      if(identical(grouped, to.rm2))
      {
        stop('No parameters selected.')
      }

      matched <- stats::na.omit(match(to.rm2, grouped))
      if (length(matched) > 0)
      {
        cols <- grouped[-matched]
      } else{
        cols <- grouped
      }

      to.rm <- which(duplicated(cols))
      if(length(to.rm) > 0)
      {
        g_filt <- cols[-to.rm]
      }else
      {
        g_filt <- cols
      }
    } else{

      to.rm <- which(duplicated(grouped))
      if(length(to.rm) > 0)
      {
        g_filt <- grouped[-to.rm]
      }else
      {
        g_filt <- grouped
      }
    }

    OUT <- temp[,g_filt]
  }

  return(OUT)
}

