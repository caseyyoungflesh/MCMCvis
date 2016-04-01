#' Extract posterior chains from MCMC output
#'
#' Extract posterior chains from MCMC output for specific parameters of interest.
#'
#' @param object Object containing MCMC output. See \code{input} argument and DETAILS below.
#' @param params Character string (or vector of character strings) denoting parameters to be
#' returned in summary output. Partial names may be used to return all parameters containing
#' that set of characters.
#'
#' Default \code{all} returns chains for all parameters.
#' @param input Indicates the nature of the \code{object} argument.
#'
#' Valid entries are \code{jags_object}, \code{mcmc_list}, and \code{chains}. See DETAILS below.
#' @section Details:
#' For \code{posummary(object, input = 'jags_object')}, input must be JAGS model object from \code{R2jags} package.
#'
#' For \code{posummary(object, input = 'mcmc_list')}, input must be of type \code{mcmc.list}.
#'
#' For \code{posummary(object, input = 'chains')}, each column of \code{object} should contain a posterior
#' chain for a single parameter. Each row represents one iteration in the chain.
#'
#' @return \code{posummary(params='all')} returns posterior chains for all parameters contained within
#' JAGS model object.
#'
#' \code{posummary(params=c('beta[1]', 'beta[2]'))} returns posterior chains for just parameters
#' \code{beta[1]} and \code{beta[2]}.
#'
#' \code{posummary(params=c('beta'))} returns posterior chains for all parameters containing \code{beta}
#'  in their name.
#'
#' @export


pochains <- function(object,
                     params = 'all')
{
    if(coda::is.mcmc.list(object) == TRUE)
    {
      temp <- object
      names <- colnames(temp[[1]])
      n_chains <- length(temp)

      ch_bind <- do.call('rbind', temp)
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
      OUT <- temp[,get.cols]
    }
  }else
  {
    grouped <- c()
    for (i in 1:length(params))
    {
      #i <- 1
      get.cols <- grep(paste(params[i]), names, fixed=TRUE)
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

