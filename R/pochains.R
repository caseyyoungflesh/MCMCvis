#' Chain extraction for jags object
#'
#' Extracts posterior MCMC chains from JAGS object output by JAGS model fit in R
#'
#' @param object JAGS model object - output from JAGS model
#' @param params Character string (or vector of character strings) denoting parameters to be
#' returned in summary output. Partial names may be used to return all parameters containing
#' that set of characters. Default \code{'all'} returns all parameters in summary output.
#' @section Details:
#' Function currently only supports JAGS model objects output by \code{R2jags} package.
#'
#' @return \code{pochains(params='all')} returns MCMC chains for all parameters contained within
#' JAGS model object.
#'
#' \code{pochains(params=c('beta[1]', 'beta[2]'))} returns MCMC chains for just parameters
#' \code{'beta[1]'} and \code{'beta[2]'}.
#'
#' \code{pochains(params=c('beta'))} returns MCMC chains for all parameters containing \code{'beta'}
#'  in their name.
#'
#' @export


pochains <- function(object, params = 'all')
{
  temp <- object$BUGSoutput$sims.matrix
  names <- colnames(temp)

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

