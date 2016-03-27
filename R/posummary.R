#' Summary function for JAGS object
#'
#' Extracts summary information for specific parameters of interests from JAGS object
#' output from JAGS model fit in R
#'
#' @param object JAGS model object - output from JAGS model
#' @param params Character string (or vector of character strings) denoting parameters to be
#' returned in summary output. Partial names may be used to return all parameters containing
#' that set of characters. Default \code{'all'} returns all parameters in summary output.
#' @param Rhat If \code{TRUE}, summary information contains Gelman-Rubin convergence statistic (Rhat)
#' and if \code{FALSE}, Rhat output is masked.
#' @section Details:
#' Function currently only supports JAGS model objects output by \code{R2jags} package.
#' @return \code{posummary(params='all')} returns summary data for all parameters contained within
#' JAGS model object.
#'
#' \code{posummary(params=c('beta[1]', 'beta[2]'))} returns summary data for just parameters
#' \code{'beta[1]'} and \code{'beta[2]'}.
#'
#' \code{posummary(params=c('beta'))} returns summary data for all parameters containing \code{'beta'}
#'  in their name.
#'
#'  Default summary information includes (posterior mean, 2.5\% quantile, median, 97.5\% quantile,
#'  and Rhat).
#'
#' @export

posummary <- function(object, params = 'all', Rhat = TRUE)
{
  temp <- object$BUGSoutput$summary
  names <- rownames(temp)

  if (length(params) == 1)
  {
    if (params == 'all')
    {
      OUT <- temp
    }else
    {
        get.rows <- grep(paste(params), names, fixed=TRUE)
        OUT <- temp[get.rows,]
    }
  }else
  {
    grouped <- c()
    for (i in 1:length(params))
    {
      get.rows <- grep(paste(params[i]), names, fixed=TRUE)
      grouped <- c(grouped, get.rows)
    }

    to.rm <- which(duplicated(grouped))
    if(length(to.rm) >0)
    {
      g_filt <- grouped[-to.rm]
    }else
    {
      g_filt <- grouped
    }

    OUT <- temp[g_filt,]
  }

  if(is.null(dim(OUT)))
  {
    if(Rhat == TRUE)
    {
      print(OUT[c(1,3,5,7,8)])
    }else
    {
      print(OUT[c(1,3,5,7)])
    }
  }else
  {
    if(Rhat == TRUE)
    {
      print(OUT[,c(1,3,5,7,8)])
    }else
    {
      print(OUT[,c(1,3,5,7)])
    }
  }
}
