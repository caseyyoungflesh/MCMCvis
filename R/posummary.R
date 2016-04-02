#' Summary function for MCMC output
#'
#' Extract summary information from MCMC output for specific parameters of interest.
#'
#' @param object Object containing MCMC output. See DETAILS below.
#' @param params Character string (or vector of character strings) denoting parameters to be
#' returned in summary output. Partial names may be used to return all parameters containing
#' that set of characters.
#'
#' Default \code{all} returns all parameters in summary output.
#' @param Rhat If \code{TRUE}, summary information contains Gelman-Rubin convergence statistic (Rhat)
#' and if \code{FALSE}, Rhat output is masked.
#' @section Details:
#' \code{object} argument can be an \code{mcmc.list} object, an \code{R2jags} model object (output from the \code{R2jags}
#' package), or a matrix containing MCMC chains (each column representing MCMC output for a single parameter, rows
#' representing iterations in the chain).
#'
#' @section Notes:
#'  Default summary information includes (parameter posterior mean, 2.5\% quantile, median, 97.5\%
#'  quantile, and Gelman-Rubin convergence statistic (Rhat)).
#'
#' For \code{mcmc.list} objects, Gelman-Rubin convergence statistic  (Rhat) is calculated using the
#' \code{gelman.diag} function in the \code{coda} package.
#'
#' @return \code{posummary(params='all')} returns summary data for all parameters.
#'
#' \code{posummary(params=c('beta[1]', 'beta[2]'))} returns summary data for just parameters
#' \code{beta[1]} and \code{beta[2]}.
#'
#' \code{posummary(params=c('beta'))} returns summary data for all parameters containing \code{beta}
#'  in their name.
#'
#' @export
#' @import coda


posummary <- function(object,
                      params = 'all',
                      Rhat = TRUE)
{
     if(coda::is.mcmc.list(object) == TRUE)
    {
      temp <- object
      names <- colnames(temp[[1]])

      ch_bind <- do.call('rbind', temp)

      bind_mn <- apply(ch_bind, 2, mean)
      bind_LCI <- apply(ch_bind, 2, quantile, probs= 0.025)
      bind_med <- apply(ch_bind,2, median)
      bind_UCI <- apply(ch_bind, 2, quantile, probs= 0.975)
      r_hat <- coda::gelman.diag(temp)$psrf[,1]

      mcmc_summary <- cbind(bind_mn, bind_LCI, bind_med, bind_UCI, r_hat)
      colnames(mcmc_summary) <- c('mean','2.5%','50%','97.5%', 'Rhat')

      if (length(params) == 1)
      {
        if (params == 'all')
        {
          OUT <- mcmc_summary
        }else
        {
          get.rows <- grep(paste(params), names, fixed=TRUE)
          OUT <- mcmc_summary[get.rows,]
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

        OUT <- mcmc_summary[g_filt,]
      }
    }



    if(typeof(object) == 'double')
    {
      temp <- object
      names <- colnames(temp)

      bind_mn <- apply(temp, 2, mean)
      bind_LCI <- apply(temp, 2, quantile, probs= 0.025)
      bind_med <- apply(temp,2, median)
      bind_UCI <- apply(temp, 2, quantile, probs= 0.975)
      if(Rhat == TRUE)
      {
        warning('Rhat statistic cannot be calculated without individaul chains. NAs inserted.')
      }
      r_hat <- rep(NA, NCOL(temp))

      mcmc_summary <- cbind(bind_mn, bind_LCI, bind_med, bind_UCI, r_hat)
      colnames(mcmc_summary) <- c('mean','2.5%','50%','97.5%', 'Rhat')

      if (length(params) == 1)
      {
        if (params == 'all')
        {
          OUT <- mcmc_summary
        }else
        {
          get.rows <- grep(paste(params), names, fixed=TRUE)
          OUT <- mcmc_summary[get.rows,]
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

        OUT <- mcmc_summary[g_filt,]
      }

    }

    if(typeof(object) == 'list' & coda::is.mcmc.list(object) == FALSE)
    {
      temp <- object$BUGSoutput$summary
      names <- rownames(temp)

      if (length(params) == 1)
      {
        if (params == 'all')
        {
          OUT <- temp[, c(1,3,5,7,8)]
        }else
        {
          get.rows <- grep(paste(params), names, fixed=TRUE)
          OUT <- temp[get.rows, c(1,3,5,7,8)]
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

        OUT <- temp[g_filt, c(1,3,5,7,8)]
      }
    }


  if(coda::is.mcmc.list(object) != TRUE &
     typeof(object) != 'double' &
     typeof(object) != 'list')
  {
      stop('Invalid object type. Input must be mcmc.list object, rjags object, or matrix with MCMC chains.')
  }


  if(is.null(dim(OUT)))
  {
    if(Rhat == TRUE)
    {
      print(OUT)
    }else
    {
      print(OUT[-5])
    }
  }else
  {
    if(Rhat == TRUE)
    {
      print(OUT)
    }else
    {
      print(OUT[,-5])
    }
  }
}


