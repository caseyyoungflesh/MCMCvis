#' Summary function for MCMC output
#'
#' Extract summary information from MCMC output for specific parameters of interest.
#'
#' @param object Object containing MCMC output. See DETAILS below.
#' @param params Character string (or vector of character strings) denoting parameters to be
#' returned in summary output. Partial names may be used to return all parameters containing
#' that set of characters.
#'
#' Default \code{'all'} returns all parameters in summary output.
#' @param Rhat If \code{TRUE}, summary information contains Gelman-Rubin convergence statistic (Rhat)
#' and if \code{FALSE}, Rhat output is masked.
#' @section Details:
#' \code{object} argument can be a \code{stanfit} object (\code{rstan} package), an \code{mcmc.list} object
#' (\code{coda} package), an \code{R2jags} model object (\code{R2jags} package), or a matrix containing MCMC
#' chains (each column representing MCMC output for a single parameter, rows representing iterations in the chain).
#'
#' @section Notes:
#'
#' For \code{mcmc.list} objects, Gelman-Rubin convergence statistic  (Rhat) is calculated using the
#' \code{gelman.diag} function in the \code{coda} package.
#'
#' @return Function returns summary information (including parameter posterior mean, 2.5\% quantile, median, 97.5\%
#'  quantile, and Gelman-Rubin convergence statistic (Rhat)) for specified parameters.
#'
#' @examples
#' #Load data
#' data(MCMC_data)
#'
#' #Summary information for MCMC output
#' posummary(MCMC_data)
#'
#' #Just 'beta' parameters
#' posummary(MCMC_data, params= 'beta')
#'
#' #Just 'beta[1]', 'gamma[4]', and 'alpha[3]'
#' posummary(MCMC_data, params= c('beta[1]', 'gamma[4]', 'alpha[3]'))
#'
#' @export
#' @import coda

posummary <- function(object,
                      params = 'all',
                      Rhat = TRUE)
{
  if(typeof(object) == 'S4')
  {
    object2 <- rstan::As.mcmc.list(object)
  } else {
    object2 <- object
  }

     if(coda::is.mcmc.list(object2) == TRUE)
    {
      temp <- object2
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
          if (length(get.rows) < 1)
          {
            stop(paste0('"', params, '"', ' not found in MCMC ouput.'))
          }
          OUT <- mcmc_summary[get.rows,]
        }
      }else
      {
        grouped <- c()
        for (i in 1:length(params))
        {
          get.rows <- grep(paste(params[i]), names, fixed=TRUE)
          if (length(get.rows) < 1)
          {
            stop(paste0('"', params[i], '"', ' not found in MCMC ouput.'))
          }
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



    if(typeof(object2) == 'double')
    {
      temp <- object2
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
          if (length(get.rows) < 1)
          {
            stop(paste0('"', params, '"', ' not found in MCMC ouput.'))
          }
          OUT <- mcmc_summary[get.rows,]
        }
      }else
      {
        grouped <- c()
        for (i in 1:length(params))
        {
          get.rows <- grep(paste(params[i]), names, fixed=TRUE)
          if (length(get.rows) < 1)
          {
            stop(paste0('"', params[i], '"', ' not found in MCMC ouput.'))
          }
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

    if(typeof(object2) == 'list' & coda::is.mcmc.list(object2) == FALSE)
    {
      temp <- object2$BUGSoutput$summary
      names <- rownames(temp)

      if (length(params) == 1)
      {
        if (params == 'all')
        {
          OUT <- temp[, c(1,3,5,7,8)]
        }else
        {
          get.rows <- grep(paste(params), names, fixed=TRUE)
          if (length(get.rows) < 1)
          {
            stop(paste0('"', params, '"', ' not found in MCMC ouput.'))
          }
          OUT <- temp[get.rows, c(1,3,5,7,8)]
        }
      }else
      {
        grouped <- c()
        for (i in 1:length(params))
        {
          get.rows <- grep(paste(params[i]), names, fixed=TRUE)
          if (length(get.rows) < 1)
          {
            stop(paste0('"', params[i], '"', ' not found in MCMC ouput.'))
          }
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
     typeof(object) != 'list' &
     typeof(object) != 'S4')
  {
      stop('Invalid object type. Input must be stanfit object, mcmc.list object, rjags object, or matrix with MCMC chains.')
  }


  if(is.null(dim(OUT)))
  {
    if(Rhat == TRUE)
    {
      return(OUT)
    }else
    {
      return(OUT[-5])
    }
  }else
  {
    if(Rhat == TRUE)
    {
      return(OUT)
    }else
    {
      return(OUT[,-5])
    }
  }
}


