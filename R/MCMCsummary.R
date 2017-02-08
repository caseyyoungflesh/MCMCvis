#' Summary function for MCMC output
#'
#' Extract summary information from MCMC output (mean, median, quantiles, and Gelman-Rubin convergence statistic)
#' for specific parameters of interest.
#'
#' @param object Object containing MCMC output. See DETAILS below.
#' @param params Character string (or vector of character strings) denoting parameters to be
#' returned in summary output. Partial names may be used to return all parameters containing
#' that set of characters.
#'
#' Default \code{'all'} returns all parameters in summary output.
#'
#' @param excl Character string (or vector of character strings) denoting parameters to exclude.
#' Partial names may be used to exclude all parameters containing that set of characters. Used in
#' conjunction with \code{params} argument to select parameters of interest.
#'
#' @param digits Number of digits to include for posterior summary. Values will be rounded to the specified
#' number of digits.
#' Default is \code{digits = 2}.
#'
#' @param Rhat If \code{TRUE}, summary information contains Gelman-Rubin convergence statistic (Rhat)
#' and if \code{FALSE}, Rhat output is masked.
#' @section Details:
#' \code{object} argument can be a \code{stanfit} object (\code{rstan} package), an \code{mcmc.list} object
#' (\code{coda} package), an \code{R2jags} model object (\code{R2jags} package), or a matrix containing MCMC
#' chains (each column representing MCMC output for a single parameter, rows representing iterations in the chain).
#' The function automatically detects the object type and proceeds accordingly.
#'
#' @section Notes:
#'
#' For \code{mcmc.list} objects, Gelman-Rubin convergence statistic (Rhat) is calculated using the
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
#' MCMCsummary(MCMC_data)
#'
#' #Just 'beta' parameters
#' MCMCsummary(MCMC_data, params= 'beta')
#'
#' #Just 'beta[1]', 'gamma[4]', and 'alpha[3]'
#' MCMCsummary(MCMC_data, params= c('beta[1]', 'gamma[4]', 'alpha[3]'))
#'
#' @export


MCMCsummary <- function(object,
                      params = 'all',
                      excl = NULL,
                      digits = 2,
                      Rhat = TRUE)
{
  if(coda::is.mcmc.list(object) != TRUE &
     typeof(object) != 'double' &
     typeof(object) != 'list' &
     typeof(object) != 'S4')
  {
    stop('Invalid object type. Input must be stanfit object, mcmc.list object, rjags object, or matrix with MCMC chains.')
  }


  if(typeof(object) == 'list' & coda::is.mcmc.list(object) == FALSE)
  {
    x <- round(object$BUGSoutput$summary[,c(1, 3, 5, 7, 8)], digits = digits)
    #to.rm <- which(rownames(x) == 'deviance')
    mcmc_summary <- x[,] #already have mcmc_summary
    names <- rownames(mcmc_summary)
  } else {
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

      bind_mn <- round(apply(ch_bind, 2, mean), digits = digits)
      bind_LCI <- round(apply(ch_bind, 2, stats::quantile, probs= 0.025), digits = digits)
      bind_med <- round(apply(ch_bind,2, stats::median), digits = digits)
      bind_UCI <- round(apply(ch_bind, 2, stats::quantile, probs= 0.975), digits = digits)
      r_hat <- round(coda::gelman.diag(temp, multivariate = FALSE)$psrf[,1], digits = digits)

      mcmc_summary <- cbind(bind_mn, bind_LCI, bind_med, bind_UCI, r_hat)
      colnames(mcmc_summary) <- c('mean','2.5%','50%','97.5%', 'Rhat')
    }

    if(typeof(object2) == 'double')
    {
      temp <- object2
      names <- colnames(temp)

      bind_mn <- round(apply(temp, 2, mean), digits = digits)
      bind_LCI <- round(apply(temp, 2, stats::quantile, probs= 0.025), digits = digits)
      bind_med <- round(apply(temp,2, stats::median), digits = digits)
      bind_UCI <- round(apply(temp, 2, stats::quantile, probs= 0.975), digits = digits)
      if(Rhat == TRUE)
      {
        warning('Rhat statistic cannot be calculated without individaul chains. NAs inserted.')
      }
      r_hat <- rep(NA, NCOL(temp))

      mcmc_summary <- cbind(bind_mn, bind_LCI, bind_med, bind_UCI, r_hat)
      colnames(mcmc_summary) <- c('mean','2.5%','50%','97.5%', 'Rhat')
    }
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
        OUT <- mcmc_summary
      }else{
        OUT <- mcmc_summary[-to.rm2,]
      }
    }else
     {
      get.rows <- grep(paste(params), names, fixed=TRUE)
      if (length(get.rows) < 1)
      {
        stop(paste0('"', params, '"', ' not found in MCMC ouput.'))
      }

      if(!is.null(excl))
      {
        if(identical(get.rows, to.rm2))
        {
          stop('No parameters selected.')
        }

        matched <- which(get.rows == to.rm2)
        if (length(matched) > 0)
        {
          g_filt <- get.rows[-matched]
        }else {
          g_filt <- get.rows
        }

      }else{
        g_filt <- get.rows
       }

      OUT <- mcmc_summary[g_filt,]
    }
  }else {
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

    if(!is.null(excl))
    {
      if(identical(grouped, to.rm2))
      {
        stop('No parameters selected.')
      }

      matched <- stats::na.omit(match(to.rm2, grouped))
      if (length(matched) > 0)
      {
        rows <- grouped[-matched]
      } else{
        rows <- grouped
      }

      to.rm <- which(duplicated(rows))
      if(length(to.rm) > 0)
      {
        g_filt <- rows[-to.rm]
      }else
      {
        g_filt <- rows
      }
      OUT <- mcmc_summary[g_filt,]
    } else{

      to.rm <- which(duplicated(grouped))
      if(length(to.rm) > 0)
      {
        g_filt <- grouped[-to.rm]
      }else
      {
        g_filt <- grouped
      }

      OUT <- mcmc_summary[g_filt,]
    }
  }

  if(Rhat == TRUE)
  {
    return(OUT)
  }else
  {
    return(OUT[,-5])
  }
}
