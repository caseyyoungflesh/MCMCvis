#' Summary function for MCMC output
#'
#' Extract summary information from MCMC output for specific parameters of interest.
#'
#' @param object Object containing MCMC output. See DETAILS below.
#' @param par Character string (or vector of character strings) denoting parameters to be
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
#' posummary(MCMC_data, par= 'beta')
#'
#' #Just 'beta[1]', 'gamma[4]', and 'alpha[3]'
#' posummary(MCMC_data, par= c('beta[1]', 'gamma[4]', 'alpha[3]'))
#'
#' @export
#' @import coda
object <- fit
object <- SD_out
posummary <- function(object,
                      par = c('all'),
                      excl = NULL,
                      Rhat = TRUE)
{
  if(typeof(object) == 'list' & coda::is.mcmc.list(object) == FALSE)
  {
    #modified from R2jags package
    x <- object$BUGSoutput
    mclist <- vector("list", x$n.chains)
    mclis <- vector("list", x$n.chains)
    strt <- x$n.burnin + 1
    end <- x$n.iter
    ord <- dimnames(x$sims.array)[[3]]
    for (i in 1:x$n.chains)
    {
      tmp1 <- x$sims.array[, i, ord]
      mclis[[i]] <- mcmc(tmp1, start = strt, end = end, thin = x$n.thin)
    }
    object2 <- as.mcmc.list(mclis)

  } else {
    if(typeof(object) == 'S4')
    {
      object2 <- rstan::As.mcmc.list(object)
    } else {
      object2 <- object
    }
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

  if (length(par) == 1)
  {
    if (par == 'all')
    {
      if(is.null(excl))
      {
        OUT <- mcmc_summary
      }else{
        OUT <- mcmc_summary[-to.rm2,]
      }
    }else
     {
      get.rows <- grep(paste(par), names, fixed=TRUE)
      if (length(get.rows) < 1)
      {
        stop(paste0('"', par, '"', ' not found in MCMC ouput.'))
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
    for (i in 1:length(par))
    {
      get.rows <- grep(paste(par[i]), names, fixed=TRUE)
      if (length(get.rows) < 1)
      {
        stop(paste0('"', par[i], '"', ' not found in MCMC ouput.'))
      }
      grouped <- c(grouped, get.rows)
    }

    if(!is.null(excl))
    {
      if(identical(grouped, to.rm2))
      {
        stop('No parameters selected.')
      }

      matched <- na.omit(match(to.rm2, grouped))
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

  if(coda::is.mcmc.list(object2) != TRUE &
     typeof(object2) != 'double' &
     typeof(object2) != 'list' &
     typeof(object2) != 'S4')
  {
      stop('Invalid object type. Input must be stanfit object, mcmc.list object, rjags object, or matrix with MCMC chains.')
  }


  if(Rhat == TRUE)
  {
    return(OUT)
  }else
  {
    return(OUT[,-5])
  }
}
