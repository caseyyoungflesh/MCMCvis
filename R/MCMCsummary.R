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

###names and temp
object <- out
params = c('beta', 'mu')
excl <- c('mu[1]', 'mu[3]')

object <- MCMC_data
params = c('alpha', 'beta')
excl <- c('beta[1]', 'beta[3]')


MCMCsummary(object, digits = 2,
            params = 'beta', excl = 'beta[1]', ISB = FALSE,
            Rhat = FALSE)


MCMCsummary <- function(object,
                      params = 'all',
                      excl = NULL,
                      ISB = TRUE,
                      digits = 2,
                      Rhat = TRUE) #ISB = TRUE ignores [], otherwise matches exact string. ISB = FALSE does not ignore [], and behaves like grep
{
  if(coda::is.mcmc.list(object) != TRUE &
     typeof(object) != 'double' &
     typeof(object) != 'list' &
     typeof(object) != 'S4')
  {
    stop('Invalid object type. Input must be stanfit object, mcmc.list object, rjags object, or matrix with MCMC chains.')
  }

  #NAME SORTING BLOCK
  #R2jags object
  if(class(object) == 'rjags')
  {
    if(ISB == TRUE)
    {
      names <- vapply(strsplit(rownames(object$BUGSoutput$summary),
                               split = "[", fixed = TRUE), `[`, 1, FUN.VALUE=character(1))
    }else{
      names <- rownames(object$BUGSoutput$summary)
    }
  }else {
    #Stan object
    if(typeof(object) == 'S4')
    {
      object2 <- rstan::As.mcmc.list(object)
    } else {
      object2 <- object
    }
    #MCMClist object
    if(coda::is.mcmc.list(object2) == TRUE)
    {
      temp <- object2
      if(ISB == TRUE)
      {
        names <- vapply(strsplit(colnames(temp[[1]]),
                                 split = "[", fixed = TRUE), `[`, 1, FUN.VALUE=character(1))
      }else{
        names <- colnames(temp[[1]])
      }
    }
    #matrix object
    if(typeof(object2) == 'double')
    {
      temp <- object2
      if(ISB == TRUE)
      {
        names <- vapply(strsplit(colnames(temp),
                                 split = "[", fixed = TRUE), `[`, 1, FUN.VALUE=character(1))
      }else{
        names <- colnames(temp)
      }
    }
    #jags.samples object (mcarray)
    if(class(object2[[1]]) == 'mcarray')
    {
      stop('Invalid object type. jags.samples objects not currently supported. Input must be stanfit object, mcmc.list object, rjags object, or matrix with MCMC chains.')
    }
  }

  #INDEX BLOCK
  #exclusions
  if(!is.null(excl))
  {
    rm_ind <- c()
    for (i in 1:length(excl))
    {
      if(ISB == TRUE)
      {
        n_excl <- vapply(strsplit(excl,
                                 split = "[", fixed = TRUE), `[`, 1, FUN.VALUE=character(1))
        rm_ind <- c(rm_ind, which(names %in% n_excl[i]))
      }else{
        n_excl <- excl
        rm_ind <- c(rm_ind, grep(n_excl[i], names, fixed = TRUE))
      }
    }
    if(length(rm_ind) < 1)
    {
      stop(paste0('"', excl, '"', ' not found in MCMC ouput.'))
    }
    dups <- which(duplicated(rm_ind))
    if(length(dups) > 0)
    {
      rm_ind2 <- rm_ind[-dups]
    }else{
      rm_ind2 <- rm_ind
    }
  }

  #selections
  if (length(params) == 1)
  {
    if (params == 'all')
    {
      if(is.null(excl))
      {
        f_ind <- 1:length(names)
      }else{
        f_ind <- (1:length(names))[-rm_ind2]
      }
    }else {
      if(ISB == TRUE)
      {
        get_ind <- which(names %in% params)
      }else{
        get_ind <- grep(paste(params), names, fixed = TRUE)
      }

      if (length(get_ind) < 1)
      {
        stop(paste0('"', params, '"', ' not found in MCMC ouput.'))
      }
      if(!is.null(excl))
      {
        if(identical(get_ind, rm_ind2))
        {
          stop('No parameters selected.')
        }
        matched <- stats::na.omit(match(rm_ind2, get_ind))
        if (length(matched) > 0)
        {
          f_ind <- get_ind[-matched]
        }else {
          f_ind <- get_ind
        }
      }else {
        f_ind <- get_ind
      }
    }
  }else {
    grouped <- c()
    for (i in 1:length(params))
    {
      if(ISB == TRUE)
      {
        get_ind <- which(names %in% params[i])
      }else{
        get_ind <- grep(paste(params[i]), names, fixed=TRUE)
      }

      if (length(get_ind) < 1)
      {
        stop(paste0('"', params[i], '"', ' not found in MCMC ouput.'))
      }
      grouped <- c(grouped, get_ind)
    }
    if(!is.null(excl))
    {
      if(identical(grouped, rm_ind2))
      {
        stop('No parameters selected.')
      }
      matched <- stats::na.omit(match(rm_ind2, grouped))
      if (length(matched) > 0)
      {
        t_ind <- grouped[-matched]
      } else{
        t_ind <- grouped
      }
      to.rm <- which(duplicated(t_ind))
      if(length(to.rm) > 0)
      {
        f_ind <- t_ind[-to.rm]
      }else
      {
        f_ind <- t_ind
      }
    } else{
      to.rm <- which(duplicated(grouped))
      if(length(to.rm) > 0)
      {
        f_ind <- grouped[-to.rm]
      }else
      {
        f_ind <- grouped
      }
    }
  }

  #PROCESSING BLOCK
  if(typeof(object) == 'list' & coda::is.mcmc.list(object) == FALSE)
  {
    if(Rhat == TRUE)
    {
      x <- round(object$BUGSoutput$summary[,c(1, 2, 3, 5, 7, 8)], digits = digits)
      mcmc_summary <- x[f_ind,]
    }else{
      x <- round(object$BUGSoutput$summary[,c(1, 2, 3, 5, 7)], digits = digits)
      mcmc_summary <- x[f_ind,]
    }
  } else {

    if(coda::is.mcmc.list(object2) == TRUE)
    {
      if(length(f_ind) > 1)
      {
        dsort <- do.call(coda::mcmc.list, temp[,f_ind])
        ch_bind <- do.call('rbind', dsort)
      }else{
        dsort <- do.call(coda::mcmc.list, temp[,c(f_ind, f_ind)])
        ch_bind <- do.call('rbind', dsort)
      }

      bind_mn <- round(apply(ch_bind, 2, mean), digits = digits)
      bind_sd <- round(apply(ch_bind, 2, sd), digits = digits)
      bind_LCI <- round(apply(ch_bind, 2, stats::quantile, probs= 0.025), digits = digits)
      bind_med <- round(apply(ch_bind,2, stats::median), digits = digits)
      bind_UCI <- round(apply(ch_bind, 2, stats::quantile, probs= 0.975), digits = digits)

      if(Rhat == TRUE)
      {
        r_hat <- round(coda::gelman.diag(dsort, multivariate = FALSE)$psrf[,1], digits = digits)
        mcmc_summary <- cbind(bind_mn, bind_sd, bind_LCI, bind_med, bind_UCI, r_hat)
        colnames(mcmc_summary) <- c('mean', 'sd', '2.5%','50%','97.5%', 'Rhat')
      }else{
        mcmc_summary <- cbind(bind_mn, bind_sd, bind_LCI, bind_med, bind_UCI)[1,]
        colnames(mcmc_summary) <- c('mean', 'sd', '2.5%','50%','97.5%')
      }
    }

    if(typeof(object2) == 'double')
    {
      if(length(f_ind) > 1)
      {
        dsort <- temp[,f_ind]
      }else{
        dsort <- temp[,c(f_ind, f_ind)]
      }

      bind_mn <- round(apply(dsort, 2, mean), digits = digits)
      bind_sd <- round(apply(ch_bind, 2, sd), digits = digits)
      bind_LCI <- round(apply(dsort, 2, stats::quantile, probs= 0.025), digits = digits)
      bind_med <- round(apply(dsort,2, stats::median), digits = digits)
      bind_UCI <- round(apply(dsort, 2, stats::quantile, probs= 0.975), digits = digits)

      if(Rhat == TRUE)
      {
        warning('Rhat statistic cannot be calculated without individaul chains. NAs inserted.')
        r_hat <- rep(NA, NCOL(dsort))
        mcmc_summary <- cbind(bind_mn, bind_sd, bind_LCI, bind_med, bind_UCI, r_hat)
        colnames(mcmc_summary) <- c('mean', 'sd', '2.5%','50%','97.5%', 'Rhat')
      }else{
        mcmc_summary <- cbind(bind_mn, bind_sd, bind_LCI, bind_med, bind_UCI)[1,]
        colnames(mcmc_summary) <- c('mean', 'sd', '2.5%','50%','97.5%')
      }
    }
  }
  return(mcmc_summary)
}
