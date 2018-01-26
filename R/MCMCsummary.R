#' Summary function for MCMC output
#'
#' Extract summary information from MCMC output (mean, median, quantiles, Gelman-Rubin convergence statistic, number of effective samples, and specified custom metrics)
#' for specific parameters of interest.
#'
#' @param object Object containing MCMC output. See DETAILS below.
#' @param params Character string (or vector of character strings) denoting parameters to be returned in summary output.
#'
#' Default \code{'all'} returns all parameters in summary output.
#'
#' @param excl Character string (or vector of character strings) denoting parameters to exclude. Used in conjunction with \code{params} argument to select parameters of interest.
#'
#' @param ISB Ignore Square Brackets (ISB). Logical specifying whether square brackets should be ignored in the \code{params} and \code{excl} arguments. If \code{TRUE}, square brackets are ignored - input from \code{params} and \code{excl} are otherwise matched exactly. If \code{FALSE}, square brackets are not ignored - input from \code{params} and \code{excl} are matched using grep, which can take arguments in regular expression format. This allows partial names to be used when specifying parameters of interest.
#'
#' @param digits Number of digits to include for posterior summary. Values will be rounded to the specified number of digits (except for Rhat which is always rounded to 2 digits).
#'
#' Default is \code{digits = 2}.
#'
#' @param Rhat Logical specifying whether to calculate and display the Gelman-Rubin convergence statistic (Rhat). Values near 1 suggest convergence (Brooks and Gelman 1998). \code{Rhat = FALSE} will prevent display of this column in summary output. Specifying \code{Rhat = FALSE}, will increase function speed, particularly with very large `mcmc.list` objects.
#'
#' @param n.eff Logical specifying whether to calculate and display the number of effective samples for each parameter. Kruschke (2014) recommends n.eff > 10,000 for reasonably stable posterior estimates from M-H/Gibbs sampling derived output. \code{n.eff = FALSE} will prevent display of this column in summary output. Specifying \code{n.eff = FALSE}, will increase function speed, particularly with very large `mcmc.list` objects.
#'
#' @param func Function to be performed on MCMC output. If a function is specified, it will be evaluated on posteriors for each specified parameter and returned as a column in the summary output (or multiple columns if the function returns more than one value).
#'
#' @param func_name Character string (or vector of character strings) specifying labels for output from \code{func} argument. If \code{func_name} is not specified, columns with \code{func} argument will be labeled 'func'.
#'
#' @section Details:
#' \code{object} argument can be a \code{stanfit} object (\code{rstan} package), an \code{mcmc.list} object (\code{coda} package), an \code{R2jags} model object (\code{R2jags} package), a \code{jagsUI} model object (\code{jagsUI} package), or a matrix containing MCMC chains (each column representing MCMC output for a single parameter, rows representing iterations in the chain). The function automatically detects the object type and proceeds accordingly.
#'
#' @section Notes:
#'
#' For \code{mcmc.list} objects, Gelman-Rubin convergence statistic (Rhat) is calculated using the
#' \code{gelman.diag} function in the \code{coda} package.
#'
#' For \code{mcmc.list} objects, the number of effective samples is calculated using the \code{effectiveSize} function in the \code{coda} package.
#'
#' @return Function returns summary information (including parameter posterior mean, posterior sd, 2.5\% quantile, median, 97.5\% quantile, Gelman-Rubin convergence statistic (Rhat), number of effective samples, and other specified metrics) for specified parameters.
#'
#'
#' @section References:
#'
#' Brooks, S. P., and A. Gelman. 1998. General methods for monitoring convergence of iterative simulations. Journal of Computational and Graphical Statistics 7:434.
#'
#' Kruschke, J. 2014. Doing Bayesian data analysis: A tutorial with R, JAGS, and Stan. Academic Press.
#'
#' @examples
#' #Load data
#' data(MCMC_data)
#'
#' #Summary information for MCMC output
#' MCMCsummary(MCMC_data)
#'
#' #Just 'beta' parameters
#' MCMCsummary(MCMC_data, params = 'beta')
#'
#' #Just 'beta[1]', 'gamma[4]', and 'alpha[3]'
#' #'params' takes regular expressions when ISB = FALSE, square brackets must be escaped with '\\'
#' MCMCsummary(MCMC_data, params = c('beta\\[1\\]', 'gamma\\[4\\]', 'alpha\\[3\\]'), ISB = FALSE)
#'
#' @export


MCMCsummary <- function(object,
                      params = 'all',
                      excl = NULL,
                      ISB = TRUE,
                      digits = 2,
                      Rhat = TRUE,
                      n.eff = FALSE,
                      func = NULL,
                      func_name = NULL)
{
  #SORTING BLOCK
  if(class(object) == 'jagsUI')
  {
    object <- object$samples
  }

  if(typeof(object) == 'double')
  {
    object2 <- MCMCchains(object, params, excl, ISB, mcmc.list = FALSE)
  }else{
    object2 <- MCMCchains(object, params, excl, ISB, mcmc.list = TRUE)
  }


  #PROCESSING BLOCK
    if(coda::is.mcmc.list(object2) == TRUE)
    {
      np <- NCOL(object2[[1]])

      if(np > 1)
      {
        ch_bind <- do.call('rbind', object2)
      }else{
        ch_bind <- as.matrix(object2)
      }

      bind_mn <- round(apply(ch_bind, 2, mean), digits = digits)
      bind_sd <- round(apply(ch_bind, 2, stats::sd), digits = digits)
      bind_LCI <- round(apply(ch_bind, 2, stats::quantile, probs= 0.025), digits = digits)
      bind_med <- round(apply(ch_bind, 2, stats::median), digits = digits)
      bind_UCI <- round(apply(ch_bind, 2, stats::quantile, probs= 0.975), digits = digits)

      if(Rhat == TRUE)
      {
        if(length(object2) > 1)
        {
          #If > 750 params use loop to calc Rhat
          if(NCOL(object2[[1]]) > 750)
          {
            r_hat <- c(rep(NA, NCOL(object2[[1]])))
            for (v in 1:length(r_hat))
            {
              r_hat[v] <- round(coda::gelman.diag(object2[,v])$psrf[,1], digits = 2)
            }

          }else{
            r_hat <- round(coda::gelman.diag(object2, multivariate = FALSE)$psrf[,1], digits = 2)
          }

        }else{
          warning('Rhat statistic cannot be calculated with one chain. NAs inserted.')
          r_hat <- rep(NA, NCOL(object2))
        }
        x <- cbind(bind_mn, bind_sd, bind_LCI, bind_med, bind_UCI, r_hat)
        colnames(x) <- c('mean', 'sd', '2.5%','50%','97.5%', 'Rhat')
      }else{
        x <- cbind(bind_mn, bind_sd, bind_LCI, bind_med, bind_UCI)
        colnames(x) <- c('mean', 'sd', '2.5%','50%','97.5%')
      }

      if(n.eff == TRUE)
      {
        bind_neff <- round(coda::effectiveSize(object2), digits = 0)
        x2 <- cbind(x, bind_neff)
        colnames(x2)[ncol(x2)] <- 'n.eff'
      }else{
        x2 <- x
      }

      if(!is.null(func))
      {
        tmp <- round(apply(ch_bind, 2, func), digits = digits)

        if(!is.null(dim(tmp)) & NROW(tmp) > 1)
        {
          x3 <- x2
          for (i in 1:NROW(tmp))
          {
            x3 <- cbind(x3, tmp[i,])
            if (!is.null(func_name))
            {
              if (length(func_name) != NROW(tmp))
              {
                stop('length(func_name) must equal number of func outputs')
              }
              colnames(x3)[ncol(x2) + i] <- func_name[i]
            }else{
              colnames(x3)[ncol(x2) + i] <- 'func'
            }
          }
        }else{
          x3 <- cbind(x2, tmp)
          if (!is.null(func_name))
          {
            if (length(func_name) > 1)
            {
              stop('length(func_name) must equal number of func outputs')
            }
            colnames(x3)[ncol(x3)] <- func_name
          }else{
            colnames(x3)[ncol(x3)] <- 'func'
          }
        }
      }else{
        x3 <- x2
      }
      mcmc_summary <- x3
    }


    if(typeof(object2) == 'double')
    {
      np <- NCOL(object2)
      ch_bind <- object2

      bind_mn <- round(apply(ch_bind, 2, mean), digits = digits)
      bind_sd <- round(apply(ch_bind, 2, stats::sd), digits = digits)
      bind_LCI <- round(apply(ch_bind, 2, stats::quantile, probs= 0.025), digits = digits)
      bind_med <- round(apply(ch_bind, 2, stats::median), digits = digits)
      bind_UCI <- round(apply(ch_bind, 2, stats::quantile, probs= 0.975), digits = digits)

      if(Rhat == TRUE)
      {
        warning('Rhat statistic cannot be calculated without individual chains. NAs inserted.')
        r_hat <- rep(NA, np)
        x <- cbind(bind_mn, bind_sd, bind_LCI, bind_med, bind_UCI, r_hat)
        colnames(x) <- c('mean', 'sd', '2.5%','50%','97.5%', 'Rhat')
      }else{
        x <- cbind(bind_mn, bind_sd, bind_LCI, bind_med, bind_UCI)
        colnames(x) <- c('mean', 'sd', '2.5%','50%','97.5%')
      }

      if(n.eff == TRUE)
      {
        warning('Number of effective samples cannot be calculated without individual chains. NAs inserted.')
        bind_neff <- rep(NA, np)
        x2 <- cbind(x, bind_neff)
        colnames(x2)[ncol(x2)] <- 'n.eff'
      }else{
        x2 <- x
      }

      if(!is.null(func))
      {
        tmp <- round(apply(ch_bind, 2, func), digits = digits)

        if(!is.null(dim(tmp)) & NROW(tmp) > 1)
        {
          x3 <- x2
          for (i in 1:NROW(tmp))
          {
            x3 <- cbind(x3, tmp[i,])
            if (!is.null(func_name))
            {
              if (length(func_name) != NROW(tmp))
              {
                stop('length(func_name) must equal number of func outputs')
              }
              colnames(x3)[ncol(x2) + i] <- func_name[i]
            }else{
              colnames(x3)[ncol(x2) + i] <- 'func'
            }
          }
        }else{
          x3 <- cbind(x2, tmp)
          if (!is.null(func_name))
          {
            if (length(func_name) > 1)
            {
              stop('length(func_name) must equal number of func outputs')
            }
            colnames(x3)[ncol(x3)] <- func_name
          }else{
            colnames(x3)[ncol(x3)] <- 'func'
          }
        }
      }else{
        x3 <- x2
      }
      mcmc_summary <- x3
    }
  return(mcmc_summary)
}
