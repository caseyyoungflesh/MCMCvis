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
#' @param ISB Ignore Square Brackets (ISB). Logical specifying whether square brackets should be ignored in the \code{params} and \code{excl} arguments. If \code{FALSE}, square brackets are ignored - input from \code{params} and \code{excl} are otherwise matched exactly. If \code{TRUE}, square brackets are not ignored - input from \code{params} and \code{excl} are matched using grep, allowing partial names to be used when specifying parameters of interest.
#'
#' @param digits Number of digits to include for posterior summary. Values will be rounded to the specified number of digits (except for Rhat which is always rounded to 2 digits).
#'
#' Default is \code{digits = 2}.
#'
#' @param Rhat Logical specifying whether to calculate and display the Gelman-Rubin convergence statistic (Rhat). Values near 1 suggest convergence (Brooks and Gelman 1998). \code{Rhat = FALSE} will prevent display of this column in summary output. Specifying \code{Rhat = FALSE}, will increase function speed, particularly with very large `mcmc.list` objects.
#'
#' @param n.eff Logical specifying whether to calculate and display the number of effective samples for each parameter. Kruschke (2014) recommends n.eff > 10,000 for reasonably stable posterior estimates. \code{n.eff = FALSE} will prevent display of this column in summary output. Specifying \code{n.eff = FALSE}, will increase function speed, particularly with very large `mcmc.list` objects.
#'
#' @param func Function to be performed on MCMC output. If a function is specified, it will be evaluated on posteriors for each specified parameter and returned as a column in the summary output (or multiple columns if the function returns more than one value).
#'
#' @param func_name Character string (or vector of character strings) specifying labels for output from \code{func} argument. If \code{func_name} is not specified, columns with \code{func} argument will be labeled 'func'.
#'
#' @section Details:
#' \code{object} argument can be a \code{stanfit} object (\code{rstan} package), an \code{mcmc.list} object (\code{coda} package), an \code{R2jags} model object (\code{R2jags} package), or a matrix containing MCMC chains (each column representing MCMC output for a single parameter, rows representing iterations in the chain). The function automatically detects the object type and proceeds accordingly.
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
#' MCMCsummary(MCMC_data, params = c('beta[1]', 'gamma[4]', 'alpha[3]'), ISB = FALSE)
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
      rh <- round(object$BUGSoutput$summary[,8], digits = 2)
      x <- cbind(round(object$BUGSoutput$summary[,c(1, 2, 3, 5, 7)], digits = digits), rh)
      colnames(x)[6] <- 'Rhat'
    }else{
      x <- cbind(round(object$BUGSoutput$summary[,c(1, 2, 3, 5, 7)], digits = digits))
    }

    if(n.eff == TRUE)
    {
      nf <- round(object$BUGSoutput$summary[,9], digits = 0)
      x2 <- cbind(x, nf)
      colnames(x2)[ncol(x2)] <- 'n.eff'
    }else{
      x2 <- x
    }

    x2 <- x2[f_ind, , drop = FALSE]

    if(!is.null(func))
    {
      if(length(f_ind) > 1)
      {
        ch_bind <- object$BUGSoutput$sims.matrix[,f_ind]
      }else{
        ch_bind <- as.matrix(object$BUGSoutput$sims.matrix[,f_ind], ncol = 1)
      }

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

  }else{

    if(coda::is.mcmc.list(object2) == TRUE)
    {
      if(length(f_ind) > 1)
      {
        dsort <- do.call(coda::mcmc.list, temp[,f_ind])
        ch_bind <- do.call('rbind', dsort)
      }else{
        dsort <- do.call(coda::mcmc.list, temp[,f_ind])
        ch_bind <- as.matrix(do.call(coda::mcmc.list, temp[,f_ind]), ncol = 1)
      }

      bind_mn <- round(apply(ch_bind, 2, mean), digits = digits)
      bind_sd <- round(apply(ch_bind, 2, stats::sd), digits = digits)
      bind_LCI <- round(apply(ch_bind, 2, stats::quantile, probs= 0.025), digits = digits)
      bind_med <- round(apply(ch_bind,2, stats::median), digits = digits)
      bind_UCI <- round(apply(ch_bind, 2, stats::quantile, probs= 0.975), digits = digits)

      if(Rhat == TRUE)
      {
        if(length(dsort) > 1)
        {
          r_hat <- round(coda::gelman.diag(dsort, multivariate = FALSE)$psrf[,1], digits = 2)
        }else{
          warning('Rhat statistic cannot be calculated with one chain. NAs inserted.')
          r_hat <- rep(NA, NCOL(dsort))
        }
        x <- cbind(bind_mn, bind_sd, bind_LCI, bind_med, bind_UCI, r_hat)
        colnames(x) <- c('mean', 'sd', '2.5%','50%','97.5%', 'Rhat')
      }else{
        x <- cbind(bind_mn, bind_sd, bind_LCI, bind_med, bind_UCI)
        colnames(x) <- c('mean', 'sd', '2.5%','50%','97.5%')
      }

      if(n.eff == TRUE)
      {
        bind_neff <- round(coda::effectiveSize(dsort), digits = 0)
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

      #needed when length(f_ind) == 1
      rownames(x3) <- colnames(temp[[1]])[f_ind]
      mcmc_summary <- x3
    }

    if(typeof(object2) == 'double')
    {
      if(length(f_ind) > 1)
      {
        dsort <- temp[,f_ind]
      }else{
        dsort <- as.matrix(temp[,f_ind], ncol = 1)
      }

      bind_mn <- round(apply(dsort, 2, mean), digits = digits)
      bind_sd <- round(apply(dsort, 2, stats::sd), digits = digits)
      bind_LCI <- round(apply(dsort, 2, stats::quantile, probs= 0.025), digits = digits)
      bind_med <- round(apply(dsort,2, stats::median), digits = digits)
      bind_UCI <- round(apply(dsort, 2, stats::quantile, probs= 0.975), digits = digits)

      if(Rhat == TRUE)
      {
        warning('Rhat statistic cannot be calculated without individual chains. NAs inserted.')
        r_hat <- rep(NA, NCOL(dsort))
        x <- cbind(bind_mn, bind_sd, bind_LCI, bind_med, bind_UCI, r_hat)
        colnames(x) <- c('mean', 'sd', '2.5%','50%','97.5%', 'Rhat')
      }else{
        x <- cbind(bind_mn, bind_sd, bind_LCI, bind_med, bind_UCI)
        colnames(x) <- c('mean', 'sd', '2.5%','50%','97.5%')
      }

      if(n.eff == TRUE)
      {
        warning('Number of effective samples cannot be calculated without individual chains. NAs inserted.')
        bind_neff <- rep(NA, NCOL(dsort))
        x2 <- cbind(x, bind_neff)
        colnames(x2)[ncol(x2)] <- 'n.eff'
      }else{
        x2 <- x
      }

      if(!is.null(func))
      {
        tmp <- round(apply(dsort, 2, func), digits = digits)

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

      #needed when length(f_ind) == 1
      rownames(x3) <- colnames(temp)[f_ind]
      mcmc_summary <- x3
    }
  }
  return(mcmc_summary)
}
