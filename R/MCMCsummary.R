#' Summarize MCMC output
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
#' @param digits Number of significant digits to include for posterior summary. All computed digits will be included by default. Note that Rhat is always rounded to 2 decimal places.
#'
#' @param round Number of decimal places to to round to for posterior summary. Cannot be used in conjunction with \code{digits} argument Note that Rhat is always rounded to 2 decimal places.
#'
#' @param Rhat Logical specifying whether to calculate and display the potential scale reduction statistic (Rhat). Values near 1 suggest convergence (Brooks and Gelman 1998). \code{Rhat = FALSE} will prevent display of this column in summary output. Specifying \code{Rhat = FALSE}, may increase function speed for very large \code{mcmc.list} objects.
#'
#' @param n.eff Logical specifying whether to calculate and display the number of effective samples for each parameter. \code{n.eff = FALSE} will prevent display of this column in summary output. Specifying \code{n.eff = FALSE}, may increase function speed for very large \code{mcmc.list} objects.
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
#' For \code{mcmc.list} objects, the potential scale reduction statistic statistic (Rhat) is calculated using the \code{gelman.diag} function in the \code{coda} package (what is typically displayed in the summary output from models fit with JAGS). For \code{stanfit} objects, Rhat is calculated using the \code{rstan} package which computes a 'split chain' Rhat, which is thought to be a more conservative diagnostic (Stan Development Team 2018).
#'
#' For \code{mcmc.list} objects, the number of effective samples is calculated using the \code{effectiveSize} function in the \code{coda} package. For \code{stanfit} objects, Rhat is calculated using the \code{rstan} package which (in a similar way to the Rhat computation noted above) employs a slightly different (and more conservative) method of computation for the number of effective samples (Stan Development Team 2018).
#'
#' @return Function returns summary information (including parameter posterior mean, posterior sd, 2.5\% quantile, median, 97.5\% quantile, potential scale reduction statistic statistic (Rhat), number of effective samples, and other specified metrics) for specified parameters.
#'
#'
#' @section References:
#'
#' Brooks, S. P., and A. Gelman. 1998. General methods for monitoring convergence of iterative simulations. Journal of Computational and Graphical Statistics 7:434.
#'
#' Stan Development Team. 2018. Stan Modeling Language Users Guide and Reference Manual, Version 2.18.0.   http://mc-stan.org
#'
#' @examples
#' #Load data
#' data(MCMC_data)
#'
#' #Summary information for MCMC output - display 2 significant digits
#' MCMCsummary(MCMC_data, digits = 2)
#'
#' #Just 'beta' parameters - round to 2 decimal places
#' MCMCsummary(MCMC_data, params = 'beta', round = 2)
#'
#' #Just 'beta[1]', 'beta[4]', and 'alpha[3]'
#' #'params' takes regular expressions when ISB = FALSE, square brackets must be escaped with '\\'
#' MCMCsummary(MCMC_data, params = c('beta\\[1\\]', 'beta\\[4\\]', 'alpha\\[3\\]'), ISB = FALSE)
#'
#' @export

MCMCsummary <- function(object,
                      params = 'all',
                      excl = NULL,
                      ISB = TRUE,
                      digits = NULL,
                      round = NULL,
                      Rhat = TRUE,
                      n.eff = FALSE,
                      func = NULL,
                      func_name = NULL)
{
  #SORTING BLOCK
  if (typeof(object) == 'double')
  {
    object2 <- MCMCchains(object, params, excl, ISB, mcmc.list = FALSE)
  } else {
    if (typeof(object) == 'S4')
    {
      object2 <- object
    } else {
      object2 <- MCMCchains(object, params, excl, ISB, mcmc.list = TRUE)
    }
  }

  #PROCESSING BLOCK
  if (coda::is.mcmc.list(object2) == TRUE)
  {
    np <- NCOL(object2[[1]])
    if (np > 1)
    {
      ch_bind <- do.call('rbind', object2)
    } else{
      ch_bind <- as.matrix(object2)
    }
    if (!is.null(digits))
    {
      if (!is.null(round))
      {
        warning("'digits' and 'round' arguments cannot be used together. Using 'digits'.")
      }
      bind_mn <- signif(apply(ch_bind, 2, mean), digits = digits)
      bind_sd <- signif(apply(ch_bind, 2, stats::sd), digits = digits)
      bind_LCI <- signif(apply(ch_bind, 2, stats::quantile, probs= 0.025), digits = digits)
      bind_med <- signif(apply(ch_bind, 2, stats::median), digits = digits)
      bind_UCI <- signif(apply(ch_bind, 2, stats::quantile, probs= 0.975), digits = digits)
    }
    if (is.null(digits) & !is.null(round))
    {
      bind_mn <- round(apply(ch_bind, 2, mean), digits = round)
      bind_sd <- round(apply(ch_bind, 2, stats::sd), digits = round)
      bind_LCI <- round(apply(ch_bind, 2, stats::quantile, probs= 0.025), digits = round)
      bind_med <- round(apply(ch_bind, 2, stats::median), digits = round)
      bind_UCI <- round(apply(ch_bind, 2, stats::quantile, probs= 0.975), digits = round)
    }
    if (is.null(digits) & is.null(round))
    {
      bind_mn <- apply(ch_bind, 2, mean)
      bind_sd <- apply(ch_bind, 2, stats::sd)
      bind_LCI <- apply(ch_bind, 2, stats::quantile, probs= 0.025)
      bind_med <- apply(ch_bind, 2, stats::median)
      bind_UCI <- apply(ch_bind, 2, stats::quantile, probs= 0.975)
    }
    if (Rhat == TRUE)
    {
      if (length(object2) > 1)
      {
        #If > 750 params use loop to calc Rhat
        if (NCOL(object2[[1]]) > 750)
        {
          r_hat <- c(rep(NA, NCOL(object2[[1]])))
          for (v in 1:length(r_hat))
          {
            r_hat[v] <- round(coda::gelman.diag(object2[,v])$psrf[,1], digits = 2)
          }
        } else{
          r_hat <- round(coda::gelman.diag(object2, multivariate = FALSE)$psrf[,1], digits = 2)
        }
      } else{
        warning('Rhat statistic cannot be calculated with one chain. NAs inserted.')
        r_hat <- rep(NA, NCOL(object2))
      }
      x <- cbind(bind_mn, bind_sd, bind_LCI, bind_med, bind_UCI, r_hat)
      colnames(x) <- c('mean', 'sd', '2.5%','50%','97.5%', 'Rhat')
    } else {
      x <- cbind(bind_mn, bind_sd, bind_LCI, bind_med, bind_UCI)
      colnames(x) <- c('mean', 'sd', '2.5%','50%','97.5%')
    }
    if (n.eff == TRUE)
    {
      bind_neff <- round(coda::effectiveSize(object2), digits = 0)
      x2 <- cbind(x, bind_neff)
      colnames(x2)[ncol(x2)] <- 'n.eff'
    } else {
      x2 <- x
    }
    if (!is.null(func))
    {
      if (!is.null(digits))
      {
        tmp <- signif(apply(ch_bind, 2, func), digits = digits)
      }
      if (is.null(digits) & !is.null(round))
      {
        tmp <- round(apply(ch_bind, 2, func), digits = round)
      }
      if (is.null(digits) & is.null(round))
      {
        tmp <- apply(ch_bind, 2, func)
      }
      if (!is.null(dim(tmp)) & NROW(tmp) > 1)
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
          } else {
            colnames(x3)[ncol(x2) + i] <- 'func'
          }
        }
      } else {
        x3 <- cbind(x2, tmp)
        if (!is.null(func_name))
        {
          if (length(func_name) > 1)
          {
            stop('length(func_name) must equal number of func outputs')
          }
          colnames(x3)[ncol(x3)] <- func_name
        } else {
          colnames(x3)[ncol(x3)] <- 'func'
        }
      }
    } else {
      x3 <- x2
    }
    mcmc_summary <- x3
  }
  
  if (typeof(object2) == 'double')
  {
    np <- NCOL(object2)
    ch_bind <- object2
    if (!is.null(digits))
    {
      if (!is.null(round))
      {
        warning("'digits' and 'round' arguments cannot be used together. Using 'digits'.")
      }
      bind_mn <- signif(apply(ch_bind, 2, mean), digits = digits)
      bind_sd <- signif(apply(ch_bind, 2, stats::sd), digits = digits)
      bind_LCI <- signif(apply(ch_bind, 2, stats::quantile, probs= 0.025), digits = digits)
      bind_med <- signif(apply(ch_bind, 2, stats::median), digits = digits)
      bind_UCI <- signif(apply(ch_bind, 2, stats::quantile, probs= 0.975), digits = digits)
    }
    if (is.null(digits) & !is.null(round))
    {
      bind_mn <- round(apply(ch_bind, 2, mean), digits = round)
      bind_sd <- round(apply(ch_bind, 2, stats::sd), digits = round)
      bind_LCI <- round(apply(ch_bind, 2, stats::quantile, probs= 0.025), digits = round)
      bind_med <- round(apply(ch_bind, 2, stats::median), digits = round)
      bind_UCI <- round(apply(ch_bind, 2, stats::quantile, probs= 0.975), digits = round)
    }
    if (is.null(digits) & is.null(round))
    {
      bind_mn <- apply(ch_bind, 2, mean)
      bind_sd <- apply(ch_bind, 2, stats::sd)
      bind_LCI <- apply(ch_bind, 2, stats::quantile, probs= 0.025)
      bind_med <- apply(ch_bind, 2, stats::median)
      bind_UCI <- apply(ch_bind, 2, stats::quantile, probs= 0.975)
    }
    if (Rhat == TRUE)
    {
      warning('Rhat statistic cannot be calculated without individual chains. NAs inserted.')
      r_hat <- rep(NA, np)
      x <- cbind(bind_mn, bind_sd, bind_LCI, bind_med, bind_UCI, r_hat)
      colnames(x) <- c('mean', 'sd', '2.5%','50%','97.5%', 'Rhat')
    } else {
      x <- cbind(bind_mn, bind_sd, bind_LCI, bind_med, bind_UCI)
      colnames(x) <- c('mean', 'sd', '2.5%','50%','97.5%')
    }
    if (n.eff == TRUE)
    {
      warning('Number of effective samples cannot be calculated without individual chains. NAs insrted.')
      bind_neff <- rep(NA, np)
      x2 <- cbind(x, bind_neff)
      colnames(x2)[ncol(x2)] <- 'n.eff'
    } else {
      x2 <- x
    }
    if (!is.null(func))
    {
      if (!is.null(digits))
      {
        tmp <- signif(apply(ch_bind, 2, func), digits = digits)
      }
      if (is.null(digits) & !is.null(round))
      {
        tmp <- round(apply(ch_bind, 2, func), digits = round)
      }
      if (is.null(digits) & is.null(round))
      {
        tmp <- apply(ch_bind, 2, func)
      }
      if (!is.null(dim(tmp)) & NROW(tmp) > 1)
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
          } else {
            colnames(x3)[ncol(x2) + i] <- 'func'
          }
        }
      } else {
        x3 <- cbind(x2, tmp)
        if (!is.null(func_name))
        {
          if (length(func_name) > 1)
          {
            stop('length(func_name) must equal number of func outputs')
          }
          colnames(x3)[ncol(x3)] <- func_name
        } else {
          colnames(x3)[ncol(x3)] <- 'func'
        }
      }
    } else {
      x3 <- x2
    }
    mcmc_summary <- x3
  }
  
  if (typeof(object2) == 'S4')
  {
    #rhat and n_eff directly from rstan output
    all_params <- row.names(data.frame(rstan::summary(object2)$summary)['n_eff'])
    rs_df <- data.frame(rstan::summary(object2)$summary)
    
    #filtering of parameters from rstan object - from MCMCchains
    if (ISB == TRUE)
    {
      names <- vapply(strsplit(all_params,
                               split = "[", fixed = TRUE), `[`, 1, FUN.VALUE=character(1))
    } else {
      names <- all_params
    }
    
    #INDEX BLOCK
    #exclusions
    if (!is.null(excl))
    {
      rm_ind <- c()
      for (i in 1:length(excl))
      {
        if (ISB == TRUE)
        {
          n_excl <- vapply(strsplit(excl,
                                    split = "[", fixed = TRUE), `[`, 1, FUN.VALUE=character(1))
          ind_excl <- which(names %in% n_excl[i])
          if (length(ind_excl) < 1)
          {
            warning(paste0('"', excl[i], '"', ' not found in MCMC output.'))
          }
          rm_ind <- c(rm_ind, ind_excl)
        } else {
          n_excl <- excl
          ind_excl <- grep(n_excl[i], names, fixed = FALSE)
          if (length(ind_excl) < 1)
          {
            warning(paste0('"', excl[i], '"', ' not found in MCMC output.'))
          }
          rm_ind <- c(rm_ind, ind_excl)
        }
      }
      if (length(rm_ind) > 0)
      {
        dups <- which(duplicated(rm_ind))
        if (length(dups) > 0)
        {
          rm_ind2 <- rm_ind[-dups]
        } else {
          rm_ind2 <- rm_ind
        }
      } else {
        excl <- NULL
      }
    }
    
    #selections
    if (length(params) == 1)
    {
      if (params == 'all')
      {
        if (is.null(excl))
        {
          f_ind <- 1:length(names)
        } else {
          f_ind <- (1:length(names))[-rm_ind2]
        }
      } else {
        if (ISB == TRUE)
        {
          get_ind <- which(names %in% params)
        } else {
          get_ind <- grep(paste(params), names, fixed = FALSE)
        }
        
        if (length(get_ind) < 1)
        {
          stop(paste0('"', params, '"', ' not found in MCMC output.'))
        }
        if (!is.null(excl))
        {
          if (identical(get_ind, rm_ind2))
          {
            stop('No parameters selected.')
          }
          matched <- stats::na.omit(match(rm_ind2, get_ind))
          if (length(matched) > 0)
          {
            f_ind <- get_ind[-matched]
          } else {
            f_ind <- get_ind
          }
        } else {
          f_ind <- get_ind
        }
      }
    } else {
      grouped <- c()
      for (i in 1:length(params))
      {
        if (ISB == TRUE)
        {
          get_ind <- which(names %in% params[i])
        } else {
          get_ind <- grep(paste(params[i]), names, fixed=FALSE)
        }
        
        if (length(get_ind) < 1)
        {
          warning(paste0('"', params[i], '"', ' not found in MCMC output.'))
          next()
        }
        grouped <- c(grouped, get_ind)
      }
      if (!is.null(excl))
      {
        if (identical(grouped, rm_ind2))
        {
          stop('No parameters selected.')
        }
        matched <- stats::na.omit(match(rm_ind2, grouped))
        if (length(matched) > 0)
        {
          t_ind <- grouped[-matched]
        } else {
          t_ind <- grouped
        }
        to.rm <- which(duplicated(t_ind))
        if (length(to.rm) > 0)
        {
          f_ind <- t_ind[-to.rm]
        }else
        {
          f_ind <- t_ind
        }
      } else {
        to.rm <- which(duplicated(grouped))
        if (length(to.rm) > 0)
        {
          f_ind <- grouped[-to.rm]
        } else {
          f_ind <- grouped
        }
      }
    }
    #end sort
    
    if (!is.null(digits))
    {
      if (!is.null(round))
      {
        warning("'digits' and 'round' arguments cannot be used together. Using 'digits'.")
      }
      bind_mn <- signif(rs_df['mean'][f_ind, 1], digits = digits)
      bind_sd <- signif(rs_df['sd'][f_ind, 1], digits = digits)
      bind_LCI <- signif(rs_df['X2.5.'][f_ind, 1], digits = digits)
      bind_med <- signif(rs_df['X50.'][f_ind, 1], digits = digits)
      bind_UCI <- signif(rs_df['X97.5.'][f_ind, 1], digits = digits)
    }
    if (is.null(digits) & !is.null(round))
    {
      bind_mn <- round(rs_df['mean'][f_ind, 1], digits = round)
      bind_sd <- round(rs_df['sd'][f_ind, 1], digits = round)
      bind_LCI <- round(rs_df['X2.5.'][f_ind, 1], digits = round)
      bind_med <- round(rs_df['X50.'][f_ind, 1], digits = round)
      bind_UCI <- round(rs_df['X97.5.'][f_ind, 1], digits = round)
    }
    if (is.null(digits) & is.null(round))
    {
      bind_mn <- rs_df['mean'][f_ind, 1]
      bind_sd <- rs_df['sd'][f_ind, 1]
      bind_LCI <- rs_df['X2.5.'][f_ind, 1]
      bind_med <- rs_df['X50.'][f_ind, 1]
      bind_UCI <- rs_df['X97.5.'][f_ind, 1]
    }
    
    if (Rhat == TRUE)
    {
      if (dim(rstan::summary(object)$c_summary)[3] > 1)
      {
        r_hat <- round(rs_df['Rhat'][f_ind,1], digits = 2)
      } else {
        warning('Rhat statistic cannot be calculated with one chain. NAs inserted.')
        r_hat <- rep(NA, length(f_ind))
      }
      x <- cbind(bind_mn, bind_sd, bind_LCI, bind_med, bind_UCI, r_hat)
      colnames(x) <- c('mean', 'sd', '2.5%','50%','97.5%', 'Rhat')
      row.names(x) <- all_params[f_ind]
    } else {
      x <- cbind(bind_mn, bind_sd, bind_LCI, bind_med, bind_UCI)
      colnames(x) <- c('mean', 'sd', '2.5%','50%','97.5%')
      row.names(x) <- all_params[f_ind]
    }
    if (n.eff == TRUE)
    {
      neff <- round(rs_df['n_eff'][f_ind, 1], digits = 0)
      x2 <- cbind(x, neff)
      colnames(x2)[ncol(x2)] <- 'n.eff'
    } else {
      x2 <- x
    }
    
    if(!is.null(func))
    {
      #convert stan object to matrix
      ch_bind <- as.matrix(object)[,f_ind]
      
      if (!is.null(digits))
      {
        tmp <- signif(apply(ch_bind, 2, func), digits = digits)
      }
      if (is.null(digits) & !is.null(round))
      {
        tmp <- round(apply(ch_bind, 2, func), digits = round)
      }
      if (is.null(digits) & is.null(round))
      {
        tmp <- apply(ch_bind, 2, func)
      }
      
      if (!is.null(dim(tmp)) & NROW(tmp) > 1)
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
      } else {
        x3 <- cbind(x2, tmp)
        if (!is.null(func_name))
        {
          if (length(func_name) > 1)
          {
            stop('length(func_name) must equal number of func outputs')
          }
          colnames(x3)[ncol(x3)] <- func_name
        } else {
          colnames(x3)[ncol(x3)] <- 'func'
        }
      }
    } else {
      x3 <- x2
    }
    mcmc_summary <- x3
  }
  return(mcmc_summary)
}
