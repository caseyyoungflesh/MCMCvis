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
#' @param ISB Ignore Square Brackets (ISB). Logical specifying whether square brackets should be ignored in the \code{params} and \code{excl} arguments. If \code{TRUE}, square brackets are ignored. If \code{FALSE}, square brackets are not ignored.  This allows partial names to be used when specifying parameters of interest. Use \code{exact} argument to specify whether input from \code{params} and \code{excl} arguments should be matched exactly.
#'
#' @param exact Logical specifying whether input from \code{params} and \code{excl} arguments should be matched exactly (after ignoring square brackets if \code{ISB = FALSE}). If \code{TRUE}, input from \code{params} and \code{excl} are matched exactly (after taking \code{ISB} argument into account). If \code{FALSE}, input from \code{params} and \code{excl} are matched using regular expression format (after taking \code{ISB} argument into account).
#'
#' @param probs Numeric vector where each element in (0,1) representing probabilities used to calculate posterior sample quantiles for the selected parameters. Default is c(0.025, 0.5, 0.975).
#' 
#' @param hpd_prob Scalar in (0,1) representing probability used to calculate highest posterior density intervals for the selected parameters. Default is 0.95.
#' 
#' @param HPD Logical specifying whether to calculate equal-tailed credible intervals (\code{HPD = FALSE}) or highest posterior density intervals (\code{HPD = TRUE}) for the selected parameters. Default is \code{HPD = FALSE}.
#' 
#' @param pg0 Logical specifying whether to calculate the proportion of the posterior that is greater than 0, rounded to 2 digits.
#' 
#' @param digits Number of significant digits to include for posterior summary. All computed digits will be included by default. Note that Rhat is always rounded to 2 decimal places.
#'
#' @param round Number of decimal places to round to for posterior summary. Cannot be used in conjunction with \code{digits} argument. Note that Rhat is always rounded to 2 decimal places.
#'
#' @param Rhat Logical specifying whether to calculate and display the potential scale reduction statistic (Rhat). Values near 1 suggest convergence (Brooks and Gelman 1998). \code{Rhat = FALSE} will prevent display of this column in summary output. Specifying \code{Rhat = FALSE}, may increase function speed for very large \code{mcmc.list} objects.
#'
#' @param n.eff Logical specifying whether to calculate and display the number of effective samples for each parameter. \code{n.eff = FALSE} will prevent display of this column in summary output. Specifying \code{n.eff = FALSE}, may increase function speed for very large \code{mcmc.list} objects. Default is \code{n.eff = TRUE}.
#'
#' @param func Function to be performed on MCMC output. If a function is specified, it will be evaluated on posteriors for each specified parameter and returned as a column in the summary output (or multiple columns if the function returns more than one value).
#'
#' @param func_name Character string (or vector of character strings) specifying labels for output from \code{func} argument. If \code{func_name} is not specified, columns with \code{func} argument will be labeled 'func'.
#'
#' @section Details:
#' \code{object} argument can be a \code{stanfit} object (\code{rstan} package), a \code{CmdStanMCMC} object (\code{cmdstanr} package), a \code{stanreg} object (\code{rstanarm} package), a \code{brmsfit} object (\code{brms} package), an \code{mcmc.list} object (\code{coda} and \code{rjags} packages), \code{mcmc} object (\code{coda} and \code{nimble} packages), \code{list} object (\code{nimble} package), an \code{R2jags} model object (\code{R2jags} package), a \code{jagsUI} model object (\code{jagsUI} package), or a matrix containing MCMC chains (each column representing MCMC output for a single parameter, rows representing iterations in the chain). The function automatically detects the object type and proceeds accordingly.

#'
#' @section Notes:
#'
#' For \code{mcmc.list}, \code{mcmc}, and \code{list} objects, the potential scale reduction statistic statistic (Rhat) is calculated using the \code{gelman.diag} function in the \code{coda} package (what is typically displayed in the summary output from models fit with JAGS). For \code{stanfit} (as well as \code{CmdStanMCMC}, \code{stanreg}, and \code{brmsfit} objects) and \code{jagsUI} objects, Rhat is calculated using a 'split chain' Rhat (in their respective packages), which is thought to be a more conservative diagnostic (Stan Development Team 2018).
#'
#' For \code{mcmc.list}, \code{mcmc}, and \code{list} objects, the number of effective samples is calculated using the \code{effectiveSize} function in the \code{coda} package. For \code{stanfit} (as well as \code{CmdStanMCMC}, \code{stanreg}, and \code{brmsfit} objects) and \code{jagsUI} objects, n.eff is calculated using a slightly different method of computation for the number of effective samples (Stan Development Team 2018). For \code{CmdStanMCMC} objects, both bulk and tail n.eff is calculated.
#'
#' @return Function returns summary information (including parameter posterior mean, posterior sd, quantiles, potential scale reduction statistic (Rhat), number of effective samples, and other specified metrics) for specified parameters.
#'
#'
#' @section References:
#'
#' Brooks, S. P., and A. Gelman. 1998. General methods for monitoring convergence of iterative simulations. Journal of Computational and Graphical Statistics 7:434.
#'
#' Stan Development Team. 2018. Stan Modeling Language Users Guide and Reference Manual, Version 2.18.0.   https://mc-stan.org
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
#' MCMCsummary(MCMC_data, params = c('beta[1]', 'beta[4]', 'alpha[3]'), ISB = FALSE, exact = TRUE)
#'
#' @export

MCMCsummary <- function(object,
                        params = 'all',
                        excl = NULL,
                        ISB = TRUE,
                        exact = TRUE,
                        probs = c(0.025, 0.5, 0.975),
                        hpd_prob = 0.95,
                        HPD = FALSE,
                        pg0 = FALSE,
                        digits = NULL,
                        round = NULL,
                        Rhat = TRUE,
                        n.eff = TRUE,
                        func = NULL,
                        func_name = NULL)
{ 
  # SORTING BLOCK

  if (methods::is(object, 'matrix'))
  {
    object2 <- MCMCchains(object, params, excl, ISB, exact = exact, mcmc.list = FALSE)
  } else {
    if (methods::is(object, 'stanfit') | methods::is(object, 'CmdStanMCMC'))
    {
      object2 <- object
    } else {
      # rstanarm
      if (methods::is(object, 'stanreg'))
      {
        object2 <- object$stanfit
      } else {
        # brms
        if (methods::is(object, 'brmsfit'))
        {
          object2 <- object$fit
        } else {
          #jagsUI
          if (methods::is(object, 'jagsUI'))
          { 
            object2 <- MCMCchains(object)
          } else {
            object2 <- MCMCchains(object, params, excl, ISB, exact = exact, mcmc.list = TRUE)
          }
        }
      }
    }
  }

#--------------------------------------------------------------------------------------------------------------           

# PROCESSING BLOCK - JAGS AND MATRIX MCMC OUTPUT

  if (coda::is.mcmc.list(object2) == TRUE | methods::is(object, 'matrix'))
  {
    if (methods::is(object, 'matrix'))
    {
      np <- NCOL(object2)
      ch_bind <- object2
    } else {
      np <- NCOL(object2[[1]])
      if (np > 1) ch_bind <- do.call("rbind", object2) else ch_bind <- as.matrix(object2)
    }

    x <- list()
  
# mean, sd, and quantiles  
    
    if (!is.null(digits)) 
    {
      if (!is.null(round)) 
      {
        warning("'digits' and 'round' arguments cannot be used together. Using 'digits'.") 
      }
      
      bind_mn <- data.frame(signif(apply(ch_bind, 2, mean), digits = digits))
      bind_sd <- data.frame(signif(apply(ch_bind, 2, stats::sd), digits = digits))
      colnames(bind_mn) <- "mean"  
      colnames(bind_sd) <- "sd"  
      
      if (HPD == FALSE)
      {
        if (length(probs)==1)
        {
            bind_q <- data.frame(signif(apply(ch_bind, 2, stats::quantile, probs = probs), digits = digits))
            colnames(bind_q) <-  paste0(signif(probs * 100, digits = 3), "%")                   
        } else { 
          bind_q <- data.frame(t(signif(apply(ch_bind, 2, stats::quantile, probs = probs), digits = digits)))
          colnames(bind_q) <-  paste0(signif(probs * 100, digits = 3), "%")                    
        }
      }
      if (HPD == TRUE)
      {
        if (length(hpd_prob) > 1)
        {
          stop('specify only a single probability for HPD interval computation.')
        }
        bind_q <- data.frame(signif(coda::HPDinterval(coda::as.mcmc(ch_bind), prob = hpd_prob), digits = digits))
        colnames(bind_q) <- c(paste0(signif(hpd_prob * 100, digits = 3), "%_HPDL"), paste0(signif(hpd_prob * 100, digits = 3), "%_HPDU"))  
      }
    }
    
    if (is.null(digits) & !is.null(round))
    {
      bind_mn <- data.frame(round(apply(ch_bind, 2, mean), digits = round))
      bind_sd <- data.frame(round(apply(ch_bind, 2, stats::sd), digits = round))
      colnames(bind_mn) <- "mean"  
      colnames(bind_sd) <- "sd"  
     
      if (HPD == FALSE)
      {
        if (length(probs)==1)
        {
          bind_q <- data.frame(round(apply(ch_bind, 2, stats::quantile, probs = probs), digits = round))
          colnames(bind_q) <-  paste0(signif(probs * 100, digits = 3), "%")                    
        } else { 
          bind_q <- data.frame(t(round(apply(ch_bind, 2, stats::quantile, probs = probs), digits = round)))
          colnames(bind_q) <-  paste0(signif(probs * 100, digits = 3), "%")
        }
      }
      if (HPD == TRUE)
      {
        if (length(hpd_prob) > 1)
        {
          stop('specify only a single probability for HPD interval computation.')
        }
        bind_q <- data.frame(round(coda::HPDinterval(coda::as.mcmc(ch_bind), prob = hpd_prob), digits = round))
        colnames(bind_q) <- c(paste0(signif(hpd_prob * 100, digits = 3), "%_HPDL"), paste0(signif(hpd_prob * 100, digits = 3), "%_HPDU"))  
      }
    }
    
    if (is.null(digits) & is.null(round))
    {
      bind_mn <- data.frame(apply(ch_bind, 2, mean))
      bind_sd <- data.frame(apply(ch_bind, 2, stats::sd))
      colnames(bind_mn) <- "mean"  
      colnames(bind_sd) <- "sd"  
      
      if (HPD == FALSE)
      {
        if (length(probs) == 1)
        {
          bind_q <- data.frame(apply(ch_bind, 2, stats::quantile, probs = probs))
          colnames(bind_q) <-  paste0(signif(probs * 100, digits =3), "%")
        } else { 
          bind_q <- data.frame(t(apply(ch_bind, 2, stats::quantile, probs = probs)))
          colnames(bind_q) <-  paste0(signif(probs * 100, digits = 3), "%")
        }
      }
      if (HPD == TRUE)
      {
        if (length(hpd_prob) > 1)
        {
          stop('specify only a single probability for HPD interval computation.')
        }
        bind_q <- data.frame(coda::HPDinterval(coda::as.mcmc(ch_bind), prob = hpd_prob))
        colnames(bind_q) <- c(paste0(signif(hpd_prob * 100, digits = 3), "%_HPDL"), paste0(signif(hpd_prob * 100, digits = 3), "%_HPDU"))  
      }
    }
    
    x[[1]] <- cbind(bind_mn, bind_sd, bind_q) 
  
# rhat 

    if (Rhat == TRUE)
    {
      if (!methods::is(object, 'matrix'))
      {
        if (length(object2) > 1)
        {
          # If > 750 params use loop to calculate Rhat
          if (NCOL(object2[[1]]) > 750)
          {
            r_hat <- c(rep(NA, NCOL(object2[[1]])))
            for (v in 1:length(r_hat)) r_hat[v] <- round(coda::gelman.diag(object2[, v])$psrf[, 1], digits = 2)
            r_hat <- data.frame(r_hat)
            colnames(r_hat) <- "Rhat"
          } else { 
            r_hat <- data.frame(round(coda::gelman.diag(object2, multivariate = FALSE)$psrf[, 1], digits = 2))
            colnames(r_hat) <- "Rhat"
          } 
        } else {
          warning("Rhat statistic cannot be calculated with one chain. NAs inserted.")
          r_hat <- data.frame(rep(NA, np))
          colnames(r_hat) <- "Rhat"
        }
      } else {
        warning("Rhat statistic cannot be calculated with one chain (matrix input). NAs inserted.")
        r_hat <- data.frame(rep(NA, np))
        colnames(r_hat) <- "Rhat"
      }
      x[[(length(x) + 1)]] <- r_hat  
    }

# neff
   
    if (n.eff == TRUE)
    {
      if (!methods::is(object, 'matrix'))
      {
        neff <- data.frame(round(coda::effectiveSize(object2), digits = 0))
        colnames(neff) <- "n.eff"
      } else {
        warning('Number of effective samples cannot be calculated without individual chains (matrix input). NAs inserted.')
        neff <- data.frame(rep(NA, np))
        colnames(neff) <- "n.eff"
      } 
      x[[(length(x) + 1)]] <- neff
    }

# p>0
    
    if (pg0 == TRUE)
    {
      tpg <- data.frame(apply(ch_bind, 2, function(x) round(sum(x > 0) / length(x), 2)))
      colnames(tpg) <- 'p>0'
      x[[(length(x) + 1)]] <- tpg
    }
    
  
# custom function
  
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
      if (!is.null(dim(tmp)))
      { 
        tmp <- data.frame(t(tmp)) 
      } else { 
        tmp <- data.frame(tmp)
      }
      if (!is.null(func_name))
      {
        if (length(func_name) != NCOL(tmp))
        {
          stop("length(func_name) must equal number of func outputs")
        }  
        colnames(tmp) <- func_name
      } else {
        colnames(tmp) <- 'func'  
      }
      
      x[[(length(x) + 1)]] <- tmp
    }
  
# bind them
    mcmc_summary <- do.call("cbind", x)
  }
  
#--------------------------------------------------------------------------------------------------------------                        
# PROCESSING BLOCK - STAN OR CMDSTAN OR JAGSUI MCMC OUTPUT
  
  if (methods::is(object2, 'stanfit') | 
      methods::is(object, 'jagsUI') | 
      methods::is(object2, 'CmdStanMCMC'))
  {
    if (methods::is(object2, 'stanfit'))
    {
      # rhat and n_eff directly from rstan output
      all_params <- row.names(rstan::summary(object2)$summary)
      rs_df <- data.frame(rstan::summary(object2)$summary)
    
      #if brms, reassign names without b_ and r_ (as in MCMCchains)
      if (methods::is(object, 'brmsfit'))
      {
        sp_names_p <- names(object2@sim$samples[[1]])
        #remove b_ and r_
        st_nm <- substr(sp_names_p, start = 1, stop = 2)
        sp_names <- rep(NA, length(sp_names_p))
        b_idx <- which(st_nm == 'b_')
        r_idx <- which(st_nm == 'r_')
        ot_idx <- which(st_nm != 'b_' & st_nm != 'r_')
        #fill names vec with b_ and r_ removed
        sp_names[b_idx] <- gsub('b_', '', sp_names_p[b_idx])
        sp_names[r_idx] <- gsub('r_', '', sp_names_p[r_idx])
        sp_names[ot_idx] <- sp_names_p[ot_idx]
      
        #assign names to df
        all_params <- sp_names
        row.names(rs_df) <- all_params
      }
    }
    
    if (methods::is(object, 'jagsUI'))
    {
      all_params <- row.names(object$summary)
      rs_df <- data.frame(object$summary)
    }
    
    if (methods::is(object2, 'CmdStanMCMC'))
    {
      rs_df <- posterior::summarize_draws(object2)
      all_params <- rs_df$variable
    }
    
    # filtering of parameters from rstan/cmdstan/jagsUI object - from MCMCchains
    if (ISB == TRUE)
    {
      names <- vapply(strsplit(all_params, split = "[", fixed = TRUE), `[`, 1, FUN.VALUE = character(1))
    } else {
      names <- all_params
    }
    
    x <- list()
  
# INDEX BLOCK exclusions
    
    if (!is.null(excl))
    {
      rm_ind <- c()
      for (i in 1:length(excl))
      {
        if (ISB == TRUE)
        {
          n_excl <- vapply(strsplit(excl,
                                    split = "[", fixed = TRUE), `[`, 1, FUN.VALUE=character(1))
        } else {
          n_excl <- excl
        }
        
        if (exact == TRUE)
        {
          ind_excl <- which(names %in% n_excl[i])
        } else {
          ind_excl <- grep(n_excl[i], names, fixed = FALSE)
        }
        
        if (length(ind_excl) < 1)
        {
          warning(paste0("\"", excl[i], "\"", " not found in MCMC output. Check 'ISB'' and 'exact' arguments to make sure the desired parsing methods are being used."))
        }
        rm_ind <- c(rm_ind, ind_excl)
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
    
# selections
    
    if (length(params) == 1)
    {
      if (params == "all")
      {
        if (is.null(excl))
        {
          f_ind <- 1:length(names)
        } else {
          f_ind <- (1:length(names))[-rm_ind2]
        }
      } else {
        if (exact == TRUE)
        {
          get_ind <- which(names %in% params)
        } else {
          get_ind <- grep(paste(params), names, fixed = FALSE)
        }
        
        if (length(get_ind) < 1)
        {
          stop(paste0("\"", params, "\"", " not found in MCMC output. Check 'ISB' and 'exact' arguments to make sure the desired parsing methods are being used."))
        }
        if (!is.null(excl))
        {
          if (identical(get_ind, rm_ind2))
          {
            stop("No parameters selected.")
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
        if (exact == TRUE)
        {
          get_ind <- which(names %in% params[i])
        } else {
          get_ind <- grep(paste(params[i]), names, fixed = FALSE)
        }
        
        if (length(get_ind) < 1)
        {
          warning(paste0("\"", params[i], "\"", " not found in MCMC output. Check 'ISB' and 'exact' arguments to make sure the desired parsing methods are being used."))
          (next)()
        }
        grouped <- c(grouped, get_ind)
      }
      if (!is.null(excl))
      {
        if (identical(grouped, rm_ind2))
        {
          stop("No parameters selected.")
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
        } else {
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

# end sort

# convert object to matrix if computing non default intervals or using custom func
    if (!is.null(func) | HPD == TRUE | 
        identical(probs, c(0.025, 0.5, 0.975)) == FALSE | pg0 == TRUE |
        methods::is(object2, 'CmdStanMCMC'))
    {
      if (methods::is(object2, 'stanfit'))
      {
        #ensure is matrix, not vector
        ch_bind <- as.matrix(as.matrix(object2)[, f_ind])
      }
      if (methods::is(object, 'jagsUI'))
      {
        ch_bind <- MCMCchains(object, params = params, excl = excl, ISB = ISB)
      }
      if (methods::is(object2, 'CmdStanMCMC'))
      {
        ch_bind <- MCMCchains(object2, params = params, ISB = ISB)
      }
    } 

# mean, sd, and quantiles  
    
    if (!is.null(digits))
    {
      if (!is.null(round))
      {
        warning("'digits' and 'round' arguments cannot be used together. Using 'digits'.")
      }
      
      bind_mn <- data.frame(signif(rs_df["mean"][f_ind, 1], digits = digits))
      bind_sd <- data.frame(signif(rs_df["sd"][f_ind, 1], digits = digits))
      colnames(bind_mn) <- "mean"  
      colnames(bind_sd) <- "sd"  
      
      if (HPD == FALSE)
      {
        if (length(probs)==1)
        {
          bind_q <- data.frame(signif(apply(ch_bind, 2, stats::quantile, probs = probs), digits = digits))
          colnames(bind_q) <-  paste0(signif(probs * 100, digits = 3), "%")
        } else {
          if (identical(probs, c(0.025, 0.5, 0.975)) == TRUE & 
              (methods::is(object2, 'stanfit') | 
               methods::is(object, 'jagsUI')))
          {
            bind_LCI <- signif(rs_df["X2.5."][f_ind, 1], digits = digits)
            bind_med <- signif(rs_df["X50."][f_ind, 1], digits = digits)
            bind_UCI <- signif(rs_df["X97.5."][f_ind, 1], digits = digits)
            bind_q <- data.frame(cbind(bind_LCI, bind_med, bind_UCI))
            colnames(bind_q) <-  paste0(signif(probs * 100, digits = 3), "%")
          } else {
            bind_q <- data.frame(t(signif(apply(ch_bind, 2, stats::quantile, probs = probs), 
                                          digits = digits)))
            colnames(bind_q) <-  paste0(signif(probs * 100, digits = 3), "%")                    
          }
        }
      } 
      if (HPD == TRUE) 
      {
        if (length(hpd_prob) > 1)
        {
          stop('Specify only a single probability for HPD interval computation.')
        }
        bind_q <- data.frame(signif(coda::HPDinterval(coda::as.mcmc(ch_bind), prob = hpd_prob), digits = digits))
        colnames(bind_q) <- c(paste0(signif(hpd_prob * 100, digits = 3), "%_HPDL"), paste0(signif(hpd_prob * 100, digits = 3), "%_HPDU"))  
      }
    }
      
    if (is.null(digits) & !is.null(round))
    {
      bind_mn <- data.frame(round(rs_df["mean"][f_ind, 1], digits = round))
      bind_sd <- data.frame(round(rs_df["sd"][f_ind, 1], digits = round))
      colnames(bind_mn) <- "mean"  
      colnames(bind_sd) <- "sd"  
      
      if (HPD == FALSE)
      {
        if (length(probs) == 1)
        {
          bind_q <- data.frame(round(apply(ch_bind, 2, stats::quantile, probs = probs), digits = round))
          colnames(bind_q) <-  paste0(signif(probs * 100, digits = 3), "%")         
        } else {
          if (identical(probs, c(0.025, 0.5, 0.975)) == TRUE & 
              (methods::is(object2, 'stanfit') | 
               methods::is(object, 'jagsUI')))
          {    
            bind_LCI <- round(rs_df["X2.5."][f_ind, 1], digits = round)
            bind_med <- round(rs_df["X50."][f_ind, 1], digits = round)
            bind_UCI <- round(rs_df["X97.5."][f_ind, 1], digits = round)
            bind_q <- data.frame(cbind(bind_LCI, bind_med, bind_UCI))
            colnames(bind_q) <-  paste0(signif(probs * 100, digits = 3), "%")
          } else {
            bind_q <- data.frame(t(round(apply(ch_bind, 2, stats::quantile, probs = probs), 
                                         digits = round)))
            colnames(bind_q) <-  paste0(signif(probs * 100, digits = 3), "%")
          }
        }
      }
      if (HPD == TRUE)
      {
        if (length(hpd_prob) > 1)
        {
          stop('Specify only a single probability for HPD interval computation.')
        }
        bind_q <- data.frame(round(coda::HPDinterval(coda::as.mcmc(ch_bind), prob = hpd_prob), digits = round))
        colnames(bind_q) <- c(paste0(signif(hpd_prob * 100, digits = 3), "%_HPDL"), paste0(signif(hpd_prob * 100, digits = 3), "%_HPDU"))  
      }
    }
    
    if (is.null(digits) & is.null(round))
    {
      bind_mn <- data.frame(rs_df["mean"][f_ind, 1])
      bind_sd <- data.frame(rs_df["sd"][f_ind, 1])
      colnames(bind_mn) <- "mean"  
      colnames(bind_sd) <- "sd"  
      
      if (HPD == FALSE)
      {
        if (length(probs)==1)
        {
          bind_q <- data.frame(apply(ch_bind, 2, stats::quantile, probs = probs))
          colnames(bind_q) <-  paste0(signif(probs * 100, digits = 3), "%")                    
        } else {
          if (identical(probs, c(0.025, 0.5, 0.975)) == TRUE & 
              (methods::is(object2, 'stanfit') | 
               methods::is(object, 'jagsUI')))
          {    
            bind_LCI <- rs_df["X2.5."][f_ind, 1]
            bind_med <- rs_df["X50."][f_ind, 1]
            bind_UCI <- rs_df["X97.5."][f_ind, 1]
            bind_q <- data.frame(cbind(bind_LCI, bind_med, bind_UCI))
            colnames(bind_q) <-  paste0(signif(probs * 100, digits = 3), "%")
          } else {
            bind_q <- data.frame(t(apply(ch_bind, 2, stats::quantile, probs = probs)))
            colnames(bind_q) <-  paste0(signif(probs * 100, digits = 3), "%")
          }
        }
      }
      if (HPD == TRUE)
      {
        if (length(hpd_prob) > 1)
        {
          stop('Specify only a single probability for HPD interval computation.')
        }
        bind_q <- data.frame(coda::HPDinterval(coda::as.mcmc(ch_bind), prob = hpd_prob))
        colnames(bind_q) <- c(paste0(signif(hpd_prob * 100, digits = 3), "%_HPDL"), paste0(signif(hpd_prob * 100, digits = 3), "%_HPDU"))  
      }
    }
    x[[1]] <- cbind(bind_mn, bind_sd, bind_q) 

# rhat - rhat in Stan calculated within chain (different than with coda package)
  
    if (Rhat == TRUE)
    {
      if (methods::is(object2, 'stanfit') | 
           methods::is(object, 'jagsUI'))
      {
        r_hat <- data.frame(round(rs_df["Rhat"][f_ind, 1], digits = 2))
      } else {
        r_hat <- data.frame(round(rs_df["rhat"][f_ind, 1], digits = 2))
      }
      colnames(r_hat) <- "Rhat"
      x[[(length(x) + 1)]] <- r_hat  
    }  
    
# neff - neff in Stan is calculated within chain (different than with coda package)    
    
    if (n.eff == TRUE & 
        (methods::is(object2, 'stanfit') | 
         methods::is(object2, 'jagsUI')))
    {
      if (methods::is(object2, 'stanfit'))
      {
        neff <- data.frame(round(rs_df["n_eff"][f_ind, 1], digits = 0))
      }
      if (methods::is(object, 'jagsUI'))
      {
        neff <- data.frame(round(rs_df["n.eff"][f_ind, 1], digits = 0))
      }
      colnames(neff) <- "n.eff"
      x[[(length(x) + 1)]] <- neff
    }
    
    #both bulk and tail for CmdStan
    if (n.eff == TRUE & 
        methods::is(object2, 'CmdStanMCMC'))
        {
          neff_bulk <- data.frame(round(rs_df["ess_bulk"][f_ind, 1], digits = 0))
          colnames(neff_bulk) <- "n.eff_bulk"
          neff_tail <- data.frame(round(rs_df["ess_tail"][f_ind, 1], digits = 0))
          colnames(neff_tail) <- "n.eff_tail"
          x[[(length(x) + 1)]] <- neff_bulk
          x[[(length(x) + 1)]] <- neff_tail
        }  
        
# p>0
    
    if (pg0 == TRUE)
    {
      tpg <- data.frame(apply(ch_bind, 2, function(x) round(sum(x > 0) / length(x), 2)))
      colnames(tpg) <- 'p>0'
      x[[(length(x) + 1)]] <- tpg
    }
    
# custom function
  
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
      if (!is.null(dim(tmp)))
      { 
        tmp <- data.frame(t(tmp)) 
      } else { 
        tmp <- data.frame(tmp)
      }
      if (!is.null(func_name))
      {
        if (length(func_name) != NCOL(tmp))
        {
          stop("length(func_name) must equal number of func outputs")
        }
        colnames(tmp) <- func_name
      } else {
        colnames(tmp) <- 'func'  
      }
      
      x[[(length(x) + 1)]] <- tmp
    }

# bind them
    
    mcmc_summary <- do.call("cbind", x)
    row.names(mcmc_summary) <- all_params[f_ind]
  }
  return(mcmc_summary)
}
