#' Diagnostics summaries for models
#'
#' Model diagnostics for models. Output includes: time elapsed, max Rhat, min number of effective samples, summary information (output from MCMCsummary - optional), as well as Stan-specific sampler diagnostics (for each chain: adapt delta, max tree depth, initial step size, number of divergences, mean tree depth, mean accept state, mean energy, step size, number of iterations exceeding max tree depth, and whether the low BFMI threshold was reached). More information on Stan-specific diagnostic can be found [here - ref manual] and [here - warnings page]. Max Rhat and min neff are calculated using only the parameters specified.
#' 
#' @param object Object containing MCMC output. See DETAILS below.
#'
#' @param txt Logical whether to print results to a .txt file
#' 
#' @param open_txt Logical - if \code{open_txt = TRUE} .txt file will open in viewer after being generated.
#'
#' @param file_name Name of .txt file to be printed. Default is 'MCMCtxt.txt'.
#'
#' @param wd Working directory for .txt file output and obj output. Default is current directory.
#' 
#' @param summary Logical specifying whether or not to output summary information from MCMCsummary (posterior mean, sd, 2.5th and 97.5th quantiles, Rhat, and n.eff).
#' 
#' @param params Character string (or vector of character strings) denoting parameters to be returned in summary output. Argument is ignored if \code{summary = FALSE}.
#' 
#' #' Default \code{'all'} returns all parameters in summary output.
#'
#' @param excl Character string (or vector of character strings) denoting parameters to exclude. Used in conjunction with \code{params} argument to select parameters of interest. Argument is ignored if \code{summary = FALSE}.
#'
#' @param ISB Ignore Square Brackets (ISB). Logical specifying whether square brackets should be ignored in the \code{params} and \code{excl} arguments. If \code{TRUE}, square brackets are ignored. If \code{FALSE}, square brackets are not ignored.  This allows partial names to be used when specifying parameters of interest. Use \code{exact} argument to specify whether input from \code{params} and \code{excl} arguments should be matched exactly.
#'
#' @param exact Logical specifying whether input from \code{params} and \code{excl} arguments should be matched exactly (after ignoring square brackets if \code{ISB = FALSE}). #' If \code{TRUE}, input from \code{params} and \code{excl} are matched exactly (after taking \code{ISB} argument into account). If \code{FALSE}, input from \code{params} and \code{excl} are matched using regular expression format (after taking \code{ISB} argument into account).
#' 
#' @param digits Number of significant digits to include for posterior summary. All computed digits will be included by default. Note that Rhat is always rounded to 2 decimal places.
#'
#' @param round Number of decimal places to round to for posterior summary. Cannot be used in conjunction with \code{digits} argument. Note that Rhat is always rounded to 2 decimal places.
#'
#' @section Details:
#' Some diagnostic information is only provided for models fit with Stan. This information provides beter understanding of for the NUTS sampler performed during the model run. Many of these diagnostics are not relevant for model fit with samplers using Gibbs or Metropolis-Hastings sampling methods (as used by JAGS).
#'
#' \code{object} argument can be a \code{stanfit} object (\code{rstan} package), a \code{stanreg} object (\code{rstanarm} package), a \code{brmsfit} object (\code{brms} package), an \code{mcmc.list} object (\code{coda} package), an \code{R2jags} model object (\code{R2jags} package), a \code{jagsUI} model object (\code{jagsUI} package), or a matrix containing MCMC chains (each column representing MCMC output for a single parameter, rows representing iterations in the chain). The function automatically detects the object type and proceeds accordingly.
#'
#'
#' @examples
#' #Load data
#' data(MCMC_data)
#'
#' @export
#'


MCMCdiag <- function(object, 
                     open_txt = TRUE,
                     model_name,
                     file_name,
                     add_field,
                     add_field_names,
                     save_obj = TRUE,
                     obj_name,
                     add_obj, #additional objects to save
                     add_obj_names, #additional object names
                     wd = getwd(),
                     summary = TRUE,
                     params = 'all',
                     excl = NULL,
                     ISB = TRUE,
                     exact = TRUE,
                     digits = NULL,
                     round = NULL)
{
  #sort object types 
  if (methods::is(object, 'brmsfit'))
  {
    #extract stanfit portion of object
    object2 <- object$fit
  } else {
    if (methods::is(object, 'stanreg'))
    {
      object2 <- object$stanfit
    } else {
      #if mcmc object (from nimble) - convert to mcmc.list
      if (methods::is(object, 'mcmc'))
      {
        object2 <- coda::mcmc.list(object)
      } else {
        #if mcmc object (from nimble) - convert to mcmc.list
        if (methods::is(object, 'list'))
        {
          object2 <- coda::mcmc.list(lapply(object, function(x) coda::mcmc(x)))
        } else {
          object2 <- object
        }
      }
    }
  }
  
  if (coda::is.mcmc.list(object2) != TRUE &
      !methods::is(object2, 'matrix') &
      !methods::is(object2, 'rjags') &
      !methods::is(object2, 'stanfit') &
      !methods::is(object2, 'jagsUI'))
  {
    stop('Invalid object type. Input must be stanfit object (rstan), stanreg object (rstanarm), brmsfit object (brms), mcmc.list object (coda/rjags), mcmc object (coda/nimble), list object (nimble), rjags object (R2jags), jagsUI object (jagsUI), or matrix with MCMC chains.')
  }
  
  if (!missing(add_field))
  {
    if (missing(add_field_names))
    {
      add_field_names <- NULL
      warning('Name(s) for additional field(s) missing. Using arbitrary names in text document.')
    } else {
      if (length(add_field) != length(add_field_names))
      {
        warning('length(add_field_names) does not equal length(add_field). Using arbitrary names in text document for additional fields.')
        add_field_names <- NULL
      }
    }
  }
  
  #filename
  if (!missing(file_name))
  {
    #add .txt if it isn't in file_name
    if (length(grep('.txt', file_name)) > 0)
    {
      fn2 <- paste0(wd, '/', file_name)
    } else {
      fn2 <- paste0(wd, '/', file_name, '.txt')
    }
  } else {
    fn2 <- paste0(wd, '/MCMCdiag.txt')
  }
  
  
  #Stan object class
  if (methods::is(object2, 'stanfit'))
  {
    #get total elapsed time - in mins
    t_elapsed_time_pch <- rstan::get_elapsed_time(object2)
    time <- round(max(apply(t_elapsed_time_pch, 1, sum)) / 60, 2)
    
    #sampler inputs
    stan_args <- object2@stan_args[[1]]
    
    iter <- stan_args$iter
    warmup <- stan_args$warmup
    thin <- stan_args$thin
    nchains <- length(object2@stan_args)
    burnin <- NULL
    
    #may or may not have this info
    adapt_delta <- stan_args$control$adapt_delta
    max_treedepth <- stan_args$control$max_treedepth
    initial_stepsize <- stan_args$control$stepsize
    
    SUMMARY <- MCMCvis::MCMCsummary(object2, params = params, excl = excl, 
                                    ISB = ISB, exact = exact, digits = digits, round = round)
    max_Rhat <- max(SUMMARY[,'Rhat'], na.rm = TRUE)
    min_n.eff <- min(SUMMARY[,'n.eff'], na.rm = TRUE)
    
    sampler_params <- rstan::get_sampler_params(object2, inc_warmup = FALSE)
    mn_stepsize <- sapply(sampler_params, 
                          function(x) mean(x[, 'stepsize__']))
    mn_treedepth <- sapply(sampler_params, 
                           function(x) mean(x[, 'treedepth__']))
    accept_stat <- sapply(sampler_params, 
                          function(x) mean(x[, 'accept_stat__']))
    num_diverge <- rstan::get_num_divergent(object2)
    num_tree <- rstan::get_num_max_treedepth(object2)
    num_BFMI <- length(rstan::get_low_bfmi_chains(object2))
  }
  
  #jagsUI object
  if (methods::is(object2, 'jagsUI'))
  {
    #jagsUI has runtime, iter, burnin, iter, thin
    time <- round(object2$mcmc.info$elapsed.mins, 2)
    iter <- object2$mcmc.info$n.iter
    burnin <- object2$mcmc.info$n.burnin
    thin <- object2$mcmc.info$n.thin
    nchains <- object2$mcmc.info$n.chains
    
    #define as NULL for non-stan models
    warmup <- NULL
    adapt_delta <- NULL
    max_treedepth <- NULL
    initial_stepsize <- NULL
    mn_stepsize <- NULL
    mn_treedepth <- NULL
    accept_stat <- NULL
    num_diverge <- NULL
    num_tree <- NULL
    num_BFMI <- NULL
    
    SUMMARY <- MCMCvis::MCMCsummary(object2, params = params, excl = excl, 
                                    ISB = ISB, exact = exact, digits = digits, round = round)
    max_Rhat <- max(SUMMARY[,'Rhat'], na.rm = TRUE)
    min_n.eff <- min(SUMMARY[,'n.eff'], na.rm = TRUE)
  }
  
  #mcmc.list objects
  if (methods::is(object2, 'mcmc.list'))
  {
    #mcmc.list
    time <- NULL #no runtime for mcmc.list
    set <- attr(object2[[1]], 'mcpar') #start, end, thin  
    iter <- set[2]
    burnin <- set[1] - 1
    thin <- set[3]
    nchains <- length(object2)
    
    if (burnin == 0)
    {
      warning("According to provided 'object', burn-in = 0. Burn-in will not be printed to the diagnostic text file. Note that burn-in information cannot be extracted from 'nimble' output.")
      burnin <- NULL
    }
    
    #define as NULL for non-stan models
    warmup <- NULL
    adapt_delta <- NULL
    max_treedepth <- NULL
    initial_stepsize <- NULL
    mn_stepsize <- NULL
    mn_treedepth <- NULL
    accept_stat <- NULL
    num_diverge <- NULL
    num_tree <- NULL
    num_BFMI <- NULL
    
    SUMMARY <- MCMCvis::MCMCsummary(object2, params = params, excl = excl, 
                                    ISB = ISB, exact = exact, digits = digits, round = round)
    max_Rhat <- max(SUMMARY[,'Rhat'], na.rm = TRUE)
    min_n.eff <- min(SUMMARY[,'n.eff'], na.rm = TRUE)
  }
  
  #rjags objects
  if (methods::is(object2, 'rjags'))
  {
    #rjags
    time <- NULL #no runtime for mcmc.list
    iter <- object2$BUGSoutput$n.iter
    burnin <- object2$BUGSoutput$n.burnin
    thin <- object2$BUGSoutput$n.thin
    nchains <- object2$BUGSoutput$n.chains
    
    #define as NULL for non-stan models
    warmup <- NULL
    adapt_delta <- NULL
    max_treedepth <- NULL
    initial_stepsize <- NULL
    mn_stepsize <- NULL
    mn_treedepth <- NULL
    accept_stat <- NULL
    num_diverge <- NULL
    num_tree <- NULL
    num_BFMI <- NULL
    
    SUMMARY <- MCMCvis::MCMCsummary(object2, params = params, excl = excl, 
                                    ISB = ISB, exact = exact, digits = digits, round = round)
    max_Rhat <- max(SUMMARY[,'Rhat'], na.rm = TRUE)
    min_n.eff <- min(SUMMARY[,'n.eff'], na.rm = TRUE)
  }
  
  #matrix objects
  if (methods::is(object2, 'matrix'))
  {
    #mcmc.list
    time <- NULL #no runtime for mcmc.list
    iter <- NROW(object2)
    burnin <- NULL
    thin <- NULL
    nchains <- 1
    
    #define as NULL for non-stan models
    warmup <- NULL
    adapt_delta <- NULL
    max_treedepth <- NULL
    initial_stepsize <- NULL
    mn_stepsize <- NULL
    mn_treedepth <- NULL
    accept_stat <- NULL
    num_diverge <- NULL
    num_tree <- NULL
    num_BFMI <- NULL
    
    SUMMARY <- MCMCvis::MCMCsummary(object2, params = params, excl = excl, 
                                    ISB = ISB, exact = exact, digits = digits, round = round)
    max_Rhat <- NULL
    min_n.eff <- NULL
  }
  
  #create .txt file
  options(max.print = 1e8)
  sink(fn2)
  cat(paste0('Information and diagnostics \n'))
  cat(paste0('=========================== \n'))
  if (!missing(model_name))
  {
  cat(paste0('Model name:                       ', model_name, ' \n'))
  }
  if (!is.null(time))
  {
  cat(paste0('Run time (min):                   ', time, ' \n'))
  }
  cat(paste0('Total iter:                       ', iter, ' \n'))
  if (!is.null(warmup))
  {
  cat(paste0('Warmup:                           ', warmup, ' \n'))
  }
  if (!is.null(burnin))
  {
  cat(paste0('Burn-in:                          ', burnin, ' \n'))
  }
  if (!is.null(thin))
  {
  cat(paste0('Thin:                             ', thin, ' \n'))
  }
  cat(paste0('Num chains:                       ', nchains, ' \n'))
  if (!is.null(adapt_delta))
  {
  cat(paste0('Adapt delta (specified):          ', adapt_delta, ' \n'))
  }
  if (!is.null(max_treedepth))
  {
  cat(paste0('Max tree depth (specified):       ', max_treedepth, ' \n'))
  }
  if (!is.null(initial_stepsize))
  {
  cat(paste0('Initial step size (specified):    ', initial_stepsize, ' \n'))
  }
  if (!is.null(accept_stat))
  {
  cat(paste0('Mean accept stat:                 ', round(mean(accept_stat), 2), ' \n'))
  }
  if (!is.null(mn_treedepth))
  {
  cat(paste0('Mean tree depth:                  ', round(mean(mn_treedepth), 1), ' \n'))
  }
  if (!is.null(mn_stepsize))
  {
  cat(paste0('Mean step size:                   ', round(mean(mn_stepsize), 4), ' \n'))
  }
  if (!is.null(num_diverge))
  {
  cat(paste0('Num divergent transitions:        ', num_diverge, ' \n'))
  }
  if (!is.null(num_tree))
  {
  cat(paste0('Num max tree depth exceeds:       ', num_tree, ' \n'))
  }
  if (!is.null(num_BFMI))
  {
  cat(paste0('Num chains with BFMI warnings:    ', num_BFMI, ' \n'))
  }
  if (!is.null(max_Rhat))
  {
  cat(paste0('Max Rhat:                         ', max_Rhat, ' \n'))
  }
  if (!is.null(min_n.eff))
  {
  cat(paste0('Min n.eff:                        ', min_n.eff, ' \n'))
  }
  if (!missing(add_field))
  {
    if (is.null(add_field_names))
    {
      add_field_names <- paste0('Additional field ', 1:length(add_field))
    }
    for (i in 1:length(add_field))
    {
      sp_rep <- 34 - nchar(add_field_names[i]) - 1
      cat(paste0(add_field_names[i], ':', paste0(rep(' ', sp_rep), collapse = ''), add_field[i], ' \n'))
    }
  }
  cat(paste0('\n'))
  cat(paste0('\n'))
  if (summary == TRUE)
  {
    cat(paste0('Model summary \n'))
    cat(paste0('============= \n'))
    print(SUMMARY)
  }
  sink()
  
  
  
  #save model object
  if (save_obj == TRUE)
  {
    if (!missing(obj_name))
    {
      #add .rds if it isn't in file_name
      if (length(grep('.rds', obj_name)) > 0)
      {
        on2 <- paste0(wd, '/', obj_name)
      } else {
        on2 <- paste0(wd, '/', obj_name, '.rds')
      }
      saveRDS(object, file = on2)
    } else {
      saveRDS(object, file = paste0(wd, '/model_fit.rds'))
    }
  }
  
  if (open_txt == TRUE)
  {
    system(paste0('open ', fn2))
  }
}