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
                     save_object = TRUE,
                     object_name,
                     wd = getwd(),
                     summary = TRUE,
                     params = 'all',
                     excl = NULL,
                     ISB = TRUE,
                     exact = TRUE,
                     digits = NULL,
                     round = NULL)
{
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
  
  # #rstanarm and brms object types -> convert to stan
  # if (methods::is(object, ''))
  # {
  # }
  
  # #other object types -> convert to mcmc.list
  # if (methods::is(object, ''))
  # {
  # }
  
  #stan object class
  if (methods::is(object, 'stanfit'))
  {
    #get total elapsed time - in mins
    t_elapsed_time_pch <- rstan::get_elapsed_time(object)
    t_elapsed_time <- round(max(apply(t_elapsed_time_pch, 1, sum)) / 60, 2)
    
    #sampler inputs
    stan_args <- object@stan_args[[1]]
    
    iter <- stan_args$iter
    warmup <- stan_args$warmup
    thin <- stan_args$thin
    nchains <- length(object@stan_args)
    
    #may or may not have this info
    adapt_delta <- stan_args$control$adapt_delta
    max_treedepth <- stan_args$control$max_treedepth
    initial_stepsize <- stan_args$control$stepsize
    
    SUMMARY <- MCMCvis::MCMCsummary(object, params = params, excl = excl, 
                                    ISB = ISB, exact = exact, digits = digits, round = round)
    max_Rhat <- max(SUMMARY[,'Rhat'], na.rm = TRUE)
    min_n.eff <- min(SUMMARY[,'n.eff'], na.rm = TRUE)
    
    sampler_params <- rstan::get_sampler_params(object, inc_warmup = FALSE)
    mn_stepsize <- sapply(sampler_params, 
                          function(x) mean(x[, 'stepsize__']))
    mn_treedepth <- sapply(sampler_params, 
                           function(x) mean(x[, 'treedepth__']))
    accept_stat <- sapply(sampler_params, 
                          function(x) mean(x[, 'accept_stat__']))
    num_diverge <- rstan::get_num_divergent(object)
    num_tree <- rstan::get_num_max_treedepth(object)
    num_BFMI <- length(rstan::get_low_bfmi_chains(object))
    
    #sampler stats for each chain
    # sampler_params <- rstan::get_sampler_params(object, inc_warmup = inc_warmup)
    # mn_stepsize <- round(sapply(sampler_params, 
    #                           function(x) mean(x[, 'stepsize__'])), round)
    # mn_treedepth <- round(sapply(sampler_params, 
    #                           function(x) mean(x[, 'treedepth__'])), round)
    # # # mn_energy <- round(sapply(sampler_params, 
    # # #                           function(x) mean(x[, 'energy__'])), round)
    # accept_stat <- round(sapply(sampler_params, 
    #                             function(x) mean(x[, 'accept_stat__'])), round)

    # #get elapsed time per chain
    # etime <- get_elapsed_time(object)
    # warmup_time <- round(etime[,1], 0)
    # sampling_time <- round(etime[,1], 0)
    # #elapsed_time <- as.character(round(get_elapsed_time(object), 0)) #in seconds
    # #colnames(elapsed_time) <- c('warmup_time', 'sampling_time')
    
    # #chain specific metrics
    # ch_metrics <- data.frame(warmup_time, sampling_time, 
    #                          n_divergent, n_exceed_max_tree, 
    #                          BFMI, mn_stepsize, mn_treedepth,
    #                          mn_energy, accept_stat)

    # colnames(ch_metrics) <- c('Total\ runtime\ (sec)',
    #                           'Total iterations',
    #                           'Warmup',
    #                           'Number\ divergent\ iter',
    #                           'Number\ iter\ exceeding\ max\ tree\ depth',
    #                           'BFMI\ (warnings\ generated\ < 0.2)',
    #                           'Mean\ stepsize',
    #                           'Mean\ treedepth',
    #                           'Mean\ energy',
    #                           'Accept\ stat')
    
    # #calc max digits for row names
    # digits_rows <- nchar(colnames(ch_metrics))
    # max_dr <- max(digits_rows) + 2
    
    # #calc max digits for spaces to separate cols
    # #NEED SEPARATE METRICS FOR EACH COLUMN (due to potential differing lengths)
    # for (i in 1:nchains)
    # {
    #   print(nchar(ch_metrics[i,]))
    # }
    # apply(paste0(ch_metrics), 2, function(x) nchar(paste0(x)))
    # digits_ch <- paste0(ch_metrics[1,])
    # max_dch <- max(digits_ch) + 2
    # 
    # 
    # #need at least 9 spaces between cols for chain-level diagnostics (for colnames)
    # vdch <- as.vector(digits_ch)
    # if (max_dch > 9)
    # {
    #   spaces <- (max_dch - vdch)
    #   ch_sp <- max_dch + 2
    # } else {
    #   spaces <- rep(9, length(digits_ch)) - vdch
    #   ch_sp <- 2
    # }
    
    #create .txt file
    options(max.print = 1e8)
    sink(fn2)
    cat(paste0('Information and diagnostics \n'))
    cat(paste0('=========================== \n'))
    if (!missing(model_name))
    {
    cat(paste0('Model name:                       ', model_name, ' \n'))
    }
    cat(paste0('Time elapsed (min):               ', t_elapsed_time, ' \n'))
    cat(paste0('Total iter:                       ', iter, ' \n'))
    cat(paste0('Warmup:                           ', warmup, ' \n'))
    cat(paste0('Thin:                             ', thin, ' \n'))
    cat(paste0('Num chains:                       ', nchains, ' \n'))
    if (!is.null(adapt_delta))
    {
    cat(paste0('Adapt delta (specified):          ', adapt_delta, ' \n'))
    }
    if (!is.null(max_treedepth))
    {
    cat(paste0('Max treedepth (specified):        ', max_treedepth, ' \n'))
    }
    if (!is.null(initial_stepsize))
    {
    cat(paste0('Initial stepsize (specified):     ', initial_stepsize, ' \n'))
    }
    cat(paste0('Mean accept stat:                 ', round(mean(accept_stat), 2), ' \n'))
    cat(paste0('Mean treedepth:                   ', round(mean(mn_treedepth), 1), ' \n'))
    cat(paste0('Mean stepsize:                    ', round(mean(mn_stepsize), 5), ' \n'))
    cat(paste0('Num divergent transitions:        ', num_diverge, ' \n'))
    cat(paste0('Num max tree depth exceeds:       ', num_tree, ' \n'))
    cat(paste0('Num chains with BFMI warnings:    ', num_BFMI, ' \n'))
    cat(paste0('Max Rhat:                         ', max_Rhat, ' \n'))
    cat(paste0('Min n.eff:                        ', min_n.eff, ' \n'))
    cat(paste0('\n'))
    cat(paste0('\n'))
    # cat(paste0('Chain-level diagnostics \n'))
    # cat(paste0('======================= \n'))
    # cat(paste0(paste(rep(' ', max_dr + 2), collapse = ''),
    #            'chain 1',
    #            paste(rep(' ', ch_sp), collapse = ''),
    #            'chain 2',
    #            paste(rep(' ', ch_sp), collapse = ''),
    #            'chain 3',
    #            ' \n'))
    # 
    # #generalize to arbitrary # of chains
    # cat(paste0('Warmup time (sec): ', 
    #            paste(rep(' ', c(max_dr - digits_rows[1])), collapse = ''),
    #            warmup_time[1], 
    #            paste(rep(' ', spaces[1]), collapse = ''),
    #            warmup_time[2],
    #            paste(rep(' ', spaces[2]), collapse = ''),
    #            warmup_time[3], ' \n'))
    # cat(paste0('Sampling time (sec): ',
    #            paste(rep(' ', c(max_dr - digits_rows[2])), collapse = ''),
    #            sampling_time[1], 
    #            paste(rep(' ', spaces[4]), collapse = ''),
    #            sampling_time[2],
    #            paste(rep(' ', spaces[5]), collapse = ''),
    #            sampling_time[3], ' \n'))
    # cat(paste0('Number iter divergent: ', 
    #            paste(rep(' ', c(max_dr - digits_rows[3])), collapse = ''),
    #            n_divergent[1], 
    #            paste(rep(' ', spaces[7]), collapse = ''),
    #            n_divergent[2],
    #            paste(rep(' ', spaces[8]), collapse = ''),
    #            n_divergent[3], ' \n'))
    # cat(paste0('Number iter exceeding max treedepth: ',
    #            paste(rep(' ', c(max_dr - digits_rows[4]) + 1), collapse = ''),
    #            n_exceed_max_tree[1], 
    #            paste(rep(' ', spaces[10]), collapse = ''),
    #            n_exceed_max_tree[2],
    #            paste(rep(' ', spaces[11]), collapse = ''),
    #            n_exceed_max_tree[3], ' \n'))
    # cat(paste0('BFMI (warnings generated < 0.2): ', 
    #            paste(rep(' ', c(max_dr - digits_rows[5])), collapse = ''),
    #            BFMI[1], 
    #            paste(rep(' ', spaces[13]), collapse = ''),
    #            BFMI[2],
    #            paste(rep(' ', spaces[14]), collapse = ''),
    #            BFMI[3], ' \n'))
    # cat(paste0('Mean stepsize: ', 
    #            paste(rep(' ', c(max_dr - digits_rows[6])), collapse = ''),
    #            mn_stepsize[1], 
    #            paste(rep(' ', spaces[16]), collapse = ''),
    #            mn_stepsize[2],
    #            paste(rep(' ', spaces[17]), collapse = ''),
    #            mn_stepsize[3], ' \n'))
    # cat(paste0('Mean treedepth: ', 
    #            paste(rep(' ', c(max_dr - digits_rows[7])), collapse = ''),
    #            mn_treedepth[1], 
    #            paste(rep(' ', spaces[19]), collapse = ''),
    #            mn_treedepth[2],
    #            paste(rep(' ', spaces[20]), collapse = ''),
    #            mn_treedepth[3], ' \n'))
    # cat(paste0('Mean energy: ', 
    #            paste(rep(' ', c(max_dr - digits_rows[8])), collapse = ''),
    #            mn_energy[1], 
    #            paste(rep(' ', spaces[22]), collapse = ''),
    #            mn_energy[2],
    #            paste(rep(' ', spaces[23]), collapse = ''),
    #            mn_energy[3], ' \n'))
    # cat(paste0('Accept stat: ', 
    #            paste(rep(' ', c(max_dr - digits_rows[9])), collapse = ''),
    #            accept_stat[1], 
    #            paste(rep(' ', spaces[25]), collapse = ''),
    #            accept_stat[2],
    #            paste(rep(' ', spaces[26]), collapse = ''),
    #            accept_stat[3], ' \n'))
    # cat(paste0('\n'))
    # cat(paste0('\n'))
    if (summary == TRUE)
    {
      cat(paste0('Model summary \n'))
      cat(paste0('============= \n'))
      print(SUMMARY)
    }
    sink()
  }
  
  #non-stan/non-jagsUI object types - convert to mcmc.list object
  if (methods::is(object, 'mcmc.list'))
  {
    #jagsUI has runtime, iter, burnin, iter, thin
    t_elapsed_time <- round(object$mcmc.info$elapsed.mins, 2)
    iter <- object$mcmc.info$n.samples #n.samples is total iter from jagsUI
    burnin <- object$mcmc.info$n.burnin
    thin <- object$mcmc.info$n.thin
    nchains <- object$mcmc.info$n.chains
    
    SUMMARY <- MCMCvis::MCMCsummary(object, params = params, excl = excl, 
                                    ISB = ISB, exact = exact, digits = digits, round = round)
    max_Rhat <- max(SUMMARY[,'Rhat'], na.rm = TRUE)
    min_n.eff <- min(SUMMARY[,'n.eff'], na.rm = TRUE)
    
    #create .txt file
    options(max.print = 1e8)
    sink(fn2)
    cat(paste0('Model-level information and diagnostics \n'))
    cat(paste0('======================================= \n'))
    if (!missing(model_name))
    {
      cat(paste0('Model name:             ', model_name, ' \n'))
    }
    cat(paste0('Time elapsed (min):     ', t_elapsed_time, ' \n'))
    cat(paste0('Total iter:             ', iter, ' \n'))
    cat(paste0('Burn-in:                ', burnin, ' \n'))
    cat(paste0('Thin:                   ', thin, ' \n'))
    cat(paste0('Num chains:             ', nchains, ' \n'))
    cat(paste0('Max Rhat:               ', max_Rhat, ' \n'))
    cat(paste0('Min n.eff:              ', min_n.eff, ' \n'))
    cat(paste0('\n'))
    cat(paste0('\n'))
    if (summary == TRUE)
    {
      cat(paste0('Model summary \n'))
      cat(paste0('============= \n'))
      print(SUMMARY)
    }
    sink()
  }
  
  
  #non-stan/non-jagsUI object types - convert to mcmc.list object
  if (methods::is(object, 'mcmc.list'))
  {
    #mcmc.list
    #no runtime for mcmc.list
    set <- attr(object[[1]], 'mcpar') #start, end, thin  
    iter <- set[2]
    burnin <- set[1] - 1
    thin <- set[3]
    nchains <- length(object)
    
    SUMMARY <- MCMCvis::MCMCsummary(object, params = params, excl = excl, 
                                    ISB = ISB, exact = exact, digits = digits, round = round)
    max_Rhat <- max(SUMMARY[,'Rhat'], na.rm = TRUE)
    min_n.eff <- min(SUMMARY[,'n.eff'], na.rm = TRUE)
    
    #create .txt file
    options(max.print = 1e8)
    sink(fn2)
    cat(paste0('Model-level information and diagnostics \n'))
    cat(paste0('======================================= \n'))
    if (!missing(model_name))
    {
    cat(paste0('Model name:             ', model_name, ' \n'))
    }
    if (!missing(t_elapsed_time))
    {
    cat(paste0('Time elapsed (min):     ', t_elapsed_time, ' \n'))
    }
    cat(paste0('Total iter:             ', iter, ' \n'))
    cat(paste0('Burn-in:                ', burnin, ' \n'))
    cat(paste0('Thin:                   ', thin, ' \n'))
    cat(paste0('Num chains:             ', nchains, ' \n'))
    cat(paste0('Max Rhat:               ', max_Rhat, ' \n'))
    cat(paste0('Min n.eff:              ', min_n.eff, ' \n'))
    cat(paste0('\n'))
    cat(paste0('\n'))
    if (summary == TRUE)
    {
      cat(paste0('Model summary \n'))
      cat(paste0('============= \n'))
      print(SUMMARY)
    }
    sink()
  }
  
  #save model object
  if (save_object == TRUE)
  {
    if (!missing(object_name))
    {
      #add .rds if it isn't in file_name
      if (length(grep('.rds', object_name)) > 0)
      {
        on2 <- paste0(wd, '/', object_name)
      } else {
        on2 <- paste0(wd, '/', object_name, '.rds')
      }
      saveRDS(object, file = on2)
    } else {
      saveRDS(object, files = paste0(wd, '/model_fit.rds'))
    }
  }
  
  if (open_txt == TRUE)
  {
    system(paste0('open ', fn2))
  }
}
