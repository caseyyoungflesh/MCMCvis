#' Diagnostics summaries for models
#'
#' Model diagnostics and summary. Function reads information embedded in model fit object. Output varies by model fit object type but includes model run inputs, diagnostic information, and parameter summary. See DETAILS below for more information.
#'
#' 
#' @param object Object containing MCMC output. See DETAILS below.
#'
#' @param file_name Character string with name of .txt file to be saved to \code{dir} (or \code{mkdir} if specified). If not specified, 'MCMCdiag.txt' will be used.
#' 
#' @param dir Character string with directory where file(s) (or directory is argument for \code{mkdir} is specified) will be created. Defaults to working directory. An absolute or relative (to the working directory) path can be used.
#' 
#' @param mkdir Character string with name of directory to be created. If specified, a directory will be created within the directory specified by \code{dir}. If directory exists, '_1' will be appended to specified name.
#' 
#' @param add_field Object (or vector of objects) to be added to the .txt file.
#' 
#' @param add_field_names Character string (or vector of character strings) specifying the name(s) of the \code{add_field} object(s).
#' 
#' @param save_object Logical specifying whether the model output provided to the function (\code{object}) should be saved as a \code{.rds} file to \code{dir} (or \code{mkdir} if specified). Note that \code{.rds} files can be opened with \code{rdsRDS()}.
#' 
#' @param obj_name Character string specifying the file name of the \code{.rds} file (created from \code{object}) if \code{save_object = TRUE}.
#' 
#' @param add_obj List with additional object(s) to be saved as \code{.rds} files to \code{dir} (or \code{mkdir} if specified). Objects can be of any types. Multiple objects can be specified. Note that \code{.rds} files can be opened with \code{rdsRDS()}.
#' 
#' @param add_obj_names Character string (or vector of character strings) specifying the name(s) of the objects to be saved as \code{.rds} files, specified with \code{add_obj}.
#' 
#' @param cp_file Character string (or vector of character strings) specifying file(s) to be copied to \code{dir} (or \code{mkdir} if specified). Absolute or relative (to the working directory) paths can be used.
#' 
#' @param cp_file_names Character string (or vector of character strings) specifying new names for files to be copied specified by \code{cp_file}. If no argument is provided, the copy names will be identical to the originals.
#'
#' @param open_txt Logical - if \code{open_txt = TRUE} .txt file will open in default .txt viewer after being generated.
#'
#' @param summary Logical specifying whether or not to output summary information from \link[MCMCvis]{MCMCsummary} at the bottom of the .txt file.
#'
#' @param ... Arguments to be passed to \link[MCMCvis]{MCMCsummary} when generating summary if \code{summary = TRUE}. 
#' 
#' @section Details:
#' Some diagnostic information is only provided for models fit with particular pieces of software. For example, \code{rstan} output includes additional diagnostics related to the NUTS sampler. Output from \code{jagsUI} includes runtime information, but output from \code{rjags} does not. Note that this information could be fed manually to the function using the \code{add_field} argument.
#'
#' \code{object} argument can be a \code{stanfit} object (\code{rstan} package), a \code{CmdStanMCMC} object (\code{cmdstanr} package), a \code{stanreg} object (\code{rstanarm} package), a \code{brmsfit} object (\code{brms} package), an \code{mcmc.list} object (\code{coda} and \code{rjags} packages), \code{mcmc} object (\code{coda} and \code{nimble} packages), \code{list} object (\code{nimble} package), an \code{R2jags} model object (\code{R2jags} package), a \code{jagsUI} model object (\code{jagsUI} package), or a matrix containing MCMC chains (each column representing MCMC output for a single parameter, rows representing iterations in the chain). The function automatically detects the object type and proceeds accordingly.
#' 
#' Output presented in .txt file varies by model fit object type but includes: model run time, number of warmup/burn-in iterations, total iterations, thinning rate, number of chains, specified adapt delta, specified max tree depth, specific initial step size, mean accept stat, mean tree depth, mean step size, number of divergent transitions, number max tree depth exceeds, number of chains with BFMI warnings, max Rhat (maximum Rhat of any parameter printed), min n.eff (minimum n.eff of any parameter printed), parameter summary information (passed from \code{MCMCsummary}), and any additional information fed to the \code{add_field} argument. See documentation for specific software used to fit model for more information on particular diagnostics.
#'
#'
#' @examples
#' #Load data
#' data(MCMC_data)
#'
#' # MCMCdiag(MCMC_data,
#' #          #name of .txt file to be saved
#' #          file_name = 'blog-model-summary-2021-01-15',
#' #          #directory within which to save .txt file
#' #          dir = '~/Desktop',
#' #          #round MCMCsummary output in .txt file to two digits
#' #          round = 2,
#' #          #add two fields to .txt file
#' #          add_field = c(50, '1.0'),
#' #          #names of two additional fields to add to .txt file
#' #          add_field_names = c('Time (min)', 'Data version'))
#'
#' @export
#'


MCMCdiag <- function(object, 
                     file_name,
                     dir = getwd(),
                     mkdir,
                     add_field,
                     add_field_names,
                     save_object = FALSE,
                     obj_name,
                     add_obj,
                     add_obj_names,
                     cp_file,
                     cp_file_names,
                     open_txt = TRUE,
                     summary = TRUE,
                     ...)
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
      !methods::is(object2, 'jagsUI') &
      !methods::is(object, 'CmdStanMCMC'))
  {
    stop('Invalid object type. Input must be stanfit object (rstan), CmdStanMCMC object (cmdstanr), stanreg object (rstanarm), brmsfit object (brms), mcmc.list object (coda/rjags), mcmc object (coda/nimble), list object (nimble), rjags object (R2jags), jagsUI object (jagsUI), or matrix with MCMC chains.')
  }
  
  #additional fields anmes
  if (!missing(add_field))
  {
    if (missing(add_field_names))
    {
      add_field_names <- NULL
      warning('Name(s) for additional field(s) missing. Using arbitrary name(s) in .txt files')
    } else {
      if (length(add_field) != length(add_field_names))
      {
        warning('length(add_field_names) does not equal length(add_field). Using arbitrary name(s) in .txt file for additional field(s).')
        add_field_names <- NULL
      }
    }
  }
  
  #additional objects names
  if (!missing(add_obj))
  {
    if (missing(add_obj_names))
    {
      add_obj_names <- NULL
      warning('Name(s) for additional object(s) missing. Using arbitrary name(s) for file(s).')
    } else {
      if (length(add_obj) != length(add_obj_names))
      {
        warning('length(add_obj_names) does not equal length(add_obj). Using arbitrary name(s) for file(s).')
        add_obj_names <- NULL
      }
    }
  }
  
  #cp file names
  if (!missing(cp_file) & !missing(cp_file_names))
  {
    if (length(cp_file) != length(cp_file_names))
    {
      warning('length(cp_file_names) does not equal length(cp_file). Using original file name(s).')
      cp_file_names <- NULL
    }
  }
  
  #mkdir
  if (!missing(mkdir))
  {
    # if specified dir to write to doesn't exist, write to home dir
    if (!dir.exists(dir))
    {
      warning("Specified `dir` does not exist. Using current working directory.")
      dir <- getwd()
    }
    dir <- paste0(dir, '/', mkdir)
    # if specified dir name (mkdir) exists, append _1 to name to avoid overwriting files
    if (dir.exists(dir))
    {
      warning("Specified `mkdir` exists. Appending '_1' to directory name.")
      dir <- paste0(dir, '_1')
    }
    dir.create(dir)
  }
  
  #filename
  if (!missing(file_name))
  {
    #add .txt if it isn't in file_name
    if (length(grep('.txt', file_name)) > 0)
    {
      fn2 <- paste0(dir, '/', file_name)
    } else {
      fn2 <- paste0(dir, '/', file_name, '.txt')
    }
  } else {
    fn2 <- paste0(dir, '/MCMCdiag.txt')
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
    
    SUMMARY <- MCMCvis::MCMCsummary(object2, ...)
    max_Rhat <- max(SUMMARY[,'Rhat'], na.rm = TRUE)
    min_n.eff <- min(SUMMARY[,'n.eff'], na.rm = TRUE)
    min_n.eff_bulk <- NULL
    min_n.eff_tail <- NULL
    
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
  
  #CmdStanMCMC object class
  if (methods::is(object2, 'CmdStanMCMC'))
  {
    #get total elapsed time - in mins
    t_elapsed_time_pch <- object2$time()$total
    time <- round(t_elapsed_time_pch / 60, 2)

    #sampler inputs
    metadata <- object2$metadata()

    iter <- metadata$iter_sampling
    warmup <- metadata$iter_warmup
    thin <- metadata$thin
    nchains <- length(metadata$id)
    burnin <- NULL

    #may or may not have this info
    adapt_delta <- metadata$adapt_delta
    max_treedepth <- metadata$max_treedepth
    initial_stepsize <- mean(metadata$step_size)

    SUMMARY <- MCMCvis::MCMCsummary(object2, ...)
    max_Rhat <- max(SUMMARY[,'Rhat'], na.rm = TRUE)
    min_n.eff_bulk <- min(SUMMARY[,'n.eff_bulk'], na.rm = TRUE)
    min_n.eff_tail <- min(SUMMARY[,'n.eff_tail'], na.rm = TRUE)
    min_n.eff <- NULL

    sampler_params <- object2$sampler_diagnostics(inc_warmup = FALSE, format = 'df')

    mn_stepsize <- mean(sampler_params$stepsize__)
    mn_treedepth <- mean(sampler_params$treedepth__)
    accept_stat <- mean(sampler_params$accept_stat__)
    num_diverge <- sum(object2$diagnostic_summary()$num_divergent)
    num_tree <- sum(object2$diagnostic_summary()$num_max_treedepth)
    num_BFMI <- sum(object2$diagnostic_summary()$ebfmi < 0.2)
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
    
    SUMMARY <- MCMCvis::MCMCsummary(object2, ...)
    max_Rhat <- max(SUMMARY[,'Rhat'], na.rm = TRUE)
    min_n.eff <- min(SUMMARY[,'n.eff'], na.rm = TRUE)
    min_n.eff_bulk <- NULL
    min_n.eff_tail <- NULL
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
      warning("According to provided 'object', burn-in = 0. Burn-in will not be printed to the diagnostic .txt file. Note that burn-in and thin information cannot be extracted from 'nimble' output; 'Total iter' reflects the number of iterations kept after thinning and burn-in, and 'Thin' is given as 1 regardless of the thinning rate used.")
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
    
    SUMMARY <- MCMCvis::MCMCsummary(object2, ...)
    max_Rhat <- max(SUMMARY[,'Rhat'], na.rm = TRUE)
    min_n.eff <- min(SUMMARY[,'n.eff'], na.rm = TRUE)
    min_n.eff_bulk <- NULL
    min_n.eff_tail <- NULL
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
    
    SUMMARY <- MCMCvis::MCMCsummary(object2, ...)
    max_Rhat <- max(SUMMARY[,'Rhat'], na.rm = TRUE)
    min_n.eff <- min(SUMMARY[,'n.eff'], na.rm = TRUE)
    min_n.eff_bulk <- NULL
    min_n.eff_tail <- NULL
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
    
    SUMMARY <- MCMCvis::MCMCsummary(object2, ...)
    max_Rhat <- NULL
    min_n.eff <- NULL
    min_n.eff_bulk <- NULL
    min_n.eff_tail <- NULL
  }
  
  #create .txt file
  options(max.print = 1e8)
  sink(fn2)
  cat(paste0('Information and diagnostics \n'))
  cat(paste0('=========================== \n'))
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
  if (!is.null(min_n.eff_bulk))
  {
    cat(paste0('Min n.eff_bulk:                   ', min_n.eff_bulk, ' \n'))
  }
  if (!is.null(min_n.eff_tail))
  {
    cat(paste0('Min n.eff_tail:                   ', min_n.eff_tail, ' \n'))
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
      if (sp_rep < 1)
      {
        sp_rep <- 1
      }
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
  
  #copy files
  if (!missing(cp_file))
  {
    if (!missing(cp_file_names))
    {
      for (i in 1:length(cp_file))
      {
        if (file.exists(cp_file[i]) == FALSE)
        {
          warning(paste0("Could not copy '", cp_file[i],"' as it does not exist. Check that the provided file path is correct."))
        } else {
          invisible(file.copy(from = cp_file[i], 
                              to = paste0(dir, '/', cp_file_names[i]))) 
        }
      }
    } else {
      for (i in 1:length(cp_file))
      {
        if (file.exists(cp_file[i]) == FALSE)
        {
          warning(paste0("Could not copy '", cp_file[i],"' as it does not exist. Check that the provided file path is correct."))
        } else {
          invisible(file.copy(from = cp_file[i], 
                              to = dir))
        }
      }
    }
  }
  
  #save model object
  if (save_object == TRUE)
  {
    if (!missing(obj_name))
    {
      #add .rds if it isn't in file_name
      if (length(grep('.rds', obj_name)) > 0)
      {
        on2 <- paste0(dir, '/', obj_name)
      } else {
        on2 <- paste0(dir, '/', obj_name, '.rds')
      }
      saveRDS(object, file = on2)
    } else {
      saveRDS(object, file = paste0(dir, '/model_fit.rds'))
    }
  }
  
  #save add objects
  if (!missing(add_obj))
  {
    if (is.null(add_obj_names))
    {
      add_obj_names <- paste0('Additional_object_', 1:length(add_obj), '.rds')
    }
    if (!is.null(add_obj_names))
    {
      ao2 <- rep(NA, length(add_obj_names))
      for (i in 1:length(add_obj_names))
      {
        #add .rds if it isn't in file_name
        if (length(grep('.rds', add_obj_names[i])) > 0)
        {
          ao2[i] <- paste0(dir, '/', add_obj_names[i])
        } else {
          ao2[i] <- paste0(dir, '/', add_obj_names[i], '.rds')
        }
      }
    }
    for (i in 1:length(add_obj))
    {
      saveRDS(add_obj[[i]], file = ao2[i])
    }
  }
  
  if (open_txt == TRUE)
  {
    system(paste0('open ', fn2))
  }
}
