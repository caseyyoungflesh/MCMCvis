#' Extract posterior chains from MCMC output
#'
#' Extract posterior chains from MCMC output for specific parameters of interest.
#'
#' @param object Object containing MCMC output. See DETAILS below.
#' @param params Character string (or vector of character strings) denoting parameters of interest.
#'
#' Default \code{'all'} returns chains for all parameters.
#'
#' @param excl Character string (or vector of character strings) denoting parameters to exclude. Used in conjunction with \code{params} argument to select parameters of interest.
#'
#' @param ISB Ignore Square Brackets (ISB). Logical specifying whether square brackets should be ignored in the \code{params} and \code{excl} arguments. If \code{TRUE}, square brackets are ignored. If \code{FALSE}, square brackets are not ignored.  This allows partial names to be used when specifying parameters of interest. Use \code{exact} argument to specify whether input from \code{params} and \code{excl} arguments should be matched exactly.
#'
#' @param exact Logical specifying whether input from \code{params} and \code{excl} arguments should be matched exactly (after ignoring square brackets if \code{ISB = FALSE}). If \code{TRUE}, input from \code{params} and \code{excl} are matched exactly (after taking \code{ISB} argument into account). If \code{FALSE}, input from \code{params} and \code{excl} are matched using regular expression format (after taking \code{ISB} argument into account).
#'
#' @param mcmc.list Logical specifying whether to return an mcmc.list. If \code{TRUE}, an \code{mcmc.list} object is returned, rather than a matrix.
#'
#' @param chain_num Numeric - specifies posterior chain number. When a value is specified, posterior for only that chain is output. Useful for determining the last iteration for each parameter, to be used as initial values in a subsequent model, to effectively 'continue' a model run.
#'
#' @section Details:
#' Function returns matrix with one parameter per column (for specified parameters). Each iteration is represented as a row. Multiple chains for each parameter are combined to one posterior chain (unless \code{chain_num} is specified, in which case only the specified chain will be returned). Parameters are arranged in columns alphabetically.
#'
#' \code{object} argument can be a \code{stanfit} object (\code{rstan} package), a \code{stanreg} object (\code{rstanarm} package), a \code{brmsfit} object (\code{brms} package), an \code{mcmc.list} object (\code{coda} and \code{rjags} packages), \code{mcmc} object (\code{coda} and \code{nimble} packages), \code{list} object (\code{nimble} package), an \code{R2jags} model object (\code{R2jags} package), a \code{jagsUI} model object (\code{jagsUI} package), or a matrix containing MCMC chains (each column representing MCMC output for a single parameter, rows representing iterations in the chain). The function automatically detects the object type and proceeds accordingly.
#'
#'
#' @examples
#' #Load data
#' data(MCMC_data)
#'
#' #Extract MCMC chains
#' ex <- MCMCchains(MCMC_data)
#' apply(ex, 2, mean)
#'
#' #Extract MCMC chains for just 'beta' parameters
#' ex2 <- MCMCchains(MCMC_data, params = 'beta')
#' apply(ex2, 2, mean)
#'
#' #Just 'beta[1]', 'beta[4]', and 'alpha[3]'
#' ex3 <- MCMCchains(MCMC_data, params = c('beta[1]', 'beta[4]', 'alpha[3]'), 
#'                  ISB = FALSE, exact = TRUE)
#' apply(ex3, 2, sd)
#'
#' @export

MCMCchains <- function(object,
                       params = 'all',
                       excl = NULL,
                       ISB = TRUE,
                       exact = TRUE,
                       mcmc.list = FALSE,
                       chain_num = NULL)
{
  #for rstanarm/brms obejcts - set to NULL by default
  sp_names <- NULL
  
  #if from R2jags::jags.parallel
  if (methods::is(object, 'rjags.parallel'))
  {
    #modified coda::as.mcmc (removing ordering of param names)
    x <- object$BUGSoutput
    mclist <- vector('list', x$n.chains)
    mclis <- vector('list', x$n.chains)
    ord <- dimnames(x$sims.array)[[3]]
    for (i in 1:x$n.chains)
    {
      tmp1 <- x$sims.array[, i, ord]
      mclis[[i]] <- coda::mcmc(tmp1, thin = x$n.thin)
    }
    object <- coda::as.mcmc.list(mclis)
    #end mod as.mcmc
  }
  
  #if mcmc object (from nimble) - convert to mcmc.list
  if (methods::is(object, 'mcmc'))
  {
    object <- coda::mcmc.list(object)
  }
  
  #if list object of matrices (from nimble) - convert to mcmc.list
  if (methods::is(object, 'list'))
  {
    object <- coda::mcmc.list(lapply(object, function(x) coda::mcmc(x)))
  }

  #if from rstanarm::stan_glm
  if (methods::is(object, 'stanreg'))
  {
    object <- object$stanfit
    sp_names <- object@sim$fnames_oi
  }
  
  if (coda::is.mcmc.list(object) != TRUE &
     !methods::is(object, 'matrix') &
     !methods::is(object, 'mcmc') &
     !methods::is(object, 'list') &
     !methods::is(object, 'rjags') &
     !methods::is(object, 'stanfit') &
     !methods::is(object, 'brmsfit') &
     !methods::is(object, 'jagsUI'))
  {
    stop('Invalid object type. Input must be stanfit object (rstan), stanreg object (rstanarm), brmsfit object (brms), mcmc.list object (coda/rjags), mcmc object (coda/nimble), list object (nimble), rjags object (R2jags), jagsUI object (jagsUI), or matrix with MCMC chains.')
  }
  
  #if from brms::brm
  if (methods::is(object, 'brmsfit'))
  {
    #extract stanfit portion of object
    object <- object$fit
    #Stan names
    sp_names_p <- names(object@sim$samples[[1]])
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
  }
  
  #NAME SORTING BLOCK
  if (methods::is(object, 'stanfit'))
  {
    #convert to mcmc.list
    temp_in <- rstan::As.mcmc.list(object)
    
    #assign new colnames for mcmc.list object if object exists (for stanreg and brms objs so parameter names are interpretable) - do not rename params for model fit directly with Stan
    if (!is.null(sp_names))
    {
      coda::varnames(temp_in) <- sp_names
    }
    
    if (ISB == TRUE)
    {
      names <- vapply(strsplit(colnames(temp_in[[1]]),
                               split = '[', fixed = TRUE), `[`, 1, FUN.VALUE=character(1))
    } else {
      names <- colnames(temp_in[[1]])
    }
  }

  if (methods::is(object, 'jagsUI'))
  {
    object <- object$samples
  }

  if (coda::is.mcmc.list(object) == TRUE)
  {
    temp_in <- object
    if (is.null(colnames(temp_in[[1]])))
    {
      warning('No parameter names provided. Assigning arbitrary names.')
      sub_cn <- paste0('Param_', 1:NCOL(temp_in[[1]]))
      colnames(temp_in[[1]]) <- sub_cn
    }
    
    if (ISB == TRUE)
    {
      names <- vapply(strsplit(colnames(temp_in[[1]]),
                               split = "[", fixed = TRUE), `[`, 1, FUN.VALUE=character(1))
    } else {
      names <- colnames(temp_in[[1]])
    }
  }

  if (methods::is(object, 'matrix'))
  {
    temp_in <- object
    if (is.null(colnames(temp_in)))
    {
      warning('No parameter names (column names) provided. Assigning arbitrary names.')
      sub_cn <- paste0('Param_', 1:NCOL(temp_in))
      colnames(temp_in) <- sub_cn
    }

    if (ISB == TRUE)
    {
      names <- vapply(strsplit(colnames(temp_in),
                               split = "[", fixed = TRUE), `[`, 1, FUN.VALUE=character(1))
    } else {
      names <- colnames(temp_in)
    }
  }

  if (methods::is(object, 'rjags'))
  {
    temp_in <- object$BUGSoutput$sims.matrix
    if (ISB == TRUE)
    {
      names <- vapply(strsplit(rownames(object$BUGSoutput$summary),
                               split = "[", fixed = TRUE), `[`, 1, FUN.VALUE=character(1))
    } else {
      names <- rownames(object$BUGSoutput$summary)
    }
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
        warning(paste0("\"", excl[i], "\"", " not found in MCMC output. Check 'ISB' and 'exact' arguments to make sure the desired parsing methods are being used."))
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
      if (exact == TRUE)
      {
        get_ind <- which(names %in% params)
      } else {
        get_ind <- grep(paste(params), names, fixed = FALSE)
      }

      if (length(get_ind) < 1)
      {
        stop(paste0("\"", params, "\"", " not found in MCMC output. Check `ISB` and `exact` arguments to make sure the desired parsing methods are being used."))
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
      if (exact == TRUE)
      {
        get_ind <- which(names %in% params[i])
      } else {
        get_ind <- grep(paste(params[i]), names, fixed = FALSE)
      }

      if (length(get_ind) < 1)
      {
        warning(paste0("\"", params[i], "\"", " not found in MCMC output. Check 'ISB' and 'exact' arguments to make sure the desired parsing methods are being used."))
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

  #PROCESSING BLOCK
  if (is.null(chain_num))
  {
    if (coda::is.mcmc.list(object) == TRUE | typeof(object) == 'S4')
    {
      if (length(f_ind) > 1)
      {
        dsort_mcmc <- do.call(coda::mcmc.list, temp_in[,f_ind])
        OUT <- do.call('rbind', dsort_mcmc)
      } else {
        dsort_mcmc <- do.call(coda::mcmc.list, temp_in[,f_ind, drop = FALSE])
        OUT <- as.matrix(do.call(coda::mcmc.list, temp_in[,f_ind, drop = FALSE]), ncol = 1)
      }
    }
    if (methods::is(object, 'matrix'))
    {
      OUT <- temp_in[,f_ind, drop = FALSE]
      if (mcmc.list == TRUE)
      {
        stop('Cannot produce mcmc.list output with matrix input')
      }
    }

    if (methods::is(object, 'rjags'))
    {
      OUT <- temp_in[,f_ind, drop = FALSE]
      if (mcmc.list == TRUE)
      {
        #modified coda::as.mcmc (removing ordering of param names)
        x <- object$BUGSoutput
        mclist <- vector("list", x$n.chains)
        mclis <- vector("list", x$n.chains)
        ord <- dimnames(x$sims.array)[[3]]
        for (i in 1:x$n.chains)
        {
          tmp1 <- x$sims.array[, i, ord]
          mclis[[i]] <- coda::mcmc(tmp1, thin = x$n.thin)
        }
        temp2 <- coda::as.mcmc.list(mclis)
        #end mod as.mcmc
        dsort_mcmc <- do.call(coda::mcmc.list, temp2[,f_ind, drop = FALSE])
      }
    }
  }

  if (!is.null(chain_num))
  {
    if (coda::is.mcmc.list(object) == TRUE | typeof(object) == 'S4')
    {
      if (length(f_ind) > 1)
      {
        dsort <- do.call(coda::mcmc.list, temp_in[,f_ind])

        if (chain_num > length(dsort))
        {
          stop('Invalid value for chain_num specified.')
        }
        dsort_mcmc <- dsort[[chain_num]]
        OUT <- as.matrix(dsort_mcmc)
      } else {
        dsort <- do.call(coda::mcmc.list, temp_in[,f_ind, drop = FALSE])

        if (chain_num > length(dsort))
        {
          stop('Invalid value for chain_num specified.')
        }
        dsort_mcmc <- dsort[[chain_num]]
        OUT <- as.matrix(dsort_mcmc)
      }
    }

    if (methods::is(object, 'matrix'))
    {
      stop('Cannot extract posterior information for individual chains from matrix input.')
    }

    if (methods::is(object, 'rjags'))
    {
      #modified coda::as.mcmc (removing ordering of param names)
      x <- object$BUGSoutput
      mclist <- vector("list", x$n.chains)
      mclis <- vector("list", x$n.chains)
      ord <- dimnames(x$sims.array)[[3]]
      for (i in 1:x$n.chains)
      {
        tmp1 <- x$sims.array[, i, ord]
        mclis[[i]] <- coda::mcmc(tmp1, thin = x$n.thin)
      }
      temp2 <- coda::as.mcmc.list(mclis)
      #end mod as.mcmc
      dsort <- do.call(coda::mcmc.list, temp2[,f_ind, drop = FALSE])
      if (chain_num > length(dsort))
      {
        stop('Invalid value for chain_num specified.')
      }
      dsort_mcmc <- dsort[[chain_num]]
      OUT <- as.matrix(dsort_mcmc)
    }
  }

  if (mcmc.list == FALSE)
  {
    return(OUT)
  }
  if (mcmc.list == TRUE)
  {
    return(dsort_mcmc)
  }
}

