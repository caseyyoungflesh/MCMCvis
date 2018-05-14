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
#' @param ISB Ignore Square Brackets (ISB). Logical specifying whether square brackets should be ignored in the \code{params} and \code{excl} arguments. If \code{TRUE}, square brackets are ignored - input from \code{params} and \code{excl} are otherwise matched exactly. If \code{FALSE}, square brackets are not ignored - input from \code{params} and \code{excl} are matched using regular expression format. This allows partial names to be used when specifying parameters of interest.
#'
#' @param mcmc.list Logical specifying whether to return an mcmc.list. If \code{TRUE}, an \code{mcmc.list} object is returned, rather than a matrix.
#'
#' @param chain_num Numeric - specifies posterior chain number. When a value is specified, posterior for only that chain is output. Useful for determining the last iteration for each parameter, to be used as initial values in a subsequent model, to effectively 'continue' a model run.
#'
#' @section Details:
#' Function returns matrix with one chain per column for specified parameters. Multiple input chains for each parameter are combined to one posterior chain. Parameters are arranged in columns alphabetically.
#'
#' \code{object} argument can be a \code{stanfit} object (\code{rstan} package), an \code{mcmc.list} object (\code{coda} package), an \code{R2jags} model object (\code{R2jags} package), a \code{jagsUI} model object (\code{jagsUI} package), or a matrix containing MCMC chains (each column representing MCMC output for a single parameter, rows representing iterations in the chain). The function automatically detects the object type and proceeds accordingly.
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
#' #'params' takes regular expressions when ISB = FALSE, square brackets must be escaped with '\\'
#' ex3 <- MCMCchains(MCMC_data, params = c('beta\\[1\\]', 'beta\\[4\\]', 'alpha\\[3\\]'), ISB = FALSE)
#' apply(ex3, 2, sd)
#'
#' @export

MCMCchains <- function(object,
                     params = 'all',
                     excl = NULL,
                     ISB = TRUE,
                     mcmc.list = FALSE,
                     chain_num = NULL)
{
  if(length(class(object)) > 1)
  {
    #modified coda::as.mcmc (removing ordering of param names)
    x <- object$BUGSoutput
    mclist <- vector("list", x$n.chains)
    mclis <- vector("list", x$n.chains)
    strt <- x$n.burnin + 1
    end <- x$n.iter
    ord <- dimnames(x$sims.array)[[3]]
    for (i in 1:x$n.chains)
    {
      tmp1 <- x$sims.array[, i, ord]
      mclis[[i]] <- coda::mcmc(tmp1, start = strt, end = end, thin = x$n.thin)
    }
    object <- coda::as.mcmc.list(mclis)
    #end mod as.mcmc
  }
  if(coda::is.mcmc.list(object) != TRUE &
     typeof(object) != 'double' &
     class(object) != 'rjags' &
     typeof(object) != 'S4' &
     class(object) != 'jagsUI')
  {
    stop('Invalid object type. Input must be stanfit object (rstan), mcmc.list object (coda), rjags object (R2jags), jagsUI object (jagsUI), or matrix with MCMC chains.')
  }

  #NAME SORTING BLOCK
  if(typeof(object) == 'S4')
  {
    temp_in <- rstan::As.mcmc.list(object)
    if(ISB == TRUE)
    {
      names <- vapply(strsplit(colnames(temp_in[[1]]),
                               split = "[", fixed = TRUE), `[`, 1, FUN.VALUE=character(1))
    }else{
      names <- colnames(temp_in[[1]])
    }
  }

  if(class(object) == 'jagsUI')
  {
    object <- object$samples
  }

  if(coda::is.mcmc.list(object) == TRUE)
  {
    temp_in <- object
    if(ISB == TRUE)
    {
      names <- vapply(strsplit(colnames(temp_in[[1]]),
                               split = "[", fixed = TRUE), `[`, 1, FUN.VALUE=character(1))
    }else{
      names <- colnames(temp_in[[1]])
    }
  }

  if(typeof(object) == 'double')
  {
    temp_in <- object
    if(ISB == TRUE)
    {
      names <- vapply(strsplit(colnames(temp_in),
                               split = "[", fixed = TRUE), `[`, 1, FUN.VALUE=character(1))
    }else{
      names <- colnames(temp_in)
    }
  }

  if(class(object) == 'rjags')
  {
    temp_in <- object$BUGSoutput$sims.matrix
    if(ISB == TRUE)
    {
      names <- vapply(strsplit(rownames(object$BUGSoutput$summary),
                               split = "[", fixed = TRUE), `[`, 1, FUN.VALUE=character(1))
    }else{
      names <- rownames(object$BUGSoutput$summary)
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
        ind_excl <- which(names %in% n_excl[i])
        if(length(ind_excl) < 1)
        {
          warning(paste0('"', excl[i], '"', ' not found in MCMC output.'))
        }
        rm_ind <- c(rm_ind, ind_excl)
      }else{
        n_excl <- excl
        ind_excl <- grep(n_excl[i], names, fixed = FALSE)
        if(length(ind_excl) < 1)
        {
          warning(paste0('"', excl[i], '"', ' not found in MCMC output.'))
        }
        rm_ind <- c(rm_ind, ind_excl)
      }
    }
    if(length(rm_ind) > 0)
    {
      dups <- which(duplicated(rm_ind))
      if(length(dups) > 0)
      {
        rm_ind2 <- rm_ind[-dups]
      }else{
        rm_ind2 <- rm_ind
      }
    }else{
      excl <- NULL
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
        get_ind <- grep(paste(params), names, fixed = FALSE)
      }

      if (length(get_ind) < 1)
      {
        stop(paste0('"', params, '"', ' not found in MCMC output.'))
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
        get_ind <- grep(paste(params[i]), names, fixed=FALSE)
      }

      if (length(get_ind) < 1)
      {
        warning(paste0('"', params[i], '"', ' not found in MCMC output.'))
        next()
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
  if (is.null(chain_num))
  {
    if(coda::is.mcmc.list(object) == TRUE | typeof(object) == 'S4')
    {
      if(length(f_ind) > 1)
      {
        dsort_mcmc <- do.call(coda::mcmc.list, temp_in[,f_ind])
        OUT <- do.call('rbind', dsort_mcmc)
      }else{
        dsort_mcmc <- do.call(coda::mcmc.list, temp_in[,f_ind, drop = FALSE])
        OUT <- as.matrix(do.call(coda::mcmc.list, temp_in[,f_ind, drop = FALSE]), ncol = 1)
      }
    }
    if(typeof(object) == 'double')
    {
      OUT <- temp_in[,f_ind, drop = FALSE]
      if(mcmc.list == TRUE)
      {
        stop('Cannot produce mcmc.list output with matrix input')
      }
    }

    if(class(object) == 'rjags')
    {
      OUT <- temp_in[,f_ind, drop = FALSE]
      if(mcmc.list == TRUE)
      {
        #modified coda::as.mcmc (removing ordering of param names)
        x <- object$BUGSoutput
        mclist <- vector("list", x$n.chains)
        mclis <- vector("list", x$n.chains)
        strt <- x$n.burnin + 1
        end <- x$n.iter
        ord <- dimnames(x$sims.array)[[3]]
        for (i in 1:x$n.chains)
        {
          tmp1 <- x$sims.array[, i, ord]
          mclis[[i]] <- coda::mcmc(tmp1, start = strt, end = end, thin = x$n.thin)
        }
        temp2 <- coda::as.mcmc.list(mclis)
        #end mod as.mcmc
        dsort_mcmc <- do.call(coda::mcmc.list, temp2[,f_ind, drop = FALSE])
      }
    }
  }

  if (!is.null(chain_num))
  {
    if(coda::is.mcmc.list(object) == TRUE | typeof(object) == 'S4')
    {
      if(length(f_ind) > 1)
      {
        dsort <- do.call(coda::mcmc.list, temp_in[,f_ind])

        if (chain_num > length(dsort))
        {
          stop('Invalid value for chain_num specified.')
        }
        dsort_mcmc <- dsort[[chain_num]]
        OUT <- as.matrix(dsort_mcmc)
      }else{
        dsort <- do.call(coda::mcmc.list, temp_in[,f_ind, drop = FALSE])

        if (chain_num > length(dsort))
        {
          stop('Invalid value for chain_num specified.')
        }
        dsort_mcmc <- dsort[[chain_num]]
        OUT <- as.matrix(dsort_mcmc)
      }
    }

    if(typeof(object) == 'double')
    {
      stop('Cannot extract posterior information for individual chains from matrix input.')
    }

    if(class(object) == 'rjags')
    {
      #modified coda::as.mcmc (removing ordering of param names)
      x <- object$BUGSoutput
      mclist <- vector("list", x$n.chains)
      mclis <- vector("list", x$n.chains)
      strt <- x$n.burnin + 1
      end <- x$n.iter
      ord <- dimnames(x$sims.array)[[3]]
      for (i in 1:x$n.chains)
      {
        tmp1 <- x$sims.array[, i, ord]
        mclis[[i]] <- coda::mcmc(tmp1, start = strt, end = end, thin = x$n.thin)
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

  if(mcmc.list == FALSE)
  {
    return(OUT)
  }
  if(mcmc.list == TRUE)
  {
    return(dsort_mcmc)
  }
}

