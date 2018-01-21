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
#' #Just 'beta[1]', 'gamma[4]', and 'alpha[3]'
#' #'params' takes regular expressions when ISB = FALSE, square brackets must be escaped with '\\'
#' ex3 <- MCMCchains(MCMC_data, params = c('beta\\[1\\]', 'gamma\\[4\\]', 'alpha\\[3\\]'), ISB = FALSE)
#' apply(ex3, 2, sd)
#'
#' @export

MCMCchains <- function(object,
                     params = 'all',
                     excl = NULL,
                     ISB = TRUE,
                     mcmc.list = FALSE)
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
    stop('Invalid object type. Input must be stanfit object (rstan), mcmc.list object (coda),
         rjags object (R2jags), jagsUI object (jagsUI), or matrix with MCMC chains.')
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

  if(class(object[[1]]) == 'mcarray')
  {
    stop('Invalid object type. jags.samples objects not currently supported. Input must be stanfit object, mcmc.list object, rjags object, jagsUI object, or matrix with MCMC chains.')
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
        rm_ind <- c(rm_ind, grep(n_excl[i], names, fixed = FALSE))
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
        get_ind <- grep(paste(params), names, fixed = FALSE)
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
        get_ind <- grep(paste(params[i]), names, fixed=FALSE)
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
  if(coda::is.mcmc.list(object) == TRUE)
  {
    if(length(f_ind) > 1)
    {
      dsort <- do.call(coda::mcmc.list, temp_in[,f_ind])
      OUT <- do.call('rbind', dsort)
    }else{
      dsort <- do.call(coda::mcmc.list, temp_in[,f_ind, drop = FALSE])
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
      dsort <- do.call(coda::mcmc.list, temp2[,f_ind, drop = FALSE])
    }
  }

  if(mcmc.list == FALSE)
  {
    return(OUT)
  }

  if(mcmc.list == TRUE)
  {
    return(dsort)
  }
}

