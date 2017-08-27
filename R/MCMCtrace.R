#' Trace and density plots from MCMC output
#'
#' Trace and density plots of MCMC chains for specific parameters of interest. Option to
#' print plots to pdf.
#'
#' @param object Object containing MCMC output. See DETAILS below.
#' @param params Character string (or vector of character strings) denoting parameters of interest.
#'
#' Default \code{'all'} returns chains for all parameters.
#'
#' @param excl Character string (or vector of character strings) denoting parameters to exclude. Used in conjunction with \code{params} argument to select parameters of interest.
#'
#' @param ISB Ignore Square Brackets (ISB). Logical specifying whether square brackets should be ignored in the \code{params} and \code{excl} arguments. If \code{FALSE}, square brackets are ignored - input from \code{params} and \code{excl} are otherwise matched exactly. If \code{TRUE}, square brackets are not ignored - input from \code{params} and \code{excl} are matched using grep, allowing partial names to be used when specifying parameters of interest.
#'
#' @param iter Number of iterations to plot for trace and density plots. The default value is 5000, meaning the last 5000 iterations of the chain will be plotted.
#'
#' @param pdf Logical - if \code{pdf = TRUE} plots will be exported to a pdf.
#' @param filename Name of pdf file to be printed. Default is 'MCMCtrace'.
#' @param wd Working directory for pdf output. Default is current directory.
#' @param type Type of plot to be output. \code{'both'} outputs both trace and density plots, \code{'trace'}
#' outputs only trace plots, and \code{'density'} outputs only density plots.
#' @param ind Logical - if \code{ind = TRUE}, separate density lines will be plotted for each chain. If
#' \code{ind= FALSE}, one density line will be plotted for all chains.
#' @section Details:
#' \code{object} argument can be a \code{stanfit} object (\code{rstan} package), an \code{mcmc.list} object
#' (\code{coda} package), or an \code{R2jags} model object (\code{R2jags} package). The function automatically
#' detects the object type and proceeds accordingly.
#'
#'
#' @examples
#' #Load data
#' data(MCMC_data)
#'
#' #Traceplots for all 'beta' parameters
#' MCMCtrace(MCMC_data, params='beta')
#'
#' #Traceplots (individual density lines for each chain) for just 'beta[1]'
#' MCMCtrace(MCMC_data, params = 'beta[1]', ISB = FALSE, filename = 'PDF_file.pdf', ind = TRUE)
#'
#' @export

MCMCtrace <- function(object,
                    params = 'all',
                    excl = NULL,
                    ISB = TRUE,
                    iter = 5000,
                    pdf = TRUE,
                    filename,
                    wd = getwd(),
                    type = 'both',
                    ind = FALSE)
{

  .pardefault <- graphics::par(no.readonly = T)

  if(typeof(object) == 'S4')
  {
    temp <- rstan::As.mcmc.list(object)
  }
  if(coda::is.mcmc.list(object) == TRUE)
  {
    temp <- object
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
    temp <- coda::as.mcmc.list(mclis)
  }
  if(coda::is.mcmc.list(object) == FALSE & typeof(object) != 'list' & typeof(object) != 'S4')
  {
    stop('Invalid input type. Object must be of type stanfit, mcmc.list, or R2jags.')
  }
  if(coda::is.mcmc.list(object) == FALSE & typeof(object) == 'list' & class(object[[1]]) == 'mcarray')
  {
    stop('Invalid object type. jags.samples objects not currently supported. Input must be stanfit object, mcmc.list object, rjags object, or matrix with MCMC chains.')
  }

  #NAME SORTING BLOCK
  if(ISB == TRUE)
  {
    names <- vapply(strsplit(colnames(temp[[1]]),
                             split = "[", fixed = TRUE), `[`, 1, FUN.VALUE=character(1))
    a_names <- colnames(temp[[1]])
  }else{
    names <- colnames(temp[[1]])
    a_names <- names
  }
  n_chains <- length(temp)
  if (nrow(temp[[1]]) > iter)
  {
    it <- (nrow(temp[[1]]) - iter) : nrow(temp[[1]])
  }else {
    it <- 1 : nrow(temp[[1]])
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
    }else
    {
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

  #OUTPUT BLOCK
  if(pdf == TRUE)
  {
    setwd(wd)
    if(missing(filename))
    {
      file_out <- 'MCMCtrace.pdf'
    }else{
      if(grepl('.pdf', filename, fixed = TRUE))
      {
        file_out <- paste0(filename)
      }else{
        file_out <- paste0(filename, '.pdf')
      }
    }
    pdf(file= file_out)
  }


  #PLOT BLOCK
  A_VAL <- 0.5 #alpha value
  graphics::layout(matrix(c(1, 2, 3, 4, 5, 6), 3, 2, byrow = TRUE))
  graphics::par(mar = c(4.1,4.1,2.1,1.1)) # bottom, left, top, right
  graphics::par(mgp = c(2.5,1,0)) #axis text distance
  gg_color_hue <- function(n)
  {
    hues = seq(15, 375, length = n+1)
    grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
  }
  colors <- gg_color_hue(n_chains)
  gg_cols <- grDevices::col2rgb(colors)/255

  if (type == 'both')
  {
    for (j in 1: length(f_ind))
    {
      #trace
      tmlt <- do.call('cbind', temp[it, f_ind[j]])
      graphics::matplot(it, tmlt, lwd = 1, lty= 1, type='l', main = paste0('Trace - ', a_names[f_ind[j]]),
              col= grDevices::rgb(red= gg_cols[1,], green= gg_cols[2,],
                       blue= gg_cols[3,], alpha = A_VAL),
              xlab= 'Iteration', ylab= 'Value')
      if (ind == TRUE & n_chains > 1)
      {
        dens <- apply(tmlt, 2, stats::density)
        max_den <- c()
        for (k in 1:NCOL(tmlt))
        {
          max_den <- c(max_den, max(dens[[k]]$y))
        }
        ylim <- c(0, max(max_den))

        graphics::plot(dens[[1]], xlab = 'Parameter estimate', ylim = ylim,
             lty = 1, lwd = 1, main = paste0('Density - ', a_names[f_ind[j]]),
             col = grDevices::rgb(red= gg_cols[1,1], green= gg_cols[2,1], blue= gg_cols[3,1]))

        for (l in 2:NCOL(tmlt))
        {
          graphics::lines(dens[[l]],
                col = grDevices::rgb(red= gg_cols[1,l], green= gg_cols[2,l],
                          blue= gg_cols[3,l]))
        }
      }else{
        #density plot
        graphics::plot(stats::density(rbind(tmlt)), xlab = 'Parameter estimate',
             lty = 1, lwd = 1, main = paste0('Density - ', a_names[f_ind[j]]))
      }
    }
  }

  if (type == 'trace')
  {
    for (j in 1: length(f_ind))
    {
      #chains
      tmlt <- do.call('cbind', temp[it,f_ind[j]])
      graphics::matplot(it, tmlt, lwd = 1, lty= 1, type='l', main = paste0('Trace - ', a_names[f_ind[j]]),
              col= grDevices::rgb(red= gg_cols[1,], green= gg_cols[2,],
                       blue= gg_cols[3,], alpha = A_VAL),
              xlab= 'Iteration', ylab= 'Value')
    }
  }

  if (type == 'density')
  {
    for (j in 1: length(f_ind))
    {
      #trace
      tmlt <- do.call('cbind', temp[it,f_ind[j]])

      if (ind == TRUE & n_chains > 1)
      {
        dens <- apply(tmlt, 2, stats::density)
        max_den <- c()
        for (k in 1:NCOL(tmlt))
        {
          max_den <- c(max_den, max(dens[[k]]$y))
        }
        ylim <- c(0, max(max_den))

        graphics::plot(dens[[1]], xlab = 'Parameter estimate', ylim = ylim,
             lty = 1, lwd = 1, main = paste0('Density - ', a_names[f_ind[j]]),
             col = grDevices::rgb(red= gg_cols[1,1], green= gg_cols[2,1], blue= gg_cols[3,1]))

        for (l in 2:NCOL(tmlt))
        {
          graphics::lines(dens[[l]],
                col = grDevices::rgb(red= gg_cols[1,l], green= gg_cols[2,l],
                          blue= gg_cols[3,l]))
        }
      }else{
        #density plot
        graphics::plot(stats::density(rbind(tmlt)), xlab = 'Parameter estimate',
             lty = 1, lwd = 1, main = paste0('Density - ', a_names[f_ind[j]]))
      }
    }
  }

  if (type != 'both' & type != 'density' & type != 'trace')
  {
    stop('Invalid argument for "type". Valid inputs are "both", "trace", and "density".')
  }

  if(pdf == TRUE)
  {
    invisible(grDevices::dev.off())
    system(paste0('open ', paste0('"', file_out, '"')))
  }else{
    graphics::par(.pardefault)
  }
}
