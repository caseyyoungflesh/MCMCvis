#' Trace and density plots from MCMC output
#'
#' Trace and density plots of MCMC chains for specific parameters of interest. Option to
#' print plots to pdf.
#'
#' @param object Object containing MCMC output. See DETAILS below.
#' @param params Character string (or vector of character strings) denoting parameters of interest.
#' Partial names may be used to return all parameters containing that set of characters.
#'
#' Default \code{'all'} returns chains for all parameters.
#'
#' @param excl Character string (or vector of character strings) denoting parameters to exclude.
#' Partial names may be used to exclude all parameters containing that set of characters. Used in
#' conjunction with \code{params} argument to select parameters of interest.
#'
#' @param iter Number of iterations to plot for trace and density plots. The default value is 2000,
#'  meaning the last 2000 iterations of the chain will be plotted.
#'
#' @param pdf Logical - if \code{pdf = TRUE} plots will be exported to a pdf.
#' @param filename Name of pdf file to be printed.
#' @param wd Working directory for pdf output. Default is current directory.
#' @param type Type of plot to be output. \code{'both'} outputs both trace and density plots, \code{'trace'}
#' outputs only trace plots, and \code{'density'} outputs only density plots.
#' @param ind Logical - if \code{ind = TRUE}, different lines will be plotted for each chain. If
#' \code{ind= FALSE}, one line will be plotted for all chains.
#' @section Details:
#' \code{object} argument can be a \code{stanfit} object (\code{rstan} package), an \code{mcmc.list} object
#' (\code{coda} package), or an \code{R2jags} model object (\code{R2jags} package). The function automatically
#' detects the object type and proceeds accordingly.
#'
#' @examples
#' #Load data
#' data(MCMC_data)
#'
#' #Traceplot for all 'beta' parameters
#' MCMCtrace(MCMC_data, params='beta')
#'
#' #Print traceplot output to pdf
#' MCMCtrace(MCMC_data, pdf= TRUE, filename = 'PDF_file.pdf')
#'
#' @export

MCMCtrace <- function(object,
                    params = 'all',
                    excl = NULL,
                    iter = 2000,
                    pdf = FALSE,
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
  if(coda::is.mcmc.list(object) == FALSE & typeof(object) == 'list')
  {
    #modified coda::as.mcmc (removing ordering of param names)
    x <- object$BUGSoutput
    mclist <- vector("list", x$n.chains)
    mclis <- vector("list", x$n.chains)
    strt <- x$n.burnin + 1
    end <- x$n.iter
    ord <- dimnames(x$sims.array)[[3]]
    for (i in 1:x$n.chains) {
      tmp1 <- x$sims.array[, i, ord]
      mclis[[i]] <- coda::mcmc(tmp1, start = strt, end = end, thin = x$n.thin)
    }
    temp <- coda::as.mcmc.list(mclis)
  }
  if(coda::is.mcmc.list(object) == FALSE & typeof(object) != 'list' & typeof(object) != 'S4')
  {
    stop('Invalid input type. Object must be of type stanfit, mcmc.list, or R2jags.')
  }

  names <- colnames(temp[[1]])
  n_chains <- length(temp)
  if (nrow(temp[[1]]) > iter)
  {
    it <- (nrow(temp[[1]]) - iter) : nrow(temp[[1]])
  }else {
    it <- 1 : nrow(temp[[1]])
  }

  if(!is.null(excl))
  {
    to.rm1 <- c()
    for (i in 1:length(excl))
    {
      to.rm1 <- c(to.rm1, grep(excl[i], names, fixed = TRUE))
    }
    dups <- -which(duplicated(to.rm1))
    if(length(dups) > 0)
    {
      to.rm2 <- to.rm1[-dups]
    }else{
      to.rm2 <- to.rm1
    }
  }


  if (length(params) == 1)
  {
    if (params == 'all')
    {
      if(is.null(excl))
      {
        g_filt <- 1:length(names)
      }else{
        g_filt <- (1:length(names))[-to.rm2]
      }
    }else
    {
      get.cols <- grep(paste(params), names, fixed=TRUE)
      if (length(get.cols) < 1)
      {
        stop(paste0('"', params, '"', ' not found in MCMC ouput.'))
      }

      if(!is.null(excl))
      {
        if(identical(get.cols, to.rm2))
        {
          stop('No parameters selected.')
        }

        matched <- which(get.cols == to.rm2)
        if (length(matched) > 0)
        {
          g_filt <- get.cols[-matched]
        }else {
          g_filt <- get.cols
        }

      }else{
        g_filt <- get.cols
      }
    }
  }else
  {
    grouped <- c()
    for (i in 1:length(params))
    {
      get.cols <- grep(paste(params[i]), names, fixed=TRUE)
      grouped <- c(grouped, get.cols)
    }

    if(!is.null(excl))
    {
      if(identical(grouped, to.rm2))
      {
        stop('No parameters selected.')
      }

      matched <- stats::na.omit(match(to.rm2, grouped))
      if (length(matched) > 0)
      {
        cols <- grouped[-matched]
      } else{
        cols <- grouped
      }

      to.rm <- which(duplicated(cols))
      if(length(to.rm) > 0)
      {
        g_filt <- cols[-to.rm]
      }else
      {
        g_filt <- cols
      }
    } else{

      to.rm <- which(duplicated(grouped))
      if(length(to.rm) > 0)
      {
        g_filt <- grouped[-to.rm]
      }else
      {
        g_filt <- grouped
      }
    }
  }


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
    for (j in 1: length(g_filt))
    {
      #trace
      tmlt <- do.call('cbind', temp[it, g_filt[j]])
      graphics::matplot(it, tmlt, lwd = 1, lty= 1, type='l', main = paste0('Trace - ', names[g_filt[j]]),
              col= grDevices::rgb(red= gg_cols[1,], green= gg_cols[2,],
                       blue= gg_cols[3,], alpha = 0.5),
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
             lty = 1, lwd = 1, main = paste0('Density - ', names[g_filt[j]]),
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
             lty = 1, lwd = 1, main = paste0('Density - ', names[g_filt[j]]))
      }
    }
  }

  if (type == 'trace')
  {
    for (j in 1: length(g_filt))
    {
      #chains
      tmlt <- do.call('cbind', temp[it,g_filt[j]])
      graphics::matplot(it, tmlt, lwd = 1, lty= 1, type='l', main = paste0('Trace - ', names[g_filt[j]]),
              col= grDevices::rgb(red= gg_cols[1,], green= gg_cols[2,],
                       blue= gg_cols[3,], alpha = 0.5),
              xlab= 'Iteration', ylab= 'Value')
    }
  }

  if (type == 'density')
  {
    for (j in 1: length(g_filt))
    {
      #trace
      tmlt <- do.call('cbind', temp[it,g_filt[j]])

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
             lty = 1, lwd = 1, main = paste0('Density - ', names[g_filt[j]]),
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
             lty = 1, lwd = 1, main = paste0('Density - ', names[g_filt[j]]))
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
  }

  graphics::par(.pardefault)

}
