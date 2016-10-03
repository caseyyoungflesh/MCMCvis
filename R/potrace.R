#' Trace and density plots from MCMC output
#'
#' Trace plots and density plots of MCMC chains for specific parameters of interest. Option to
#' print plots to pdf.
#'
#' @param object Object containing MCMC output. See DETAILS below.
#' @param par Character string (or vector of character strings) denoting parameters of interest.
#' Partial names may be used to return all parameters containing that set of characters.
#'
#' Default \code{'all'} returns chains for all parameters.
#'
#' @param pdf Logical - if \code{pdf = TRUE} plots will be exported to a pdf.
#' @param filename Name of pdf file to be printed.
#' @param wd Working directory for pdf output. Default is current directory.
#' @param type Type of plot to be output. \code{'both'} outputs both trace and density plots, \code{'trace'}
#' outputs only trace plots, and \code{'density'} outputs only density plots.
#' @param ind Logical - if \code{ind = TRUE}, different lines will be plotted for each chain. If
#' \code{ind= FALSE}, one line will be plotted for all chains.
#' @section Details:
#' Plots created similar to that of \code{traceplot} from \code{coda} package.
#'
#' \code{object} argument can be a \code{stanfit} object (\code{rstan} package), an \code{mcmc.list} object
#' (\code{coda} package), or an \code{R2jags} model object (\code{R2jags} package).
#'
#' @examples
#' #Load data
#' data(MCMC_data)
#'
#' #Traceplot for all 'beta' parameters
#' potrace(MCMC_data, par='beta')
#'
#' #Print traceplot output to pdf
#' potrace(MCMC_data, pdf= TRUE, filename = 'PDF_file.pdf')
#'
#' @export

potrace <- function(object,
                    par = 'all',
                    pdf = FALSE,
                    filename,
                    wd = getwd(),
                    type = 'both',
                    ind = FALSE)
{

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
    x <- object$BUGSoutput
    mclist <- vector("list", x$n.chains)
    mclis <- vector("list", x$n.chains)
    strt <- x$n.burnin + 1
    end <- x$n.iter
    ord <- dimnames(x$sims.array)[[3]]
    for (i in 1:x$n.chains) {
      tmp1 <- x$sims.array[, i, ord]
      mclis[[i]] <- mcmc(tmp1, start = strt, end = end, thin = x$n.thin)
    }
    temp <- as.mcmc.list(mclis)
  }
  if(coda::is.mcmc.list(object) == FALSE & typeof(object) != 'list' & typeof(object) != 'S4')
  {
    stop('Invalid input type. Object must be of type stanfit, mcmc.list, or R2jags.')
  }

  names <- colnames(temp[[1]])
  n_chains <- length(temp)
  it <- 1:nrow(temp[[1]])

  if (length(par) == 1)
  {
    if (par == 'all')
    {
      g_filt <- 1:length(names)
    }else
    {
      g_filt <- grep(paste(par), names, fixed=TRUE)
    }
  }else
  {
    grouped <- c()
    for (i in 1:length(par))
    {
      get.cols <- grep(paste(par[i]), names, fixed=TRUE)
      grouped <- c(grouped, get.cols)
    }

    to.rm <- which(duplicated(grouped))
    if(length(to.rm) >0)
    {
      g_filt <- grouped[-to.rm]
    }else
    {
      g_filt <- grouped
    }
  }

  if(pdf == TRUE)
  {
    setwd(wd)
    if(missing(filename))
    {
      file_out <- 'potrace.pdf'
    }else{
      if(grepl('pdf', filename))
      {
        file_out <- paste0(filename)
      }else{
        file_out <- paste0(filename, '.pdf')
      }
    pdf(file= file_out)
    }
  }

  .pardefault <- par(no.readonly = T)

  layout(matrix(c(1, 2, 3, 4, 5, 6), 3, 2, byrow = TRUE))
  par(mar=c(4.1,4.1,2.1,1.1)) # bottom, left, top, right
  par(mgp=c(2.5,1,0)) #axis text distance
  gg_color_hue <- function(n)
  {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
  }
  colors <- gg_color_hue(n_chains)
  gg_cols <- col2rgb(colors)/255

  if (type == 'both')
  {
    for (j in 1: length(g_filt))
    {
      #trace
      tmlt <- do.call('cbind', temp[,g_filt[j]])
      matplot(it, tmlt, lwd = 1, lty= 1, type='l', main = paste0('Trace - ', names[g_filt[j]]),
              col= rgb(red= gg_cols[1,], green= gg_cols[2,],
                       blue= gg_cols[3,], alpha = 0.6),
              xlab= 'Iteration', ylab= 'Value')
      if (ind == TRUE & n_chains > 1)
      {
        dens <- apply(tmlt, 2, density)
        max_den <- c()
        for (k in 1:NCOL(tmlt))
        {
          max_den <- c(max_den, max(dens[[k]]$y))
        }
        ylim <- c(0, max(max_den))

        plot(dens[[1]], xlab = 'Parameter estimate', ylim = ylim,
             lty = 1, lwd = 1, main = paste0('Density - ', names[g_filt[j]]),
             col = rgb(red= gg_cols[1,1], green= gg_cols[2,1], blue= gg_cols[3,1]))

        for (l in 2:NCOL(tmlt))
        {
          lines(dens[[l]],
                col = rgb(red= gg_cols[1,l], green= gg_cols[2,l],
                          blue= gg_cols[3,l]))
        }
      }else{
        #density plot
        plot(density(rbind(tmlt)), xlab = 'Parameter estimate',
             lty = 1, lwd = 1, main = paste0('Density - ', names[g_filt[j]]))
      }
    }
  }

  if (type == 'trace')
  {
    for (j in 1: length(g_filt))
    {
      #chains
      tmlt <- do.call('cbind', temp[,g_filt[j]])
      matplot(it, tmlt, lwd = 1, lty= 1, type='l', main = paste0('Trace - ', names[g_filt[j]]),
              col= rgb(red= gg_cols[1,], green= gg_cols[2,],
                       blue= gg_cols[3,], alpha = 0.6),
              xlab= 'Iteration', ylab= 'Value')
    }
  }

  if (type == 'density')
  {
    for (j in 1: length(g_filt))
    {
      #trace
      tmlt <- do.call('cbind', temp[,g_filt[j]])

      if (ind == TRUE & n_chains > 1)
      {
        dens <- apply(tmlt, 2, density)
        max_den <- c()
        for (k in 1:NCOL(tmlt))
        {
          max_den <- c(max_den, max(dens[[k]]$y))
        }
        ylim <- c(0, max(max_den))

        plot(dens[[1]], xlab = 'Parameter estimate', ylim = ylim,
             lty = 1, lwd = 1, main = paste0('Density - ', names[g_filt[j]]),
             col = rgb(red= gg_cols[1,1], green= gg_cols[2,1], blue= gg_cols[3,1]))

        for (l in 2:NCOL(tmlt))
        {
          lines(dens[[l]],
                col = rgb(red= gg_cols[1,l], green= gg_cols[2,l],
                          blue= gg_cols[3,l]))
        }
      }else{
        #density plot
        plot(density(rbind(tmlt)), xlab = 'Parameter estimate',
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
    invisible(dev.off())
    system(paste0('open ', paste0('"', file_out, '"')))
  }
  par(.pardefault)
}
