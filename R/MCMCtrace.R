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

object <- out_R2jags
params = c('alpha', 'beta')
excl = NULL
ISB = TRUE
iter = 5000
pdf = FALSE
wd = getwd()
type = 'both'
ind = FALSE

#priors used in R2jags model
priors <- cbind(rnorm(15000, 0, 31.6), rnorm(15000, 0, 31.6))
ab <- MCMCchains(object, params = c('alpha', 'beta'))

#calc overlap - see paper that notes for mark recapture using uniform prior overlap < 0.3 indicate robust identifiability

tt <- list(priors[,1], ab[,1])

f1 <- rnorm(1000, 0, 1)
f2 <- rnorm(1000, 0, 5)
tt <- list(f1, f2)

t1 <- proc.time()
overlapping::overlap(tt, nbins= 1000, plot = FALSE)
proc.time() - t1


#CASES:
#1 prior, 2 parms
  #more iters
  #less iters
  #equal iters
#2 priors, 2 params
#2 priors, 3 params
#3 priors, 2 params
#2 priors, 1 param
#1 prior, 1 param

PP <- cbind(rnorm(20000, 0, 31.6))
MCMCtrace(object, params = 'alpha', priors = PP, pdf = FALSE, type = 'both')

priors = (PP)
#If nothing is specified for priors, no prior is plotted and overlap is not calculated

#priors need to be in a matrix format (with each parameter in a different column)
#the number of iterations for the prior must match that of the number of iterations plotted with MCMCtrace - default is 5000, though this can be changed using the `iter` argument
#these matrices can be generate using rnorm, rgamma, runif, etc. in R. Distirbutions not supported by base R can be used by loaded the appropriate packages. It is important to note any discrepencies in the parameterization of the distribution between JAGS and R - the most obvious of this is the use of precision in JAGS and standard deviation in base R (`rnorm`).

#list packages that can distirbutions of interest (cauchy, negative binomial, etc.)


MCMCtrace <- function(object,
                    params = 'all',
                    excl = NULL,
                    ISB = TRUE,
                    iter = 5000,
                    priors = NULL,
                    pdf = TRUE,
                    filename,
                    wd = getwd(),
                    type = 'both',
                    ind = FALSE)
{

  .pardefault <- graphics::par(no.readonly = T)

  #SORTING BLOCK
  if(typeof(object) == 'double')
  {
    stop('Invalid input type. Object must be of type stanfit, mcmc.list, or R2jags.')
  }else{
    object2 <- MCMCchains(object, params, excl, ISB, mcmc.list = TRUE)
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
  n_chains <- length(object2)
  colors <- gg_color_hue(n_chains)
  gg_cols <- grDevices::col2rgb(colors)/255

  #see how many iter is in output and if it is at least specified iter (5k by default)
  if (nrow(object2[[1]]) > iter)
  {
    it <- (nrow(object2[[1]]) - iter+1) : nrow(object2[[1]])
  }else {
    it <- 1 : nrow(object2[[1]])
  }

  #number of parameters
  np <- colnames(object2[[1]])

  if (!is.null(priors))
  {
    if (NCOL(priors) == 1 & length(np) > 1)
    {
      warning('Using a single prior for all parameters.')
    }
    if ((NCOL(priors) > 1 & NCOL(priors) != length(np)))
    {
      stop('Number of priors does not equal number of specified parameters.')
    }
  }


  if (type == 'both')
  {
    for (j in 1:length(np))
    {
      #j <- 1
      #trace
      tmlt <- do.call('cbind', object2[it, np[j]]) #make into matrix with three chains in columns
      graphics::matplot(it, tmlt, lwd = 1, lty= 1, type='l', main = paste0('Trace - ', np[j]),
              col= grDevices::rgb(red= gg_cols[1,], green= gg_cols[2,],
                       blue= gg_cols[3,], alpha = A_VAL),
              xlab= 'Iteration', ylab= 'Value')

      #density
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
             lty = 1, lwd = 1, main = paste0('Density - ', np[j]),
             col = grDevices::rgb(red = gg_cols[1,1], green = gg_cols[2,1], blue = gg_cols[3,1]))

        for (l in 2:NCOL(tmlt))
        {
          graphics::lines(dens[[l]],
                col = grDevices::rgb(red = gg_cols[1,l], green = gg_cols[2,l],
                          blue = gg_cols[3,l]))
        }
      }else{
        #density plot
        graphics::plot(stats::density(rbind(tmlt)), xlab = 'Parameter estimate',
             lty = 1, lwd = 1, main = paste0('Density - ', np[j]))
      }

      if (!is.null(priors))
      {
        if (NCOL(priors) == 1 & length(np) > 1)
        {
          wp <- priors
        }else{
          wp <- priors[,j]
        }
        lwp <- length(wp)
        if (lwp > length(it)*n_chains)
        {
          warning(paste0('Number of samples in prior is greater than number of total or specified iterations (for all chains) for specified parameter. Only last ', length(it)*n_chains, ' iterations will be used.'))

          pit <- (lwp - (length(it)*n_chains)+1) : lwp
          wp2 <- wp[pit]
        }
        if (lwp < length(it)*n_chains)
        {
          warning(paste0('Number of samples in prior is less than number of total or specified iterations (for all chains) for specified parameter. Resampling from prior to generate ', length(it)*n_chains, ' total iterations.'))

          samps <- sample(wp, size = ((length(it)*n_chains)-lwp), replace = TRUE)
          wp2 <- c(wp, samps)
        }
        if (lwp == length(it)*n_chains)
        {
          wp2 <- wp
        }

        #calculate percent ovelap
        tmlt_1c <- matrix(tmlt, ncol = 1)
        pp <- list(wp2, tmlt_1c)
        ovrlap <- paste0(round((overlapping::overlap(pp)$OV[[1]])*100, digits = 1), '% overlap')

        #plot prior and overlap text
        dpr <- stats::density(wp2)
        graphics::lines(dpr, col = 'red')
        graphics::legend('topright', legend = ovrlap, bty = 'n', pch = NA, text.col = 'red')
      }
    }
  }

  if (type == 'trace')
  {
    for (j in 1: length(np))
    {
      #trace
      tmlt <- do.call('cbind', object2[it, np[j]])
      graphics::matplot(it, tmlt, lwd = 1, lty= 1, type='l', main = paste0('Trace - ', np[j]),
              col= grDevices::rgb(red= gg_cols[1,], green= gg_cols[2,],
                       blue= gg_cols[3,], alpha = A_VAL),
              xlab= 'Iteration', ylab= 'Value')
    }
  }

  if (type == 'density')
  {
    for (j in 1: length(np))
    {
      #trace
      tmlt <- do.call('cbind', object2[it, np[j]])

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
             lty = 1, lwd = 1, main = paste0('Density - ', np[j]),
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
             lty = 1, lwd = 1, main = paste0('Density - ', np[j]))
      }

      if (!is.null(priors))
      {
        if (NCOL(priors) == 1 & length(np) > 1)
        {
          wp <- priors
        }else{
          wp <- priors[,j]
        }
        lwp <- length(wp)
        if (lwp > length(it)*n_chains)
        {
          warning(paste0('Number of samples in prior is greater than number of total or specified iterations (for all chains) for specified parameter. Only last ', length(it)*n_chains, ' iterations will be used.'))

          pit <- (lwp - (length(it)*n_chains)+1) : lwp
          wp2 <- wp[pit]
        }
        if (lwp < length(it)*n_chains)
        {
          warning(paste0('Number of samples in prior is less than number of total or specified iterations (for all chains) for specified parameter. Resampling from prior to generate ', length(it)*n_chains, ' total iterations.'))

          samps <- sample(wp, size = ((length(it)*n_chains)-lwp), replace = TRUE)
          wp2 <- c(wp, samps)
        }
        if (lwp == length(it)*n_chains)
        {
          wp2 <- wp
        }

        #calculate percent ovelap
        tmlt_1c <- matrix(tmlt, ncol = 1)
        pp <- list(wp2, tmlt_1c)
        ovrlap <- paste0(round((overlapping::overlap(pp)$OV[[1]])*100, digits = 1), '% overlap')

        #plot prior and overlap text
        dpr <- stats::density(wp2)
        graphics::lines(dpr, col = 'red')
        graphics::legend('topright', legend = ovrlap, bty = 'n', pch = NA, text.col = 'red')
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
