#' Trace and density plots from MCMC output
#'
#' Trace and density plots of MCMC chains for specific parameters of interest. Print plots to pdf by default.
#'
#' @param object Object containing MCMC output. See DETAILS below.
#'
#' @param params Character string (or vector of character strings) denoting parameters of interest.
#'
#' Default \code{'all'} returns chains for all parameters.
#'
#' @param excl Character string (or vector of character strings) denoting parameters to exclude. Used in conjunction with \code{params} argument to select parameters of interest.
#'
#' @param ISB Ignore Square Brackets (ISB). Logical specifying whether square brackets should be ignored in the \code{params} and \code{excl} arguments. If \code{TRUE}, square brackets are ignored - input from \code{params} and \code{excl} are otherwise matched exactly. If \code{FALSE}, square brackets are not ignored - input from \code{params} and \code{excl} are matched using grep, which can take arguments in regular expression format. This allows partial names to be used when specifying parameters of interest.
#'
#' @param iter Number of iterations to plot for trace and density plots. The default value is 5000, meaning the last 5000 iterations of the chain will be plotted.
#'
#' @param gvals Vector containing generating values if simulated data was used to fit model. These values will be plotted as vertical lines on the density plots to compare posterior distributions with the true parameter values used to generate the data. No line will be apparent if the generating value is outside the plotted range of the posterior distribution.
#'
#' @param priors Matrix containing random draws from prior distributions corresponding to parameters of interest. If specified, priors are plotted along with posterior density plots. Percent overlap between prior and posterior (PPO) is also calculated and displayed on each plot. Each column of the matrix represents a prior for a different parameter. Parameters are plotted alphabetically - priors should be sorted accordingly. If \code{priors} contains only one prior and more than one parameter is specified for the \code{params} argument, this prior will be used for all parameters. The number of draws for each prior should equal the number of iterations specified by \code{iter} (or total draws if less than \code{iter}) times the number of chains, though the function will automatically adjust if more or fewer iterations are specified. See DETAILS below.
#'
#' @param post_zm Logical - if \code{post_zm = FALSE} x- and y-limits of density plots are scaled so that both the prior and posterior can be visualized on a single density plot (rather than zoomed on the posterior).
#'
#' @param PPO_out Logical - if \code{PPO_out = TRUE} percent overlap between prior and posterior (PPO) will be output to a dataframe.
#'
#' @param Rhat Logical - if \code{Rhat = TRUE} potential scale reduction factor (Rhat) for each parameter is plotted on the trace plots.
#'
#' @param n.eff Logical - if \code{n.eff = TRUE} number of effective samples for each parameter is plotted on the trace plots.
#' 
#' @param ind Logical - if \code{ind = TRUE}, separate density lines will be plotted for each chain. If \code{ind= FALSE}, one density line will be plotted for all chains.
#' 
#' @param pdf Logical - if \code{pdf = TRUE} plots will be exported to a pdf.
#'
#' @param plot Logical - if \code{plot = FALSE} no plot will be output. Designed to be used in conjunction with \code{PPO_out = TRUE}, to calculate PPO without displaying plot output.
#'
#' @param open_pdf Logical - if \code{open_pdf = TRUE} pdf will open in viewer after being generated.
#'
#' @param filename Name of pdf file to be printed. Default is 'MCMCtrace'.
#'
#' @param wd Working directory for pdf output. Default is current directory.
#'
#' @param type Type of plot to be output. \code{'both'} outputs both trace and density plots, \code{'trace'}
#' outputs only trace plots, and \code{'density'} outputs only density plots.
#'
#' @param xlim Vector of length two specifying limits for x-axis on density plots only. If specified, overrides argument \code{post_zm}.
#'
#' @param ylim Vector of length two specifying limits for y-axis on density plots only. If specified, overrides argument \code{post_zm}.
#'
#' @param xlab_den Character string specifying label for x-axis on density plots.
#'
#' @param ylab_den Character string specifying label for x-axis on density plots.
#'
#' @param xlab_tr Character string specifying label for x-axis on trace plots.
#'
#' @param ylab_tr Character string specifying label for x-axis on trace plots.
#'
#' @param main_den Character string (or vector of character strings if plotting > 1 parameter) specifying title(s) of density plot(s).
#'
#' @param main_tr Character string (or vector of character strings if plotting > 1 parameter) specifying title(s) of trace plot(s).
#'
#' @param lwd_den Number specifying thickness of density line on density plots.
#'
#' @param lwd_pr Number specifying thickness of prior line on density plots.
#'
#' @param lty_den Number specifying the line type for the density line on density plots.
#'
#' @param lty_pr Number specifying the line type for the prior line on density plots.
#'
#' @param col_den Character string specifying color of density line on density plots. Does not specify color if \code{ind = TRUE}.
#'
#' @param col_pr Character string specifying color of prior line on density plots.
#'
#' @param col_txt Character string specifying color of text (denoting PPO) on plot when value specified for \code{priors}. If \code{NULL} is specified, no text will be plot.
#' 
#' @param sz_txt Number specifying size of text (denoting PPO) when value specified for \code{priors}. If \code{NULL} is specified, no text will be plot.
#'
#' @param sz_ax Number specifying thickness of axes and ticks.
#'
#' @param sz_ax_txt Number specifying size of text for axes labels.
#'
#' @param sz_tick_txt Number specifying size of text for tick labels on axis.
#'
#' @param sz_main_txt Number specifying size of text for main title.
#' 
#' @param pos_tick_x_tr Numeric vector specifying where ticks on x-axis should be placed for trace plots.
#'
#' @param pos_tick_y_tr Numeric vector specifying where ticks on y-axis should be placed for trace plots.
#'
#' @param pos_tick_x_den Numeric vector specifying where ticks on x-axis should be placed for density plots.
#'
#' @param pos_tick_y_den Numeric vector specifying where ticks on y-axis should be placed for density plots.
#'
#' @section Details:
#' \code{object} argument can be a \code{stanfit} object (\code{rstan} package), a \code{stanreg} object (\code{rstanarm} package), a \code{brmsfit} object (\code{brms} package), an \code{mcmc.list} object (\code{coda} package), an \code{R2jags} model object (\code{R2jags} package), a \code{jagsUI} model object (\code{jagsUI} package), or a matrix containing MCMC chains (each column representing MCMC output for a single parameter, rows representing iterations in the chain). The function automatically detects the object type and proceeds accordingly.
#'
#' Matrices for the \code{priors} argument can be generated using commands such as rnorm, rgamma, runif, etc. Distributions not supported by base R can be generated by using the appropriate packages. It is important to note that some discrepancies between MCMC samplers and R may exist regarding the parameterization of distributions - one example of this is the use of precision in JAGS but standard deviation in R for the 'second parameter' of the normal distribution. If the number of draws for each prior distribution is greater than the total number used for the density plot (\code{iter} times the number of chains), the function will use a subset of the prior draws. If the number of draws for each prior distribution is less than the total number used for the density plot, the function will resample (with replacement) from the prior to obtain the appropriate number of draws.
#'
#' @examples
#' #Load data
#' data(MCMC_data)
#'
#' #Traceplots for all 'beta' parameters - pdf is generated by default
#' MCMCtrace(MCMC_data, params = 'beta', pdf = FALSE)
#'
#' #Traceplots (individual density lines for each chain) just for 'beta[1]'
#' #'params' takes regular expressions when ISB = FALSE, square brackets 
#' #must be escaped with '\\'
#' MCMCtrace(MCMC_data, params = 'beta\\[1\\]', 
#'           ISB = FALSE, ind = TRUE, pdf = FALSE)
#'
#' #Plot prior on top of posterior, calculate prior/posterior overlap (PPO) 
#' #just for 'beta[1]'
#' #Add Rhat and n.eff values to density plots
#' #'params' takes regular expressions when ISB = FALSE, square brackets must 
#' #be escaped with '\\'
#' PR <- rnorm(15000, 0, 32)
#' MCMCtrace(MCMC_data, params = 'beta\\[1\\]', ISB = FALSE, priors = PR, 
#'           pdf = FALSE, Rhat = TRUE, n.eff = TRUE)
#'
#' #Output PPO to R object without plotting trace plots
#' PR <- rnorm(15000, 0, 32)
#' PPO <- MCMCtrace(MCMC_data, params = 'beta\\[1\\]', ISB = FALSE, 
#'                  priors = PR, plot = FALSE, PPO_out = TRUE)
#' @export


MCMCtrace <- function(object,
                      params = 'all',
                      excl = NULL,
                      ISB = TRUE,
                      iter = 5000,
                      gvals = NULL,
                      priors = NULL,
                      post_zm = TRUE,
                      PPO_out = FALSE,
                      Rhat = FALSE,
                      n.eff = FALSE,
                      ind = FALSE,
                      pdf = TRUE,
                      plot = TRUE,
                      open_pdf = TRUE,
                      filename,
                      wd = getwd(),
                      type = 'both',
                      ylim = NULL,
                      xlim = NULL,
                      xlab_tr,
                      ylab_tr,
                      xlab_den,
                      ylab_den,
                      main_den = NULL,
                      main_tr = NULL,
                      lwd_den = 1,
                      lwd_pr = 1,
                      lty_den = 1,
                      lty_pr = 1,
                      col_den,
                      col_pr,
                      col_txt,
                      sz_txt = 1.2,
                      sz_ax = 1, 
                      sz_ax_txt = 1, 
                      sz_tick_txt = 1, 
                      sz_main_txt = 1.2, 
                      pos_tick_x_tr = NULL,
                      pos_tick_y_tr = NULL,
                      pos_tick_x_den = NULL,
                      pos_tick_y_den = NULL)
{
  .pardefault <- graphics::par(no.readonly = T)
  
  #SORTING BLOCK
  if (methods::is(object, 'matrix'))
  {
    warning('Input type matrix - assuming only one chain for each parameter.')
    object1 <- coda::as.mcmc.list(coda::as.mcmc(object))
    object2 <- MCMCchains(object1, params, excl, ISB, mcmc.list = TRUE)
  } else {
    object2 <- MCMCchains(object, params, excl, ISB, mcmc.list = TRUE)
  }
  
  
  #PREP BLOCK
  #parameter names
  np <- colnames(object2[[1]])
  
  #number of chains
  n_chains <- length(object2)
  
  #see how many iter is in output and if it is at least specified iter (5k by default)
  if (nrow(object2[[1]]) > iter)
  {
    it <- (nrow(object2[[1]]) - iter+1) : nrow(object2[[1]])
  } else {
    it <- 1 : nrow(object2[[1]])
  }
  
  
  #WARNING/ERROR BLOCK
  if (!is.null(priors))
  {
    if (NCOL(priors) == 1 & length(np) > 1)
    {
      warning('Only one prior specified for > 1 parameter. Using a single prior for all parameters.')
    }
    if ((NCOL(priors) > 1 & NCOL(priors) != length(np)))
    {
      stop('Number of priors does not equal number of specified parameters.')
    }
    if (NROW(priors) > length(it)*n_chains)
    {
      warning(paste0('Number of samples in prior is greater than number of total or specified iterations (for all chains) for specified parameter. Only last ', length(it)*n_chains, ' iterations will be used.'))
    }
    if (NROW(priors) < length(it)*n_chains)
    {
      warning(paste0('Number of samples in prior is less than number of total or specified iterations (for all chains) for specified parameter. Resampling from prior to generate ', length(it)*n_chains, ' total iterations.'))
    }
    if (type == 'trace')
    {
      warning("Prior posterior overlap (PPO) cannot be plotting without density plots. Use type = 'both' or type = 'density'.")
    }
  }
  
  if (!is.null(gvals))
  {
    if (length(gvals) == 1 & length(np) > 1)
    {
      warning('Only one generating value specified for > 1 parameter. Using a single generating value for all parameters.')
    }
    if (length(gvals) > 1 & length(gvals) != length(np))
    {
      stop('Number of generating values does not equal number of specified parameters.')
    }
  }
  
  if (PPO_out == TRUE)
  {
    PPO_df <- data.frame(param = rep(NA, length(np)), percent_PPO = rep(NA, length(np)))
  }
  
  
  
  #PLOT BLOCK
  if (plot == TRUE)
  {
    if (pdf == TRUE)
    {
      if (missing(filename))
      {
        file_out <- paste0(wd, '/MCMCtrace.pdf')
      } else {
        if (grepl('.pdf', filename, fixed = TRUE))
        {
          file_out <- paste0(wd, '/', filename)
        } else {
          file_out <- paste0(wd, '/', filename, '.pdf')
        }
      }
      pdf(file = file_out)
    }
    
    #plotting parameters
    ref_col <- 'red'
    A_VAL <- 0.5 #alpha value
    
    #adjust layout based on # params - differs based on den, tr, both
    if (type == 'both')
    {
      if (length(np) >= 3)
      {
        graphics::layout(matrix(c(1, 2, 3, 4, 5, 6), 3, 2, byrow = TRUE))
        graphics::par(mar = c(4.1,4.1,2.1,1.1)) # bottom, left, top, right
        MN_LINE <- NULL
      }
      if (length(np) == 2)
      {
        graphics::layout(matrix(c(1, 2, 3, 4, 5, 6), 2, 2, byrow = TRUE))
        graphics::par(mar = c(4.1,4.1,2.1,1.1)) # bottom, left, top, right
        MN_LINE <- NULL
      }
      if (length(np) == 1)
      {
        graphics::layout(matrix(c(1, 2, 3, 4, 5, 6), 1, 2, byrow = TRUE))
        graphics::par(mar = c(8.1,4.1,7.1,1.1)) # bottom, left, top, right
        MN_LINE <- 1.1   
      }
    } else {
      if (length(np) >= 5)
      {
        graphics::layout(matrix(c(1, 2, 3, 4, 5, 6), 3, 2, byrow = TRUE))
        graphics::par(mar = c(4.1,4.1,2.1,1.1)) # bottom, left, top, right
        MN_LINE <- NULL
      }
      if (length(np) == 3 | length(np) == 4)
      {
        graphics::layout(matrix(c(1, 2, 3, 4, 5, 6), 2, 2, byrow = TRUE))
        graphics::par(mar = c(4.1,4.1,2.1,1.1)) # bottom, left, top, right
        MN_LINE <- NULL
      }
      if (length(np) == 2)
      {
        graphics::layout(matrix(c(1, 2, 3, 4, 5, 6), 1, 2, byrow = TRUE))
        graphics::par(mar = c(8.1,4.1,7.1,1.1)) # bottom, left, top, right
        MN_LINE <- 1.1
      }
      if (length(np) == 1)
      {
        graphics::layout(matrix(c(1, 2, 3, 4, 5, 6), 1, 1, byrow = TRUE))
        graphics::par(mar = c(5.1, 4.1, 4.1, 2.1)) # bottom, left, top, right
        MN_LINE <- NULL
      }
    }
    
    graphics::par(mgp = c(2.5,1,0)) #axis text distance
    gg_color_hue <- function(n)
    {
      hues = seq(15, 375, length = n+1)
      grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
    }
    
    colors <- gg_color_hue(n_chains)
    gg_cols <- grDevices::col2rgb(colors)/255
    
    
    #from input - for lim
    YLIM <- ylim
    XLIM <- xlim
    
    #from input - for lab
    if (!missing(xlab_tr))
    {
      xlab_tr <- xlab_tr
    } else {
      xlab_tr <- 'Iteration'
    }
    
    if (!missing(ylab_tr))
    {
      ylab_tr <- ylab_tr
    } else {
      ylab_tr <- 'Value'
    }
    
    if (!missing(xlab_den))
    {
      xlab_den <- xlab_den
    } else {
      xlab_den <- 'Parameter estimate'
    }
    
    if (!missing(ylab_den))
    {
      ylab_den <- ylab_den
    } else {
      ylab_den <- 'Density'
    }
    
    #from input - for main - if length one and given, print - if length > 1, print that element 
    if (!missing(main_den))
    {
      if (length(main_den) != 1 & length(main_den) != length(np))
      {
        stop("Number of elements for 'main_den' does not equal number of specified parameters.")
      }
    }
    if (!missing(main_tr))
    {
      if (length(main_tr) !=1 & length(main_tr) != length(np))
      {
        stop("Number of elements for 'main_tr' does not equal number of specified parameters.")
      }
    }
    
    MAIN_DEN <- function(x, md = main_den, idx)
    {
      if (!is.null(md))
      {
        if (length(md) == 1)
        {
          return(paste0(md))
        } else {
          return(paste0(md[idx]))
        }
      } else {
        return(paste0('Density - ', x))
      }
    }
    
    MAIN_TR <- function(x, mtr = main_tr, idx)
    {
      if (!is.null(mtr))
      {
        if (length(mtr) == 1)
        {
          return(paste0(mtr))
        } else {
          return(paste0(mtr[idx]))
        }
      } else {
        return(paste0('Trace - ', x))
      }
    }
    
    #color for density and prior lines
    if (missing(col_den))
    {
      COL_DEN <- 'black'
    } else {
      COL_DEN <- col_den
    }
    
    if (missing(col_pr))
    {
      COL_PR <- 'red'
    } else {
      COL_PR <- col_pr
    }
    
    if (missing(col_txt))
    {
      COL_TXT <- 'red'
    } else {
      COL_TXT <- col_txt
    }
    
    #tick position
    if (is.null(pos_tick_x_den))
    {
      XAXT_DEN <- 's'
    } else {
      XAXT_DEN <- 'n'
    }
    if (is.null(pos_tick_y_den))
    {
      YAXT_DEN <- 's'
    } else {
      YAXT_DEN <- 'n'
    }
    if (is.null(pos_tick_x_tr))
    {
      XAXT_TR <- 's'
    } else {
      XAXT_TR <- 'n'
    }
    if (is.null(pos_tick_y_tr))
    {
      YAXT_TR <- 's'
    } else {
      YAXT_TR <- 'n'
    }
    
    if (Rhat == TRUE | n.eff == TRUE)
    {
      summ <- MCMCsummary(object, params = params, excl = excl, ISB = ISB, 
                          Rhat = Rhat, n.eff = n.eff)
      if (Rhat == TRUE)
      {
        rhat <- summ[,grep('Rhat', colnames(summ))]
      }
      if (n.eff == TRUE)
      {
        neff <- summ[,grep('n.eff', colnames(summ))]
      }
    }
    
    if (type == 'both')
    {
      for (j in 1:length(np))
      {
        #j <- 1
        #trace
        tmlt <- do.call('cbind', object2[it, np[j]]) #make into matrix with three chains in columns
        graphics::matplot(it, tmlt, lwd = 1, lty = 1, type = 'l', main = NULL,
                          col = grDevices::rgb(red = gg_cols[1,], green = gg_cols[2,],
                                              blue = gg_cols[3,], alpha = A_VAL),
                          xlab = xlab_tr, ylab = ylab_tr, cex.axis = sz_tick_txt, 
                          cex.lab = sz_ax_txt, xaxt = XAXT_TR, yaxt = YAXT_TR)
        
        graphics::title(main = MAIN_TR(np[j], main_tr, j), line = MN_LINE, cex.main = sz_main_txt)
        
        graphics::axis(side = 1, at = pos_tick_x_tr, cex.axis = sz_tick_txt)
        graphics::axis(side = 2, at = pos_tick_y_tr, cex.axis = sz_tick_txt)
        
        #if sz_ax is specified
        if (!missing(sz_ax))
        {
          graphics::box(lwd = sz_ax)
          graphics::axis(1, lwd.ticks = sz_ax, labels = FALSE, at = pos_tick_x_tr)
          graphics::axis(2, lwd.ticks = sz_ax, labels = FALSE, at = pos_tick_y_tr)
        }
        
        #PPO
        #if priors are specified
        if (!is.null(priors))
        {
          if (NCOL(priors) == 1)
          {
            wp <- priors
          } else {
            wp <- priors[,j]
          }
          lwp <- length(wp)
          if (lwp > length(it)*n_chains)
          {
            #warnings in block above
            pit <- (lwp - (length(it) * n_chains) + 1) : lwp
            wp2 <- wp[pit]
          }
          if (lwp < length(it)*n_chains)
          {
            #warnings in block above
            samps <- sample(wp, size = ((length(it) * n_chains) - lwp),
                            replace = TRUE)
            wp2 <- c(wp, samps)
          }
          if (lwp == length(it) * n_chains)
          {
            wp2 <- wp
          }
          
          #for PPO plotting below
          dpr <- stats::density(wp2)
          PPO_x_rng <- range(dpr$x)
          PPO_y_rng <- range(dpr$y)
          
          #calculate percent ovelap
          tmlt_1c <- matrix(tmlt, ncol = 1)
          pp <- list(wp2, tmlt_1c)
          ovr_v <- round((overlapping::overlap(pp)$OV[[1]]) * 100, digits = 1)
          ovrlap <- paste0(ovr_v, '% overlap')
          
          if (PPO_out == TRUE)
          {
            PPO_df$param[j] <- np[j]
            PPO_df$percent_PPO[j] <- ovr_v
          }
        }
        
        #density
        if (ind == TRUE & n_chains > 1)
        {
          dens <- apply(tmlt, 2, stats::density)
          max_den_y <- c()
          rng_den_x <- c()
          for (k in 1:NCOL(tmlt))
          {
            max_den_y <- c(max_den_y, max(dens[[k]]$y))
            rng_den_x <- c(rng_den_x, range(dens[[k]]$x))
          }
          
          #set axes limits according to inputs
          if (!is.null(priors) & post_zm == FALSE & is.null(ylim) & is.null(xlim))
          {
            ylim <- range(c(0, max(max_den_y), PPO_y_rng))
            xlim <- range(c(range(rng_den_x), PPO_x_rng))
          } else {
            if (!is.null(ylim))
            {
              ylim <- YLIM
              xlim <- XLIM
            }
            if (is.null(ylim) & is.null(xlim))
            {
              ylim <- c(0, max(max_den_y))
              xlim <- NULL
            }
          }
          
          graphics::plot(dens[[1]], xlab = xlab_den, ylab = ylab_den, ylim = ylim, xlim = xlim,
                         lty = lty_den, lwd = lwd_den, main = '',
                         col = grDevices::rgb(red = gg_cols[1,1], green = gg_cols[2,1], 
                                              blue = gg_cols[3,1]),
                         cex.axis = sz_tick_txt, cex.lab = sz_ax_txt, xaxt = XAXT_DEN, 
                         yaxt = YAXT_DEN)
          
          graphics::title(main = MAIN_DEN(np[j], main_den, j), line = MN_LINE, 
                          cex.main = sz_main_txt)
          
          graphics::axis(side = 1, at = pos_tick_x_den, cex.axis = sz_tick_txt)
          graphics::axis(side = 2, at = pos_tick_y_den, cex.axis = sz_tick_txt)
          
          for (l in 2:NCOL(tmlt))
          {
            graphics::lines(dens[[l]], lty = lty_den, lwd = lwd_den,
                            col = grDevices::rgb(red = gg_cols[1,l], green = gg_cols[2,l],
                                                 blue = gg_cols[3,l]))
          }
          
          #if sz_ax is specified
          if (!missing(sz_ax))
          {
            graphics::box(lwd = sz_ax)
            graphics::axis(1, lwd.ticks = sz_ax, labels = FALSE, at = pos_tick_x_den)
            graphics::axis(2, lwd.ticks = sz_ax, labels = FALSE, at = pos_tick_y_den)
          }
          
        } else {
          dens <- stats::density(rbind(tmlt))
          rng_den_x <- range(dens$x)
          #set axes limits according to inputs
          if (!is.null(priors) & post_zm == FALSE & is.null(ylim) & is.null(xlim))
          {
            ylim <- range(c(range(dens$y), PPO_y_rng))
            xlim <- range(c(range(dens$x), PPO_x_rng))
          } else {
            if (!is.null(ylim))
            {
              ylim <- YLIM
              xlim <- XLIM
            }
            if (is.null(ylim) & is.null(xlim))
            {
              ylim <- NULL
              xlim <- NULL
            }
          }
          
          #density plot
          graphics::plot(dens, xlab = xlab_den, ylab = ylab_den, ylim = ylim, 
                         main = '', col = COL_DEN, xlim = xlim, lty = lty_den, 
                         lwd = lwd_den, cex.axis = sz_tick_txt,
                         cex.lab = sz_ax_txt, xaxt = XAXT_DEN, yaxt = YAXT_DEN)
          
          graphics::title(main = MAIN_DEN(np[j], main_den, j), line = MN_LINE, 
                          cex.main = sz_main_txt)
          
          graphics::axis(side = 1, at = pos_tick_x_den, cex.axis = sz_tick_txt)
          graphics::axis(side = 2, at = pos_tick_y_den, cex.axis = sz_tick_txt)
          
          #if sz_ax is specified
          if (!missing(sz_ax))
          {
            graphics::box(lwd = sz_ax)
            graphics::axis(1, lwd.ticks = sz_ax, labels = FALSE, at = pos_tick_x_den)
            graphics::axis(2, lwd.ticks = sz_ax, labels = FALSE, at = pos_tick_y_den)
          }
        }
        
        #plotting PPO
        if (!is.null(priors))
        {
          #plot prior and overlap text
          graphics::lines(dpr, col = COL_PR, lwd = lwd_pr, lty = lty_pr)
          
          #don't plot text if NULL specified for SZ or COL
          if (!is.null(sz_txt) & !is.null(COL_TXT))
          {
            graphics::legend('topright', legend = ovrlap, bty = 'n', pch = NA, 
                             text.col = COL_TXT, cex = sz_txt)
          }
        }
        
        #diagnostics plotted on density plots
        if (Rhat == TRUE & n.eff == TRUE)
        {
          diag_txt <- list(paste0('Rhat: ', rhat[j]),
                           paste0('n.eff: ', neff[j]))
        } 
        if (Rhat == TRUE & n.eff == FALSE)
        {
          diag_txt <- paste0('Rhat: ', rhat[j])
        }
        if (Rhat == FALSE & n.eff == TRUE)
        {
          diag_txt <- paste0('n.eff: ', neff[j])
        }
        
        #don't plot text if NULL specified for SZ or COL
        if (Rhat == TRUE | n.eff == TRUE)
        {
          if (!is.null(sz_txt) & !is.null(COL_TXT))
          {
            #inset needs to be smaller when sz_txt is larger
            graphics::legend('topleft', x.intersp = -0.5,  
                             legend = diag_txt, 
                             bty = 'n', pch = NA, text.col = COL_TXT, cex = sz_txt)
          }
        }
        
        #if generating values are specified - warnings in block above
        if (!is.null(gvals))
        {
          if (length(gvals) == 1)
          {
            gv <- gvals
          } else {
            gv <- gvals[j]
          }
          graphics::abline(v = gv, lty = 2, lwd = 3, col = ref_col)
        }
      }
    }
    
    if (type == 'trace')
    {
      for (j in 1: length(np))
      {
        #trace
        tmlt <- do.call('cbind', object2[it, np[j]])
        graphics::matplot(it, tmlt, lwd = 1, lty = 1, type='l', main = NULL,
                          col = grDevices::rgb(red = gg_cols[1,], green = gg_cols[2,],
                                              blue = gg_cols[3,], alpha = A_VAL),
                          xlab = xlab_tr, ylab = ylab_tr, cex.axis = sz_tick_txt, 
                          cex.lab = sz_ax_txt, xaxt = XAXT_TR, yaxt = YAXT_TR)
        
        graphics::title(main = MAIN_TR(np[j], main_tr, j), line = MN_LINE, 
                        cex.main = sz_main_txt)
        
        graphics::axis(side = 1, at = pos_tick_x_tr, cex.axis = sz_tick_txt)
        graphics::axis(side = 2, at = pos_tick_y_tr, cex.axis = sz_tick_txt)
        
        #if sz_ax is specified
        if (!missing(sz_ax))
        {
          graphics::box(lwd = sz_ax)
          graphics::axis(1, lwd.ticks = sz_ax, labels = FALSE, at = pos_tick_x_tr)
          graphics::axis(2, lwd.ticks = sz_ax, labels = FALSE, at = pos_tick_y_tr)
        }
      }
    }
    
    if (type == 'density')
    {
      for (j in 1: length(np))
      {
        tmlt <- do.call('cbind', object2[it, np[j]])
        
        #PPO
        #if priors are specified
        if (!is.null(priors))
        {
          if (NCOL(priors) == 1)
          {
            wp <- priors
          } else {
            wp <- priors[,j]
          }
          lwp <- length(wp)
          if (lwp > length(it)*n_chains)
          {
            #warnings in block above
            pit <- (lwp - (length(it)*n_chains)+1) : lwp
            wp2 <- wp[pit]
          }
          if (lwp < length(it)*n_chains)
          {
            #warnings in block above
            samps <- sample(wp, size = ((length(it)*n_chains)-lwp), replace = TRUE)
            wp2 <- c(wp, samps)
          }
          if (lwp == length(it)*n_chains)
          {
            wp2 <- wp
          }
          
          #for PPO plotting below
          dpr <- stats::density(wp2)
          PPO_x_rng <- range(dpr$x)
          PPO_y_rng <- range(dpr$y)
          
          #calculate percent ovelap
          tmlt_1c <- matrix(tmlt, ncol = 1)
          pp <- list(wp2, tmlt_1c)
          ovr_v <- round((overlapping::overlap(pp)$OV[[1]])*100, digits = 1)
          ovrlap <- paste0(ovr_v, '% overlap')
          
          if (PPO_out == TRUE)
          {
            PPO_df$param[j] <- np[j]
            PPO_df$percent_PPO[j] <- ovr_v
          }
        }
        
        
        if (ind == TRUE & n_chains > 1)
        {
          dens <- apply(tmlt, 2, stats::density)
          max_den_y <- c()
          rng_den_x <- c()
          for (k in 1:NCOL(tmlt))
          {
            max_den_y <- c(max_den_y, max(dens[[k]]$y))
            rng_den_x <- c(rng_den_x, range(dens[[k]]$x))
          }
          
          #set axes limits according to inputs
          if (!is.null(priors) & post_zm == FALSE & is.null(ylim) & is.null(xlim))
          {
            ylim <- range(c(0, max(max_den_y), PPO_y_rng))
            xlim <- range(c(range(rng_den_x), PPO_x_rng))
          } else {
            if (!is.null(ylim))
            {
              ylim <- YLIM
              xlim <- XLIM
            }
            if (is.null(ylim) & is.null(xlim))
            {
              ylim <- c(0, max(max_den_y))
              xlim <- NULL
            }
          }
          
          graphics::plot(dens[[1]], xlab = xlab_den, ylab = ylab_den, ylim = ylim, xlim = xlim,
                         lty = lty_den, lwd = lwd_den, main = '',
                         col = grDevices::rgb(red = gg_cols[1,1], green = gg_cols[2,1], 
                                              blue = gg_cols[3,1]),
                         cex.axis = sz_tick_txt, cex.lab = sz_ax_txt, xaxt = XAXT_DEN, 
                         yaxt = YAXT_DEN)
          
          graphics::title(main = MAIN_DEN(np[j], main_den, j), line = MN_LINE, 
                          cex.main = sz_main_txt)
          
          graphics::axis(side = 1, at = pos_tick_x_den, cex.axis = sz_tick_txt)
          graphics::axis(side = 2, at = pos_tick_y_den, cex.axis = sz_tick_txt)
          
          for (l in 2:NCOL(tmlt))
          {
            graphics::lines(dens[[l]], lty = lty_den, lwd = lwd_den,
                            col = grDevices::rgb(red = gg_cols[1,l], green = gg_cols[2,l],
                                                 blue = gg_cols[3,l]))
          }
          
          #if sz_ax is specified
          if (!missing(sz_ax))
          {
            graphics::box(lwd = sz_ax)
            graphics::axis(1, lwd.ticks = sz_ax, labels = FALSE, at = pos_tick_x_den)
            graphics::axis(2, lwd.ticks = sz_ax, labels = FALSE, at = pos_tick_y_den)
          }
        } else {
          
          dens <- stats::density(rbind(tmlt))
          #set axes limits according to inputs
          if (!is.null(priors) & post_zm == FALSE & is.null(ylim) & is.null(xlim))
          {
            ylim <- range(c(range(dens$y), PPO_y_rng))
            xlim <- range(c(range(dens$x), PPO_x_rng))
          } else {
            if (!is.null(ylim))
            {
              ylim <- YLIM
              xlim <- XLIM
            }
            if (is.null(ylim) & is.null(xlim))
            {
              ylim <- NULL
              xlim <- NULL
            }
          }
          
          #density plot
          graphics::plot(stats::density(rbind(tmlt)), xlab = xlab_den, ylab = ylab_den, ylim = ylim,
                         col = COL_DEN, xlim = xlim, lty = lty_den, lwd = lwd_den, main = '',
                         cex.axis = sz_tick_txt, cex.lab = sz_ax_txt, xaxt = XAXT_DEN, 
                         yaxt = YAXT_DEN)
          
          graphics::title(main = MAIN_DEN(np[j], main_den, j), line = MN_LINE, 
                          cex.main = sz_main_txt)
          
          graphics::axis(side = 1, at = pos_tick_x_den, cex.axis = sz_tick_txt)
          graphics::axis(side = 2, at = pos_tick_y_den, cex.axis = sz_tick_txt)
          
          #if sz_ax is specified
          if (!missing(sz_ax))
          {
            graphics::box(lwd = sz_ax)
            graphics::axis(1, lwd.ticks = sz_ax, labels = FALSE, at = pos_tick_x_den)
            graphics::axis(2, lwd.ticks = sz_ax, labels = FALSE, at = pos_tick_y_den)
          }
        }
        
        #plotting PPO
        if (!is.null(priors))
        {
          #plot prior and overlap text
          graphics::lines(dpr, col = COL_PR, lwd = lwd_pr, lty = lty_pr)
          
          #don't plot text if NULL specified for SZ or COL
          if (!is.null(sz_txt) & !is.null(COL_TXT))
          {
            graphics::legend('topright', legend = ovrlap, bty = 'n', pch = NA, 
                             text.col = COL_TXT, cex = sz_txt)
          }
        }
        
        #diagnostics plotted on density plots
        if (Rhat == TRUE & n.eff == TRUE)
        {
          diag_txt <- list(paste0('Rhat: ', rhat[j]),
                           paste0('n.eff: ', neff[j]))
        } 
        if (Rhat == TRUE & n.eff == FALSE)
        {
          diag_txt <- paste0('Rhat: ', rhat[j])
        }
        if (Rhat == FALSE & n.eff == TRUE)
        {
          diag_txt <- paste0('n.eff: ', neff[j])
        }
        
        #don't plot text if NULL specified for SZ or COL
        if (Rhat == TRUE | n.eff == TRUE)
        {
          if (!is.null(sz_txt) & !is.null(COL_TXT))
          {
            #inset needs to be smaller when sz_txt is larger
            graphics::legend('topleft', x.intersp = -0.5, 
                             legend = diag_txt, 
                             bty = 'n', pch = NA, text.col = COL_TXT, cex = sz_txt)
          }
        }
        
        #if generating values are specified - warnings in block above
        if (!is.null(gvals))
        {
          if (length(gvals) == 1)
          {
            gv <- gvals
          } else {
            gv <- gvals[j]
          }
          graphics::abline(v = gv, lty = 2, lwd = 3, col = ref_col)
        }
      }
    }
    
    if (type != 'both' & type != 'density' & type != 'trace')
    {
      stop('Invalid argument for "type". Valid inputs are "both", "trace", and "density".')
    }
    
    if (pdf == TRUE)
    {
      invisible(grDevices::dev.off())
      if (open_pdf == TRUE)
      {
        system(paste0('open ', file_out))
      }
    } else {
      graphics::par(.pardefault)
    }
  }
  
  if (plot == FALSE)
  {
    for (j in 1: length(np))
    {
      tmlt <- do.call('cbind', object2[it, np[j]])
      
      #PPO
      #if priors are specified
      if (!is.null(priors))
      {
        if (NCOL(priors) == 1)
        {
          wp <- priors
        } else {
          wp <- priors[,j]
        }
        lwp <- length(wp)
        if (lwp > length(it)*n_chains)
        {
          #warnings in block above
          pit <- (lwp - (length(it)*n_chains)+1) : lwp
          wp2 <- wp[pit]
        }
        if (lwp < length(it)*n_chains)
        {
          #warnings in block above
          samps <- sample(wp, size = ((length(it)*n_chains)-lwp), replace = TRUE)
          wp2 <- c(wp, samps)
        }
        if (lwp == length(it)*n_chains)
        {
          wp2 <- wp
        }
        
        #for PPO plotting below
        dpr <- stats::density(wp2)
        PPO_x_rng <- range(dpr$x)
        PPO_y_rng <- range(dpr$y)
        
        #calculate percent ovelap
        tmlt_1c <- matrix(tmlt, ncol = 1)
        pp <- list(wp2, tmlt_1c)
        ovr_v <- round((overlapping::overlap(pp)$OV[[1]])*100, digits = 1)
        ovrlap <- paste0(ovr_v, '% overlap')
        
        if (PPO_out == TRUE)
        {
          PPO_df$param[j] <- np[j]
          PPO_df$percent_PPO[j] <- ovr_v
        }
      }
    }
  }
  
  if (PPO_out == TRUE)
  {
    if (is.null(priors))
    {
      warning('NAs produced for PPO dataframe as priors not specified.')
    }
    return(PPO_df)
  }
}
