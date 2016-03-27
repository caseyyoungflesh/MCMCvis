#' Plots posterior distributions
#'
#' Function plots the posterior distributions obtained from MCMC
#' chains output from a Bayesian model.
#'
#' @param object JAGS model object or matrix of MCMC chains.
#' @param params Character string (or vector of character strings) denoting parameters to be
#' returned in summary output. Partial names may be used to return all parameters containing
#' that set of characters. Default \code{all} returns all parameters in summary output.
#' @param input Indicates the nature of the \code{object} argument. Valid entries are
#' \code{jags_object} (for JAGS model object) and \code{chains} (for matrix of MCMC
#' chains).
#' @param quantiles Numerical vecor of length 2, indicating which quantiles to plot. Default plots 95\%
#' credible intervals.
#' @param centrality Indicates which measure of centrality to plot. Valid options are \code{mean}
#' and \code{median}.
#' @param xlim Numerical vector of length 2, indicating range of x-axis.
#' @param xlab Character string labeling x-axis.
#' @param ylab Character string (or vector of character strings if plotting > 1 parameter) labeling
#' y-axis. Specifying labels in the argument will use these to label axis. Default option will use
#' parameter names from model (or chains) object. Option \code{NULL} will return plot with no labels
#' on y-axis.
#' @param main Character string indicating title of plot.
#' @param dbar_height Height of density bar in plot.
#' @param tick_height Height of ticks in plot.
#' @param tick_width Width of ticks in plot.
#' @section Details:
#' For \code{input = 'chains'}, each column of \code{object} should contain MCMC output for
#' a single parameter. Each row represents one iteration in the chain.
#'
#' If JAGS model object is used, output must be from the \code{R2jags} package.
#'
#' @section Notes:
#' Plot code uses \code{denstrip} package, as highlighted in Jackson (2008) - generalized from code
#' for Zipkin et al. 2014, figure 3.
#'
#' @section References:
#' Jackson, C. H. 2008. Displaying Uncertainty With Shading. The American Statistician 62:340-347.
#'
#' Zipkin, E. F., T. S. Sillett, E. H. C. Grant, R. B. Chandler,
#' and J. A. Royle. 2014. Inferences about population dynamics
#' from count data using multistate models: a comparison to
#' capture-recapture approaches. Ecology and Evolution 4:417-426.
#'
#' @examples
#' x1 <- rnorm(1000, mean=0.5)
#' x2 <- rnorm(1000, mean=0)
#' data <- cbind(x1, x2)
#' poplot(data, input = 'chains')
#'
#' @export

#TO DO
#figure out how to make it not plot twice




poplot <- function(object,
                   params= 'all',
                   input = 'jags_object',
                   quantiles = c(0.025, 0.975),
                   centrality = 'mean',
                   xlim,
                   xlab = 'Parameter probability values',
                   ylab,
                   main,
                   dbar_height = 0.2,
                   tick_height = 0.5,
                   tick_width = 3)
{

  # Plotting parameters -----------------------------------------------------

  WID <- dbar_height #height of bar
  H <- ((tick_height-dbar_height)/2) #height of mean and CI ticks
  W <- tick_width #thickness of mean tick
  W2 <- tick_width #thickness of CI tick
  MN_col <- 'black' #color of centrality tick
  CI_col <- 'grey87' #color of CI tick

  if(missing(main))
  {
    Tmain <- ''
  }else
  {
    Tmain <- main
  }

  # Input data --------------------------------------------------------------

  if(input == 'jags_object')
  {
    if(typeof(object) == 'list')
    {
      data <- pochains(object, params= params)
    }else
    {
      stop('object type and input do not match')
    }
  }
  if(input == 'chains')
  {
   if(typeof(object) == 'double')
   {
     names <- colnames(object)

    if (length(params) == 1)
    {
      if (params == 'all')
      {
        data <- object
      }else
      {
        get.cols <- grep(paste(params), names, fixed=TRUE)
        data <- object[,get.cols]
      }
    }else
    {
      grouped <- c()
      for (i in 1:length(params))
      {
        get.cols <- grep(paste(params[i]), names, fixed=TRUE)
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

      data <- object[,g_filt]
    }
   }else
   {
     stop('object type and input do not match')
   }
  }
  if(input != 'jags_object' & input != 'chains')
  {
    stop(paste0(input,' is not a valid entry for the argument "input"'))
  }

  # Process data ------------------------------------------------------------

  if (NCOL(data) > 1)
  {
    if (length(quantiles)==2 & typeof(quantiles) == 'double')
    {
      chains <- as.data.frame(data)
      X <- NCOL(data)
      idx <- X:1

      if (missing(ylab))
      {
        labs <- colnames(data)[idx] #apply and sort labels
      }
      if (!missing(ylab))
      {
        if (is.null(ylab))
        {
          labs <- rep('', X)
        }
        if (!is.null(ylab))
        {
          if (length(ylab) == X)
          {
            labs <- ylab
          }else
          {
            stop('ylab length not equal to number of parameters')
          }
        }
      }

      mp <- suppressMessages(reshape2::melt(chains[,idx], value.name='value')) #melt
      qdata <- apply(chains, 2, quantile, probs=quantiles)
    }else
    {
      stop('quantiles must be a numerical vector of length 2')
    }
  }

  if (NCOL(data) == 1)
  {
    if (length(quantiles)==2 & typeof(quantiles) == 'double')
    {
      chains <- as.data.frame(data)
      X <- NCOL(data)
      idx <- X:1

      if (missing(ylab))
      {
        labs <- colnames(data)[idx] #apply and sort labels
      }
      if (!missing(ylab))
      {
        if (is.null(ylab))
        {
          labs <- rep('', X)
        }
        if (!is.null(ylab))
        {
          if (length(ylab) == X)
          {
            labs <- ylab
          }else
          {
            stop('ylab length not equal to number of parameters')
          }
        }
      }

      n_mp <- suppressMessages(reshape2::melt(chains[,idx], value.name='value')) #melt
      mp <- data.frame(variable = rep(paste0(params), nrow(n_mp)), value = n_mp)
      qdata <- apply(chains, 2, quantile, probs=quantiles)
    }else
    {
      stop('quantiles must be a numerical vector of length 2')
    }
  }

  # create plot object ------------------------------------------------------
?bwplot
    if (missing(xlim))
    {
      rpp <- bwplot(variable ~ value, data = mp,
               xlab = list(label = xlab,cex = 1.3),
               main = Tmain,
               panel = function(x, y)
               {
                 #grid.segments(1,0,0,0)
                 xlist <- split(mp$value, factor(mp$variable))
                 xlist <- split(x, factor(y))

                 for (i in seq(along = xlist))
                 {
                   #panel.grid(h=c(0), col='grey')
                   denstrip::panel.denstrip(x = xlist[[i]], at = i,
                                            width = WID, colmax = 'black', colmin = 'white')
                 }
               }, par.settings = list(axis.line = list(col = NA)),
               scales=list(col = 1,cex = 1, x = list(col = 1),
                           y = list(draw = T,labels = labs)))
    }else
    {
      if (length(xlim)==2 & typeof(xlim) == 'double')
      {
        rpp <- bwplot(variable~value,data=mp,
                  xlab=list(label= xlab, cex=1.3),
                  main = Tmain,
                  xlim= xlim,
                  panel = function(x, y)
                  {
                    #grid.segments(1,0,0,0)
                    xlist <- split(mp$value, factor(mp$variable))
                    xlist <- split(x, factor(y))

                    for (i in seq(along=xlist))
                    {
                      #panel.grid(h=c(0), col='grey')
                      denstrip::panel.denstrip(x=xlist[[i]], at=i,
                                               width= WID, colmax='black', colmin= 'white')
                    }
                  }, par.settings = list(axis.line = list(col=NA)),
                  scales=list(col=1,cex=1,x=list(col=1),
                              y=list(draw=T,labels=labs)))
      }else
      {
        stop('"xlim" must be a numerical vector of length 2')
      }
    }

  print(rpp)

  # Top and bottom lines ----------------------------------------------------

  lattice::trellis.focus()

  x_limits <- rpp$x.limits
  lattice::panel.lines(c(x_limits[1],x_limits[2]),c(0.4,0.4), col="black")
  lattice::panel.lines(c(x_limits[1],x_limits[2]),c(X+.6,X+.6), col="black")

  # Hashes - centrality -----------------------------------------------------------

  if (centrality == 'median')
  {
    for (i in 1:X)
    {
      TMP<- X-i+1
      lattice::panel.lines(c(median(chains[,i])), c(TMP-H, TMP+H), lwd= W, col= MN_col)
    }
  }
  if (centrality == 'mean')
  {
    for (i in 1:X)
    {
      TMP<- X-i+1
      lattice::panel.lines(c(mean(chains[,i])), c(TMP-H, TMP+H), lwd= W, col= MN_col)
    }
  }
  if (centrality != 'median' & centrality != 'mean')
  {
    stop(paste0(centrality,' is not a valid entry for the argument "centrality"'))
  }

  # Lower CI ----------------------------------------------------------------

  for(i in 1:X)
  {
    TMP<- X-i+1
    lattice::panel.lines(c(rep(qdata[1,i],2)), c(TMP-H, TMP+H), col= CI_col, lwd= W2)
  }

  # Upper CI ----------------------------------------------------------------

  for(i in 1:X)
  {
    TMP<- X-i+1
    lattice::panel.lines(c(rep(qdata[2,i],2)), c(TMP-H, TMP+H), col= CI_col, lwd= W2)
  }


  lattice::trellis.unfocus()
}
