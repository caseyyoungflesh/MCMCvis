#' Density strip plots of posterior distributions from MCMC output
#'
#' Visualize posterior distributions from MCMC output for specific parameters of interest using density
#' strips. Medians and credible intervals are indicated by ticks on density strips.
#'
#'
#' @param object Object containing MCMC output. See DETAILS below.
#' @param params Character string (or vector of character strings) denoting parameters to be
#' plotted. Partial names may be used to plot all parameters containing that set of characters.
#'
#' Default \code{'all'} plots posteriors for all parameters. See VALUE below.
#'
#' Valid entries are \code{jags_object}, \code{mcmc_list}, and \code{chains}. See DETAILS below.
#'
#'
#' @param ref_line Numerical vector indicating where vertical reference line should be created.
#'
#' Default is \code{ref_line = 0}.
#'
#' Argument \code{NULL} will plot no reference line.
#'
#'
#' @param ref_line_width Width of vertical reference lines.
#'
#'
#' @param quantiles Numerical vecor of length 2, indicating which quantiles to plot.
#'
#' Default plots 95\% credible intervals.
#' @param rank If \code{TRUE} posteriors will ranked in decreasing order (based on
#' median posterior estimate) from top down.
#'
#' @param xlim Numerical vector of length 2, indicating range of x-axis.
#' @param xlab Character string labeling x-axis.
#' @param main Character string indicating title of plot.
#' @param labels Character string (or vector of character strings if plotting > 1 parameter) labeling
#' parameter estimates along y-axis.
#'
#' Specifying labels in the argument will use these to label density strips on y-axis.
#'
#' Default option will use parameter names from \code{object}.
#'
#' Option \code{NULL} will return plot with no labels on y-axis.
#' @param colors Vector of colors, indicating which color should should be used for each parameter. Default is 'black' for
#' all parameters.
#' @param dbar_height Height of density bar in plot.
#' @param dbar_t_height Height of ticks on density bar in plot.
#' @param dbar_t_width Width of ticks on density bar in plot.
#' @param CI_t_color Color of credible interval ticks.
#' @param ax_sz Width of x-axis and ticks.
#'
#' @section Details:
#' \code{object} argument can be an \code{mcmc.list} object, an \code{R2jags} model object (output from the \code{R2jags}
#' package), or a matrix containing MCMC chains (each column representing MCMC output for a single parameter, rows
#' representing iterations in the chain).
#'
#' @section Notes:
#' Plot code uses \code{denstrip} package, as highlighted in Jackson (2008) - generalized from code
#' for Zipkin et al. 2014, figure 3.
#'
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
#' #Load data
#' data(MCMC_data)
#'
#' #Plot MCMC output
#' denplot(MCMC_data, labels=NULL)
#'
#' #Just 'beta' parameters
#' denplot(MCMC_data, params= 'beta')
#'
#' #Just 'beta[1]', 'gamma[4]', and 'alpha[3]'
#'denplot(MCMC_data, params= c('beta[1]', 'gamma[4]', 'alpha[3]'))
#'
#' #Rank parameters by posterior median
#'denplot(MCMC_data, params= 'beta', rank=TRUE)
#'
#' @export
#' @import lattice

denplot <- function(object,
                   params= 'all',
                   ref_line = 0,
                   ref_line_width = 1,
                   quantiles = c(0.025, 0.975),
                   rank = FALSE,
                   xlim,
                   xlab,
                   main,
                   labels,
                   colors,
                   dbar_height = 0.25,
                   dbar_t_height = 0.35,
                   dbar_t_width = 3,
                   ax_sz = 2,
                   CI_t_color = 'grey75')
{

  # Input data --------------------------------------------------------------

  data <- pochains(object, params= params)

  if(coda::is.mcmc.list(object) != TRUE &
     typeof(object) != 'double' &
     typeof(object) != 'list')
  {
    stop('Invalid object type. Input must be mcmc.list object, rjags object, or matrix with MCMC chains.')
  }


  # Process data ------------------------------------------------------------

  if (NCOL(data) > 1)
  {
    if (length(quantiles) == 2 & typeof(quantiles) == 'double')
    {
      chains <- as.data.frame(data)
      X <- NCOL(data)

      if (rank == TRUE)
      {
          tsrt <- apply(chains, 2, median)
          idx <- order(tsrt, decreasing = TRUE)
      }
      if (rank == FALSE)
      {
        idx <- X:1
      }

      if (missing(labels))
      {
        labs <- colnames(data)[idx]
      }
      if (!missing(labels))
      {
        if (is.null(labels))
        {
          labs <- rep('', X)
        }
        if (!is.null(labels))
        {
          if (length(labels) == X)
          {
            labs <- labels[idx]
          }else
          {
            stop('labels length not equal to number of parameters')
          }
        }
      }

      mp <- suppressMessages(reshape2::melt(chains[,idx], value.name='value')) #melt
      qdata <- apply(chains, 2, quantile, probs=quantiles)
      medians <- apply(chains, 2, quantile, probs=0.5)
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

      if (missing(labels))
      {
        labs <- colnames(data)[idx]
      }
      if (!missing(labels))
      {
        if (is.null(labels))
        {
          labs <- rep('', X)
        }
        if (!is.null(labels))
        {
          if (length(labels) == X)
          {
            labs <- labels[idx]
          }else
          {
            stop('labels length not equal to number of parameters')
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

  # Plotting parameters -----------------------------------------------------

  WID <- dbar_height #height of bar
  H <- (dbar_t_height/2) #height of median and CI ticks
  W <- dbar_t_width #thickness of median tick
  W2 <- dbar_t_width #thickness of CI tick
  MN_col <- 'black' #color of centrality tick
  CI_col <- CI_t_color #color of CI tick
  GCOL <- rgb(0,0,0,alpha = .5)
  VTHICK <- ref_line_width

  if (missing(xlab))
  {xlab = 'Parameter Estimate'}

  if(missing(main))
  {
    Tmain <- ''
  }else
  {
    Tmain <- main
  }

  if (missing(colors))
  {
    COLORS <- rep('black', X)
  }else
  {
    if (is.null(colors))
    {
      COLORS <- rep('black', X)
    }else
    {
      if(length(colors) == X)
      {
        COLORS <- colors[idx]
      }else
        stop('length(colors) does not equal number of parameters to be plotted.')
    }
  }


  # create plot object ------------------------------------------------------

    if (missing(xlim))
    {

      rpp <- lattice::bwplot(variable ~ value, data = mp,
               xlab = list(label = xlab,cex = 1.3),
               main = Tmain,
               panel = function(x, y)
               {
                 xlist <- split(mp$value, factor(mp$variable))
                 xlist <- split(x, factor(y))

                 for (i in seq(along = xlist))
                 {
                   denstrip::panel.denstrip(x = xlist[[i]], at = i,
                                            width = WID, colmax = COLORS[i], colmin = 'white')
                 }

                 if(!is.null(ref_line))
                 {
                   for (k in 1: length(ref_line))
                   {
                     panel.abline(v=ref_line[k], lty = 2, col = GCOL, lwd = VTHICK)
                   }
                 }

               }, par.settings = list(axis.line = list(col = NA),
                                      layout.heights = list(xlab.axis.padding = 0.5)),
               scales=list(col = 1,cex = 1, x = list(col = 1),
                           y = list(draw = T,labels = labs), lwd=ax_sz))
    }else
    {
      if (length(xlim)==2 & typeof(xlim) == 'double')
      {
        rpp <- lattice::bwplot(variable~value,data=mp,
                  xlab=list(label= xlab, cex=1.3),
                  main = Tmain,
                  xlim= xlim,
                  panel = function(x, y)
                  {
                    xlist <- split(mp$value, factor(mp$variable))
                    xlist <- split(x, factor(y))

                    for (i in seq(along=xlist))
                    {
                      #panel.grid(h=c(0), col='grey')
                      denstrip::panel.denstrip(x=xlist[[i]], at=i,
                                               width= WID, colmax= COLORS[i], colmin= 'white')
                    }

                    if(!is.null(ref_line))
                    {
                      for (k in 1: length(ref_line))
                      {
                        panel.abline(v=ref_line[k], lty = 2, col = GCOL, lwd = VTHICK)
                      }
                    }

                  }, par.settings = list(axis.line = list(col=NA)),
                  scales=list(col=1,cex=1,x=list(col=1),
                              y=list(draw=T,labels=labs), lwd=ax_sz))
      }else
      {
        stop('"xlim" must be a numerical vector of length 2')
      }
    }


  print(rpp)

  # Top and bottom lines ----------------------------------------------------

  lattice::trellis.focus(highlight=FALSE)

  x_limits <- rpp$x.limits
  lattice::panel.lines(c(x_limits[1],x_limits[2]),c(0.4,0.4), col="black", lwd = ax_sz+1)
  lattice::panel.lines(c(x_limits[1],x_limits[2]),c(X+0.6,X+0.6), col="black", lwd = ax_sz+1)




  # Hashes - centrality and CI -----------------------------------------------------------


    for (i in 1:X)
    {
      j <- idx[i]
      #TMP <- X-j+1
      TMP <- i

      lattice::panel.lines(c(median(chains[,j])), c(TMP-H, TMP+H), lwd= W, col= MN_col)
      lattice::panel.lines(c(rep(qdata[1,j],2)), c(TMP-H, TMP+H), col= CI_col, lwd= W2)
      lattice::panel.lines(c(rep(qdata[2,j],2)), c(TMP-H, TMP+H), col= CI_col, lwd= W2)
    }


  lattice::trellis.unfocus()
}
