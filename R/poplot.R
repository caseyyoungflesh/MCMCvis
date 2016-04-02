#' Plot posterior distributions from MCMC output
#'
#' Plot the posterior distributions from MCMC output for specific parameters of interest. All posterior
#' parameter estimates are plotted on the same plot using density strips, similar to a caterpillar plot.
#'
#'
#' @param object Object containing MCMC output. See DETAILS below.
#' @param params Character string (or vector of character strings) denoting parameters to be
#' plotted. Partial names may be used to plot all parameters containing that set of characters.
#'
#' Default \code{all} plots posteriors for all parameters. See VALUE below.
#'
#' Valid entries are \code{jags_object}, \code{mcmc_list}, and \code{chains}. See DETAILS below.
#' @param g_lines Numerical vector indicating where vertical reference lines should be created.
#'
#' Default is \code{g_lines = 0}.
#'
#' Argument \code{NULL} will plot no guidelines.
#'
#' @param quantiles Numerical vecor of length 2, indicating which quantiles to plot.
#'
#' Default plots 95\% credible intervals.
#' @param rank If \code{TRUE} posteriors will ranked in decreasing order (based on
#' specified measure of centrality) from top down.
#' @param centrality Indicates which measure of centrality to plot.
#'
#' Valid options are \code{mean} and \code{median}.
#'
#' @param xlim Numerical vector of length 2, indicating range of x-axis.
#' @param xlab Character string labeling x-axis.
#' @param ylab Character string (or vector of character strings if plotting > 1 parameter) labeling
#' y-axis.
#'
#' Specifying labels in the argument will use these to label axis.
#'
#' Default option will use parameter names from \code{object}.
#'
#' Option \code{NULL} will return plot with no labels on y-axis.
#' @param main Character string indicating title of plot.
#' @param dbar_height Height of density bar in plot.
#' @param dbar_tick_height Height of ticks on density bar in plot.
#' @param dbar_tick_width Width of ticks on density bar in plot.
#' @section Details:
#' \code{object} argument can be an \code{mcmc.list} object, an \code{R2jags} model object (output from the \code{R2jags}
#' package), or a matrix containing MCMC chains (each column representing MCMC output for a single parameter, rows
#' representing iterations in the chain).
#'
#' @section Notes:
#' Plot code uses \code{denstrip} package, as highlighted in Jackson (2008) - generalized from code
#' for Zipkin et al. 2014, figure 3.
#'
#' @return \code{poplot(params = 'all')} plots posteriors for all parameters.
#'
#' \code{poplot(params = c('beta[1]', 'beta[2]'))} plots posteriors for just parameters \code{beta[1]} and \code{beta[2]}.
#'
#' \code{poplot(params = 'beta')} plots posteriors for all parameters containing \code{beta} in their name.
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
#' poplot(data)
#'
#' @export
#' @import lattice

poplot <- function(object,
                   params= 'all',
                   g_lines = 0,
                   quantiles = c(0.025, 0.975),
                   rank = FALSE,
                   centrality = 'mean',
                   xlim,
                   xlab = 'Parameter probability values',
                   ylab,
                   main,
                   dbar_height = 0.2,
                   dbar_tick_height = 0.5,
                   dbar_tick_width = 3)
{

  # Plotting parameters -----------------------------------------------------

  WID <- dbar_height #height of bar
  H <- ((dbar_tick_height-dbar_height)/2) #height of mean and CI ticks
  W <- dbar_tick_width #thickness of mean tick
  W2 <- dbar_tick_width #thickness of CI tick
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
        if (centrality == 'median')
        {
          tsrt <- apply(chains, 2, median)
          idx <- order(tsrt, decreasing = TRUE)
        }
        if (centrality == 'mean')
        {
          tsrt <- apply(chains, 2, mean)
          idx <- order(tsrt, decreasing = FALSE)
        }
        if (centrality != 'median' & centrality != 'mean')
        {
          stop(paste0(centrality,' is not a valid entry for the argument "centrality"'))
        }
      }
      if (rank == FALSE)
      {
        idx <- X:1
      }

      if (missing(ylab))
      {
        labs <- colnames(data)[idx]
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
            labs <- sort(ylab, decreasing =FALSE)
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
        labs <- colnames(data)[idx]
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

    if (missing(xlim))
    {

      rpp <- lattice::bwplot(variable ~ value, data = mp,
               xlab = list(label = xlab,cex = 1.3),
               main = Tmain,
               panel = function(x, y)
               {
                 xlist <- split(mp$value, factor(mp$variable))
                 xlist <- split(x, factor(y))

                 if(!is.null(g_lines))
                 {
                  for (k in 1: length(g_lines))
                  {
                    panel.abline(v=g_lines[k], lty = "dotted", col = "black")
                  }
                 }

                 for (i in seq(along = xlist))
                 {
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
        rpp <- lattice::bwplot(variable~value,data=mp,
                  xlab=list(label= xlab, cex=1.3),
                  main = Tmain,
                  xlim= xlim,
                  panel = function(x, y)
                  {
                    xlist <- split(mp$value, factor(mp$variable))
                    xlist <- split(x, factor(y))

                    if(!is.null(g_lines))
                    {
                      for (k in 1: length(g_lines))
                      {
                        panel.abline(v=g_lines[k], lty = "dotted", col = "black")
                      }
                    }

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

  lattice::trellis.focus(highlight=FALSE)

  x_limits <- rpp$x.limits
  lattice::panel.lines(c(x_limits[1],x_limits[2]),c(0.4,0.4), col="black")
  lattice::panel.lines(c(x_limits[1],x_limits[2]),c(X+.6,X+.6), col="black")

  # Hashes - centrality and CI -----------------------------------------------------------



  if (centrality == 'median')
  {
    for (i in 1:X)
    {
      j <- idx[i]
      #TMP<- X-j+1
      TMP <- i

      lattice::panel.lines(c(median(chains[,i])), c(TMP-H, TMP+H), lwd= W, col= MN_col)
      lattice::panel.lines(c(rep(qdata[1,i],2)), c(TMP-H, TMP+H), col= CI_col, lwd= W2)
      lattice::panel.lines(c(rep(qdata[2,i],2)), c(TMP-H, TMP+H), col= CI_col, lwd= W2)
    }
  }
  if (centrality == 'mean')
  {
    for (i in 1:X)
    {
      j <- idx[i]
      #TMP <- X-i+1
      TMP <- i

      lattice::panel.lines(c(mean(chains[,j])), c(TMP-H, TMP+H), lwd= W, col= MN_col)
      lattice::panel.lines(c(rep(qdata[1,j],2)), c(TMP-H, TMP+H), col= CI_col, lwd= W2)
      lattice::panel.lines(c(rep(qdata[2,j],2)), c(TMP-H, TMP+H), col= CI_col, lwd= W2)
    }
  }
  if (centrality != 'median' & centrality != 'mean')
  {
    stop(paste0(centrality,' is not a valid entry for the argument "centrality"'))
  }


  lattice::trellis.unfocus()
}
