#' Plot posterior distributions from MCMC output - VERSION 2
#'
#' Plot the posterior distributions from MCMC output for specific parameters of interest. Similar to a caterpillar plot.
#'
#'
#' @param object Object containing MCMC output. See DETAILS below.
#' @param params Character string (or vector of character strings) denoting parameters to be
#' plotted. Partial names may be used to plot all parameters containing that set of characters.
#'
#' Default \code{all} plots posteriors for all parameters. See VALUE below.
#'
#' Valid entries are \code{jags_object}, \code{mcmc_list}, and \code{chains}. See DETAILS below.
#' @param g_line Numerical vector indicating where vertical reference lines should be created.
#'
#' Default is \code{g_line = 0}.
#'
#' Argument \code{NULL} will plot no guidelines.
#'
#' @param g_line_width Width of vertical reference lines.
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
#' @param colors Vector of colors, indicating which color should should be used for each parameter. Default is 'black' for
#' all parameters.
#' @param dbar_height Height of density bar in plot.
#' @param dbar_t_height Height of ticks on density bar in plot.
#' @param dbar_t_width Width of ticks on density bar in plot.
#' @param CI_t_color Color of credible interval ticks.
#'
#' @section Details:
#' \code{object} argument can be an \code{mcmc.list} object, an \code{R2jags} model object (output from the \code{R2jags}
#' package), or a matrix containing MCMC chains (each column representing MCMC output for a single parameter, rows
#' representing iterations in the chain).
#'
#' @section Notes:
#'
#' When specifying \code{rank = TRUE} and specifying labels for \code{ylab}, labels will be applied to parameters before
#' they are ranked.
#'
#' Plot code uses \code{denstrip} package, as highlighted in Jackson (2008) - generalized from code
#' for Zipkin et al. 2014, figure 3.
#'
#' @return Function returns density strip plot, similar to caterpillar plot, for all specified parameters. Plotted output
#' can be sorted by mean estimate.
#'
#'
#' @examples
#' #Load data
#' data(MCMC_data)
#'
#' #Plot MCMC output
#' poplot(MCMC_data, ylab=NULL)
#'
#' #Just 'beta' parameters
#' poplot(MCMC_data, params= 'beta')
#'
#' #Just 'beta[1]', 'gamma[4]', and 'alpha[3]'
#' poplot(MCMC_data, params= c('beta[1]', 'gamma[4]', 'alpha[3]'))
#'
#' #Rank parameters by posterior mean
#' poplot(MCMC_data, params= 'beta', rank=TRUE)
#'
#' @export

data(MCMC_data)




plot(MCMC_data)

object <- MCMC_data
params <- 'all'
g_line = 0
g_line_width = 1
thin = 95
thick = 50
rank = FALSE
xlim
xlab = 'Parameter probability values'
ylab
main
colors
dbar_height = 0.25
dbar_t_height = 0.35
dbar_t_width = 3
CI_t_color = 'grey75'



function(object)
{




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
    if (length(thin) == 1 & typeof(thin) == 'double' & length(thick) == 1 & typeof(thick) == 'double')
    {
      chains <- as.data.frame(data)
      X <- NCOL(data)

      if (rank == TRUE)
      {
          tsrt <- apply(chains, 2, median) #used to rank positions of parameter estimates
          idx <- order(tsrt, decreasing = TRUE)
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
            labs <- ylab[idx]
          }else
          {
            stop('ylab length not equal to number of parameters')
          }
        }
      }


      #Not sure what to do with mp object here
      #
      #
      #thick_q (and thin) stand in for qdata object

      mp <- suppressMessages(reshape2::melt(chains[,idx], value.name='value')) #melt

      thick_ci <- c((100-((100-thick)/2)), ((100-thick)/2))*0.01
      thin_ci <- c((100-((100-thin)/2)), ((100-thin)/2))*0.01

      thick_q <- apply(chains, 2, quantile, probs= thick_ci)
      thin_q <- apply(chains, 2, quantile, probs= thin_ci)

      medians <- apply(chains, 2, quantile, probs = 0.5)

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
            labs <- ylab[idx]
          }else
          {
            stop('ylab length not equal to number of parameters')
          }
        }
      }



      #NOT sure about n_mp/mp objects here
      #
      #
      #


      n_mp <- suppressMessages(reshape2::melt(chains[,idx], value.name='value')) #melt
      mp <- data.frame(variable = rep(paste0(params), nrow(n_mp)), value = n_mp)



      thick_ci <- c((100-((100-thick)/2)), ((100-thick)/2))*0.01
      thin_ci <- c((100-((100-thin)/2)), ((100-thin)/2))*0.01

      thick_q <- apply(chains, 2, quantile, probs= thick_ci)
      thin_q <- apply(chains, 2, quantile, probs= thin_ci)

      medians <- apply(chains, 2, quantile, probs = 0.5)


    }else
    {
      stop('quantiles must be a numerical vector of length 2')
    }
  }
}



# base --------------------------------------------------------------------

med_sz <- 1.3 #size of median circles
thick_sz <- 4 #thick CI width
thin_sz <- 2 #thin CI width
gr_col <- 'gray' #which color gray is used for CI and medians
horizontal = TRUE


len <- length(medians)

black_cl <- c() #95% CI (default) does not overlap 0
gray_cl <- c() #50% CI (default) does not overlap 0
white_cl <- c() #Both 50% and 95% CI (default) overlap 0
for (i in 1:len)
{
  #i <- 1
  if ((thin_q[1,i] > 0 & thin_q[2,i] > 0) |
      (thin_q[1,i] < 0 & thin_q[2,i] < 0))
  {
    black_cl <- c(black_cl, i)
  } else {
  if ((thick_q[1,i] > 0 & thick_q[2,i] > 0) |
      (thick_q[1,i] < 0 & thick_q[2,i] < 0))
  {
    gray_cl <- c(gray_cl, i)
  }else {
    white_cl <- c(white_cl, i)
  }
  }
}


bnd_q <- rbind(1:len, 1:len) # to bind with quantiles to plot them

if (horizontal)
{
  plot(medians, 1:len, xlim = range(thin_q))
  #abline(h=1:len, col = gray(0.8)) #horizontal lines at parameters
  abline(v=0, lty = 3, lwd = 3) #0 line

  #CI - thick
  matlines(thick_q[,black_cl], bnd_q[,black_cl],
           type = 'l', lty = 1, lwd = thick_sz, col = 'black') #black
  matlines(thick_q[,gray_cl], bnd_q[,gray_cl],
           type = 'l', lty = 1, lwd = thick_sz, col = gr_col) #gray
  matlines(thick_q[,white_cl], bnd_q[,white_cl],
           type = 'l', lty = 1, lwd = thick_sz, col = gr_col) #white (gray)

  #CI - thin
  matlines(thin_q[,black_cl], bnd_q[,black_cl],
           type = 'l', lty = 1, lwd = thin_sz, col = 'black') #black
  matlines(thin_q[,gray_cl], bnd_q[,gray_cl],
           type = 'l', lty = 1, lwd = thin_sz, col = gr_col) #gray
  matlines(thin_q[,white_cl], bnd_q[,white_cl],
           type = 'l', lty = 1, lwd = thin_sz, col = gr_col) #white (gray)


  #Medians
  points(medians, 1:len, pch = 16, col = 'white', cex = med_sz) #plot points over other plot features
  points(medians[black_cl], black_cl, pch = 16, col = 'black', cex = med_sz) #95% CI doesn't overlap 0
  points(medians[gray_cl], gray_cl, pch = 16, col = gr_col, cex = med_sz) #50% CI doesn't overlap 0
  points(medians[white_cl], white_cl, pch = 21, col = gr_col, cex = med_sz) #Both CI overlap 0
}

grid()


