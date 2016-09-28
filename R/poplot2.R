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
#' Thanks to Cinner et al. 2016, whose Fig. 1 inspired this plot.
#'
#' @return Function returns density strip plot, similar to caterpillar plot, for all specified parameters. Plotted output
#' can be sorted by mean estimate.
#'
#' @section References:
#'
#' Cinner, J.E., .... (2016) Bright spots among the world's coral reefs. Nature.
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


require(potools)



object <- MCMC_data
params <- 'gamma'
thin = 95
thick = 50
rank = FALSE


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

 # if (NCOL(data) > 1)
  #{
    #if (length(thin) == 1 & typeof(thin) == 'double' & length(thick) == 1 & typeof(thick) == 'double')
    #{
      chains <- as.data.frame(data)
      len <- NCOL(data)

      if (rank == TRUE)
      {
          tsrt <- apply(chains, 2, median) #used to rank positions of parameter estimates
          idx <- order(tsrt, decreasing = TRUE)
      }
      if (rank == FALSE)
      {
        idx <- 1:len
      }

      #if (missing(ylab))
      #{
      #  labs <- colnames(data)[idx]
      #}
      #if (!missing(ylab))
      #{
      #  if (is.null(ylab))
      #  {
      #    labs <- rep('', X)
      #  }
      #  if (!is.null(ylab))
      #  {
      #    if (length(ylab) == X)
      #    {
      #      labs <- ylab[idx]
      #    }else
      #    {
      #      stop('ylab length not equal to number of parameters')
      #    }
      #  }
      #}


      thick_ci <- c((100-((100-thick)/2)), ((100-thick)/2))*0.01
      thin_ci <- c((100-((100-thin)/2)), ((100-thin)/2))*0.01

      thick_q <- apply(chains, 2, quantile, probs= thick_ci)[,idx]
      thin_q <- apply(chains, 2, quantile, probs= thin_ci)[,idx]

      medians <- apply(chains, 2, quantile, probs = 0.5)[idx]

    #}else
    #{
    #  stop('quantiles must be a numerical vector of length 2')
    #}
  #}

  #if (NCOL(data) == 1)
  #{
  #  if (length(quantiles)==2 & typeof(quantiles) == 'double')
  #  {
  #    chains <- as.data.frame(data)
  #    len <- NCOL(data)
  #    idx <- len:1
#
 #     if (missing(ylab))
  #    {
   #     labs <- colnames(data)[idx]
    #  }
     # if (!missing(ylab))
     # {
      #   if (is.null(ylab))
      #  {
      #    labs <- rep('', len)
      #  }
      #  if (!is.null(ylab))
      #  {
      #    if (length(ylab) == len)
      #    {
      #      labs <- ylab[idx]
      #    }else
      #    {
      #      stop('ylab length not equal to number of parameters')
      #    }
      #  }
      #}


      #thick_ci <- c((100-((100-thick)/2)), ((100-thick)/2))*0.01
      #thin_ci <- c((100-((100-thin)/2)), ((100-thin)/2))*0.01

      #thick_q <- apply(chains, 2, quantile, probs= thick_ci)[,idx]
      #thin_q <- apply(chains, 2, quantile, probs= thin_ci)[,idx]

      #medians <- apply(chains, 2, quantile, probs = 0.5)[idx]


      #}else
      #{
      #stop('quantiles must be a numerical vector of length 2')
      #}
      #}
#}




# base --------------------------------------------------------------------

#plotting parameters
sp = 2 #spacing
med_sz = 2.4 #size of median circles
thick_sz = 5 #thick CI width
thin_sz = 2 #thin CI width
ax_th = 3#axis and tick thickness
tick_text_sz = 1.2 #tick label size
axis_text_sz = 1.3 #axis label size


gr_col = 'gray60' #color used for CI and medians
zero_col = 'gray60' #color used for 0 line
horizontal = TRUE


xlab = 'Parameter Estimate' #should be changed to: if (missing(xlab)){xlab <- 'Parameter Estimate'}
ylab = NULL #should be changed to: if (missing(ylab)){ylab <- NULL}
main = '' #should be changed to : if (missing(ylab)){ylab <- ''}
labels = names(medians) #y-axis labels - should be changed to : if (missing(labels)){labels <- names(medians)}

tick_pos = NULL #where ticks should be placed - should be changed to: if (missing(tick_pos)){tick_pos <- NULL}
xlim = range(thin_q)*1.1 #should be changed to: if (missing(xlim)){xlim <- range(thin_q)*1.25}
ylim = c(0.5,(len) + 0.5) #should be changed to: if (missing(ylim)){ylim <- c(0.5,(len)+0.5)}
mar = c(5,4,4,2) #should be changed to: if (missing(mar)){mar <- c(5,4,4,2)}



#Determine which params have CI that overlap 0
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




#positions bound together to plot CI
blk_bnd <- rbind(black_cl, black_cl)
gry_bnd <- rbind(gray_cl, gray_cl)
wht_bnd <- rbind(white_cl, white_cl)




#plot for horizontal
if (horizontal)
{

  m_char <- max(sapply(labels, nchar))
  #variable at LEFT position to account for differing label sizes - can be altered manually
  par(mar=c(mar[1], (1.5 + (m_char/2)) + (4 - mar[2]), mar[3], mar[4]) + 0.1)



  #plot blank plot
  plot(medians, 1:len, xlim = xlim, ylim = ylim, type = "n",
       ann = TRUE, xaxt = 'n', yaxt = "n", bty = "n", ylab = NA,
       xlab = xlab, main = main,
       #lab #number of ticks to plot on each axis
       cex.lab = axis_text_sz) #cex.axis is tick labels, lab is axis label



  #bottom x-axis line
  #abline(h = 0.1, lwd = ax_th)
  #top x-axis line
  #abline(h = len + 0.9, lwd =ax_th)
  #bottom axis params
  axis(3, lwd.tick = ax_th, labels = FALSE,
       at = tick_pos, lwd = ax_th)
  #bottom axis params
  axis(1, lwd.tick = ax_th, labels = TRUE,
       at = tick_pos, lwd = ax_th,
       cex.axis = tick_text_sz) #bottom axis
  #left axis params (labels)
  axis(2, at = (1:len)+0.1, tick = FALSE,
       labels = labels, las = 1, adj = 0, #las - 0 parallel to axis, 1 horiz, 2 perp to axis, 3 vert
       line = -1.5, cex.axis = tick_text_sz)


  #lheight par for making space!


  #zero line
  abline(v=0, lty = 2, lwd = 3, col = zero_col)

  #Black CI
  if (!is.null(black_cl))
  {
    #Thick
    matlines(thick_q[,black_cl], blk_bnd,
             type = 'l', lty = 1, lwd = thick_sz, col = 'black')
    #Thin
    matlines(thin_q[,black_cl], blk_bnd,
             type = 'l', lty = 1, lwd = thin_sz, col = 'black')
  }

  #Gray CI
  if (!is.null(gray_cl))
  {
    matlines(thick_q[,gray_cl], gry_bnd,
             type = 'l', lty = 1, lwd = thick_sz, col = gr_col) #gray
    matlines(thin_q[,gray_cl], gry_bnd,
             type = 'l', lty = 1, lwd = thin_sz, col = gr_col) #gray
  }

  #White CI
  if (!is.null(white_cl))
  {
    matlines(thick_q[,white_cl], wht_bnd,
             type = 'l', lty = 1, lwd = thick_sz, col = gr_col) #white (gray)
    matlines(thin_q[,white_cl], wht_bnd,
             type = 'l', lty = 1, lwd = thin_sz, col = gr_col) #white (gray)
  }


  #Medians
  points(medians, 1:len, pch = 16, col = 'white', cex = med_sz*.9) #plot points over other plot features
  points(medians[black_cl], black_cl, pch = 16, col = 'black', cex = med_sz) #95% CI doesn't overlap 0
  points(medians[gray_cl], gray_cl, pch = 16, col = gr_col, cex = med_sz) #50% CI doesn't overlap 0
  points(medians[white_cl], white_cl, pch = 21, col = gr_col, cex = med_sz) #Both CI overlap 0
}



#par(mar=c(5,4,4,2) + 0.1) #only needed if the ylab option is enabled

#end of function
}


#make all arguments same as base graphics
#feature to add space between lines
#clean up into function
#look at integrating other 'plot' features (look at caterplot)
#rename function(s) - perhaps name poplot, other can be dstplot
#add vertical argument

#add to potrace

#change names of functions

#look at adding stan object compatibility
