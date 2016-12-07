#' Caterpillar plots of posterior distributions from MCMC output
#'
#' Visualize posterior distributions from MCMC output for specific parameters of interest using caterpillar plots. Color of median dot represents the overlap of the posterior distrubtion with 0 (or other specified value).
#'
#'
#' @param object Object containing MCMC output. See DETAILS below.
#' @param params Character string (or vector of character strings) denoting parameters to be plotted. Partial names may be used to plot all parameters containing that set of characters.
#'
#' Default \code{'all'} plots posteriors for all parameters. See VALUE below.
#'
#'
#' @param excl Character string (or vector of character strings) denoting parameters to exclude. Partial names may be used to exclude all parameters containing that set of characters. Used in conjunction with \code{params} argument to select parameters of interest.
#'
#' @param ref Value indicating where vertical reference line should be created and what value to use a reference for caterpillar median coloration.
#'
#' Default is \code{ref = 0}.
#'
#' Argument \code{NULL} will plot no reference line.
#'
#' @param ref_ovl Logical specifying whether the style/color of plotted median dots and CI should be changed based on whether the 50 \% and 95 \% credible intervals overlap the reference line. See DETAILS for more information.
#'
#' @param rank If \code{TRUE} posteriors will be ranked in decreasing order (based on specified measure of centrality) from top down.
#'
#' @param horiz If \code{TRUE} posteriors will be plotted running horizontally (parallel to the x-axis). If \code{FALSE} posteriors will be plotted running vertically (perpendicular to the x-axis).
#' @param xlim Numerical vector of length 2, indicating range of x-axis. Only applicable if \code{horiz = TRUE}.
#' @param ylim Numerical vector of length 2, indicating range of y-axis. Only applicable if \code{horiz = FALSE}.
#' @param xlab Character string labeling x-axis. Only applicable if \code{horiz = TRUE}.
#'
#' Default label is 'Parameter Estimate'. Option \code{NULL} will return plot with no label on x-axis.
#' @param ylab Character string labeling y-axis. Only applicable if \code{horiz = FALSE}.
#'
#' Default label is 'Parameter Estimate'. Option \code{NULL} will return plot with no label on y-axis.
#' @param main Character string indicating title of plot.
#'
#' @param labels Character string (or vector of character strings if plotting > 1 parameter) labeling parameter estimates along y-axis (if \code{horiz = FALSE}) or x-axis (if \code{horiz = TRUE}).
#'
#' Default option will use parameter names from \code{object}.
#'
#' Option \code{NULL} will return plot with no labels on axis.
#'
#' @param labels_sz Number specifying size of text for parameter labels on axis.
#'
#' @param med_sz Number specifying size of points represents posterior medians.
#'
#' @param thick_sz Number specifying thickness of 50 percent CI line (thicker line).
#'
#' @param thin_sz Number specifying thickness of 95 percent CI line (thinner line).
#'
#' @param ax_sz Number specifying thickness of axis and ticks.
#'
#' @param x_axis_text_sz Number specifying size of text for axis label.
#'
#' @param x_tick_text_sz Number specifying size of text for tick labels on axis.
#'
#' @param main_text_sz Number specifying size of text for main title.
#'
#' @param tick_pos Numeric vector specifying where ticks on axis should be placed.
#'
#' @param mar Numerical vector of length 4 specifying plot margins - (BOTTOM, LEFT, TOP, RIGHT). Changes to the margin should be made within the function rather than using the \code{par} call.
#'
#' Default is c(5.1, 4.1, 4.1, 2.1) - the R plot default.
#'
#' @section Details:
#' Points represent posterior medians. Parameters where 50\% credible intervals overlap 0 (or other specified value) are indicated by 'open' circles. Parameters where 50 percent credible intervals DO NOT overlap 0 AND 95 percent credible intervals DO overlap 0 (or other specified value) are indicated by 'closed' grey circles. Parameters where 95 percent credible intervals DO NOT overlap 0 (or other specified value) are indicated by 'closed' black circles. Thick lines represent 50 percent credible intervals while thin lines represent 95 \% credible intervals. \code{ref_ovl = FALSE} can be used to disable this feature. All median dots will be represented as 'closed' black circles. \code{object} argument can be a \code{stanfit} object (\code{rstan} package), an \code{mcmc.list} object (\code{coda} package), an \code{R2jags} model object (\code{R2jags} package), or a matrix containing MCMC chains (each column representing MCMC output for a single parameter, rows representing iterations in the chain). The function automatically detects the object type and proceeds accordingly.
#'
#' @section Notes:
#'
#' When specifying \code{rank = TRUE} and specifying labels for \code{labels}, labels will be applied to parameters before they are ranked.
#'
#' Thanks to Cinner et al. 2016, whose Fig. 1 inspired this plot.
#'
#'
#' @section References:
#'
#' Cinner, J. E., C. Huchery, M. A. MacNeil, N. A. J. Graham, T. R. McClanahan, J. Maina, E. Maire, J. N. Kittinger, C. C. Hicks, C. Mora, E. H. Allison, S. D'Agata, A. Hoey, D. A. Feary, L. Crowder, I. D. Williams, M. Kulbicki, L. Vigliola, L. Wantiez, G. Edgar, R. D. Stuart-Smith, S. A. Sandin, A. L. Green, M. J. Hardt, M. Beger, A. Friedlander, S. J. Campbell, K. E. Holmes, S. K. Wilson, E. Brokovich, A. J. Brooks, J. J. Cruz-Motta, D. J. Booth, P. Chabanet, C. Gough, M. Tupper, S. C. A. Ferse, U. R. Sumaila, and D. Mouillot. 2016. Bright spots among the world's coral reefs. Nature 535:416-419.
#'
#'
#' @examples
#' #Load data
#' data(MCMC_data)
#'
#' #Plot MCMC output
#' MCMCplot(MCMC_data, labels = NULL)
#'
#' #Just 'beta' parameters
#' MCMCplot(MCMC_data, params = 'beta')
#'
#' #Just 'beta[1]', 'gamma[4]', and 'alpha[3]'
#' MCMCplot(MCMC_data, params = c('beta[1]', 'gamma[4]', 'alpha[3]'))
#'
#' #Rank parameters by posterior mean
#' MCMCplot(MCMC_data, params = 'beta', rank = TRUE)
#'
#' #Create vertical plot
#' MCMCplot(MCMC_data, params = 'beta', horiz = FALSE)
#'
#' @export
#'

MCMCplot <- function(object,
                   params = 'all',
                   excl = NULL,
                   ref = 0,
                   ref_ovl = TRUE,
                   rank = FALSE,
                   horiz = TRUE,
                   xlim,
                   ylim,
                   xlab,
                   ylab,
                   main,
                   labels,
                   labels_sz = 1.2,
                   med_sz = 1.5,
                   thick_sz = 5,
                   thin_sz = 2,
                   ax_sz = 3,
                   x_axis_text_sz = 1.3,
                   x_tick_text_sz = 1.2,
                   main_text_sz = 1.2,
                   tick_pos,
                   mar = c(5.1, 4.1, 4.1, 2.1))
{


  data <- MCMCchains(object, params= params, excl = excl)

  # Plotting parameters -----------------------------------------------------

  if (missing(tick_pos))
  {tick_pos = NULL}

  gr_col = 'gray60' #color used for CI and medians
  ref_col = 'gray60' #color used for 0 line
  thin = 95 #CI for thin line
  thick = 50 #CI for thick line

  # Process data ------------------------------------------------------------

  if (NCOL(data) > 1)
  {
    if (length(thin) == 1 & typeof(thin) == 'double' & length(thick) == 1 & typeof(thick) == 'double')
    {
      chains <- as.data.frame(data)
      len <- NCOL(data)

      if (rank == TRUE)
      {
        tsrt <- apply(chains, 2, stats::median)
        idx <- order(tsrt, decreasing = TRUE)
      }
      if (rank == FALSE)
      {
        idx <- len:1
      }


      thick_ci <- c((100-((100-thick)/2)), ((100-thick)/2))*0.01
      thin_ci <- c((100-((100-thin)/2)), ((100-thin)/2))*0.01

      thick_q <- apply(chains, 2, stats::quantile, probs= thick_ci)[,idx]
      thin_q <- apply(chains, 2, stats::quantile, probs= thin_ci)[,idx]

      medians <- apply(chains, 2, stats::quantile, probs = 0.5)[idx]

    }else
    {
      stop("'thick' and 'thin' must be single numbers")
    }
  }

  if (NCOL(data) == 1)
  {
    if (length(thin) == 1 & typeof(thin) == 'double' & length(thick) == 1 & typeof(thick) == 'double')
    {
      chains <- as.data.frame(data)
      len <- 1
      idx <- 1

      thick_ci <- c((100-((100-thick)/2)), ((100-thick)/2))*0.01
      thin_ci <- c((100-((100-thin)/2)), ((100-thin)/2))*0.01

      thick_q <- as.matrix(apply(chains, 2, stats::quantile, probs= thick_ci)[,idx])
      thin_q <- as.matrix(apply(chains, 2, stats::quantile, probs= thin_ci)[,idx])

      medians <- apply(chains, 2, stats::quantile, probs = 0.5)[idx]


    }else
    {
      stop("'thick' and 'thin' must be single numbers")
    }
  }


  # plotting ----------------------------------------------------------------

  #Determine which params have CI that overlap 0 (or ref line more technically)
  black_cl <- c() #95% CI (default) does not overlap 0
  gray_cl <- c() #50% CI (default) does not overlap 0
  white_cl <- c() #Both 50% and 95% CI (default) overlap 0

  if(!is.null(ref))
  {
    marker <- ref
  }else {
    marker <- 0
  }

  for (i in 1:len)
  {
    #i <- 1
    if ((thin_q[1,i] > marker & thin_q[2,i] > marker) |
        (thin_q[1,i] < marker & thin_q[2,i] < marker))
    {
      black_cl <- c(black_cl, i)
    } else {
      if ((thick_q[1,i] > marker & thick_q[2,i] > marker) |
          (thick_q[1,i] < marker & thick_q[2,i] < marker))
      {
        gray_cl <- c(gray_cl, i)
      }else {
        white_cl <- c(white_cl, i)
      }
    }
  }

  #positions to plot CI
  blk_bnd <- rbind(black_cl, black_cl)
  gry_bnd <- rbind(gray_cl, gray_cl)
  wht_bnd <- rbind(white_cl, white_cl)

  #plot for horizontal - CI lines parallel to x-axis
  if (horiz == TRUE)
  {

    if (missing(xlim))
    {xlim = range(thin_q)*1.2}
    ylim = c(0.5,(len)+0.5)
    if (missing(xlab))
    {xlab = 'Parameter Estimate'}
    if (is.null(xlab))
    {xlab = ''}
    if (missing(main))
    {main = ''}
    if (missing(labels))
    {
      labs = names(medians)
    }else{
      if (!missing(labels))
      {
        if (is.null(labels))
        {
          labs <- rep('', len)
        }
        if (!is.null(labels))
        {
          if (length(labels) == len)
          {
            labs <- labels[idx]
          }else
          {
            stop('Labels length not equal to number of parameters')
          }
        }
      }
    }

    #0.2 inches per line - mar measured in lines
    m_char <- (max(sapply(labs, function(x){graphics::strwidth(x, cex = labels_sz, units = 'in')}))/0.2)

    graphics::par(mar=c(mar[1], (m_char + (mar[2] - 3)), mar[3], mar[4]-1))

    #plot blank plot
    graphics::plot(medians, (1:len), xlim = xlim, ylim = ylim, type = "n",
         ann = TRUE, xaxt = 'n', yaxt = "n", bty = "n", ylab = NA,
         xlab = xlab, cex.lab = x_axis_text_sz) #cex.lab is axis label
    #lab #number of ticks to plot on each axis

    #title
    graphics::title(main, cex.main = main_text_sz)
    #bottom axis params
    graphics::axis(3, lwd.ticks = ax_sz, labels = FALSE,
         at = tick_pos, lwd = ax_sz)
    #bottom axis params
    graphics::axis(1, lwd.ticks = ax_sz, labels = TRUE,
         at = tick_pos, lwd = ax_sz,
         cex.axis = x_tick_text_sz) #bottom axis
    #left axis params (labels)
    graphics::axis(2, at = ((1:len)+(0.0015*len)), tick = FALSE,
         labels = labs, las = 1, adj = 0, #las - 0 parallel to axis, 1 horiz, 2 perp to axis, 3 vert
         line = -1, cex.axis = labels_sz)

    #ref line
    if(!is.null(ref))
    {
      graphics::abline(v=ref, lty = 2, lwd = 3, col = ref_col)
    }

    if (ref_ovl == TRUE)
    {
      #Black CI
      if (!is.null(black_cl))
      {
        #Thick
        graphics::matlines(thick_q[,black_cl], blk_bnd,
                 type = 'l', lty = 1, lwd = thick_sz, col = 'black')
        #Thin
        graphics::matlines(thin_q[,black_cl], blk_bnd,
                 type = 'l', lty = 1, lwd = thin_sz, col = 'black')
      }

      #Gray CI
      if (!is.null(gray_cl))
      {
        #Thick
        graphics::matlines(thick_q[,gray_cl], gry_bnd,
                 type = 'l', lty = 1, lwd = thick_sz, col = gr_col)
        #Thin
        graphics::matlines(thin_q[,gray_cl], gry_bnd,
                 type = 'l', lty = 1, lwd = thin_sz, col = gr_col)
      }

      #White CI
      if (!is.null(white_cl))
      {
        graphics::matlines(thick_q[,white_cl], wht_bnd,
                 type = 'l', lty = 1, lwd = thick_sz, col = gr_col) #white (gray)
        graphics::matlines(thin_q[,white_cl], wht_bnd,
                 type = 'l', lty = 1, lwd = thin_sz, col = gr_col) #white (gray)
      }

      #Medians
      graphics::points(medians, 1:len, pch = 16, col = 'white', cex = med_sz)
      graphics::points(medians[black_cl], black_cl, pch = 16, col = 'black', cex = med_sz)
      graphics::points(medians[gray_cl], gray_cl, pch = 16, col = gr_col, cex = med_sz)
      graphics::points(medians[white_cl], white_cl, pch = 21, col = gr_col, cex = med_sz, lwd = 2)
    } else{
      graphics::matlines(thick_q[,1:len], rbind(1:len, 1:len),
               type = 'l', lty = 1, lwd = thick_sz, col = 'black')
      graphics::matlines(thin_q[,1:len], rbind(1:len, 1:len),
               type = 'l', lty = 1, lwd = thin_sz, col = 'black')
      #medians
      graphics::points(medians[1:len], 1:len, pch = 16,
             col = 'black', cex = med_sz)
    }
  }

  #vertical plot - CI lines perpendicular to x-axis
  if (horiz == FALSE)
  {

    if (missing(ylim))
    {ylim = range(thin_q)*1.2}
    xlim = c(0.5,(len)+0.5)
    if (missing(ylab))
    {ylab = 'Parameter Estimate'}
    if (is.null(ylab))
    {ylab = ''}
    if (missing(main))
    {main = ''}
    if (missing(labels))
    {
      labs = names(medians)
    }else{
      if (!missing(labels))
      {
        if (is.null(labels))
        {
          labs <- rep('', len)
        }
        if (!is.null(labels))
        {
          if (length(labels) == len)
          {
            labs <- labels[idx]
          }else
          {
            stop('Labels length not equal to number of parameters')
          }
        }
      }
    }

    #0.2 inches per line - mar measured in lines
    m_char <- (max(sapply(labs, function(x){graphics::strwidth(x, cex = labels_sz, units = 'in')}))/0.2)
    graphics::par(mar = c((m_char + (mar[1] - 4)), mar[2]+1, mar[3] - 1, mar[4]))

    #plot blank plot
    graphics::plot((len:1), medians, xlim = xlim, ylim = ylim, type = "n",
                   ann = TRUE, xaxt = 'n', yaxt = "n", bty = "n", ylab = ylab,
                   xlab = NA, cex.lab = x_axis_text_sz)

    #title
    graphics::title(main, cex.main = main_text_sz)
    #right axis params
    graphics::axis(4, lwd.ticks = ax_sz, labels = FALSE,
                   at = tick_pos, lwd = ax_sz)
    #left axis params
    graphics::axis(2, lwd.ticks = ax_sz, labels = TRUE,
                   at = tick_pos, lwd = ax_sz,
                   cex.axis = x_tick_text_sz) #bottom axis
    #bottom axis params (labels)
    graphics::axis(1, at = (len:1) + 0.013, tick = FALSE,
                   labels = labs, las = 2, adj = 0, #las - 0 parallel to axis, 1 horiz, 2 perp to axis, 3 vert
                   line = -1, cex.axis = labels_sz)

    #ref line
    if(!is.null(ref))
    {
      graphics::abline(h=ref, lty = 2, lwd = 3, col = ref_col)
    }

    v_black_cl <- len+1 - black_cl
    v_gray_cl <- len+1 - gray_cl
    v_white_cl <- len+1 - white_cl

    if (ref_ovl == TRUE)
    {
      #Black CI
      if (!is.null(black_cl))
      {
        #Thick
        graphics::matlines(blk_bnd, thick_q[,v_black_cl],
                           type = 'l', lty = 1, lwd = thick_sz, col = 'black')
        #Thin
        graphics::matlines(blk_bnd, thin_q[,v_black_cl],
                           type = 'l', lty = 1, lwd = thin_sz, col = 'black')
      }

      #Gray CI
      if (!is.null(gray_cl))
      {
        #Thick
        graphics::matlines(gry_bnd, thick_q[,v_gray_cl],
                           type = 'l', lty = 1, lwd = thick_sz, col = gr_col)
        #Thin
        graphics::matlines(gry_bnd, thin_q[,v_gray_cl],
                           type = 'l', lty = 1, lwd = thin_sz, col = gr_col)
      }

      #White CI
      if (!is.null(white_cl))
      {
        graphics::matlines(wht_bnd, thick_q[,v_white_cl],
                           type = 'l', lty = 1, lwd = thick_sz, col = gr_col) #white (gray)
        graphics::matlines(wht_bnd, thin_q[,v_white_cl],
                           type = 'l', lty = 1, lwd = thin_sz, col = gr_col) #white (gray)
      }

      #Medians
      graphics::points(len:1, medians, pch = 16, col = 'white', cex = med_sz)
      graphics::points(black_cl, medians[v_black_cl], pch = 16, col = 'black', cex = med_sz)
      graphics::points(gray_cl, medians[v_gray_cl], pch = 16, col = gr_col, cex = med_sz)
      graphics::points(white_cl, medians[v_white_cl], pch = 21, col = gr_col, cex = med_sz, lwd = 2)
    } else{
      graphics::matlines(rbind(1:len, 1:len), thick_q[,len:1],
                         type = 'l', lty = 1, lwd = thick_sz, col = 'black')
      graphics::matlines(rbind(1:len, 1:len), thin_q[,len:1],
                         type = 'l', lty = 1, lwd = thin_sz, col = 'black')
      #medians
      graphics::points(1:len, medians[len:1], pch = 16,
                       col = 'black', cex = med_sz)
    }
  }


  graphics::par(mar=c(5,4,4,2) + 0.1)

}
