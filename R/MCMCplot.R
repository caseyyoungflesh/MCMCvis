#' Caterpillar plots of posterior distributions from MCMC output
#'
#' Visualize posterior distributions from MCMC output for specific parameters of interest using caterpillar plots. Color of median dot represents the overlap of the posterior distribution with 0 (or other specified value).
#'
#'
#' @param object Object containing MCMC output. See DETAILS below.
#' 
#' @param object2 Optional second object containing MCMC output. If specified, parameter estimates from each model will be displayed in a paired manner. Parameter names for \code{'object'} and \code{'object2'} must be identical. See DETAILS below.
#' 
#' @param params Character string (or vector of character strings) denoting parameters to be plotted.
#'
#' Default \code{'all'} plots posteriors for all parameters. See VALUE below.
#'
#' @param HPD Logical specifying whether to calculate equal-tailed credible intervals (\code{HPD = FALSE}) or highest posterior density intervals (\code{HPD = TRUE}) for the selected parameters. Default is \code{HPD = FALSE}.
#' 
#' @param ci Numeric vector of length 2, where each element is (0,100] and represents the width of an equal-tailed (\code{HPD = FALSE}) or highest posterior density (\code{HPD = TRUE}) credible interval. The first element of this vector corresponds to the thicker (narrower) credible interval displayed on the plot (default is 0.5) and the second element of this vector corresponds to the thinner (wider) credible interval (default is 0.95). The first credible interval width (\code{ci[1]}) must be less than or equal to the width of the second credible interval (\code{ci[2]}).
#'
#' @param excl Character string (or vector of character strings) denoting parameters to exclude. Used in conjunction with \code{params} argument to select parameters of interest.
#'
#' @param ISB Ignore Square Brackets (ISB). Logical specifying whether square brackets should be ignored in the \code{params} and \code{excl} arguments. If \code{TRUE}, square brackets are ignored. If \code{FALSE}, square brackets are not ignored.  This allows partial names to be used when specifying parameters of interest. Use \code{exact} argument to specify whether input from \code{params} and \code{excl} arguments should be matched exactly.
#'
#' @param exact Logical specifying whether input from \code{params} and \code{excl} arguments should be matched exactly (after ignoring square brackets if \code{ISB = FALSE}). If \code{TRUE}, input from \code{params} and \code{excl} are matched exactly (after taking \code{ISB} argument into account). If \code{FALSE}, input from \code{params} and \code{excl} are matched using regular expression format (after taking \code{ISB} argument into account).
#'
#' @param ref Value indicating where vertical reference line should be created and what value to use a reference for caterpillar median coloration.
#'
#' Default is \code{ref = 0}.
#'
#' Argument \code{NULL} will plot no reference line.
#'
#' @param ref_ovl Logical specifying whether the style/color of plotted median dots and CI should be changed based on whether the specified credible intervals (50 \% and 95 \% by default) overlap the reference line. See DETAILS for more information.
#'
#' @param col Character string (or vector of character strings) specifying which color to render estimates on plot. When \code{ref_ovl = TRUE}, this argument has no effect and colors plotted will be based on the credible intervals and reference line. Number of specified colors must equal the number of specified parameters or one.
#' 
#' @param col2 Character string (or vector of character strings) specifying which color to render estimates on plot for \code{object2} (if specified). Number of specified colors must equal the number of specified parameters or one. Red by default.
#'
#' @param offset Value indicating how much to offset plotted posteriors when \code{object2} is specified (i.e., control the amount of space between the two  plotted posteriors for each parameter). The distance from one set of parameters to another corresponds to a value of 1.
#'
#' @param rank Logical specifying whether output should be ranked. If \code{TRUE} posteriors will be ranked in decreasing order (based on specified measure of centrality) from top down.
#'
#' @param horiz Logical specifying orientation of plot. If \code{TRUE} posteriors will be plotted running horizontally (parallel to the x-axis). If \code{FALSE} posteriors will be plotted running vertically (perpendicular to the x-axis).
#' 
#' @param xlim Numerical vector of length 2, indicating range of x-axis. Only applicable if \code{horiz = TRUE}.
#' 
#' @param ylim Numerical vector of length 2, indicating range of y-axis. Only applicable if \code{horiz = FALSE}.
#' 
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
#' @param guide_lines Logical specifying whether to plot reference lines for each parameter in order to better visualize which parameter names correspond to each posterior.
#'
#' @param guide_axis Logical specifying whether a second axis should be plotted (x-axis if \code{HORIZ = TRUE}, y-axis if \code{HORIZ = FALSE}) to help interpret values on plot.
#'
#' @param sz_labels Number specifying size of text for parameter labels on axis.
#'
#' @param sz_med Number specifying size of points represents posterior medians.
#' 
#' @param sz_thick Number specifying thickness of 50 percent CI line (thicker line).
#'
#' @param sz_thin Number specifying thickness of 95 percent CI line (thinner line).
#'
#' @param sz_ax Number specifying thickness of axis and ticks.
#'
#' @param sz_ax_txt Number specifying size of text for axis label.
#'
#' @param sz_tick_txt Number specifying size of text for tick labels on axis.
#'
#' @param sz_main_txt Number specifying size of text for main title.
#'
#' @param pos_tick Numeric vector specifying where ticks on axis should be placed.
#'
#' @param mar Numerical vector of length 4 specifying plot margins - (BOTTOM, LEFT, TOP, RIGHT). Changes to the margin should be made within the function rather than using the \code{par} call.
#'
#' Default is c(5.1, 4.1, 4.1, 2.1) - the R plot default.
#'
#' @section Details:
#' Points represent posterior medians. Parameters where the smaller specified credible intervals (first value in vector provided to \code{ci} argument; 50\% by default) overlap 0 (or other value specified by the \code{ref} argument) are indicated by 'open' circles and a lighter color than that specified (e.g., gray will be displayed with the default black is specified for \code{col}). Parameters where the specific smaller credible intervals DO NOT overlap 0 (or other specified value) AND the larger specified credible intervals (second value in vector provided to \code{ci} argument; 95\% by default) DO overlap 0 (or other specified value) are indicated by 'closed' circles and a lighter color that that specified. Parameters where the larger specified credible intervals DO NOT overlap 0 (or other specified value) are indicated by 'closed' circles and the color specified (black by default). Thick lines represent the smaller specified credible intervals percent credible intervals (50\% by default) while thin lines represent the larger specified credible intervals (95 \% by default). \code{ref_ovl = TRUE} can be used to enable this feature. When two model objects are supplied to the function (with \code{object} and \code{object2}) and no argument is supplied to \code{col} or \code{col2}, light red lines can be interpreted as analogous to gray lines. When \code{ref_ovl = TRUE} and a color (or colors) other than the default are specified (with the \code{col} and \code{col2} arguments) lighter versions of the color specified are used in black of the light gray and/or light red lines.
#' 
#' When \code{object2} is specified, paired caterpillar plots of each parameter are produced. For this reason, parameter names of \code{object} and \code{object2} specified with the \code{params} argument must be identical (to be used for comparing posterior estimates of similar models). \code{col} and \code{col2} arguments can be specified to change the color of output from \code{object} and \code{object2}, respectively. By default, output from \code{object} is plotted in black and \code{object2} is plotted in red. The \code{ref_ovl} argument can also be specified.
#'
#' \code{object} argument can be a \code{stanfit} object (\code{rstan} package), a \code{stanreg} object (\code{rstanarm} package), a \code{brmsfit} object (\code{brms} package), an \code{mcmc.list} object (\code{coda} and \code{rjags} packages), \code{mcmc} object (\code{coda} and \code{nimble} packages), \code{list} object (\code{nimble} package), an \code{R2jags} model object (\code{R2jags} package), a \code{jagsUI} model object (\code{jagsUI} package), or a matrix containing MCMC chains (each column representing MCMC output for a single parameter, rows representing iterations in the chain). The function automatically detects the object type and proceeds accordingly.
#'
#' @section Notes:
#'
#' When specifying \code{rank = TRUE} and specifying labels for \code{labels}, labels will be applied to parameters before they are ranked.
#'
#' Thanks to Cinner et al. 2016, whose Fig. 1 inspired this plot.
#'
#' @section References:
#'
#' Cinner, J. E., C. Huchery, M. A. MacNeil, N. A. J. Graham, T. R. McClanahan, J. Maina, E. Maire, J. N. Kittinger, C. C. Hicks, C. Mora, E. H. Allison, S. D'Agata, A. Hoey, D. A. Feary, L. Crowder, I. D. Williams, M. Kulbicki, L. Vigliola, L. Wantiez, G. Edgar, R. D. Stuart-Smith, S. A. Sandin, A. L. Green, M. J. Hardt, M. Beger, A. Friedlander, S. J. Campbell, K. E. Holmes, S. K. Wilson, E. Brokovich, A. J. Brooks, J. J. Cruz-Motta, D. J. Booth, P. Chabanet, C. Gough, M. Tupper, S. C. A. Ferse, U. R. Sumaila, and D. Mouillot. 2016. Bright spots among the world's coral reefs. Nature 535:416-419.
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
#' #Just 'beta' parameters using highest posterior density intervals
#' MCMCplot(MCMC_data, params = 'beta', HPD = TRUE)
#
#' #Just 'beta[1]', 'beta[4]', and 'alpha[3]'
#' MCMCplot(MCMC_data, params = c('beta[1]', 'beta[4]', 'alpha[3]'), ISB = FALSE, exact = TRUE)
#'
#' #Just 'beta[1]', 'beta[4]', and 'alpha[3]' and change the credible interval widths
#' MCMCplot(MCMC_data, ci = c(50, 89), params = c('beta[1]', 'beta[4]', 'alpha[3]'), 
#'   ISB = FALSE, exact = TRUE)
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
                     object2 = NULL,
                     params = 'all',
                     HPD = FALSE,
                     ci = c(50, 95),
                     excl = NULL,
                     ISB = TRUE,
                     exact = TRUE,
                     ref = 0,
                     ref_ovl = FALSE,
                     col = 'black',
                     col2 = 'red',
                     offset = 0.1,
                     rank = FALSE,
                     horiz = TRUE,
                     xlim,
                     ylim,
                     xlab,
                     ylab,
                     main,
                     labels,
                     guide_lines = FALSE,
                     guide_axis = TRUE,
                     sz_labels = 1.2,
                     sz_med = 1.5,
                     sz_thick = 5,
                     sz_thin = 2,
                     sz_ax = 3,
                     sz_ax_txt = 1.3,
                     sz_tick_txt = 1.2,
                     sz_main_txt = 1.2,
                     pos_tick,
                     mar = c(5.1, 4.1, 4.1, 2.1)) {

####################################################################################################  
# load object(s)
  
  data <- MCMCchains(object, params = params, excl = excl, ISB = ISB, exact = exact)
  if (!is.null(object2)) {
    data2 <- MCMCchains(object2, params = params, excl = excl, ISB = ISB, exact = exact)
    if (!identical(colnames(data), colnames(data2))) { 
      stop('`object2` must have identical column names to `object`.')
    }
  }

####################################################################################################  
# user and non-user defined plotting parameters  

  if (missing(pos_tick)) pos_tick = NULL
  LA <- 0.4 #amount to lighten colors when ref_ovl = TRUE
  ref_col <- 'gray60' # color used for 0 line
  guide_col <- 'gray80'
  PL_SC <- 0.3 # how much whitespace flanks plotted estimates
  
  if (length(ci) == 2 & typeof(ci) == 'double' & ci[1] <= ci[2] & ci[1] > 0 & 
    ci[1] <= 100 & ci[2] > 0 & ci[2] <= 100) {
      thin <- ci[2] # CI for thin line
      thick <- ci[1]  # CI for thick line
  } else {
    stop("ci must be a numeric vector of length 2 where each element is (0, 100] and ci[1] <= ci[2]")  
  }

  # Colors 
  # number of parameters
  np <- NCOL(data)
  if (np > 1) {
    if (length(col) == 1) {
      COL <- rep(col, np)
    } else {
      if (length(col) != np) {
        stop('Number of specified colors must equal number of plotted parameters (or one).')
      } else {
        COL <- col
      }
    } # closes if (length(col) == 1)
    if (!is.null(object2)) {
      if (length(col2) == 1) {
        COL2 <- rep(col2, np)
      } else {
        if (length(col2) != np) {
          stop('Number of specified colors must equal number of plotted parameters (or one).')
        } else {
          COL2 <- col2
        }
      } # closes if (length(col2) == 1)
    } # closes if (!is.null(object2))
  } # closes f (np > 1)
  if (np == 1) {
    if (length(col) == 1) {
      COL <- col
    } else {
      stop('Number of specified colors must equal number of plotted parameters (or one).')
    }
    if (!is.null(object2)) {
      if (length(col2) == 1) {
        COL2 <- col2
      } else {
        stop('Number of specified colors must equal number of plotted parameters (or one).')
      }
    } # closes if (!is.null(object2))
  } # closes if (np == 1)

####################################################################################################  
# Process data 

  pro_fun <- function(input, ...) {
    chains <- as.data.frame(input)
    tsrt <- apply(chains, 2, stats::median)
    if (np > 1) {
      len <- np
      if (rank == TRUE) idx <- order(tsrt, decreasing = TRUE)
      if (rank == FALSE) idx <- len:1
    } 
    if (np == 1) {
      len <- 1
      idx <- 1
    } 
    medians <- tsrt[idx]
    if (HPD == FALSE) {
      thick_ci <- c((100 - ((100 - thick) / 2)), ((100 - thick) / 2)) * 0.01
      thin_ci <- c((100 - ((100 - thin) / 2)), ((100 - thin) / 2)) * 0.01

      if (np == 1)
      {
        thick_q <- as.matrix(apply(chains, 2, function(x) stats::quantile(x, probs = thick_ci, na.rm = TRUE))[, idx], nrow = 2)
        thin_q <- as.matrix(apply(chains, 2, function(x) stats::quantile(x, probs = thin_ci, na.rm = TRUE))[, idx], nrow = 2)
      } else {
        thick_q <- apply(chains, 2, function(x) stats::quantile(x, probs = thick_ci, na.rm = TRUE))[, idx]
        thin_q <- apply(chains, 2, function(x) stats::quantile(x, probs = thin_ci, na.rm = TRUE))[, idx]
      }
    } else {
      if (np == 1)
      {
        thick_q <- as.matrix(coda::HPDinterval(coda::as.mcmc(chains), prob = ci[1] / 100)[idx, ])
        thin_q <- as.matrix(coda::HPDinterval(coda::as.mcmc(chains), prob = ci[2] / 100)[idx, ])
      }
       else {
        thick_q <- t(coda::HPDinterval(coda::as.mcmc(chains), prob = ci[1] / 100)[idx, ])
        thin_q <- t(coda::HPDinterval(coda::as.mcmc(chains), prob = ci[2] / 100)[idx, ])
      }  
    }

    return(list(len, idx, thick_q, thin_q, medians))
  } # closes pro_fun
  
  pf_out <- pro_fun(input = data)
  # assign list elements to objects
  len <- pf_out[[1]]
  idx <- pf_out[[2]]
  thick_q <- pf_out[[3]]
  thin_q <- pf_out[[4]]
  medians <- pf_out[[5]]
  if (!is.null(object2)) {
    pf_out2 <- pro_fun(input = data2)
    thick_q2 <- pf_out2[[3]]
    thin_q2 <- pf_out2[[4]]
    medians2 <- pf_out2[[5]]
  } else {
    thin_q2 <- NULL
  }
 
####################################################################################################  
# Plot prep 

  # reference line
  if(!is.null(ref)) {
    marker <- ref
  } else {
    marker <- 0
  }

  # position for coloring with regards to reference line for credible intervals
  if (ref_ovl == TRUE) {
    # Determine which params have CI that overlap 0 (or ref line if specified)
    # kept object names to avoid changing variable names - SHOULD CHANGE
    black_cl <- c() #95% CI (default) does not overlap 0
    gray_cl <- c() #50% CI (default) does not overlap 0
    white_cl <- c() #Both 50% and 95% CI overlap 0
    for (i in 1:len) {
      if ((thin_q[1, i] > marker & thin_q[2, i] > marker) | 
          (thin_q[1, i] < marker & thin_q[2, i] < marker)) {
        black_cl <- c(black_cl, i)
      } else {
        if ((thick_q[1, i] > marker & thick_q[2, i] > marker) |
            (thick_q[1, i] < marker & thick_q[2, i] < marker)) {
          gray_cl <- c(gray_cl, i)
        } else {
          white_cl <- c(white_cl, i)
        }
      }
    } # closes for (i in 1:len)
    if (horiz == FALSE) {
      v_black_cl <- rev(black_cl)
      v_gray_cl <- rev(gray_cl)
      v_white_cl <- rev(white_cl)
    }
    if (!is.null(object2)) {
      black_cl2 <- c()
      gray_cl2 <- c()
      white_cl2 <- c()
      for (i in 1:len) {
        if ((thin_q2[1, i] > marker & thin_q2[2, i] > marker) |
            (thin_q2[1, i] < marker & thin_q2[2, i] < marker)) {
          black_cl2 <- c(black_cl2, i)
        } else {
          if ((thick_q2[1, i] > marker & thick_q2[2, i] > marker) |
              (thick_q2[1, i] < marker & thick_q2[2, i] < marker)) {
            gray_cl2 <- c(gray_cl2, i)
          } else {
            white_cl2 <- c(white_cl2, i)
          }
        }
      } # closes for (i in 1:len)
      if (horiz == FALSE) {
        v_black_cl2 <- rev(black_cl2)
        v_gray_cl2 <- rev(gray_cl2)
        v_white_cl2 <- rev(white_cl2)
      }
    } # closes (!is.null(object2))
  } # closes (ref_ovl  == TRUE)

####################################################################################################  
# Horizontal plot - CI lines parallel to x-axis 
  
  if (horiz == TRUE) {
    if (missing(xlim)) {
      thin_comb <- c(thin_q, thin_q2)
      rn <- diff(range(thin_comb, na.rm = TRUE)) * PL_SC
      XLIM <- c((min(thin_comb, na.rm = TRUE) - rn), (max(thin_comb, na.rm = TRUE) + rn))
    } else {
      XLIM <- xlim
    }
    if (len > 20) {
      tt <- 1.5
    } else {
      alpha <- 0.667
      beta <- 0.0412
      tt <- alpha + beta * len
    } # closes if (len > 20)
    YLIM <- c((1 - tt), (len + tt))
    if (missing(xlab)) xlab <- 'Parameter Estimate'
    if (is.null(xlab)) xlab <- ''
    if (missing(main)) main <- ''
    if (missing(labels)) {
      labs <- names(medians)
    } else {
      if (!missing(labels)) {
        if (is.null(labels)) {
          labs <- rep('', len)
        }
        if (!is.null(labels)) {
          if (length(labels) == len) {
            labs <- labels[idx]
          } else {
            stop('Labels length not equal to number of parameters')
          }
        }
      }
    } # closes if (missing(labels))
    
    # 0.2 inches per line - mar measured in lines
    m_char <- max(sapply(labs, function(x){graphics::strwidth(x, cex = sz_labels, units = 'in')})) / 0.2
    graphics::par(mar = c(mar[1], (m_char + (mar[2] - 3)), mar[3], mar[4] - 1))

    #plot blank plot
    graphics::plot(
      medians, 
      (1:len), 
      xlim = XLIM, 
      ylim = YLIM, 
      type = "n", 
      ann = TRUE, 
      xaxt = 'n', 
      yaxt = "n", 
      bty = "n", 
      ylab = NA, 
      xlab = xlab, 
      cex.lab = sz_ax_txt, 
      yaxs = 'i') #cex.lab is axis label

    # title
    graphics::title(main, cex.main = sz_main_txt)

    if (guide_axis == TRUE) {
      # top axis params
      graphics::axis(3, lwd.ticks = sz_ax, labels = FALSE, at = pos_tick, lwd = sz_ax)
    }
    # bottom axis params
    graphics::axis(1, lwd.ticks = sz_ax, labels = TRUE, at = pos_tick, lwd = sz_ax, cex.axis = sz_tick_txt)
     
     # left axis params (labels)
     # las - 0 parallel to axis, 1 horiz, 2 perp to axis, 3 vert
     graphics::axis(
       2, 
       at = ((1:len)), 
       tick = FALSE,
       labels = labs, 
       las = 1, 
       adj = 0,
       line = -1, 
       cex.axis = sz_labels)

     if (guide_lines == TRUE) {
       # limits
       # determine where ticks were placed - from axTicks source code:
       xaxp <- graphics::par("xaxp")
       XLIM2 <- c(xaxp[1], xaxp[2])
       xs <- matrix(rep(XLIM2, len), nrow = 2)
       ys <- rbind(1:len, 1:len)
       graphics::matlines(xs, ys, lty = 1, col = guide_col)
     }
     
     # ref line
     if(!is.null(ref)) graphics::abline(v=ref, lty = 2, lwd = 3, col = ref_col)
     if (ref_ovl == TRUE) {
       # positions to plot CI
       blk_bnd <- rbind(black_cl, black_cl)
       gry_bnd <- rbind(gray_cl, gray_cl)
       wht_bnd <- rbind(white_cl, white_cl)
       
       #lightened colors
       lcol <- colorspace::lighten(COL, LA)

       bl_col <- COL[black_cl]
       gr_col <- lcol[gray_cl]
       wh_col <- lcol[white_cl]

       #if any colors are 'black', change them to 'grey60' (for consistency with older pkg versions)
       bl_idx <- which(COL == 'black')
       gr_col[which(gray_cl %in% bl_idx)] <- 'grey60'
       wh_col[which(white_cl %in% bl_idx)] <- 'grey60'
       
       # Black CI
       if (!is.null(black_cl)) {
         #only one object
         if (is.null(object2)) {
           # Thick
           graphics::matlines(
             thick_q[,black_cl], 
             blk_bnd, 
             type = 'l', 
             lty = 1, 
             lwd = sz_thick, 
             col = bl_col)
           # Thin
           graphics::matlines(
             thin_q[,black_cl], 
             blk_bnd, 
             type = 'l', 
             lty = 1, 
             lwd = sz_thin, 
             col = bl_col)
         } else {
           #two objects
           blk_bnd2 <- rbind(black_cl2, black_cl2)
           
           bl_col2 <- COL2[black_cl2]
           
           # object
           graphics::matlines(
             thick_q[,black_cl], 
             (blk_bnd + offset), 
             type = 'l', 
             lty = 1, 
             lwd = sz_thick, 
             col = bl_col)
           graphics::matlines(
             thin_q[,black_cl], 
             (blk_bnd + offset),
             type = 'l', 
             lty = 1, 
             lwd = sz_thin, 
             col = bl_col)
           # object2
           graphics::matlines(
             thick_q2[, black_cl2], 
             (blk_bnd2 - offset),
             type = 'l', 
             lty = 1, 
             lwd = sz_thick, 
             col = bl_col2)
           graphics::matlines(
             thin_q2[, black_cl2], 
             (blk_bnd2 - offset),
             type = 'l', 
             lty = 1, 
             lwd = sz_thin, 
             col = bl_col2)
         } # closes if (is.null(object2))
       } # closes if (!is.null(black_cl))

       # Gray CI
       if (!is.null(gray_cl)) {
         #only one object
         if (is.null(object2)) {
          # Thick
          graphics::matlines(
            thick_q[,gray_cl], 
            gry_bnd,
            type = 'l', 
            lty = 1, 
            lwd = sz_thick, 
            col = gr_col)
          # Thin
          graphics::matlines(
            thin_q[,gray_cl], 
            gry_bnd,
            type = 'l', 
            lty = 1, 
            lwd = sz_thin, 
            col = gr_col)
         } else {
           #two objects
           gry_bnd2 <- rbind(gray_cl2, gray_cl2)
           
           #lightened colors for obj 2
           lcol2 <- colorspace::lighten(COL2, LA)
           gr_col2 <- lcol2[gray_cl2]
           
           #if any colors are 'black', change them to 'grey60' (for consistency with older pkg versions)
           bl_idx2 <- which(COL2 == 'black')
           gr_col2[which(gray_cl2 %in% bl_idx2)] <- 'grey60'
           
           # object
           graphics::matlines(
             thick_q[,gray_cl], 
             (gry_bnd + offset),
             type = 'l', 
             lty = 1, 
             lwd = sz_thick,
             col = gr_col)
           graphics::matlines(
             thin_q[,gray_cl],
             (gry_bnd + offset),
             type = 'l', 
             lty = 1, 
             lwd = sz_thin,
             col = gr_col)
           #object2
           graphics::matlines(
             thick_q2[,gray_cl2], 
             (gry_bnd2 - offset),
             type = 'l', 
             lty = 1, 
             lwd = sz_thick,
             col = gr_col2)
           graphics::matlines(
             thin_q2[,gray_cl2],
             (gry_bnd2 - offset),
             type = 'l', 
             lty = 1, 
             lwd = sz_thin,
             col = gr_col2)
         } # closes if (is.null(object2))
       } # closes if (!is.null(gray_cl))

       # White CI
       if (!is.null(white_cl)) {
         #only one object
         if (is.null(object2)) {
          graphics::matlines(
            thick_q[,white_cl], 
            wht_bnd,
            type = 'l', 
            lty = 1, 
            lwd = sz_thick,
            col = wh_col) # white (gray)
          graphics::matlines(
            thin_q[,white_cl], 
            wht_bnd,
            type = 'l',
            lty = 1,
            lwd = sz_thin, 
            col = wh_col) # white (gray)
         } else {
           #two objects
           wht_bnd2 <- rbind(white_cl2, white_cl2)
           
           #lightened colors for obj 2
           lcol2 <- colorspace::lighten(COL2, LA)
           wh_col2 <- lcol2[white_cl2]
           
           #if any colors are 'black', change them to 'grey60' (for consistency with older pkg versions)
           bl_idx2 <- which(COL2 == 'black')
           wh_col2[which(white_cl2 %in% bl_idx2)] <- 'grey60'
           
           # object
           graphics::matlines(
             thick_q[,white_cl], 
             (wht_bnd + offset),
             type = 'l', 
             lty = 1, 
             lwd = sz_thick,
             col = wh_col) # white (gray)
           graphics::matlines(
             thin_q[,white_cl], 
             (wht_bnd + offset),
             type = 'l', 
             lty = 1,
             lwd = sz_thin, 
             col = wh_col) # white (gray)
           # object2
           graphics::matlines(
             thick_q2[,white_cl2], 
             (wht_bnd2 - offset),
             type = 'l', 
             lty = 1, 
             lwd = sz_thick,
             col = wh_col2) # white (gray)
           graphics::matlines(
             thin_q2[,white_cl2], 
             (wht_bnd2 - offset),
             type = 'l',
             lty = 1, 
             lwd = sz_thin,
             col = wh_col2) #white (gray)
         } # closes if (is.null(object2))
       } # closes  if (!is.null(white_cl))

       # Medians
       if (is.null(object2)) {
         #only one object
        graphics::points(medians, 1:len, pch = 16, col = 'white', cex = sz_med)
        graphics::points(medians[black_cl], black_cl, pch = 16, col = bl_col, cex = sz_med)
        graphics::points(medians[gray_cl], gray_cl, pch = 16, col = gr_col, cex = sz_med)
        graphics::points(medians[white_cl], white_cl, pch = 21, col = wh_col, cex = sz_med, lwd = 2) 
       } else {
         # two objects
         graphics::points(medians, (1:len + offset), pch = 16, col = 'white', cex = sz_med)
         graphics::points(medians[black_cl], (black_cl + offset), pch = 16, col = bl_col, cex = sz_med)
         graphics::points(medians[gray_cl], (gray_cl + offset), pch = 16, col = gr_col, cex = sz_med)
         graphics::points(medians[white_cl], (white_cl + offset), pch = 21, col = wh_col, cex = sz_med, lwd = 2) 
         
         # object2
         graphics::points(medians2, (1:len - offset), pch = 16, col = 'white', cex = sz_med)
         graphics::points(medians2[black_cl2], (black_cl2 - offset), pch = 16, col = bl_col2, cex = sz_med)
         graphics::points(medians2[gray_cl2], (gray_cl2 - offset), pch = 16, col = gr_col2, cex = sz_med)
         graphics::points(medians2[white_cl2], (white_cl2 - offset), pch = 21, col = wh_col2, cex = sz_med, lwd = 2) 
       }
     } else {
       if (is.null(object2))
       {
         graphics::matlines(thick_q[,1:len], rbind(1:len, 1:len), type = 'l', lty = 1, lwd = sz_thick, col = COL[idx])
         graphics::matlines(thin_q[,1:len], rbind(1:len, 1:len), type = 'l', lty = 1, lwd = sz_thin, col = COL[idx])
         graphics::points(medians[1:len], 1:len, pch = 16, col = COL[idx], cex = sz_med)
       } else {
         # object
         graphics::matlines(thick_q[,1:len], (rbind(1:len, 1:len) + offset), type = 'l', 
           lty = 1, lwd = sz_thick, col = COL[idx])
         graphics::matlines(thin_q[,1:len], (rbind(1:len, 1:len) + offset), type = 'l', 
           lty = 1, lwd = sz_thin, col = COL[idx])
         graphics::points(medians[1:len], (1:len + offset), pch = 16, col = COL[idx], cex = sz_med)
         # object2
         graphics::matlines(thick_q2[,1:len], (rbind(1:len, 1:len) - offset), type = 'l', 
           lty = 1, lwd = sz_thick, col = COL2[idx])
         graphics::matlines(thin_q2[,1:len], (rbind(1:len, 1:len) - offset), type = 'l', 
           lty = 1, lwd = sz_thin, col = COL2[idx])
         graphics::points(medians2[1:len], (1:len - offset), pch = 16, col = COL2[idx], cex = sz_med)
       }
     }
   } # closes if (horiz == TRUE)
 
####################################################################################################  
# Vertical plot - CI lines parallel to x-axis 

  if (horiz == FALSE) {
    if (missing(ylim)) {
      thin_comb <- c(thin_q, thin_q2)
      rn <- diff(range(thin_comb, na.rm = TRUE)) * PL_SC
      YLIM <- c((min(thin_comb, na.rm = TRUE) - rn), (max(thin_comb, na.rm = TRUE) + rn))
    } else {
      YLIM <- ylim
    }
    if (len > 20) {
      tt <- 1.5
    } else {
      alpha <- 0.667
      beta <- 0.0412
      tt <- alpha + beta * len
    }
    XLIM = c((1 - tt), (len + tt))
    if (missing(ylab)) ylab = 'Parameter Estimate'
    if (is.null(ylab)) ylab = ''
    if (missing(main)) main = ''
    if (missing(labels)) {
      labs = names(medians)
    } else {
      if (!missing(labels)) {
        if (is.null(labels)) {
          labs <- rep('', len)
        }
        if (!is.null(labels)) {
          if (length(labels) == len) {
            labs <- labels[idx]
          } else {
           stop('Labels length not equal to number of parameters')
          }
        }
      }
    } # closes if (missing(labels))
  
    # 0.2 inches per line - mar measured in lines
    m_char <- (max(sapply(labs, function(x){graphics::strwidth(x, cex = sz_labels, units = 'i')}))  / 0.2)
     
    # blank plot - do not display
    grDevices::pdf(file = NULL)
    graphics::plot(
      (len:1), 
      medians, 
      xlim = XLIM, 
      ylim = YLIM, 
      type = "n",
      ann = TRUE, 
      xaxt = 'n', 
      yaxt = "n", 
      bty = "n", 
      ylab = NA,
      xlab = NA, 
      cex.lab = sz_ax_txt, 
      xaxs = 'i')
     
    # create invisible ticks to determine where to put y-axis label
    tickp <- graphics::axis(
      2, 
      lwd.ticks = sz_ax, 
      labels = FALSE,
      at = pos_tick,
      lwd = sz_ax,
      cex.axis = sz_tick_txt,
      las = 1,
      col = 'white',
      col.ticks = 'white')
    invisible(grDevices::dev.off())
     
    # determine how long labels are
    ml_tickp <- max(graphics::strwidth(tickp, cex = sz_tick_txt, units = 'in'))
    # 5 lines/inch
    ll <- 1.8 + 5 * ml_tickp
    # set plot margins according to labels
    graphics::par(mar = c((m_char + (mar[1] - 4)), mar[2] + 1 + ll - 3, mar[3] - 1.5, mar[4]))

    # new blank plot
    graphics::plot(
      (len:1), 
      medians, 
      xlim = XLIM, 
      ylim = YLIM, 
      type = "n",
      ann = TRUE, 
      xaxt = 'n', 
      yaxt = "n", 
      bty = "n", 
      ylab = NA,
      xlab = NA, 
      cex.lab = sz_ax_txt, 
      xaxs = 'i')
     
    # ticks
    graphics::axis(
      2, 
      lwd.ticks = sz_ax, 
      labels = TRUE, 
      at = pos_tick, 
      lwd = sz_ax, 
      cex.axis = sz_tick_txt, 
      las = 1)
     
    # y-axis label
    graphics::title(ylab = ylab, cex.lab = sz_ax_txt, line = ll)
    
    # title
    graphics::title(main, cex.main = sz_main_txt, line = 0.5)

    if (guide_axis == TRUE) {
      # right axis params
      graphics::axis(4, lwd.ticks = sz_ax, labels = FALSE, at = pos_tick, lwd = sz_ax)
     }

    # bottom axis params (labels)
    # las - 0 parallel to axis, 1 horizontal, 2 perpendicular to axis, 3 vertical
    graphics::axis(
      1, 
      at = (len:1) + 0.013, 
      tick = FALSE, 
      labels = labs, 
      las = 2, 
      adj = 0,
      line = -1, 
      cex.axis = sz_labels)

    if (guide_lines == TRUE) {
       # limits
       yaxp <- graphics::par("yaxp")
       YLIM2 <- c(yaxp[1], yaxp[2])
       ys <- matrix(rep(YLIM2, len), nrow = 2)
       xs <- rbind(1:len, 1:len)
       graphics::matlines(xs, ys, lty = 1, col = guide_col)
    }
 
    # ref line
    if (!is.null(ref)) {
      graphics::abline(h=ref, lty = 2, lwd = 3, col = ref_col)
    }
    if (ref_ovl == TRUE) {
      # positions to plot CI
      v_blk_bnd <- matrix(rep(rev(len + 1 - black_cl), 2), nrow = 2, byrow = TRUE)
      v_gry_bnd <- matrix(rep(rev(len + 1 - gray_cl), 2), nrow = 2, byrow = TRUE)
      v_wht_bnd <- matrix(rep(rev(len + 1 - white_cl), 2), nrow = 2, byrow = TRUE)
      # Black CI
      if (!is.null(black_cl)) {
        if (is.null(object2)) {
          # Thick
          graphics::matlines(v_blk_bnd, thick_q[,v_black_cl], type = 'l', lty = 1, lwd = sz_thick, col = 'black')
          # Thin
          graphics::matlines(v_blk_bnd, thin_q[,v_black_cl], type = 'l', lty = 1, lwd = sz_thin, col = 'black')
         } else {
           v_blk_bnd2 <- matrix(rep(rev(len + 1 - black_cl2), 2), nrow = 2, byrow = TRUE)
           # object
           graphics::matlines((v_blk_bnd - offset), thick_q[,v_black_cl], type = 'l', lty = 1, lwd = sz_thick, col = 'black')
           graphics::matlines((v_blk_bnd - offset), thin_q[,v_black_cl], type = 'l', lty = 1, lwd = sz_thin, col = 'black')
           # object2
           graphics::matlines((v_blk_bnd2 + offset), thick_q2[,v_black_cl2], type = 'l', lty = 1, lwd = sz_thick, col = 'black')
           graphics::matlines((v_blk_bnd2 + offset), thin_q2[,v_black_cl2], type = 'l', lty = 1, lwd = sz_thin, col = 'black')
         } # closes if (is.null(object2))
       } # closes if (!is.null(black_cl))

       # Gray CI
       if (!is.null(gray_cl)) {
         if (is.null(object2)) {
          # Thick
          graphics::matlines(v_gry_bnd, thick_q[,v_gray_cl], type = 'l', lty = 1, lwd = sz_thick, col = gr_col)
          # Thin
          graphics::matlines(v_gry_bnd, thin_q[,v_gray_cl], type = 'l', lty = 1, lwd = sz_thin, col = gr_col)
         } else {
           v_gry_bnd2 <- matrix(rep(rev(len + 1 - gray_cl2), 2), nrow = 2, byrow = TRUE)
           # object
           graphics::matlines((v_gry_bnd - offset), thick_q[,v_gray_cl], type = 'l', lty = 1, lwd = sz_thick, col = gr_col)
           graphics::matlines((v_gry_bnd - offset), thin_q[,v_gray_cl], type = 'l', lty = 1, lwd = sz_thin, col = gr_col)
           #object2
           graphics::matlines((v_gry_bnd2 + offset), thick_q2[,v_gray_cl2], type = 'l', lty = 1, lwd = sz_thick, col = gr_col)
           graphics::matlines((v_gry_bnd2 + offset), thin_q2[,v_gray_cl2], type = 'l', lty = 1, lwd = sz_thin, col = gr_col)
         }
       }

       # White CI
       if (!is.null(white_cl)) {
         if (is.null(object2)) {
           graphics::matlines(v_wht_bnd, thick_q[,v_white_cl], type = 'l', lty = 1, lwd = sz_thick, col = gr_col) #white (gray)
           graphics::matlines(v_wht_bnd, thin_q[,v_white_cl], type = 'l', lty = 1, lwd = sz_thin, col = gr_col) #white (gray)
         } else {
           v_wht_bnd2 <- matrix(rep(rev(len + 1 - white_cl2), 2), nrow = 2, byrow = TRUE)
           # object
           graphics::matlines((v_wht_bnd - offset), thick_q[,v_white_cl], type = 'l', lty = 1, lwd = sz_thick, col = gr_col)
           graphics::matlines((v_wht_bnd - offset), thin_q[,v_white_cl], type = 'l', lty = 1, lwd = sz_thin, col = gr_col)
           # object2
           graphics::matlines((v_wht_bnd2 + offset), thick_q2[,v_white_cl2], type = 'l', lty = 1, lwd = sz_thick, col = gr_col)
           graphics::matlines((v_wht_bnd2 + offset), thin_q2[,v_white_cl2], type = 'l', lty = 1, lwd = sz_thin, col = gr_col)
         }
       }

       # Medians
      if (is.null(object2)) {
        graphics::points(len:1, medians, pch = 16, col = 'white', cex = sz_med)
        graphics::points(v_blk_bnd[1,], medians[v_black_cl], pch = 16, col = 'black', cex = sz_med)
        graphics::points(v_gry_bnd[1,], medians[v_gray_cl], pch = 16, col = gr_col, cex = sz_med)
        graphics::points(v_wht_bnd[1,], medians[v_white_cl], pch = 21, col = gr_col, cex = sz_med, lwd = 2)
       } else {
         # object
         graphics::points((len:1 - offset), medians, pch = 16, col = 'white', cex = sz_med)
         graphics::points((v_blk_bnd[1,] - offset), medians[v_black_cl], pch = 16, col = 'black', cex = sz_med)
         graphics::points((v_gry_bnd[1,] - offset), medians[v_gray_cl], pch = 16, col = gr_col, cex = sz_med)
         graphics::points((v_wht_bnd[1,] - offset), medians[v_white_cl], pch = 21, col = gr_col, cex = sz_med, lwd = 2)
         
         # object2
         graphics::points((len:1 + offset), medians2, pch = 16, col = 'white', cex = sz_med)
         graphics::points((v_blk_bnd2[1,] + offset), medians2[v_black_cl2], pch = 16, col = 'black', cex = sz_med)
         graphics::points((v_gry_bnd2[1,] + offset), medians2[v_gray_cl2], pch = 16, col = gr_col, cex = sz_med)
         graphics::points((v_wht_bnd2[1,] + offset), medians2[v_white_cl2], pch = 21, col = gr_col, cex = sz_med, lwd = 2)
       }
    } else {
      if (is.null(object2)) {
        graphics::matlines(rbind(1:len, 1:len), thick_q[,len:1], type = 'l', lty = 1, lwd = sz_thick, col = COL[idx])
        graphics::matlines(rbind(1:len, 1:len), thin_q[,len:1], type = 'l', lty = 1, lwd = sz_thin, col = COL[idx])
        graphics::points(1:len, medians[len:1], pch = 16, col = COL[idx], cex = sz_med)
       } else {
         # object
         graphics::matlines((rbind(1:len, 1:len) - offset), thick_q[,len:1], type = 'l', lty = 1, lwd = sz_thick, col = COL[idx])
         graphics::matlines((rbind(1:len, 1:len) - offset), thin_q[,len:1], type = 'l', lty = 1, lwd = sz_thin, col = COL[idx])
         graphics::points((1:len - offset), medians[len:1], pch = 16, col = COL[idx], cex = sz_med)
         # object2
         graphics::matlines((rbind(1:len, 1:len) + offset), thick_q2[,len:1], type = 'l', lty = 1, lwd = sz_thick, col = COL2[idx])
         graphics::matlines((rbind(1:len, 1:len) + offset), thin_q2[,len:1], type = 'l', lty = 1, lwd = sz_thin, col = COL2[idx])
         graphics::points((1:len + offset), medians2[len:1], pch = 16, col = COL2[idx], cex = sz_med)
       } # closes if (is.null(object2))
     } # closes if (ref_ovl == TRUE)
   } # closes if (horiz == FALSE)
  
  # restore graphics
  graphics::par(mar = c(5, 4, 4, 2) + 0.1)
}
