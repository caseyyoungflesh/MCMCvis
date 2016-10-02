#' Plot MCMC chains to check for convergence
#'
#' Plot MCMC chains for specific parameters of interest. Option to print trace plots
#' to pdf output.
#'
#' @param object Object containing MCMC output. See DETAILS below.
#' @param params Character string (or vector of character strings) denoting parameters of interest.
#' Partial names may be used to return all parameters containing that set of characters.
#'
#' Default \code{all} returns chains for all parameters.
#'
#' @param pdf Logical - if \code{pdf = TRUE} plots will be exported to a pdf.
#' @param wd Working directory for pdf output. Default is current directory.
#' @section Details:
#'
#' Plots created similar to that of \code{traceplot} from the \code{coda} package.
#'
#' \code{object} argument can be an \code{mcmc.list} object or an \code{R2jags} model object (output from the \code{R2jags}
#' package).
#'
#'
#' @examples
#' #Load data
#' data(MCMC_data)
#'
#' #Traceplot for all 'beta' parameters
#' potrace(MCMC_data, params='beta')
#'
#' #Print traceplot output to pdf
#' potrace(MCMC_data, pdf= TRUE)
#'
#' @export


potrace(MCMC_data,
        pdf = TRUE,
        filename = 'TEST')


potrace <- function(object,
                    params = 'all',
                    pdf = FALSE,
                    filename,
                    wd = getwd())
{
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
  if(coda::is.mcmc.list(object) == FALSE & typeof(object) != 'list')
  {
    stop('Invalid input type. Object must be of type mcmc.list or R2jags output.')
  }

  names <- colnames(temp[[1]])
  n_chains <- length(temp)
  it <- 1:nrow(temp[[1]])

  if (length(params) == 1)
  {
    if (params == 'all')
    {
      g_filt <- 1:length(names)
    }else
    {
      g_filt <- grep(paste(params), names, fixed=TRUE)
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
  }

  if(pdf == TRUE)
  {
    setwd(wd)
    if(missing(filename))
    {
      file_out <- 'potrace.pdf'
    }else{
      file_out <- paste0(filename,'.pdf')
    }
    pdf(file= file_out)
  }

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


  for (j in 1: length(g_filt))
  {
    #chains
    tmlt <- do.call('cbind', temp[,g_filt[j]])
    matplot(it, tmlt, lwd = 1, lty= 1, type='l', main = paste0('Trace - ', names[g_filt[j]]),
            col= rgb(red= gg_cols[1,], green= gg_cols[2,],
                     blue= gg_cols[3,], alpha = 0.6),
            xlab= 'Iteration', ylab= 'Value')
    #density plot
    plot(density(rbind(tmlt)), xlab = 'Parameter estimate',
         lty = 1, lwd = 1, main = paste0('Density - ', names[g_filt[j]]))
  }

  if(pdf == TRUE)
  {
    invisible(dev.off())
    system(paste0('open ', paste0('"', file_out, '"')))
  }
}
