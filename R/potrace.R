#' Plot MCMC chains to check for convergence
#'
#' Plot MCMC chains for specific parameters of interest.
#' @param object Object containing MCMC output. See DETAILS below.
#' @param params Character string (or vector of character strings) denoting parameters of interest.
#' Partial names may be used to return all parameters containing that set of characters.
#'
#' Default \code{all} returns chains for all parameters.

#' @section Details:
#' \code{object} argument can be an \code{mcmc.list} object, an \code{R2jags} model object (output from the \code{R2jags}
#' package), or a matrix containing MCMC chains (each column representing MCMC output for a single parameter, rows
#' representing iterations in the chain).
#'
#' @return \code{potrace(params='all')} returns chains for all parameters.
#'
#' \code{potrace(params=c('beta[1]', 'beta[2]'))} returns chains for just parameters \code{beta[1]} and \code{beta[2]}.
#'
#' \code{potrace(params=c('beta'))} returns chains for all parameters containing \code{beta} in their name.
#'
#' @import ggplot2
#' @export

potrace <- function(object,
                    params = 'all')
{
  if(coda::is.mcmc.list(object) == TRUE)
  {
    temp <- object
  }
  if(coda::is.mcmc.list(object) == FALSE & typeof(object) == 'list')
  {
    temp <- as.mcmc(object)
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

    g_out <- NULL
    for (j in 1: length(g_filt))
    {
      tmlt <- do.call('cbind', temp[,g_filt[j]])

      tplt <- reshape2::melt(data.frame(Iteration= it,tmlt), id='Iteration')
      g <- with(tplt, ggplot2::ggplot(tplt, aes(Iteration, value, color= variable)) +
        geom_line(alpha=.6) +
        theme_bw() +
        ggtitle(paste0(names[g_filt[j]])) +
        guides(color=FALSE))

      g_out[[j]] <- g
    }

  lgf <- length(g_filt)
  for (i in seq(1,lgf, by=6))
  {
    if (lgf-i >= 5)
    {
      gridExtra::grid.arrange(g_out[[i]], g_out[[i+1]],
                              g_out[[i+2]], g_out[[i+3]],
                              g_out[[i+4]], g_out[[i+5]],
                              nrow=3, ncol=2)
    }
    if (lgf-i >=4 & lgf-i < 5)
    {
      gridExtra::grid.arrange(g_out[[i]], g_out[[i+1]],
                              g_out[[i+2]], g_out[[i+3]],
                              g_out[[i+4]], nrow=3, ncol=2)
    }
    if (lgf-i >=3 & lgf-i < 4)
    {
      gridExtra::grid.arrange(g_out[[i]], g_out[[i+1]],
                              g_out[[i+2]], g_out[[i+3]],
                              nrow=3, ncol=2)
    }
    if (lgf-i >=2 & lgf-i < 3)
    {
      gridExtra::grid.arrange(g_out[[i]], g_out[[i+1]],
                              g_out[[i+2]], nrow=3, ncol=2)
    }
    if (lgf-i >=1 & lgf-i < 2)
    {
      gridExtra::grid.arrange(g_out[[i]], g_out[[i+1]],
                              nrow=3, ncol=2)
    }
    if (lgf-i < 1)
    {
      gridExtra::grid.arrange(g_out[[i]], nrow=3, ncol=2)
    }
  }
}
