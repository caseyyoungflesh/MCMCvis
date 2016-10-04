#' The 'potools' package
#'
#' The potools package assists with manipulating and summarizing MCMC output
#' (including parameter estimate extraction, chain extraction, and
#' visualizing posterior distributions). MCMC output may be derived from
#' Bayesian model ouput fit with Stan, JAGS, or other MCMC samplers.
#'
#' @section Details:
#' The following functions are currently available:
#'
#'      -\code{posummary} (summarize MCMC output for particular parameters of interest)
#'
#'      -\code{potrace} (create trace and density plots from MCMC chains for particular parameters of interest)
#'
#'      -\code{pochains} (easily extract posterior chains from MCMC output for particular parameters of interest)
#'
#'      -\code{poplot} (create caterpillar plots from MCMC output for particular parameters of interest)
#'
#' Example data for a \code{mcmc.list} object can be loaded using \code{data(MCMC_data)}.
#'
#' These tools make it possible to input the MCMC output object directly into the function, without significant
#' prior manipulation. Summary information, trace plots, posterior chains, and posterior disritbution plots can be
#' easily obtained for any parameter of interest within each function. Each function is designed to operate independently
#' of one another.
#'
#' The vignette can be run using \code{vignette('potools')} if vignette is built when installing package.
#'
#' @section Author(s):
#' Casey Youngflesh <caseyyoungflesh@gmail.com>
#'
#' @docType package
#' @name potools

NULL

