#' The 'potools' package
#'
#' The potools package assists with manipulating and summarizing MCMC output
#' (including parameter estimate extraction, chain extraction, and
#' visualizing posterior distributions). MCMC output may be derived from
#' Bayesian model ouput fit wit JAGS or other MCMC samplers.
#'
#' @section Details:
#' The following functions are currently available:
#'
#'      -\code{posummary} (summarize MCMC output for particular parameters of interest)
#'
#'      -\code{potrace} (plots chains for particular parameters of interest to check for convergence)
#'
#'      -\code{pochains} (extract posterior chains from MCMC output for particular parameters of interest)
#'
#'      -\code{poplot} (plot posterior distributions from MCMC output for paticular parameters of interest)
#'
#' Example data for a \code{mcmc.list} object can be loaded using \code{data(MCMC_data)}.
#'
#' These tools make it possible to input the MCMC output object directly into the function, without significant
#' prior manipulation. Summary information, trace plots, posterior chains, and posterior disritbution plots can be
#' easily obtained for any parameter of interest within each function. Each function is designed to operate independently
#' of one another.
#'
#' The vignette can be run using \code{vignette('potools')} if vignette is built.
#'
#' @section Author(s):
#' Casey Youngflesh <casey.youngflesh@stonybrook.edu>
#'
#' @docType package
#' @name potools

NULL

