#' The 'post' package
#'
#' The post package assists with summarizing posterior MCMC chains
#' (including parameter estimate extraction, chain extraction, and
#' visualizing posterior distributions), obtained from Bayesian model
#' ouput.
#'
#' @section Details:
#' The following functions are currently available:
#'
#'      -posummary (summarizes jags model output)
#'
#'      -pochains (extracts chains from jags model output)
#'
#'      -poplot (plots posterior chains from jags model output)
#'
#' @section Workflow:
#' Workflow of package as follows:
#'
#'      -posummary(jags_object, params)
#'
#'      -chains <- pochains(jags_object, params)
#'
#'      -poplot(chains)
#'
#'
#' @section Author(s):
#' Casey Youngflesh <casey.youngflesh@stonybrook.edu>
#'
#' @docType package
#' @name post

NULL

