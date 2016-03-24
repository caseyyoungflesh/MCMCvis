#' Summary function for jags object
#'
#' Simplifies extracting info for particular parameters from jags object
#'
#' @param jags_object Output from jags model
#' @param params Parameters of interest from jags object. Partial names
#' may be used to return all object containing that name.
#' @return Summary data (mean, 95\% quantile, R_hat) for particular parameters of interest
#' @section Details:
#' Output from R2jags only right now.
#'
#' @export


posummary <- function(jags_object, params = 'all')
{
  temp <- jags_object$BUGSoutput$summary

  if (params == 'all')
  {
    print(temp[,c(1,3,7,8)])
  }else
  {
    get.rows <- grep(paste(params), rownames(temp))
    print(temp[get.rows,c(1,3,7,8)])
  }
}
