#' Chain extraction for jags object
#'
#' Simplifies extracting particular chains from jags object
#'
#' @param jags_object Output from jags model
#' @param params Parameters of interest from jags object. Partial names
#' may be used to return all object containing that name.
#' @return Chains for particular parameters of interest
#' @section Details:
#' Output from R2jags only right now.
#'
#' @export


pochains <- function(jags_object, params = 'all')
{
  temp <- jags_object$BUGSoutput$sims.matrix

  if (params == 'all')
  {
    return(temp)
  }else
  {
    get.cols <- grep(paste(params), colnames(temp))
    return(temp[,get.cols])
  }
}
