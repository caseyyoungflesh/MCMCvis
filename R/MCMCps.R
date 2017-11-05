#' Summary function for MCMC output preserving parameter structure
#'
#' Extract summary information from MCMC output (specific function specified) for specific parameters of interest while preserving original parameter structure (e.g., scalar, vector, matrix).
#'
#' @param object Object containing MCMC output. See DETAILS below.
#' @param params Character string (or vector of character strings) denoting parameters to be returned in summary output.
#'
#' Default \code{'all'} returns all parameters in summary output.
#'
#' @param excl Character string (or vector of character strings) denoting parameters to exclude. Used in conjunction with \code{params} argument to select parameters of interest.
#'
#' @param ISB Ignore Square Brackets (ISB). Logical specifying whether square brackets should be ignored in the \code{params} and \code{excl} arguments. If \code{FALSE}, square brackets are ignored - input from \code{params} and \code{excl} are otherwise matched exactly. If \code{TRUE}, square brackets are not ignored - input from \code{params} and \code{excl} are matched using grep, allowing partial names to be used when specifying parameters of interest.
#'
#' @param digits Number of digits to include for posterior summary. Values will be rounded to the specified number of digits (except for Rhat which is always rounded to 2 digits).
#'
#' Default is \code{digits = 2}.
#'
#' @param func Function to be performed on MCMC output.
#'
#' @section Details:
#' \code{object} argument can be a \code{stanfit} object (\code{rstan} package), an \code{mcmc.list} object (\code{coda} package), an \code{R2jags} model object (\code{R2jags} package), or a matrix containing MCMC chains (each column representing MCMC output for a single parameter, rows representing iterations in the chain). The function automatically detects the object type and proceeds accordingly.
#'
#' @examples
#' #Load data
#' data(MCMC_data)
#'
#' XXXX
#'
#' @export


MCMCps <- function(object,
                   params = 'all',
                   excl = NULL,
                   ISB = TRUE,
                   digits = 2,
                   func = mean)
{
  #SORTING BLOCK
  if(typeof(object) == 'double')
  {
    object2 <- MCMCchains(object, params, excl, ISB, mcmc.list = FALSE)
  }else{
    object2 <- MCMCchains(object, params, excl, ISB, mcmc.list = TRUE)
  }

  if(coda::is.mcmc.list(object2) == TRUE)
  {
    temp_in <- object2
    if(ISB == TRUE)
    {
      names <- vapply(strsplit(colnames(temp_in[[1]]),
                               split = "[", fixed = TRUE), `[`, 1, FUN.VALUE=character(1))
    }else{
      names <- colnames(temp_in[[1]])
    }
  }

  if(typeof(object2) == 'double')
  {
    temp_in <- object2
    if(ISB == TRUE)
    {
      names <- vapply(strsplit(colnames(temp_in),
                               split = "[", fixed = TRUE), `[`, 1, FUN.VALUE=character(1))
    }else{
      names <- colnames(temp_in)
    }
  }

  np <- NCOL(object2[[1]])

  if(np > 1)
  {
    ch_bind <- do.call('rbind', object2)
  }else{
    ch_bind <- as.matrix(object2)
  }

  #how many elements will be in the list
  un <- unique(names)
  onames <- colnames(temp_in[[1]])

  #create empty list
  out_list <- vector('list', length(un))

  #iterate through each element
  for (i in 1:length(un))
  {
    #i <- 3
    ind <- which(un[i] == names)

    #if there's a , then it's a matrix
    if(length(grep(',', onames[ind[1]])) == 0)
    {
      #1 dimension
      temp_obj <- rep(NA, length(ind))
      for (j in 1:length(ind))
      {
        #j <- 1
        temp_obj[j] <- func(ch_bind[,ind[j]])
      }

      #fill list
      out_list[[i]] <- temp_obj
    }else{
      #2 dimensions
      pnames <- vapply(strsplit(onames[ind],
                                '[', fixed = TRUE), `[`, 2, FUN.VALUE=character(1))
      pnames2 <- vapply(strsplit(pnames,
                                ']', fixed = TRUE), `[`, 1, FUN.VALUE=character(1))

      #determine how large matrix should be
      RI <- c()
      CI <- c()
      for (j in 1:length(ind))
      {
        #j <- 1
        tt <- strsplit(pnames2[j], ',')
        ri <- as.numeric(tt[[1]][1])
        ci <- as.numeric(tt[[1]][2])
        RI <- c(RI, ri)
        CI <- c(CI, ci)
      }

      #create blank matrix
      temp_obj <- matrix(NA, nrow = max(RI), ncol = max(CI))

      for (j in 1:length(ind))
      {
        temp_obj[RI[j], CI[j]] <- round(func(ch_bind[,ind[j]]), digits = digits)
      }

      #fill list
      out_list[[i]] <- temp_obj
    }
  }

  #name elements in list
  names(out_list) <- un

  return(out_list)
}


