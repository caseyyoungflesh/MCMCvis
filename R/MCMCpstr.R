#' Summary function for MCMC output that preserves parameter structure
#'
#' Extract summary information from MCMC output (specific function specified) for specific parameters of interest while preserving original parameter structure (i.e., scalar, vector, matrix, array). Function outputs a \code{list} with calculated values for each specified parameter.
#'
#' @param object Object containing MCMC output. See DETAILS below.
#' @param params Character string (or vector of character strings) denoting parameters to be returned in summary output.
#'
#' Default \code{'all'} returns all parameters in summary output.
#'
#' @param excl Character string (or vector of character strings) denoting parameters to exclude. Used in conjunction with \code{params} argument to select parameters of interest.
#'
#' @param ISB Ignore Square Brackets (ISB). Logical specifying whether square brackets should be ignored in the \code{params} and \code{excl} arguments. If \code{TRUE}, square brackets are ignored - input from \code{params} and \code{excl} are otherwise matched exactly. If \code{FALSE}, square brackets are not ignored - input from \code{params} and \code{excl} are matched using grep, which can take arguments in regular expression format. This allows partial names to be used when specifying parameters of interest.
#'
#' @param digits Number of digits to include for posterior summary. Values will be rounded to the specified number of digits (except for Rhat which is always rounded to 2 digits).
#'
#' Default is \code{digits = 2}.
#'
#' @param func Function to be performed on MCMC output. Output of specified function must be of length 1.
#'
#' @section Details:
#' \code{object} argument can be a \code{stanfit} object (\code{rstan} package), an \code{mcmc.list} object (\code{coda} package), an \code{R2jags} model object (\code{R2jags} package), a \code{jagsUI} model object (\code{jagsUI} package), or a matrix containing MCMC chains (each column representing MCMC output for a single parameter, rows representing iterations in the chain). The function automatically detects the object type and proceeds accordingly.
#'
#' @examples
#' #Load data
#' data(MCMC_data)
#'
#' MCMCpstr(MCMC_data, func = function(x) quantile(x, probs = 0.01))
#'
#' @export

MCMCpstr <- function(object,
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
    np <- NCOL(object2)
    ch_bind <- object2

    #how many elements will be in the list
    un <- unique(names)
    onames <- colnames(temp_in)
  }


  #create empty list
  out_list <- vector('list', length(un))

  #iterate through each element
  for (i in 1:length(un))
  {
    #i <- 1
    ind <- which(un[i] == names)

    #determine how many ',' and therefore how many dimensions for parameter
    dimensions <- length(strsplit(onames[ind[1]], split = ',', fixed = TRUE)[[1]])

    if (length(func(ch_bind[,ind[1]])) > 1)
    {
      stop("Output from 'func' argument must be of length 1.")
    }

    #scalar or vector
    if(dimensions == 1)
    {
      #1 dimension
      temp_obj <- rep(NA, length(ind))
      for (j in 1:length(ind))
      {
        #j <- 1
        temp_obj[j] <- round(func(ch_bind[,ind[j]]), digits = digits)
      }

      #fill list
      out_list[[i]] <- temp_obj
    }
    if(dimensions == 2)
    {
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
    if(dimensions == 3)
    {
      #3 dimensions
      pnames <- vapply(strsplit(onames[ind],
                                '[', fixed = TRUE), `[`, 2, FUN.VALUE=character(1))
      pnames2 <- vapply(strsplit(pnames,
                                 ']', fixed = TRUE), `[`, 1, FUN.VALUE=character(1))

      #determine how large array should be
      RI <- c()
      CI <- c()
      TI <- c()
      for (j in 1:length(ind))
      {
        #j <- 1
        tt <- strsplit(pnames2[j], ',')
        ri <- as.numeric(tt[[1]][1])
        ci <- as.numeric(tt[[1]][2])
        ti <- as.numeric(tt[[1]][3])
        RI <- c(RI, ri)
        CI <- c(CI, ci)
        TI <- c(TI, ti)
      }

      #create blank array
      temp_obj <- array(NA, dim = c(max(RI), max(CI), max(TI)))

      for (j in 1:length(ind))
      {
        temp_obj[RI[j], CI[j], TI[j]] <- round(func(ch_bind[,ind[j]]), digits = digits)
      }

      #fill list
      out_list[[i]] <- temp_obj
    }
    if(dimensions == 4)
    {
      #4 dimensions
      pnames <- vapply(strsplit(onames[ind],
                                '[', fixed = TRUE), `[`, 2, FUN.VALUE=character(1))
      pnames2 <- vapply(strsplit(pnames,
                                 ']', fixed = TRUE), `[`, 1, FUN.VALUE=character(1))

      #determine how large matrix should be
      RI <- c()
      CI <- c()
      TI <- c()
      FI <- c()
      for (j in 1:length(ind))
      {
        #j <- 1
        tt <- strsplit(pnames2[j], ',')
        ri <- as.numeric(tt[[1]][1])
        ci <- as.numeric(tt[[1]][2])
        ti <- as.numeric(tt[[1]][3])
        fi <- as.numeric(tt[[1]][4])
        RI <- c(RI, ri)
        CI <- c(CI, ci)
        TI <- c(TI, ti)
        FI <- c(FI, fi)
      }

      #create blank array
      temp_obj <- array(NA, dim = c(max(RI), max(CI), max(TI), max(FI)))

      for (j in 1:length(ind))
      {
        temp_obj[RI[j], CI[j], TI[j], FI[j]] <- round(func(ch_bind[,ind[j]]), digits = digits)
      }

      #fill list
      out_list[[i]] <- temp_obj
    }
    if(dimensions > 4)
    {
      stop('This function does not currently support parameters with > 4 dimensions. If you have a need for this functionality, please put in a bug report on the package Github page.')
    }
  }

  #name elements in list
  names(out_list) <- un

  return(out_list)
}

