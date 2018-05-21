#' Summarize and extract posterior chains from MCMC output while preserving parameter structure
#'
#' Extract summary information and posterior chains from MCMC output (specific function specified) for specific parameters of interest while preserving original parameter structure (i.e., scalar, vector, matrix, array). Function outputs a \code{list} with calculated values or posterior chains for each specified parameter.
#'
#' @param object Object containing MCMC output. See DETAILS below.
#' @param params Character string (or vector of character strings) denoting parameters to be returned in output.
#'
#' Default \code{'all'} returns all parameters in output.
#'
#' @param excl Character string (or vector of character strings) denoting parameters to exclude. Used in conjunction with \code{params} argument to select parameters of interest.
#'
#' @param ISB Ignore Square Brackets (ISB). Logical specifying whether square brackets should be ignored in the \code{params} and \code{excl} arguments. If \code{TRUE}, square brackets are ignored - input from \code{params} and \code{excl} are otherwise matched exactly. If \code{FALSE}, square brackets are not ignored - input from \code{params} and \code{excl} are matched using grep, which can take arguments in regular expression format. This allows partial names to be used when specifying parameters of interest.
#'
#' @param func Function to be performed on MCMC output. When output of specified function is greater than length 1, an extra dimension is added. For instance, output of length 3 for a parameter with dimensions 2x2 results in a 2x2x3 output. Functions that produce output with dimensionality greater than 1 are not permitted. \code{func} is ignored when \code{type = 'chains'}.
#'
#' @param type Character string specifying whether to return summary information (calculated based on \code{func} argument) or posterior chains. Valid options are \code{'summary'} and \code{'chains'}. When \code{type = 'chains'}, the \code{'func'} argument is ignored. When \code{type = 'chains'}, posterior chains are concatenated and stored in the last dimension in the array for each element (parameter) of the list.
#'
#' @section Details:
#' \code{object} argument can be a \code{stanfit} object (\code{rstan} package), an \code{mcmc.list} object (\code{coda} package), an \code{R2jags} model object (\code{R2jags} package), a \code{jagsUI} model object (\code{jagsUI} package), or a matrix containing MCMC chains (each column representing MCMC output for a single parameter, rows representing iterations in the chain). The function automatically detects the object type and proceeds accordingly.
#'
#' @examples
#' #Load data
#' data(MCMC_data)
#'
#' MCMCpstr(MCMC_data, func = function(x) quantile(x, probs = c(0.01, 0.99)))
#'
#' @export

MCMCpstr <- function(object,
                   params = 'all',
                   excl = NULL,
                   ISB = TRUE,
                   func = mean,
                   type = 'summary')
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
    dims <- length(strsplit(onames[ind[1]], split = ',', fixed = TRUE)[[1]])

    #scalar or vector
    if (dims == 1)
    {
      if (type == 'summary')
      {
        func_out <- func(ch_bind[,ind[1]])
        if (length(func_out) > 1)
        {
          if (!is.null(dim(func_out)))
          {
            stop("Output from 'func' argument must be of dimension 1.")
          }
          temp_obj <- matrix(NA, nrow = length(ind), ncol = length(func_out))
          dimnames(temp_obj)[[2]] <- names(func_out)
          dimnames(temp_obj)[[1]] <- onames[ind]

          for (j in 1:length(ind))
          {
            temp_obj[j,] <- func(ch_bind[,ind[j]])
          }
        } else{
          temp_obj <- rep(NA, length(ind))
          for (j in 1:length(ind))
          {
            temp_obj[j] <- func(ch_bind[,ind[j]])
          }
        }
      }
      if (type == 'chains')
      {
        temp_obj <- matrix(NA, nrow = length(ind), ncol = NROW(ch_bind))
        dimnames(temp_obj)[[1]] <- onames[ind]

        for (j in 1:length(ind))
        {
          temp_obj[j,] <- ch_bind[,ind[j]]
        }
      }
      if (type != 'summary' & type != 'chains')
      {
        stop("Invalid input for argument 'type'. Valid options are 'summary' and 'chains'.")
      }

      #fill list
      out_list[[i]] <- temp_obj
    }
    #2 dimensions
    if (dims == 2)
    {
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

      if (type == 'summary')
      {
        func_out <- func(ch_bind[,ind[1]])
        if (length(func_out) > 1)
        {
          if (!is.null(dim(func_out)))
          {
            stop("Output from 'func' argument must be of dimension 1.")
          }
          #create blank array
          temp_obj <- array(NA, dim = c(max(RI), max(CI), length(func_out)))
          dimnames(temp_obj)[[3]] <- names(func_out)

          for (j in 1:length(ind))
          {
            #j <- 1
            temp_obj[RI[j], CI[j], ] <- func(ch_bind[,ind[j]])
          }
        } else {
          #create blank matrix
          temp_obj <- matrix(NA, nrow = max(RI), ncol = max(CI))

          for (j in 1:length(ind))
          {
            temp_obj[RI[j], CI[j]] <- func(ch_bind[,ind[j]])
          }
        }
      }
      if (type == 'chains')
      {
        #create blank array
        temp_obj <- array(NA, dim = c(max(RI), max(CI), NROW(ch_bind)))

        for (j in 1:length(ind))
        {
          temp_obj[RI[j], CI[j], ] <- ch_bind[,ind[j]]
        }
      }
      if (type != 'summary' & type != 'chains')
      {
        stop("Invalid input for argument 'type'. Valid options are 'summary' and 'chains'.")
      }

      #fill list
      out_list[[i]] <- temp_obj
    }
    #3 dimensions
    if (dims == 3)
    {
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

      if (type == 'summary')
      {
        func_out <- func(ch_bind[,ind[1]])
        if (length(func_out) > 1)
        {
          if (!is.null(dim(func_out)))
          {
            stop("Output from 'func' argument must be of dimension 1.")
          }
          #create blank array
          temp_obj <- array(NA, dim = c(max(RI), max(CI), max(TI), length(func_out)))
          dimnames(temp_obj)[[4]] <- names(func_out)

          for (j in 1:length(ind))
          {
            #j <- 1
            temp_obj[RI[j], CI[j], TI[j], ] <- func(ch_bind[,ind[j]])
          }
        } else {
          #create blank array
          temp_obj <- array(NA, dim = c(max(RI), max(CI), max(TI)))

          for (j in 1:length(ind))
          {
            temp_obj[RI[j], CI[j], TI[j]] <- func(ch_bind[,ind[j]])
          }
        }
      }
      if (type == 'chains')
      {
        #create blank array
        temp_obj <- array(NA, dim = c(max(RI), max(CI), max(TI), NROW(ch_bind)))

        for (j in 1:length(ind))
        {
          temp_obj[RI[j], CI[j], TI[j], ] <- ch_bind[,ind[j]]
        }
      }
      if (type != 'summary' & type != 'chains')
      {
        stop("Invalid input for argument 'type'. Valid options are 'summary' and 'chains'.")
      }

      #fill list
      out_list[[i]] <- temp_obj
    }
    #4 dimensions
    if (dims == 4)
    {
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

      if (type == 'summary')
      {
        func_out <- func(ch_bind[,ind[1]])
        if (length(func_out) > 1)
        {
          if (!is.null(dim(func_out)))
          {
            stop("Output from 'func' argument must be of dimension 1.")
          }
          #create blank array
          temp_obj <- array(NA, dim = c(max(RI), max(CI), max(TI), max(FI), length(func_out)))
          dimnames(temp_obj)[[5]] <- names(func_out)

          for (j in 1:length(ind))
          {
            #j <- 1
            temp_obj[RI[j], CI[j], TI[j], FI[j], ] <- func(ch_bind[,ind[j]])
          }
        } else {
          #create blank array
          temp_obj <- array(NA, dim = c(max(RI), max(CI), max(TI), max(FI)))

          for (j in 1:length(ind))
          {
            temp_obj[RI[j], CI[j], TI[j], FI[j]] <- func(ch_bind[,ind[j]])
          }
        }
      }
      if (type == 'chains')
      {
        #create blank array
        temp_obj <- array(NA, dim = c(max(RI), max(CI), max(TI), max(FI), NROW(ch_bind)))

        for (j in 1:length(ind))
        {
          temp_obj[RI[j], CI[j], TI[j], FI[j], ] <- ch_bind[,ind[j]]
        }
      }
      if (type != 'summary' & type != 'chains')
      {
        stop("Invalid input for argument 'type'. Valid options are 'summary' and 'chains'.")
      }

      #fill list
      out_list[[i]] <- temp_obj
    }
    if(dims > 4)
    {
      stop('This function does not currently support parameters with > 4 dimensions. If you have a need for this functionality, please put in a bug report on the package Github page.')
    }
  }

  #name elements in list
  names(out_list) <- un

  return(out_list)
}

