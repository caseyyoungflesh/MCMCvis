#' Plots posterior distributions
#'
#' Function plots the posterior distributions obtained from MCMC
#' chains output from a Bayesian model.
#'
#' @param data Matrix of MCMC chains
#' @param LML Lower bounds x-axis
#' @param LMH upper bounds x-axis
#' @section Details:
#' Chains should be in columns of input matrix
#' @section Notes:
#' Wrapper for denstrip package, as highlighted by Jackson (2008). Example can be seen in
#' Zipkin et al. 2014, figure 3.
#' @section References:
#' Jackson, C. H. 2008. Displaying Uncertainty With Shading. The American Statistician 62:340-347.
#'
#' Zipkin, E. F., T. S. Sillett, E. H. C. Grant, R. B. Chandler,
#' and J. A. Royle. 2014. Inferences about population dynamics
#' from count data using multistate models: a comparison to
#' capture-recapture approaches. Ecology and Evolution 4:417-426.
#' @examples
#' x1 <- rnorm(1000, mean=0.5)
#' x2 <- rnorm(1000, mean=0)
#' data <- cbind(x1, x2)
#' poplot(data)
#'
#' @export


poplot <- function(data, LML = -1, LMH= 1) #WID= height of density strips, LML, LMH = x-axis bounds
{
  #for debug
  #-------------#

  #LML = -1
  #LMH = 1
  #x1 <- rnorm(1000, mean=.5)
  #x2 <- rnorm(1000, mean=0)
  #data <- cbind(x1,x2)
  #-------------#


  #plotting parameters
  WID <- .2 #height of bar (don't know if this matters)
  W <- 3 #height of mean tick
  W2 <- 3 #height of 95% CI tick
  CI_col <- "grey87" #color of CI tick


  #data
  chains <- as.data.frame(data)
  X <- NCOL(data)
  idx <- X:1
  labs <- colnames(data)[idx] #apply and sort labels
  mp <- suppressMessages(reshape2::melt(chains[,idx], value.name='value')) #melt
  qdata <- apply(chains, 2, quantile, probs=c(.025,.975)) #quantiles of chains


  #create plot object
  rpp <- bwplot(variable~value,data=mp,
               xlab=list(label="Parameter probability values",cex=1.3),
               xlim=c(LML, LMH),
               panel = function(x, y)
               {
                 #grid.segments(1,0,0,0)
                 xlist <- split(mp$value, factor(mp$variable))
                 xlist <- split(x, factor(y))

                 for (i in seq(along=xlist))
                 {
                   #panel.grid(h=c(0), col='grey')
                   denstrip::panel.denstrip(x=xlist[[i]], at=i, width= WID, colmax='black', colmin= 'white')

                 }
               }, par.settings = list(axis.line = list(col=NA)),
               scales=list(col=1,cex=1,x=list(col=1),
                           y=list(draw=T,labels=labs)))

  #print plot object
  print(rpp)

  #top and bottom lines
  lattice::trellis.focus()
  lattice::panel.lines(c(LML,LMH),c(0.4,0.4), col="black")
  lattice::panel.lines(c(LML,LMH),c(X+.6,X+.6), col="black")


  #mean hashes
  for (i in 1:X)
  {
    TMP<- X-i+1
    lattice::panel.lines(c(median(chains[,i])), c(TMP-.25, TMP+.25), lwd=W, col='black')
  }


  #Lower CI
  for(i in 1:X)
  {
    TMP<- X-i+1
    #i <- 6
    lattice::panel.lines(c(rep(qdata[1,i],2)), c(TMP-.25, TMP+.25), col= CI_col, lwd=W2)
  }


  #Upper CI
  for(i in 1:X)
  {
    TMP<- X-i+1
    #i <- 6
    lattice::panel.lines(c(rep(qdata[2,i],2)), c(TMP-.25, TMP+.25), col= CI_col, lwd=W2)
  }


  lattice::trellis.unfocus()
}
