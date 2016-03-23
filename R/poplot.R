#' Plots posterior distributions
#'
#' @param data Matrix of MCMC chains
#' @param WID Height?
#' @param LML Lower bounds x-axis
#' @param LMH upper bounds x-axis
#' @section Details:
#' Chains should be in columns of input matrix
#'
#' @export


poplot <- function(data, WID=.2, LML = -1, LMH= 1) #WID= height of density strips, LML, LMH = x-axis bounds
{
  #for debug
  #WID=.2
  #LML = -1
  #LMH = 1
  #data <- OR_chains_beta

  chains <- as.data.frame(data)

  X <- NCOL(data)
  idx<-X:1

  # apply and sort labels
  labs = colnames(data)[idx]

  mp= reshape2::melt(chains[,idx])

  #x11()
  rpp = lattice::bwplot(variable~value,data=mp,
               xlab=list(label="Parameter probability values",cex=1.3),
               xlim=c(LML, LMH),
               panel = function(x, y)
               {
                 #grid.segments(1,0,0,0)
                 xlist <- split(x, factor(y))
                 for (i in seq(along=xlist))
                 {
                   #panel.grid(h=c(0), col='grey')
                   denstrip::panel.denstrip(x=xlist[[i]], at=i, width= WID)#, colmax='white', colmin= 'white')

                 }
               }, par.settings = list(axis.line = list(col=NA)),
               scales=list(col=1,cex=1,x=list(col=1),
                           y=list(draw=T,labels=labs)))

  print(rpp)

  lattice::trellis.focus()
  #top and bottom lines
  lattice::panel.lines(c(LML,LMH),c(0.4,0.4), col="black")
  lattice::panel.lines(c(LML,LMH),c(X+.6,X+.6), col="black")
  #trellis.unfocus()

  #width mean line
  W <- 3
  #width 95% CI lines
  W2 <- 3

  #mean hashes
  for (i in 1:X)
  {
    TMP<- X-i+1
    lattice::panel.lines(c(median(chains[,i])), c(TMP-.25, TMP+.25), lwd=W, col='black')
  }

  qdata <- apply(chains, 2, quantile, probs=c(.025,.975))

  #Lower CI
  for(i in 1:X)
  {
    TMP<- X-i+1
    #i <- 6
    lattice::panel.lines(c(rep(qdata[1,i],2)), c(TMP-.25, TMP+.25), col="grey87", lwd=W2)
  }

  #Upper CI
  for(i in 1:X)
  {
    TMP<- X-i+1
    #i <- 6
    lattice::panel.lines(c(rep(qdata[2,i],2)), c(TMP-.25, TMP+.25), col="grey87", lwd=W2)
  }


  lattice::trellis.unfocus()
}
