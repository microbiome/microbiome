#' @title plot_bagged_rda
#' @description Bagged RDA visualization
#'
#' @param x Output from bagged_rda
#' @param which.bac TBA
#' @param ptype TBA
#' @param comp TBA
#' @param cex.bac TBA
#' @param plot.names Plot names
#' @param group.cols TBA
#' @param ... Other arguments to be passed
#'   
#' @return TBA
#'
#' @examples \dontrun{
#'   library(microbiome)
#'   data(peerj32)
#'   x <- t(peerj32$microbes)
#'   y <- factor(peerj32$meta$time); names(y) <- rownames(peerj32$meta)
#'   res <- bagged_rda(x, y, sig.thresh=0.05, nboot=100)
#'   plot_bagged_rda(res)
#'  }
#' @export
#' @importFrom ade4 s.class
#'
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
plot_bagged_rda <- function(x, which.bac = 1:nrow(Bag.res$loadings),
	           ptype="spider", comp=1:2, cex.bac=0.5, plot.names=T,
		   group.cols = as.numeric(unique(Y)),...){

  Bag.res <- x$bagged.rda
  Y <- x$variable

  scaled.loadings <- (Bag.res$loadings/max(abs(Bag.res$loadings)))[,comp]
  scaled.scores <- (Bag.res$scores/max(abs(Bag.res$scores)))[,comp]

  plot(rbind(scaled.scores,scaled.loadings),type="n",xlab=paste(names(Bag.res$R2)[1]," (",format(100*Bag.res$R2[1],digits=2),"%)",sep=""),ylab=paste(names(Bag.res$R2)[2]," (",format(100*Bag.res$R2[2],digits=2),"%)",sep=""))
  if (ptype=="spider")
    s.class(scaled.scores,factor(Y),grid=F,col=group.cols,cellipse=0.5,cpoint=0,add.plot=T)
  if (ptype=="hull"){

    ll=split(rownames(scaled.scores),Y)
    hulls=lapply(ll,function(ii) ii[chull(scaled.scores[ii,])])
    for (i in 1:length(hulls))
      polygon(scaled.scores[hulls[[i]],],border=group.cols[i])
  }   
  if (plot.names){
     text(scaled.scores,rownames(scaled.scores),cex=0.5,...)
  }else{
    points(scaled.scores,...)
  }
  text(scaled.loadings[which.bac,],rownames(scaled.loadings)[which.bac],cex=cex.bac)
}


