#' @title Bagged RDA
#' @description Bagged (or Bootstrap Aggregated) RDA feature selection
#' @param x a matrix, samples on columns, variables (bacteria) on rows. 
#'        Or a \code{\link{phyloseq-class}} object
#' @param y vector with names(y)=rownames(X). 
#'            Or name of phyloseq sample data variable name.
#' @param sig.thresh signal p-value threshold, default 0.1
#' @param nboot Number of bootstrap iterations
#' @param verbose verbose
#' @param plot Also show a diagnostic plot of the result
#' @return List with items:
#'   \itemize{
#'     \item{loadings}{bagged loadings}
#'     \item{scores}{bagged scores}
#'     \item{significance}{significances of X variables}
#'     \item{group.centers}{group centers on latent space}
#'     \item{bootstrapped}{bootstrapped loadings}
#'     \item{data}{data set with non-significant components dropped out}
#'   }
#' @examples \dontrun{
#'   library(microbiome)
#'
#'   # Example with abundance matrix
#'   data(peerj32)
#'   phy <- peerj32$phyloseq
#'   x <- taxa_abundances(phy) 
#'   y <- factor(sample_data(phy)$gender);
#'   names(y) <- rownames(sample_data(phy))
#'   res <- bagged_rda(x, y, sig.thresh=0.05, nboot=20)
#'   plot_bagged_rda(res, y)
#'
#'   # Example with phyloseq object
#'   res <- bagged_rda(phy, "gender", sig.thresh=0.05, nboot=20)
#'   plot_bagged_rda(res, y)
#'
#'  }
#' @export
#' @details Bootstrap aggregation (bagging) is expected to improve the stability of the results. Aggregating results over several modeling runs with different boostrap samples of the data are averaged to produce the final summary.
#' @references See citation("microbiome") 
#' @author Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
bagged_rda <- function(x, y, sig.thresh = 0.1, nboot = 1000, verbose = T, plot = FALSE){

  if (class(x) == "phyloseq") {
    # Pick OTU matrix and the indicated annotation field
    if (!y %in% names(sample_data(x))) {
      stop(paste("The variable y ('", y, "') is not available in the phyloseq object i.e. sample_data(x). Only use variables listed in names(sample_data(x)) ie. one of the following: ", paste(names(sample_data(x)), collapse = " / "), sep = ""))
    }

    y <- factor(sample_data(x)[[y]])
    names(y) <- sample_data(x)$sample
    x <- taxa_abundances(x) 
  }

  stop.run=F
  class.split=split(names(y),y)
  dropped=vector("character",nrow(x))
  x.all=x
  mean.err=rep(1,nrow(x))
  while(stop.run==F){
    boot=replicate(nboot,unlist(sapply(class.split,function(x) sample(x,length(x),replace=T))),simplify=F)
    Bag.res=Bagged.RDA(x,y,boot=boot)
    min.prob=Bag.res$significance[[1]]
    if (length(levels(y))>2){
      for (i in 1:nrow(x))
         min.prob[i]=min(sapply(Bag.res$significance,function(x) x[i]))
    }
    mean.err[nrow(x)]=Bag.res$error
    dropped[nrow(x)]=rownames(x)[which.max(min.prob)]
    if (verbose) {message(c(nrow(x), Bag.res$error))}
    if (nrow(x)>max(length(class.split),2))
      x=x[-which.max(min.prob),]
    else
      stop.run=T
  }
  dropped[1:length(class.split)]=rownames(x)[order(min.prob)[1:length(class.split)]]
  best.res=which.min(mean.err)

  Bag.res=Bagged.RDA(x.all[dropped[1:best.res],],y,boot=boot)
  Bag.res$data=x.all[dropped[1:best.res],]
  Bag.res$Err.selection=mean.err
  Bag.res$dropped=dropped

  if (plot) {
    plot(mean.err,xlab="x dimension")
    points(best.res,mean.err[best.res],col="red")
  }
  
  list(bagged.rda = Bag.res, variable = y)

}

