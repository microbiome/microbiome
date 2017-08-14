#' @title Bagged RDA
#' @description Bagged (or Bootstrap Aggregated) RDA feature selection
#' @param x a matrix, samples on columns, variables (bacteria) on rows. 
#'        Or a \code{\link{phyloseq-class}} object
#' @param y vector or factor with names(y)=rownames(X). 
#'            Or name of phyloseq sample data variable name
#'            (one of sample_variables(x)).
#' @param bs.iter Number of bootstrap iterations
#' @param verbose verbose
#' @return List with items:
#'   \itemize{
#'     \item{loadings}{bagged loadings}
#'     \item{significance}{significances of X variables}
#'     \item{scores}{bagged scores}
#'     \item{group.centers}{group centers on latent space}
#'     \item{bootstrapped}{bootstrapped loadings}
#'     \item{data}{data set with non-significant components dropped out}
#'   }
#' @examples \dontrun{
#'   # Example with abundance matrix
#'   data(peerj32)
#'   pseq <- peerj32$phyloseq
#' 
#'   # Run Bagged RDA
#'   # For any real study, use bs.iter = 100 or higher
#'   res <- rda_bagged(pseq, "gender", bs.iter=50) 
#'
#'   # Visualize Bagged RDA
#'   # plot_rda_bagged(res, y)
#'
#'   # Example with phyloseq object
#'   # res <- rda_bagged(phy, "gender", bs.iter=50)
#'   # plot_rda_bagged(res, y)
#'
#'
#'  }
#' @export
#' @seealso vegan::rda and phyloseq::ordinate
#' @details Bootstrap aggregation (bagging) is expected to improve the stability of the results. Aggregating results over several modeling runs with different boostrap samples of the data are averaged to produce the final summary.
#' @references See citation("microbiome") 
#' @author Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
rda_bagged <- function(x, y, bs.iter = 100, verbose = T){

  if (class(x) == "phyloseq") {
    # Pick OTU matrix and the indicated annotation field
    if (!y %in% sample_variables(x)) {
    
      stop(paste("The variable y ('", y, "') is not available in the phyloseq object i.e. sample_data(x). Only use variables listed in sample_variables(x) ie. one of the following: ", paste(names(sample_data(x)), collapse = " / "), sep = ""))
    }

    if (!"sample" %in% sample_variables(x)) {
      warning("The sample_data(x) does not contain 'sample' field; using the rownames(sample_data(x)) as sample ID.")
      sample_data(x)$sample <- rownames(sample_data(x))
    }

    y <- factor(sample_data(x)[[y]])
    names(y) <- sample_data(x)$sample
    
    x <- abundances(x)
    
  }

  stop.run=F
  class.split=split(names(y),y)
  dropped=vector("character",nrow(x))
  x.all=x
  mean.err=rep(1,nrow(x))

  while(stop.run==F){
  
    boot = replicate(bs.iter,unlist(sapply(class.split,function(x) sample(x,length(x),replace=T))),simplify=F)
    
    Bag.res = Bagged.RDA(x,y,boot=boot)
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
  Bag.res$variable=y

  #if (plot) {
  #  plot(mean.err,xlab="x dimension")
  #  points(best.res,mean.err[best.res],col="red")
  #}
  
  Bag.res

}




#' @title Bagged RDA
#' @description Bootstrap solutions that follows the Jack-knife estimation of
#'              PLS by Martens and Martens, 2000.  Solves rotational
#'              invariance of latent space by orthogonal procrustes rotations.
#' @param X a matrix, samples on columns, variables (bacteria) on rows.
#' @param Y vector with names(Y)=rownames(X), for example
#' @param boot Number of bootstrap iterations
#' @return List with elements:
#'   \itemize{
#'     \item{loadings}{bagged loadings}
#'     \item{scores}{bagged scores}
#'     \item{significance}{significances of X variables}
#'   }
#' @examples
#'  \dontrun{
#'    data(peerj32)
#'    x <- as.matrix(peerj32$microbes)[1:20, 1:6]
#'    y <- rnorm(nrow(x))
#'    names(y) <- rownames(x)
#'    res <- Bagged.RDA(x, y , boot = 50)
#' }
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
Bagged.RDA <- function(X, Y, boot = 100){

  ## Jarkko Salojarvi 7.8.2012
  ##  #17.8.2012 fixed problem with multiclass RDA  

   if (is.numeric(boot)){
      class.split <- split(names(Y),Y)

      boot <- replicate(boot,unlist(sapply(class.split,function(x) sample(x,length(x),replace=T))),simplify=F)
   }

   boot <- length(boot)
   n.lev <- length(levels(Y))
   TT <- scores(rda(t(X)~Y),choices=1:max(n.lev-1,2),display="sites")
   nRDA <- ncol(TT)
   
   # Rotate 
   rotateMat <- function(M,TT,x){
      M.rot <- procrustes(TT[x,],M)
      return(M.rot$Yrot)
   }

   # Auxiliary function
   nearest.centers <- function(xx,cc){
     nC=nrow(cc)
     nS=nrow(xx)
     MM=as.matrix(dist(rbind(cc,xx)))[1:nC,(nC+1):(nC+nS)]
     if (nS>1)
       apply(MM,2,which.min)
     else
       which.min(MM)
   }
   
   # bootstrap loadings
   Tx <- lapply(boot,function(x){
   
      nC = length(levels(Y[x]))

      M =  rda(t(X[,x])~Y[x])

      # get scores
      TT.m=scores(M,choices=1:max(nC-1,2),display="sites")
      
      # bootstrap error
      testset=setdiff(colnames(X),x)
      err=NA
      
      if (length(testset)>0){
        Pr=predict(M,t(as.data.frame(X[,testset])),type="wa",model="CCA")
        centers=apply(TT.m,2,function(z) sapply(split(z,Y[x]),mean))
        if (nC==2)
           y.pred=apply(sapply(Pr,function(x) (x-centers[,1])^2),2,which.min)
        else
           y.pred=nearest.centers(Pr,centers) 
        err=mean(y.pred-as.numeric(Y[testset])!=0)
      }
      
      # procrustes rotation of scores
      TT.m=rotateMat(TT.m,TT,x)
      
      # solve loadings
      a=t(TT.m) %*% TT.m
      b=as.matrix(X[,x]-rowMeans(X[,x]))
      loadingsX=t(solve(a,t(b %*% TT.m)))
      list(loadingsX=loadingsX,err=err)
   })

   # significances   
   sig.prob=list()
   for (i in 1:nRDA){
     tmp=sapply(Tx,function(x) x$loadingsX[,i])
     sig.prob[[i]]=apply(tmp,1,function(x){ x1=sum(x>0)/length(x); x2=1-x1; min(x1,x2)})
    }
    names(sig.prob)=colnames(TT)
    # bagged estimates
    bagged.loadings=Tx[[1]]$loadingsX
    for (i in 2:boot)
      bagged.loadings=bagged.loadings+Tx[[i]]$loadingsX
    bagged.loadings=bagged.loadings/boot
    colnames(bagged.loadings)=colnames(TT)

    # solve scores
    a=t(bagged.loadings) %*% bagged.loadings
    b=as.matrix(X-rowMeans(X))
    bagged.scores=t(solve(a,t(bagged.loadings) %*% b))
    colnames(bagged.scores)=colnames(TT)
    
    # Group centers 
    Group.center=apply(bagged.scores,2,function(x) sapply(split(x,Y),mean))
    err.t=mean((nearest.centers(bagged.scores,Group.center)-as.numeric(Y))!=0) 

    # bagged error
    err.random=replicate(boot,mean((as.numeric(Y)-sample(as.numeric(Y)))!=0))
    bagged.error=mean(sapply(Tx, function(x) x$err),na.rm=T)
    R=(bagged.error-err.t)/(mean(err.random)-err.t)
    R=max(min(R,1),0)
    w=.632/(1-.368*R)
    bagged.R2=(1-w)*bagged.error+w*err.t

    can.cor.R=apply(X,1,function(x) cor(x,bagged.scores))^2
    Rsquare=rowSums(can.cor.R)/sum(diag(var(t(X))))
    names(Rsquare)=colnames(bagged.scores)
    Rsquare.variable=t(can.cor.R/apply(X,1,var))
    colnames(Rsquare.variable)=colnames(bagged.scores)

    # The CRAN/BioC recommendations do not allow lines over 100 chars
    list(loadings = bagged.loadings, scores = bagged.scores,
    	 significance = sig.prob,
	 error=bagged.R2,group.centers=Group.center,
	 bootstrapped=Tx,err.random=mean(err.random),
	 err.significance=sum(err.random>bagged.R2)/boot,
	 R2=Rsquare,R2.variables=Rsquare.variable)
}

