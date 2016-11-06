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
#' @export
#' @examples
#'  \dontrun{
#'    data(peerj32)
#'    x <- as.matrix(peerj32$microbes)[1:20, 1:6]
#'    y <- rnorm(nrow(x))
#'    names(y) <- rownames(x)
#'    res <- Bagged.RDA(x, y , boot = 5)
#' }
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
Bagged.RDA <- function(X, Y, boot = 1000){

  ## Jarkko Salojarvi 7.8.2012
  ##  #17.8.2012 fixed problem with multiclass RDA  

   if (is.numeric(boot)){
      class.split=split(names(Y),Y)

      boot=replicate(boot,unlist(sapply(class.split,function(x) sample(x,length(x),replace=T))),simplify=F)
   }
   nboot <- length(boot)
   n.lev <- length(levels(Y))
   TT <- scores(rda(t(X)~Y),choices=1:max(n.lev-1,2),display="sites")
   nRDA <- ncol(TT) 
   # rotate 
   rotateMat <- function(M,TT,x){
      M.rot <- procrustes(TT[x,],M)
      return(M.rot$Yrot)
   }
   nearest.centers=function(xx,cc){
     nC=nrow(cc)
     nS=nrow(xx)
     MM=as.matrix(dist(rbind(cc,xx)))[1:nC,(nC+1):(nC+nS)]
     if (nS>1)
       apply(MM,2,which.min)
     else
       which.min(MM)
   }
   # bootstrap loadings
   Tx=lapply(boot,function(x){ 
      nC=length(levels(Y[x]))
      M=rda(t(X[,x])~Y[x])
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
    for (i in 2:nboot)
      bagged.loadings=bagged.loadings+Tx[[i]]$loadingsX
    bagged.loadings=bagged.loadings/nboot
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
    #err.t=mean((Y-prediction(X,Y,bagged.loadings,bagged.loadings.Y,bagged.proj,mean(Y),rowMeans(X)))^2)
    #random.pred=Y-rmultinom(length(Y),1,
    err.random=replicate(nboot,mean((as.numeric(Y)-sample(as.numeric(Y)))!=0))
    bagged.error=mean(sapply(Tx, function(x) x$err),na.rm=T)
    R=(bagged.error-err.t)/(mean(err.random)-err.t)
    R=max(min(R,1),0)
    w=.632/(1-.368*R)
    bagged.R2=(1-w)*bagged.error+w*err.t
    #bagged.R2=1-((1-w)*err.t+w*err.b)/mean(Y^2)
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
	 err.significance=sum(err.random>bagged.R2)/nboot,
	 R2=Rsquare,R2.variables=Rsquare.variable)
}

