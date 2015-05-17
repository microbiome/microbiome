# Convert data to JSON format for corr_w_scatter visualization
# of correlation matrix linked to scatterplots
# Based on a similar example at 
# http://www.biostat.wisc.edu/~kbroman/D3/corr_w_scatter

#' Interactive scatter plot
#'
#' @param dat data matrix
#' @param group sample groups (vector or factor)
#' @param reorder reorder the data
#'
#' @return JSON
#'
#' @export
#'
#' @import rjson df2json MASS
#'
#' @examples
#' print("Check 
#' https://github.com/microbiome/d3/tree/master/corr_w_scatter 
#' for an example of the convert4corrwscatter function")
#'
#' @references Based on a similar example originally suggested by K. Broman: 
#'           https://www.biostat.wisc.edu/~kbroman/D3/corr_w_scatter
#'    To cite microbiome R package, see citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

convert4corrwscatter <- function(dat, group, reorder=TRUE)
{
  ind <- rownames(dat)
  variables <- colnames(dat)

  if(nrow(dat) != length(group))
    stop("nrow(dat) != length(group)")
  if(!is.null(names(group)) && !all(names(group) == ind))
    stop("names(group) != rownames(dat)")

  if(reorder) {
    ord <- hclust(dist(t(dat)), method="ward")$order
    variables <- variables[ord]
    dat <- dat[,ord]
  }

  # correlation matrix
  corr <- cor(dat, use="pairwise.complete.obs")

  # get rid of names
  dimnames(corr) <- dimnames(dat) <- NULL
  names(group) <- NULL

  # data structure for JSON
  output <- list("ind" = toJSON(ind),
                 "var" = toJSON(variables),
                 "corr" = matrix2json(corr),
                 "dat" =  matrix2json(t(dat)), # columns as rows
                 "group" = toJSON(group))
  paste0("{", paste0("\"", names(output), "\" :", output, collapse=","), "}")
}







#' core.which
#'
#' @param data data matrix; phylotypes vs. samples
#' @param intTr intTr
#' @param prevalenceTr prevalenceTr
#'
#' @return Number of OTUs.
#'
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'
#'   To cite the microbiome R package, see citation('microbiome') 
#'
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
core.which <- function(data, intTr, prevalenceTr) {
    d.bin <- data >= intTr
    prevalences <- rowSums(d.bin)
    nOTUs <- as.numeric(prevalences >= prevalenceTr)
    return(nOTUs)
}

#' plot_cumulative
#'
#' Plot cumulative core microbiota.
#'
#' @param d.sub d.sub
#' @param i.set i.set
#' @param type plot type 
#' @param ylim y axis limits
#' @param tax.table tax.table matrix
#'
#' @return Used for side-effects (plot)
#'
#' @examples 
#' \dontrun{
#'   data(peerj32)
#'   bs <- bootstrap.microbes(t(peerj32$microbes), Nboot = 5);
#'   tax.table <- GetPhylogeny("HITChip")
#'   plot_cumulative(bs, tax.table = tax.table)
#' }
#' @export 
#'
#' @references 
#' 
#' The core microbiota bootstrap method implemented with this function:
#' Salonen A, Salojarvi J, Lahti L, de Vos WM. The adult intestinal
#' core microbiota is determined by analysis depth and health
#' status. Clinical Microbiology and Infection 18(S4):16-20, 2012
#'  
#' To cite this R package, see citation("microbiome") 
#'
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
plot_cumulative <- function(d.sub, i.set = NULL, type = "cumulative", 
		            ylim = NULL, tax.table){

   PH.i <- unique(tax.table[,1:2])
   PH.i <- PH.i[order(PH.i[,2]),]
   rownames(PH.i)=PH.i[,2]
   d.sub$Microbe <- PH.i[d.sub[,1],1]
   d.sub <- d.sub[order(d.sub[,2], decreasing = TRUE),]

   if (is.null(i.set)) {
      i.set <- 1:length(levels(d.sub[,1]))
   }
   colmap <- colorRampPalette(c("Red", "Green","Blue"),
   	     			 space="rgb")(length(levels(d.sub[,1])))
   cnt <- 1;
   i.accept <- vector("logical",length(levels(d.sub[,1])))

   for (i in i.set){
     l.res <- as.numeric(d.sub[,1]==levels(d.sub[,1])[i])
     if (type=="cumulative"){

        out=cumsum(l.res)

     }
     if (type=="gsea"){

        l.res[l.res==0]=-1
        out=l.res
        for (j in 2:length(l.res))
           out[j]=max(l.res[j]+out[j-1],0)
     }
     t1=seq(max(d.sub[,2]),min(d.sub[,2]),-0.01)
     out=vector("numeric",length(t1))
     null.cum=matrix(NA,length(t1),3)
     for (j in 1:length(t1)){
       out[j]=sum(l.res*(d.sub[,2]>t1[j]))
       null.cum[j,]=quantile(replicate(1000,sum(sample(l.res,length(l.res))
		*(d.sub[,2]>=t1[j]))),probs=c(0.025,0.5,0.975))
     }
     yplot <- (out-null.cum[,2])/max(abs(out-null.cum[,2]))
     if (sum(out<null.cum[,1])>0 | sum(out>null.cum[,3])>0){
      if (cnt==1){
        if (is.null(ylim))
           plot(t1, yplot, type="l", 
	     main=paste(type,"prevalence of L1 taxa"), 
	     xlim=c(max(t1),min(t1)), col=colmap[i], 
	     ylab="proportion of total",
	     xlab="Frequency")
        else
           plot(t1,yplot,type="l",
	     main=paste(type,"prevalence of L1 taxa"),
	     ylim=ylim,xlim=c(max(t1),min(t1)),col=colmap[i],
	     ylab="proportion of total",xlab="Frequency")
      } else
         lines(t1,yplot,col=colmap[i])
      cnt <- cnt+1;
      i.accept[i] <- TRUE
    }
   }
   legend(max(t1),1,levels(d.sub[,1])[which(i.accept == TRUE)],fill=colmap[which(i.accept==TRUE)],cex=0.5)

  NULL

}



#' bootstrap.microbes
#'
#' Bootstrap method for core microbiota estimation from Salonen et al. (2012)
#'
#' @param D data (phylotypes x samples)
#' @param Nsample bootstrap sample size, default is the same size as data
#' @param minPrev Lower limit for number of samples where microbe needs 
#'   	    to exceed the intensity threshold for a 'present' call. 
#' @param Nboot bootstrap sample size
#' @param I.thr Lower limit for intensity threshold
#' @param ncore number of nodes for parallelization - default 1
#' @param I.max Upper limit for intensity threshold. Later addition.
#'              set to NULL (default) to replicate Salonen et al.
#'
#' @return data frame with microbes and their frequency of presence in 
#'   	     the core
#'
#' @examples data(peerj32); 
#' 	     bs <- bootstrap.microbes(t(peerj32$microbes), Nboot = 5)
#'
#' @export 
#' @import parallel
#' 
#' @references 
#' 
#' The core microbiota bootstrap method implemented with this function:
#' Salonen A, Salojarvi J, Lahti L, de Vos WM. The adult intestinal
#' core microbiota is determined by analysis depth and health
#' status. Clinical Microbiology and Infection 18(S4):16-20, 2012
#'  
#' To cite this R package, see citation("microbiome") 
#' 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
bootstrap.microbes <- function(D, Nsample = NULL, minPrev = 2, Nboot = 1000, 
		      	       I.thr = 1.8, ncore = 1, I.max=NULL){
   # added threshold for maximum intensity I.max
   if (is.null(Nsample)) {Nsample <- ncol(D)}

   boot <- replicate(Nboot, sample(ncol(D), Nsample, replace = TRUE), 
   	   		    simplify = FALSE)

   # choose intensity such that there is at least one bacteria 
   # fulfilling prevalence criterion
   if (ncore > 1) {
     boot.which <- mclapply(boot, function(x){ 
       Prev = round(runif(1, minPrev, length(x)));
       if (is.null(I.max))
         I.max = max(apply(D[,x], 1, 
        	    function(xx) quantile(xx, probs = (1 - Prev/length(x)))));
       I.max = max(I.thr, I.max);  # Ensure I.max > I.thr, otherwise Insty gives NA's / LL 13.8.2012
       Insty = runif(1, I.thr, I.max);
       return(core.which(D[,x], Insty, Prev))
    }, mc.cores = ncore)
   } else {
     boot.which <- lapply(boot, function(x){ 
       Prev = round(runif(1, minPrev, length(x)));
       if (is.null(I.max))
          I.max = max(apply(D[,x], 1, 
                      function(xx) quantile(xx, probs = (1 - Prev/length(x)))));
       I.max = max(I.thr, I.max); # Ensure I.max > I.thr, otherwise Insty gives NA's / LL 13.8.2012
       Insty = runif(1, I.thr, I.max);
       return(core.which(D[,x], Insty, Prev))
    })
   }

   boot.prob <- rowSums(as.data.frame(boot.which))/Nboot

   df <- data.frame(Microbe = rownames(D), Frequency = boot.prob)
   df <- df[order(df$Frequency,decreasing = TRUE),]

   mm <- bootstrap.microbecount(D,Nsample = Nsample, minprev = minPrev, 
      	 	Nboot = Nboot, I.thr = I.thr, ncore = ncore, I.max=I.max)

   df$suggested.core <- as.numeric(sapply(1:nrow(df),function(x) x<= mm))

   message("Bootstrapped median core size: ", mm)

   return(df)

}


#' bootstrap.microbecount
#'
#' Auxiliary function for bootstrap.microbes
#'
#' @param D data
#' @param Nsample sample size
#' @param minprev minimum prevalence
#' @param Nboot bootstrap sample size
#' @param I.thr intensity threshold
#' @param ncore number of nodes for parallelization
#' @param I.max max intensity threshold
#'
#' @return median microbe count in bootstrapped cores
#'
#' @examples 
#'   data(peerj32)
#'   tmp <- bootstrap.microbecount(t(peerj32$microbes),	Nboot = 5)
#'
#' @export 
#' 
#' @references 
#' 
#' The core microbiota bootstrap method implemented with this function:
#' Salonen A, Salojarvi J, Lahti L, de Vos WM. The adult intestinal
#' core microbiota is determined by analysis depth and health
#' status. Clinical Microbiology and Infection 18(S4):16-20, 2012
#'  
#' To cite this R package, see citation("microbiome") 
#'
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

bootstrap.microbecount <- function(D, Nsample = NULL, minprev = 1, 
		       	  	      Nboot = 1000, I.thr = 1.8, ncore = 1,I.max=NULL){

  if (is.null(Nsample)) {Nsample <- ncol(D)}

   boot <- replicate(Nboot,sample(ncol(D),Nsample, replace = TRUE),
   	   					   simplify = FALSE)

   # below: choose intensity such that there is at least one bacteria 
   # fulfilling prevalence criterion
   if (Nsample>1 && ncore > 1) {
     boot.which=mclapply(boot,function(x){
        if (is.null(I.max)) 
           I.max=max(apply(D[,x],1,min))
        Insty=runif(1,I.thr,I.max)
        sum(rowSums(D[,x]>Insty) >= minprev)
     }, mc.cores = ncore)
   } else if (Nsample>1 && ncore == 1) {
     boot.which=lapply(boot,function(x){ 
        if (is.null(I.max)) 
          I.max=max(apply(D[,x],1,min))
        Insty=runif(1,I.thr,I.max)
        sum(rowSums(D[,x]>Insty) >= minprev)
     }) 
   } else {
     boot.which=lapply(boot,function(x){ 
        if (is.null(I.max)) 
           I.max=max(D[,x])
        Insty=runif(1,I.thr,I.max)
        return(sum(D[,x]>=Insty))
     })
   }

   boot.prob <- as.matrix(as.data.frame(boot.which, check.names = FALSE))
   t1 <- quantile(boot.prob, probs = c(0.05, 0.5, 0.95))
   #t1[2] <- mean(boot.prob)

   #print(t1)
   return(t1[2])
}





