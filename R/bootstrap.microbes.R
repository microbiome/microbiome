#' @title Bootstrap microbes
#' @description Bootstrap method for core microbiota estimation from Salonen et al. (2012)
#' @param D data (phylotypes x samples)
#' @param Nsample bootstrap sample size, default is the same size as data
#' @param minPrev Lower limit for number of samples where microbe needs 
#'   	    to exceed the intensity threshold for a 'present' call. 
#' @param Nboot bootstrap sample size
#' @param I.thr Lower limit for intensity threshold
#' @param ncore number of nodes for parallelization - default 1
#' @param I.max Upper limit for intensity threshold. Later addition.
#'              set to NULL (default) to replicate Salonen et al.
#' @return data frame with microbes and their frequency of presence in 
#'   	     the core
#' @examples data(peerj32); 
#' 	     bs <- bootstrap.microbes(t(peerj32$microbes), Nboot = 5)
#' @export 
#' @importFrom parallel mclapply
#' @references 
#' The core microbiota bootstrap method implemented with this function:
#' Salonen A, Salojarvi J, Lahti L, de Vos WM. The adult intestinal
#' core microbiota is determined by analysis depth and health
#' status. Clinical Microbiology and Infection 18(S4):16-20, 2012
#' To cite this R package, see citation("microbiome")  
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




#' @title Core which
#' @description Auxiliary function 
#' @param data phylotypes vs. samples data matrix
#' @param intTr intTr
#' @param prevalenceTr prevalenceTr
#' @return Number of OTUs.
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
core.which <- function(data, intTr, prevalenceTr) {
    d.bin <- data >= intTr
    prevalences <- rowSums(d.bin)
    nOTUs <- as.numeric(prevalences >= prevalenceTr)
    return(nOTUs)
}


#' @title Bootstrap microbecount
#' @description Auxiliary function for bootstrap.microbes.
#' @param D data
#' @param Nsample sample size
#' @param minprev minimum prevalence
#' @param Nboot bootstrap sample size
#' @param I.thr intensity threshold
#' @param ncore number of nodes for parallelization
#' @param I.max max intensity threshold
#' @return median microbe count in bootstrapped cores
#' @examples 
#'   data(peerj32)
#'   tmp <- bootstrap.microbecount(t(peerj32$microbes),	Nboot = 5)
#' @export
#' @importFrom parallel mclapply
#' @references 
#' The core microbiota bootstrap method implemented with this function:
#' Salonen A, Salojarvi J, Lahti L, de Vos WM. The adult intestinal
#' core microbiota is determined by analysis depth and health
#' status. Clinical Microbiology and Infection 18(S4):16-20, 2012
#' To cite this R package, see citation("microbiome") 
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
        I.max = max(I.thr, I.max); 
        Insty=runif(1,I.thr,I.max)
        sum(rowSums(D[,x]>Insty) >= minprev)
     }, mc.cores = ncore)
   } else if (Nsample>1 && ncore == 1) {
     boot.which=lapply(boot,function(x){ 
        if (is.null(I.max)) 
          I.max=max(apply(D[,x],1,min))
        I.max = max(I.thr, I.max); 
        Insty=runif(1,I.thr,I.max)
        sum(rowSums(D[,x]>Insty) >= minprev)
     }) 
   } else {
     boot.which=lapply(boot,function(x){ 
        if (is.null(I.max)) 
           I.max=max(D[,x])
        I.max = max(I.thr, I.max); 
        Insty=runif(1,I.thr,I.max)
        return(sum(D[,x]>=Insty))
     })
   }

   boot.prob <- as.matrix(as.data.frame(boot.which, check.names = FALSE))
   t1 <- quantile(boot.prob, probs = c(0.05, 0.5, 0.95))
   return(t1[2])
}