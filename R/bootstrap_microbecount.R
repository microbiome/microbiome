#' @title Bootstrap microbecount
#' @description Auxiliary function for bootstrap_microbes.
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
#'   tmp <- bootstrap_microbecount(t(peerj32$microbes),	Nboot = 5)
#' @references 
#' The core microbiota bootstrap method implemented with this function:
#' Salonen A, Salojarvi J, Lahti L, de Vos WM. The adult intestinal
#' core microbiota is determined by analysis depth and health
#' status. Clinical Microbiology and Infection 18(S4):16-20, 2012
#' To cite this R package, see citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords internal
bootstrap_microbecount <- function(D, Nsample = NULL, minprev = 1, 
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

