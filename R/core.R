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
#' @param phylogeny.info phylogeny.info matrix
#'
#' @return Used for side-effects (plot)
#'
#' @examples 
#' \dontrun{
#'   data(peerj32)
#'   bs <- bootstrap.microbes(t(peerj32$microbes), Nboot = 5);
#'   phylogeny.info <- GetPhylogeny("HITChip")
#'   plot_cumulative(bs, phylogeny.info = phylogeny.info)
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
		   		   ylim = NULL, phylogeny.info){

   PH.i <- unique(phylogeny.info[,1:2])
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
		      	       I.thr = 1.8, ncore = 1){

   if (is.null(Nsample)) {Nsample <- ncol(D)}

   boot <- replicate(Nboot, sample(ncol(D), Nsample, replace = TRUE), 
   	   		    simplify = FALSE)

   # choose intensity such that there is at least one bacteria 
   # fulfilling prevalence criterion
   if (ncore > 1) {
     boot.which <- mclapply(boot, function(x){ 
       Prev = round(runif(1, minPrev, length(x)));
       Imax = max(apply(D[,x], 1, 
       	    function(xx) quantile(xx, probs = (1 - Prev/length(x))))); # Ensure Imax > I.thr, otherwise Insty gives NA's / LL 13.8.2012
       Imax = max(I.thr, Imax); 
       Insty = runif(1, I.thr, Imax);
       return(core.which(D[,x], Insty, Prev))
    }, mc.cores = ncore)
   } else {
     boot.which <- lapply(boot, function(x){ 
       Prev = round(runif(1, minPrev, length(x)));
       Imax = max(apply(D[,x], 1, 
       	    function(xx) quantile(xx, probs = (1 - Prev/length(x)))));
       Imax = max(I.thr, Imax); # Ensure Imax > I.thr, otherwise Insty gives NA's / LL 13.8.2012
       Insty = runif(1, I.thr, Imax);
       return(core.which(D[,x], Insty, Prev))
    })
   }

   boot.prob <- rowSums(as.data.frame(boot.which))/Nboot

   df <- data.frame(Microbe = rownames(D), Frequency = boot.prob)
   df <- df[order(df$Frequency,decreasing = TRUE),]

   mm <- bootstrap.microbecount(D,Nsample = Nsample, minprev = minPrev, 
      	 	Nboot = Nboot, I.thr = I.thr, ncore = ncore)

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
#' @param I.thr threshold
#' @param ncore number of nodes for parallelization
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
		       	  	      Nboot = 1000, I.thr = 1.8, ncore = 1){

  if (is.null(Nsample)) {Nsample <- ncol(D)}

   boot <- replicate(Nboot,sample(ncol(D),Nsample, replace = TRUE),
   	   					   simplify = FALSE)

   # below: choose intensity such that there is at least one bacteria 
   # fulfilling prevalence criterion
   if (Nsample>1 && ncore > 1) {
     boot.which=mclapply(boot,function(x){ 
        Imax=max(apply(D[,x],1,min))
        Insty=runif(1,I.thr,Imax)
        sum(rowSums(D[,x]>Insty) >= minprev)
     }, mc.cores = ncore)
   } else if (Nsample>1 && ncore == 1) {
     boot.which=lapply(boot,function(x){ 
        Imax=max(apply(D[,x],1,min))
        Insty=runif(1,I.thr,Imax)
        sum(rowSums(D[,x]>Insty) >= minprev)
     }) 
   } else {
     boot.which=lapply(boot,function(x){ 
        Imax=max(D[,x])
        Insty=runif(1,I.thr,Imax)
        return(sum(D[,x]>=Insty))
     })
   }

   boot.prob <- as.matrix(as.data.frame(boot.which, check.names = FALSE))
   t1 <- quantile(boot.prob, probs = c(0.05, 0.5, 0.95))
   #t1[2] <- mean(boot.prob)

   #print(t1)
   return(t1[2])
}


#' core.sum
#'
#' @param data data matrix; phylotypes vs. samples
#' @param intTr intTr
#' @param prevalenceTr prevalenceTr
#'
#' @return Number of OTUs
#'
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

core.sum <- function(data, intTr, prevalenceTr) {
    d.bin <- data > intTr
    prevalences <- rowSums(d.bin)
    # jos haluat tietaa lajit, ala summaa!
    nOTUs <- sum(prevalences >= prevalenceTr)
    return(nOTUs)
}




#' createCore
#'
#' create coreMatrix 
#'
#' @param data data matrix; phylotypes vs. samples
#' @param verbose verbose
#' @param prevalence.intervals a vector of prevalence percentages in [0,100]
#' @param intensity.intervals a vector of intensities around the data range
#'
#' Returns:
#' @return Estimated core microbiota
#'
#' @examples 
#'   data(peerj32)
#'   core <- createCore(t(peerj32$microbes))
#'
#' @export 
#' 
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

createCore <- function(data, verbose = FALSE, 
          prevalence.intervals = seq(20, 100, 20), 
          intensity.intervals = NULL) {
    
    ## Prevalence vector
    if (is.null(prevalence.intervals)) {
        prevalence.intervals <- seq(0, 100, 10)
    }
    # Convert prevalences from percentages to numerics
    p.seq <- 0.01 * prevalence.intervals * ncol(data)
    
    ## Intensity vector
    if (is.null(intensity.intervals)) {
        i.seq <- seq(min(data), max(data), length = 10)
    } else {
        i.seq <- intensity.intervals
    }
    
    coreMat <- matrix(NA, nrow = length(i.seq), ncol = length(p.seq), 
                      dimnames = list(i.seq, p.seq))
    
    n <- length(i.seq) * length(p.seq)
    cnt <- 0
    for (i in i.seq) {
        for (p in p.seq) {
            if (verbose) {
                cnt <- cnt + 1
                message(cnt/n)
            }
            coreMat[as.character(i), as.character(p)] <- core.sum(data, i, p)
        }
    }
    
    # Convert Prevalences to percentages
    colnames(coreMat) <- 100 * as.numeric(colnames(coreMat))/ncol(data)
    
    return(coreMat)
}






#' Core2D
#'
#' Core visualization 2D
#'
#' @param coreMat core matrix
#' @param title title
#' @param plot plot the figure 
#' @param xlabel X axis label
#' @param ylabel Y axis label
#'  
#' @return Used for its side effects
#'
#' @examples 
#'   data(peerj32)
#'   c2d <- Core2D(createCore(t(peerj32$microbes)))
#' @export 
#' 
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

Core2D <- function(coreMat, title = "Common core", plot = TRUE, 
                   xlabel = "Abundance", 
                   ylabel = "Core size (number of taxa)") {
    
    Abundance <- Prevalence <- Count <- NULL
    
    df <- melt(coreMat)
    names(df) <- c("Abundance", "Prevalence", "Count")
    theme_set(theme_bw(20))
    p <- ggplot(df, aes(x = Abundance, y = Count, color = Prevalence, 
                        group = Prevalence))
    p <- p + geom_line()
    p <- p + geom_point()
    p <- p + xlab(xlabel)
    p <- p + ylab(ylabel)
    p <- p + ggtitle("Core microbiota")
    
    if (plot) {
        print(p)
    }
    
    return(p)
    
}


#' core_heatmap
#'
#' Heatmap of core microbiota
#'
#' @param data data matrix: phylotypes vs. samples
#' @param detection.thresholds Vector of detection thresholds
#' @param plot plot the figure
#' @param palette 'bw' (grayscale) or 'spectral' (colourscale)
#'
#' @return List with the following elements: 
#'         plot: ggplot figure
#'         data: prevalence data with the varying thresholds
#'
#' @examples 
#'   data(peerj32)
#'   core <- core_heatmap(t(peerj32$microbes))
#'
#' @export 
#' @importFrom reshape melt
#' @import ggplot2
#' @import RColorBrewer
#' 
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

core_heatmap <- function(data, detection.thresholds = NULL, plot = TRUE, 
                         palette = "bw") {
    
    DetectionThreshold <- Taxa <- Prevalence <- NULL
    
    if (is.null(detection.thresholds)) {
        detection.thresholds <- seq(min(data), max(data), length = 10)
    }
    
    # Prevalences with varying detection thresholds
    taxa <- rownames(data)
    prevalences <- matrix(NA, nrow = length(taxa), 
                          ncol = length(detection.thresholds))
    rownames(prevalences) <- taxa
    colnames(prevalences) <- as.character(detection.thresholds)
    for (det.th in detection.thresholds) {
        prevalence <- 100 * sort(rowMeans(data > det.th))
        prevalences[taxa, as.character(det.th)] <- prevalence[taxa]
    }
    
    df <- melt(prevalences)
    names(df) <- c("Taxa", "DetectionThreshold", "Prevalence")
    o <- names(sort(rowSums(prevalences)))
    df$Taxa <- factor(df$Taxa, levels = o)
    theme_set(theme_bw(10))
    p <- ggplot(df, aes(x = DetectionThreshold, y = Taxa, fill = Prevalence))
    p <- p + geom_tile()
    
    if (palette == "bw") {
        colours <- c("black", "darkgray", "gray", "lightgray", "white")
    } else if (palette == "spectral") {
        myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
        colours <- myPalette(5)
    }
    
    p <- p + scale_fill_gradientn("Prevalence", 
        breaks = seq(from = 0, to = 100, 
        by = 10), colours = colours, limits = c(0, 100))
    
    p <- p + ggtitle("Core microbiota")
    
    if (plot) {
        print(p)
    }
    
    return(list(plot = p, data = prevalences))
    
} 
