#' @title ANCOM Test for Abundance Differences
#' @description Perform ANCOM test on phyloseq data.
#' @param x \code{\link{phyloseq-class}}
#' @param varname Name of the metadata variable to be used for grouping
#' @param multcorr_type takes value 1 (stringent correction on the full pairwise matrix), OR 2 (within OTU/taxa multiple correction), OR 3 (no correction).
#' @param sig significance level (default 0.05)
#' @return Vector of adjusted p-values, one for each taxon.
#' @examples
#'   # Taxa that show significant differences between nationalities
#'   data(dietswap)
#'   s <- ancom(dietswap, "nationality")
#'   print(names(which(s < 0.05)))
#' @details This function is directly modified from Weiss et al. (2017). The R code was obtained from the first author and included with permission in the microbiome package. For reference to the original ANCOM method by Mandal et al. (2015).
#' @export
#' @references 
#' Weiss et al. Normalization and microbial differential abundance strategies depend upon data characteristics. Microbiome 5:27, 2017.
#'
#' Mandal S et al. Analysis of composition of microbiomes: a novel method for studying microbial composition. Microbial Ecology in Health and Disease, 26:1–7, 2015.
#'
#' To cite the microbiome R package, see citation('microbiome')
#'
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
ancom <- function (x, varname, multcorr_type = 2, sig = 0.05) {

  physeq <- x  

  # Enforce orientation
  if (!taxa_are_rows(physeq)) {
      physeq <- t(physeq)
  }
  real.data = as(otu_table(physeq), "matrix")
  Group = as.factor(get_variable(physeq, varname))
  colnames(real.data) <- Group
  
  # Convert sample_data to AnnotatedDataFrame
  par1_new <- nrow(real.data)
  
  W.detected <- ancom.detect(real.data,par1_new,sig,multcorr_type)
  
  ### Detected using a stepwise mode detection
  if(max(W.detected)/par1_new >=0.10){
    c.start=max(W.detected)/par1_new
    cutoff=c.start-c(0.05,0.10,0.15,0.20,0.25)
    prop_cut=rep(0,length(cutoff))
    for(cut in 1:length(cutoff)){
      prop_cut[cut]=length(which(W.detected>=par1_new*cutoff[cut]))/length(W.detected)
    } 
    del=rep(0,length(cutoff)-1)
    for(i in 1:(length(cutoff)-1)){
      del[i]=abs(prop_cut[i]-prop_cut[i+1])
    }
    
    if(del[1]<0.02&del[2]<0.02&del[3]<0.02){nu=cutoff[1]
    }else if(del[1]>=0.02&del[2]<0.02&del[3]<0.02){nu=cutoff[2]
    }else if(del[2]>=0.02&del[3]<0.02&del[4]<0.02){nu=cutoff[3]                                
    }else{nu=cutoff[4]}
    
    up_point=min(W.detected[which(W.detected>=nu*par1_new)])
    
    W.detected[W.detected>=up_point]=99999
    W.detected[W.detected<up_point]=0
    W.detected[W.detected==99999]=1
  }else{W.detected=0}
  
  detected_stepwise=rownames(real.data)[which(W.detected==1)]
  
  results_list <- list(detected_stepwise)
  padj <- rep(1, length(rownames(real.data)))
  ind <- which(rownames(real.data) %in% unlist(results_list))
  padj[ind] <- 0
  
  res <- as.numeric(as.character(padj))
  names(res) <- rownames(real.data)
  res <- sort(res)

  return(res)

}






ancom.detect <- function(otu_data,n_otu,alpha,multcorr){

  ## arbit: arbitrary cutoff for the detection algorithm
  ### Values returned by function 
  ### $Arbitrary: A list of OTUs/taxa detected using an arbitrary cutoff
  ### $Stepwise: A list of OTUs/taxa detected using the stepwise mode detection

  logratio.mat=matrix(NA,nrow=n_otu,ncol=n_otu)
  Group = as.factor(colnames(otu_data))
  for(i in 1:(n_otu-1)){

    for(j in (i+1):n_otu){
      data.pair=otu_data[c(i,j),]    #,n_otu+1)]
      lr=log((0.001+as.numeric(data.pair[1,]))/(0.001+as.numeric(data.pair[2,])))
      logratio.mat[i,j]=wilcox.test(lr[colnames(data.pair)==unique(colnames(data.pair))[1]],lr[colnames(data.pair)==unique(colnames(data.pair))[2]])$p.value
    }
  }
  ind <- lower.tri(logratio.mat)
  logratio.mat[ind] <- t(logratio.mat)[ind]
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)]=1
  
  mc.pval=t(apply(logratio.mat,1,function(x){
    s=p.adjust(x, method = "BH")
    return(s)
  }))
  
  a=logratio.mat[upper.tri(logratio.mat,diag=F)==T]
  
  b=matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T]=p.adjust(a, method = "BH")
  diag(b)=NA
  ind.1 <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  if(multcorr==2){
    W=apply(mc.pval,1,function(x){
      subp=length(which(x<alpha))
    })
  }else if(multcorr==1){
    W=apply(b,1,function(x){
      subp=length(which(x<alpha))
    })
  }else if(multcorr==3){
    W=apply(logratio.mat,1,function(x){
      subp=length(which(x<alpha))
    })
  }
  return(W)
}






