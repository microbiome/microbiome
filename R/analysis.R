# Copyright (C) 2011-2012 Leo Lahti and Jarkko Salojarvi 
# Contact: <microbiome-admin@googlegroups.com>. All rights reserved.

# This file is a part of the microbiome R package
# http://microbiome.github.com/

# This program is open source software; you can redistribute it and/or
# modify it under the terms of the FreeBSD License (keep this notice):
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.




#' Description: Calculate distance matrix between the _columns_ of the 
#' input matrix. Can prduce correlation-based distance matrices, otherwise
#' uses the standard 'dist' function.
#'
#' Arguments:
#'   @param x data matrix
#'   @param method distance method
#'   @param ... other arguments to be passed
#'
#' Returns:
#'   @return distance object
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

distance.matrix <- function (x, method = "pearson", ...) {
  if (method %in% c("pearson", "spearman")) {
    cmat <- as.dist((1-cor(x, use = "pairwise.complete", method = method)))
  } else {
    cmat <- dist(x, method = method, ...) 
  }
  cmat
}




#' Description: Calculate Wilcoxon test (unpaired; BH correction) for the specified sample groups. 
#' Either provide the input data as matrix, file path, or select the file through GUI.
#'             
#' Arguments:
#'   @param dat data matrix (features x samples)
#'   @param file data file  (if data matrix not provided) 
#'   @param p.adjust.method p-value correction method for p.adjust function (default "BH"). If NULL, no correction will be performed.
#'   @param sort sort the results
#' Returns:
#'   @return (Corrected) p-values for two-group comparison.
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

check.wilcoxon <- function (dat = NULL, fnam = NULL, p.adjust.method = "BH", sort = FALSE) {

  require(svDialogs)

  ## Open your tab fnam, Level 1 & 2 Sum_BGsub_Rel.contribution

  if (is.null(dat) && is.null(fnam)) { fnam <- tk_choose.files(multi = F) }
  if (is.null(dat)) {
    dat <- read.table(fnam, sep = "\t", header = T, row.names = 1)
  } 

  samples <- colnames(dat)
  levels <- rownames(dat)

  ## To select samples you can do that in 2 ways: select G1 and those samples not in G1 are G2, for which you would use:
  ## G2 <- samples[!(samples %in% G1)]
  ## or select G1 and select G2 (this is useful when you have multiple groups) (standard below)

  G1 <- tk_select.list(samples, multiple = T, title = "Select samples for 1st group")
  G2 <- tk_select.list(samples, multiple = T, title = "Select samples for 2nd group")

  M <- matrix(data = NA, length(levels), 1)
  rownames(M) <- levels

  for (i in 1:length(levels)) {
	
    lvl  <- levels[i]
    l.g1 <- dat[lvl,G1]
    l.g2 <- dat[lvl,G2]
	
    p <- wilcox.test(as.numeric(l.g1), as.numeric(l.g2))$p.value

    message(lvl, " p-value: ", p, "\n")

    M[i, 1] <- p

  }

  ## To Adjust P-values for Multiple Comparisons with Benjamini & Hochberg (1995) 
  ## ("BH" or its alias "fdr")
  if (!is.null(p.adjust.method)) {

    cor.p <- p.adjust(M, method = p.adjust.method) 
    names(cor.p) <- rownames(M)

  } else {

    # Skip p-value correction
    cor.p <- as.vector(M)
    names(cor.p) <- rownames(M)

  }

  # Sort the values
  if (sort) { cor.p <- sort(cor.p) }

  cor.p

}




#' Description: Cross-correlate input variables
#'              
#' Arguments:
#'   @param annot annotation matrix: samples x features
#'   @param dat numeric data matrix: samples x features
#'   @param method association method (pearson, spearman for continuous; categorical for discrete
#'   @param qth q-value threshold for included features 
#'   @param cth correlation threshold
#'   @param order order the results
#'   @param n.signif mininum number of significant correlations for each element
#'   @param mode Specify the output format ("table" or "matrix")
#'
#' Returns:
#'   @return List with cor, pval, qval
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

cross.correlate <- function(annot, dat, method = "pearson", qth = NULL, cth = NULL, order = FALSE, n.signif = 0, mode = "table"){

  # annot <- metadata.df; dat <- t(genus.matrix); method = "pearson"; qth = NULL; cth = NULL; order = FALSE; n.signif = 0; verbose = TRUE; mode = "matrix"

  # annot <- meta[sample.set,]; dat <- t(ds[, sample.set]); method = cor.method; qth = NULL; cth = NULL; order = FALSE; n.signif = 0; mode = "table"

  x <- as.data.frame(annot) # numeric or discrete
  y <- dat # numeric

  if (is.null(colnames(y))) { colnames(y) <- paste("column-", 1:ncol(y), sep = "") }

  xnames <- colnames(x)
  ynames <- colnames(y)
  qv <- NULL

  numeric.methods <- c("spearman", "pearson", "bicor", "mi")
  categorical.methods <- c("categorical")

  # Rows paired.
  if (method %in% numeric.methods) {
    inds <- sapply(x, is.numeric) 
    if (any(!inds)) {
      warning("Considering only numeric annotations for pearson/spearman/bicor/mi")
    }
    inds <- names(which(inds))
  } else if (method %in% categorical.methods) {
    inds <- sapply(x, is.factor) 
    if (any(!inds)) {
      warning("Considering only categorical annotations for factors")
    }
    inds <- names(which(inds))
  }

  xnames <- inds
  x <- as.matrix(x[inds], ncol = length(inds))
  colnames(x) <- xnames

  Pc <- matrix(NA, ncol(x), ncol(y))
  Cc <- matrix(NA, ncol(x), ncol(y))
  rownames(Cc) <- colnames(x)
  colnames(Cc) <- colnames(y)
  rownames(Pc) <- colnames(x)
  colnames(Pc) <- colnames(y)

  # ----------------------------------------------------------------

  if (method %in% c("pearson", "spearman")) {

    for (j in 1:ncol(y)){
      jc <- apply(x, 2, function (xi) { 
        if (sum(!is.na(xi)) > 10) {
          res <- cor.test(xi, y[, j], method = method, use = "pairwise.complete.obs"); 
	  res <- c(res$estimate, res$p.value)	   
	} else {
	  warning("Not enough observations; skipping correlation estimation")
	  res <- c(NA, NA)
        }
	res
      })
  
      Cc[,j] <- jc[1,]        
      Pc[,j] <- jc[2,]        

    } 

  } else if (method == "bicor") {

    require(WGCNA)
    t1 <- bicorAndPvalue(x, y, use = "pairwise.complete.obs")
    Pc <- t1$p
    Cc <- t1$bicor

  } else if (method == "categorical") {  

      Cc <- matrix(NA, nrow = ncol(x), ncol = ncol(y))
      rownames(Cc) <- colnames(x)
      colnames(Cc) <- colnames(y)

      for (varname in colnames(x)) {

        for (lev in colnames(y)) {

          xvec <- x[, varname]
	  yvec <- y[, lev]
	  keep <- rowSums(is.na(cbind(xvec, yvec))) == 0	       
	  xvec <- xvec[keep]
	  yvec <- yvec[keep]
	
	  # Number of data-annotation samples for calculating the correlations
    	  n <- sum(keep)

    	  Cc[varname, lev] <- GKtau(xvec, yvec)

        }
      }
   } else if (method == "mi") {

      Cc <- matrix(NA, nrow = ncol(x), ncol = ncol(y))
      rownames(Cc) <- colnames(x)
      colnames(Cc) <- colnames(y)

      for (i in 1:ncol(x)) {
        for (j in 1:ncol(y)) {

          Cc[i,j] <- build.mim(cbind(x[,i], y[,j]), estimator = "spearman")[1, 2]

        }
      }
   }

   if (!all(is.na(Pc))) {

     rownames(Pc) <- xnames
     colnames(Pc) <- ynames

     rownames(Cc) <- xnames
     colnames(Cc) <- ynames

     # Corrected p-values
     if ((prod(dim(Pc)) - sum(is.na(Pc))) >= 100) {
       qv <- matrix.qvalue(Pc)
     } else {
       warning("Too few p-values available, q-value calculation skipped-")
       qv <- array(NA, dim = dim(Pc))
     }

  }

   # Filter
   if (!is.null(qth) || !is.null(cth)) {

     # Replace NAs with extreme values for filtering purposes
     qv[is.na(qv)] <- 1
     Pc[is.na(qv)] <- 1
     Cc[is.na(Cc)] <- 0

     # Filter by qvalues and correlations
     inds1.q <- inds2.q <- inds1.c <- inds2.c <- NULL

     if (!is.null(qth)) {
       inds1.q <- apply(abs(qv), 1, function(x) {sum(x < qth) >= n.signif}) 
       inds2.q <- apply(abs(qv), 2, function(x) {sum(x < qth) >= n.signif}) 
     }

     if (!is.null(cth)) {
       inds1.c <- apply(abs(Cc), 1, function(x) {sum(x > cth) >= n.signif})
       inds2.c <- apply(abs(Cc), 2, function(x) {sum(x > cth) >= n.signif})
     }

     if (!is.null(qth) && !is.null(cth)) {
       inds1 <- inds1.q & inds1.c
       inds2 <- inds2.q & inds2.c
     } else if (is.null(qth) && !is.null(cth)) {
       inds1 <- inds1.c
       inds2 <- inds2.c
     } else if (!is.null(qth) && is.null(cth)) {
       inds1 <- inds1.q
       inds2 <- inds2.q
     }

     Cmat <- as.matrix(0)
     # TODO: add also correlation filter, not only significance
     # Require each has at least n.signif. correlations
     if (sum(inds1) >= n.signif && sum(inds2) >=n.signif) {

       rnams <- rownames(Cc)[inds1]
       cnams <- colnames(Cc)[inds2]

       Cc <- matrix(Cc[inds1,inds2, drop = FALSE], nrow = sum(inds1))
       Pc <- matrix(Pc[inds1,inds2, drop = FALSE], nrow = sum(inds1))
       qv <- matrix(qv[inds1,inds2, drop = FALSE], nrow = sum(inds1))

       rownames(qv) <- rownames(Pc) <- rownames(Cc) <- rnams
       colnames(qv) <- colnames(Pc) <- colnames(Cc) <- cnams

       if ( order && sum(inds1) >=2 && sum(inds2)>=2 ) { # Order in visually appealing order

         tmp <- Cc
         rownames(tmp) <- NULL
         colnames(tmp) <- NULL
         h <- heatmap(tmp, xlab = NULL, ylab = NULL, xaxt = 'n', yaxt = 'n'); dev.off()

         rnams <- rownames(Cc)[h$rowInd]
         cnams <- colnames(Cc)[h$colInd]
         Cc <- Cc[h$rowInd,h$colInd]
         Pc <- Pc[h$rowInd,h$colInd]
         qv <- qv[h$rowInd,h$colInd]

         rownames(qv) <- rownames(Pc) <- rownames(Cc) <- rnams
         colnames(qv) <- colnames(Pc) <- colnames(Cc) <- cnams

       }

     } else {
       message("No significant correlations with the given criteria\n")
       Cc <- Pc <- qv <- NULL
    }
   }

   res <- list(cor = Cc, pval = Pc, qval = qv)

   if (mode == "matrix") {
     return(res)     
   } else if (mode == "table") {
     ctab <- NULL

     if (!is.null(res$cor)) {
       ctab <- melt(res$cor)
       colnames(ctab) <- c("X1", "X2", "correlation")
     }

     if (!is.null(res$qval)) {
       ctab <- cbind(ctab, melt(res$qval)$value)
       colnames(ctab) <- c("X1", "X2", "correlation", "qvalue")
       ctab <- esort(ctab, qvalue, -abs(correlation))
       colnames(ctab) <- c("X1", "X2", method, "qvalue")
     } else {
       ctab <- esort(ctab, -abs(correlation))
       colnames(ctab) <- c("X1", "X2", method)
     }

     return(ctab)
   }
}

#' Description: Stability analysis. Calculates average Pearson '
#  correlation between samples in the input data and picks the lower '
#  triangular matrix to avoid duplicating the correlations. Returns 
#  correlations and stability estimate (average of the correlations).
#'
#' Arguments:
#'   @param dat data matrix phylotypes vs. samples
#'
#' Returns:
#'   @return List with correlations and astability estimate
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

calculate.stability <- function (dat) {
  cors <- lower.triangle(cor(dat))
  list(correlations = cors, stability = mean(cors))
}


