# Copyright (C) 2011-2012 Leo Lahti and Jarkko Salojarvi 
# Contact: <leo.lahti@iki.fi>. All rights reserved.

# This file is a part of the microbiome R package

# This program is open source software; you can redistribute it and/or
# modify it under the terms of the FreeBSD License (keep this notice):
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' Description: Correlation distance. Computes the correlation 
#' distance between the _columns_ of the given matrix
#'
#' Arguments:
#'   @param x data matrix
#'
#' Returns:
#'   @return distance object
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

corDist <- function(x){
  as.dist((1-cor(x, use = "pairwise.complete")))
}


#' Description: pairwise comparison between factor levels
#'              
#' Arguments:
#'   @param x data matrix: samples x features
#'   @param y factor levels for the samples
#'   @param qth q-value threshold for included features 
#'   @param resdir Optional. If given, output files will be written in the result directory
#'
#' Returns:
#'   @return List with top findings from pairwise comparisons and their q-values
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

pairwise.comparisons <- function (x, y, qth = 0.05, resdir = NULL) {

  y <- factor(y)

  library(limma)

  # Remove categories with NA annotations and very small sample size
  keep <- !is.na(y) & (y %in% names(which(table(y) > 2)))
  y <- y[keep]
  x <- x[keep, ]

  top.findings <- list()
  top.findings.qvals <- list()

  # Compare each group to other groups
  for (i in 1:(length(unique(y))-1)) {

    qval.list <- list()

    for (j in (i+1):(length(unique(y)))) {
   
      ai <- as.character(unique(y)[[i]])
      aj <- as.character(unique(y)[[j]])

      inds1 <- which(y == ai)
      inds2 <- which(y == aj)

      lab <- factor(y[c(inds1, inds2)])
      ann <- data.frame(list(varname = lab))
      dat <- x[c(inds1, inds2),]

      # Here y refers to internal lm.matrix response variable(!)
      qvals <- lm.matrix(y ~ varname,  dat, ann, type = "qval")
      coefs <- lm.matrix(y ~ varname - 1,  dat, ann, type = "coef") # remove intercept

      nams <- names(sort(qvals))
      
      fc <- colMeans(x[inds1, nams]) - colMeans(x[inds2, nams])
    
      tab <- cbind(nams, coefs[nams, paste("varname", ai, sep = "")], coefs[nams, paste("varname", aj, sep = "")])
      colnames(tab) <- c("phylotype", ai, aj)
      tab <- as.data.frame(tab)
      tab$fold.change <- fc[nams]
      tab$qvalue <- qvals[nams]

      top.findings[[paste(ai, aj, sep = "-")]] <- subset(tab, qvalue < qth)
      qval.list[[aj]] <- qvals
      top.findings.qvals[[ai]] <- qval.list

      if (!is.null(resdir)) {

        resfile <- paste(resdir,"/limma-", ai, "-", aj, "-results.tab", sep = "")
        message(paste("Writing comparison statistics in: ",  resfile))
        write.table(tab, file = resfile, sep = "\t", quote = FALSE, row.names = FALSE)

      if (nrow(tab) > 0) {
    
        df <- as.data.frame(t(rbind(group = as.character(lab), t(dat)[rownames(tab)[1:min(nrow(tab), 9)],])))
        df$group <- factor(df$group)
        dfm <- reshape::melt.data.frame(df, id.vars = c("group"))
        dfm$value <- as.numeric(as.character(dfm$value))
        dfm$title <- apply(cbind(as.character(dfm$variable), round(as.numeric(as.character(tab$qvalue[match(dfm$variable, names(tab$qvalue))])),3)),1,function(x){paste(x, collapse = "/ qval: ")})
        p <- ggplot(dfm, aes(group, value)) + geom_boxplot() + facet_wrap(~title) + opts(title = "Top hits")
        pdfname <- paste(resdir,"/limma-", ai, "-", aj, "-tophit.boxplots.pdf", sep = "")
        print(pdfname)
        pdf(pdfname)
        print(p)
        dev.off()
  
      }
     }
    }
  }

  list(qval = top.findings.qvals, top = top.findings)

}


check.wilcoxon <- function () {

  require(svDialogs)

  ## Open your tab file, Level 1&2 Sum_BGsub_Rel.contribution

  file <- choose.files(multi=F)
  oligo <- read.table(file,sep="\t", header=T, row.names=1)
  samples <- colnames(oligo)

  ## To select samples you can do that in 2 ways: select G1 and those samples not in G1 are G2, for which you would use:
  ## G2 <- samples[!(samples %in% G1)]
  ## or select G1 and select G2 (this is useful when you have multiple groups) (standard below)

  G1 <- tk_select.list(samples, multiple=T, title="Select samples for 1st group")
  G2 <- tk_select.list(samples, multiple=T, title="Select samples for 2nd group")

  levels <- rownames(oligo)

  M <- matrix(data=NA,length(levels),1)
  rownames(M) <- levels

  for (i in 1:length(levels)) {
	
	lvl <- levels[i]
	l.g1 <- oligo[lvl,G1]
	l.g2 <- oligo[lvl,G2]
	
	p<-wilcox.test(as.numeric(l.g1),as.numeric(l.g2))$p.value

	

	cat(lvl," p-value: ",p, "\n")

	M[i,1]<-p
	}

  M

  ## To Adjust P-values for Multiple Comparisons with Benjamini & Hochberg (1995) ("BH" or its alias "fdr")
  cor.p <- p.adjust(M,method="BH") # Other methods can be applied here
  cor.p
}




#' Description: Cross-correlate input variables
#'              
#' Arguments:
#'   @param annot annotation matrix: samples x features
#'   @param dat data matrix: samples x features
#'   @param method association method (pearson, spearman for continuous; categorical for discrete
#'   @param qth q-value threshold for included features 
#'   @param cth correlation threshold
#'   @param order order the results
#'   @param n.signif mininum number of significant correlations for each element
#'   @param verbose verbose
#'
#' Returns:
#'   @return List with cor, pval, qval, N
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

cross.correlate <- function(annot, dat, method = "pearson", qth = NULL, cth = NULL, order = FALSE, n.signif = 0, verbose = TRUE, mode = "list"){

  x <- annot # numeric or discrete
  y <- dat # numeric

  if (is.null(colnames(y))) {colnames(y) <- paste("column-", 1:ncol(y), sep = "")}

  xnames <- colnames(x)
  ynames <- colnames(y)
  qv <- nmat <- NULL

  # Rows paired.

  Pc <- matrix(NA, ncol(x), ncol(y))
  Cc <- matrix(NA, ncol(x), ncol(y))
  rownames(Cc) <- colnames(x)
  colnames(Cc) <- colnames(y)
  rownames(Pc) <- colnames(x)
  colnames(Pc) <- colnames(y)


  # Ensure the correct formats

  if (ncol(x) == 1) {
    x <- matrix(t(apply(x, 1, as.numeric)), nrow = length(x))
  } else {
    if (!method == "categorical") {
      x <- matrix(t(apply(x, 1, as.numeric)), nrow(x))
    } else {
      x <- matrix(t(apply(x, 1, as.factor)), nrow(x))
    }
  }
  colnames(x) <- xnames

  if (ncol(y) == 1) {
    y <- matrix(t(apply(y, 1, as.numeric)), nrow = length(y))
  } else {
      y <- matrix(t(apply(y, 1, as.numeric)), nrow(y))
  }
  colnames(y) <- ynames
  
  # ----------------------------------------------------------------

  if (method %in% c("pearson","spearman")) {

    for (j in 1:ncol(y)){
      if (verbose) {print(j/ncol(y))}
      jc <- apply(x, 2, function (xi) {res <- cor.test(xi, y[, j], method = method, use = "pairwise.complete.obs"); c(res$estimate, res$p.value)})

      Pc[,j] <- jc[2,]        
      Cc[,j] <- jc[1,]        

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
      nmat <- Pc <- Cc

      for (varname in colnames(x)) {

        for (lev in colnames(y)) {

          xvec <- x[, varname]
	  yvec <- y[, lev]
	  keep <- rowSums(is.na(cbind(xvec,yvec))) == 0	       
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
      nmat <- Pc <- Cc

      for (i in 1:ncol(x)) {
        for (j in 1:ncol(y)) {

          Cc[i,j] <- build.mim(cbind(x[,i], y[,j]), estimator = "spearman")[1, 2]
        }
      }
   }

   if (!all(is.na(Pc))) {

     # Corrected p-values
     require(qvalue)
     if (prod(dim(Pc)) >= 100) {
       qv <- qvalue(Pc)
       if (length(qv) == 1) { 
         qv <- NULL 
       } else {
         qv <- qv$qvalue
         rownames(qv) <- xnames
         colnames(qv) <- ynames
       }
     } else {
       cat("Too few p-values available, q-value calculation skipped-")
       qv <- NULL
     }

     rownames(Pc) <- xnames
     colnames(Pc) <- ynames

     rownames(Cc) <- xnames
     colnames(Cc) <- ynames
  }

   # Filter
   if (!is.null(qth) || !is.null(cth)) {

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
       cat("No significant correlations with the given criteria\n")
       Cc <- Pc <- qv <- NULL
    }
   }

   res <- list(cor = Cc, pval = Pc, qval = qv, N = nmat)

   if (mode == "list") {
     return(res)     
   } else if (mode == "table") {
     ctab <- cbind(melt(res$cor), melt(res$qval)$value)
     colnames(ctab) <- c("X1", "X2", "correlation", "qvalue")
     o <- order(ctab$qvalue)
     ctab <- ctab[o, ]
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
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

calculate.stability <- function (dat) {
  cors <- lower.triangle(cor(dat))
  list(correlations = cors, stability = mean(cors))
}

