# Copyright (C) 2011-2014 Leo Lahti and Jarkko Salojarvi 
# Contact: <microbiome-admin@googlegroups.com>. All rights reserved.

# This file is a part of the microbiome R package
# http://microbiome.github.com/

# This program is open source software; you can redistribute it and/or
# modify it under the terms of the FreeBSD License (keep this notice):
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' Description: Cross-hybridization table between multimodal taxa as percentages of shared probes. 
#' The number indicates how many percent of oligos for the row taxon are also hybridizing 
#' the corresponding column taxon.
#'
#' Arguments:
#'   @param tax.level Taxonomic level to investigate
#'   @param chip Chip type (e.g. "HITChip")
#'   @param selected.taxa Restrict cross-hyb analysis to the selected groups.
#'   @param phylogeny.info phylogeny.info 
#'
#' Returns:
#'   @return A list containing cross-hybridization table 
#'
#' @examples ch <- CrosshybTable(tax.level = "L1")
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

CrosshybTable <- function (tax.level = "L1", chip = "HITChip", selected.taxa = NULL, phylogeny.info = NULL) {

  # Get hylogeny info
  if (is.null(phylogeny.info)) {
    phylogeny.info <- GetPhylogeny(chip, phylogeny.version = "filtered")
  }

  # Pick necessary columns
  phi <- phylogeny.info[, c(tax.level, "oligoID")]

  # Include only selected groups (if any)
  if (!is.null(selected.taxa)) {
    phi <- phi[phi[[tax.level]] %in% selected.taxa,]
  }

  # Create taxon-oligo mapping matrix
  tax.oligos <- sapply(split(phi, phi[[tax.level]]), function (x) {x$oligoID})
  tax2oligo <- matrix(0, nrow = length(unique(phi[[tax.level]])), ncol = length(unique(phi$oligoID)))
  rownames(tax2oligo) <- unique(phi[[tax.level]])
  colnames(tax2oligo) <- unique(phi$oligoID)
  for (tax in names(tax.oligos)) {
    oligos <- tax.oligos[[tax]] 
    tax2oligo[tax, oligos] <- 1
  }

  # Confusion matrix: how many overlapping oligos between two taxa
  confusion.matrix <- matrix(NA, nrow = length(tax.oligos), ncol = length(tax.oligos))
  rownames(confusion.matrix) <- colnames(confusion.matrix) <- names(tax.oligos)
  for (tax1 in rownames(tax2oligo)) {
    for (tax2 in rownames(tax2oligo)) {
      to <- tax2oligo[c(tax1, tax2),] 

      confusion.matrix[tax1, tax2] <- mean(to[tax2, to[tax1, ] == 1])
      confusion.matrix[tax2, tax1] <- mean(to[tax1, to[tax2, ] == 1])

    }
  }

  confusion.matrix

}



#' Description: Cross-hybridization between multimodal taxa as percentages of shared probes. 
#' The number indicates how many percent of oligos for the row taxon are also hybridizing 
#' the corresponding column taxon.
#'
#' Arguments:
#'   @param tax.level Taxonomic level to investigate
#'   @param chip Chip type (e.g. "HITChip")
#'   @param selected.taxa Restrict cross-hyb analysis to the selected groups.
#'   @param show.plot Produce the plot
#'   @param order.rows Order table rows
#'   @param order.cols Order table columns
#'   @param keep.empty Keep taxa that do not show any cross-hybridization
#'   @param rounding Rounding of the cell contents
#'   @param phylogeny.info phylogeny.info 
#'   @param self.correlations Show self-correlations (always 100%); or remove (indicate as 0%; default)
#'
#' Returns:
#'   @return A list containing cross-hybridization table and plot
#'
#' @examples res <- PlotCrosshyb(tax.level = "L1", rounding = 1, show.plot = FALSE)
#' @export
#' @import ggplot2
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

PlotCrosshyb <- function (tax.level = "L1", chip = "HITChip", selected.taxa = NULL, show.plot = TRUE, order.rows = TRUE, order.cols = TRUE, keep.empty = FALSE, rounding = 1, phylogeny.info = NULL, self.correlations = FALSE) {

  # tax.level = "L1"; chip = "HITChip"; selected.taxa = NULL; show.plot = TRUE; order.rows = TRUE; order.cols = TRUE; keep.empty = FALSE; rounding = 1

  # Get crosshyb matrix
  confusion.matrix <- CrosshybTable(tax.level = tax.level, chip = "HITChip", selected.taxa = NULL, phylogeny.info = NULL)

  # Remove self-correlations
  if (!self.correlations) {
    diag(confusion.matrix) <- 0
  }

  # Focus on selected taxa
  if (!is.null(selected.taxa)) {
    confusion.matrix <- confusion.matrix[selected.taxa, ]
  }

  # Remove the taxa that do not have any crosshyb
  if (!keep.empty) {
    confusion.matrix <- confusion.matrix[rowSums(confusion.matrix) > 0, colSums(confusion.matrix) > 0]
  }

  # Avoid warnings
  Taxon1 <- Taxon2 <- crosshyb <- NULL

  # Organize into data frame
  df <- reshape::melt(confusion.matrix)
  names(df) <- c("Taxon1", "Taxon2", "crosshyb")

  # Switch to percentages
  df[["crosshyb"]] <- 100*df[["crosshyb"]]

  # Order rows and cols
  if (order.rows || order.cols) {

    hc <- hclust(as.dist(1-cor(confusion.matrix)), "ward")
    colord <- hc$ord

    hc <- hclust(as.dist(1-cor(t(confusion.matrix))), "ward")
    roword <- hc$ord

    if (order.rows) {
      df[["Taxon1"]] <- factor(df[["Taxon1"]], levels = rownames(confusion.matrix)[roword])
    }

    if (order.cols) {  
      df[["Taxon2"]] <- factor(df[["Taxon2"]], levels = colnames(confusion.matrix)[colord])
    }
  }

  # Visualize
  theme_set(theme_bw(15))
  df$labels <- round(df[["crosshyb"]], rounding)
    
  aes <- NULL
  p <- ggplot(df, aes(Taxon2, Taxon1, group = Taxon1)) 
  p <- p + ggplot2::geom_tile(aes(fill = crosshyb)) 
  p <- p + ggplot2::geom_text(aes(fill = crosshyb, label = labels, size = 4))
  p <- p + scale_fill_gradient(low = "white", high = "red") 
  p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)) + ggplot2::xlab("") + ggplot2::ylab("") 
  p <- p + ggplot2::theme(legend.position = "none")

  if (show.plot) {
    print(p)
  }

  list(data = df[, 1:3], plot = p)

}




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
#' @examples data(peerj32); d <- distance.matrix(peerj32$microbes[1:10, 1:3])
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




#' Description: Calculate ANOVA test (BH correction) for multi-group comparison. 
#' Arguments:
#'   @param dat data matrix (features x samples; eg. HITChip taxa vs. samples)
#'   @param group Vector with specifying the groups
#'   @param p.adjust.method p-value correction method for p.adjust function (default "BH"). If NULL, no correction will be performed.
#'   @param sort sort the results
#'
#' Returns:
#'   @return (Corrected) p-values for multi-group comparison.
#'
#' @examples data(peerj32); pval <- check.anova(t(peerj32$microbes), peerj32$meta$time)
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
check.anova <- function (dat, group, p.adjust.method = "BH", sort = FALSE) {

  if (is.null(p.adjust.method)) {p.adjust.method <- "none"}

  pvals <- c()
  for (tax in rownames(dat)) {
    pvals[[tax]] <- anova(lm(dat[tax, ] ~ group))["group", "Pr(>F)"]
  }

  pvals <- p.adjust(pvals, method = p.adjust.method)

  if (sort) {pvals <- sort(pvals)}

  pvals
}


#' Description: Calculate Wilcoxon test (unpaired; BH correction) for the specified sample groups. 
#' Either provide the input data as matrix, file path, or select the file through GUI.
#'             
#' Arguments:
#'   @param dat data matrix (features x samples)
#'   @param fnam data file  (if data matrix not provided) 
#'   @param G1 Sample group 1 (for comparison) 
#'   @param G2 Sample group 2 (for comparison)
#'   @param p.adjust.method p-value correction method for p.adjust function (default "BH"). If NULL, no correction will be performed.
#'   @param sort sort the results
#'
#' Returns:
#'   @return (Corrected) p-values for two-group comparison.
#'
#' @examples data(peerj32); pval <- check.wilcoxon(t(peerj32$microbes), G1 = 1:22, G2 = 23:44)
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

check.wilcoxon <- function (dat = NULL, fnam = NULL, G1, G2, p.adjust.method = "BH", sort = FALSE) {

  ## Open your tab fnam, Level 1 & 2 Sum_BGsub_Rel.contribution

  if (is.null(dat) && is.null(fnam)) { stop("Provide dat or fnam in function arguments!")}
  if (is.null(dat)) {
    dat <- read.table(fnam, sep = "\t", header = T, row.names = 1)
  } 

  samples <- colnames(dat)
  levels <- rownames(dat)

  M <- matrix(data = NA, length(levels), 1)
  rownames(M) <- levels

  for (i in 1:length(levels)) {
	
    lvl  <- levels[i]
    l.g1 <- dat[lvl,G1]
    l.g2 <- dat[lvl,G2]
	
    p <- wilcox.test(as.numeric(l.g1), as.numeric(l.g2))$p.value

    # message(lvl, " p-value: ", p, "\n")

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


#' Description: Cross-correlate columns of the input matrices
#'              
#' Arguments:
#'   @param x matrix (samples x features if annotation matrix)
#'   @param y matrix (samples x features if cross-correlated with annotations)
#'   @param method association method ('pearson', 'spearman', or 'bicor' for continuous; categorical for discrete)
#'   @param p.adj.threshold q-value threshold to include features 
#'   @param cth correlation threshold to include features 
#'   @param order order the results
#'   @param n.signif mininum number of significant correlations for each element
#'   @param mode Specify output format ("table" or "matrix")
#'   @param p.adj.method p-value multiple testing correction method. Either "qvalue" or one of the methods in p.adjust function ("BH" and others; see help(p.adjust)). Default: "fdr"
#'   @param verbose verbose
#'   @param filter.self.correlations Filter out correlations between identical items.
#'
#' Returns:
#'   @return List with cor, pval, pval.adjusted
#'
#' @examples data(peerj32); cc <- cross.correlate(peerj32$microbes[1:20, 1:10], peerj32$lipids[1:20,1:10])
#' @export
#'
#' @details As the method=categorical (discrete) association measure
#'          for nominal (no order for levels) variables ' we using Goodman and
#'          Kruskal tau based on
#'          http://www.r-bloggers.com/measuring-associations-between-non-numeric-variables/
#'  	    The 'bicor' method is from the WGCNA package.
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

cross.correlate <- function(x, y = NULL, method = "pearson", p.adj.threshold = Inf, cth = NULL, order = FALSE, n.signif = 0, mode = "table", p.adj.method = "fdr", verbose = F, filter.self.correlations = F) {

  if (is.null(y)) {
    message("Cross-correlating the data with itself")
    y <- x

    if (filter.self.correlations) {
      # Ignore self-correlations in filtering
      n.signif <- n.signif + 1 
    }
  }		

  if (verbose) {message("Polishing the data")}

  x <- as.data.frame(x) # numeric or discrete
  y <- y # numeric

  if ( is.null(colnames(y)) ) { colnames(y) <- paste("column-", 1:ncol(y), sep = "") }

  xnames <- colnames(x)
  ynames <- colnames(y)
  qv <- NULL

  numeric.methods <- c("spearman", "pearson", "bicor", "mi")
  categorical.methods <- c("categorical")

  if ( verbose ) { message("Methods") }

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

    if (verbose) {message(method)}

    for (j in 1:ncol(y)){
      jc <- apply(x, 2, function (xi) { 
        if (sum(!is.na(xi)) >= 8) {
          res <- cor.test(xi, y[, j], method = method, use = "pairwise.complete.obs"); 
	  res <- c(res$estimate, res$p.value)	   
	} else {
	  warning(paste("Not enough observations; (", sum(!is.na(xi)), ") - skipping correlation estimation"))
	  res <- c(NA, NA)
        }
	res
      })
  
      Cc[,j] <- jc[1,]        
      Pc[,j] <- jc[2,]        

    } 

  } else if (method == "bicor") {

    if (verbose) {message(method)}

    InstallMarginal("WGCNA")

    t1 <- WGCNA::bicorAndPvalue(x, y, use = "pairwise.complete.obs")
    Pc <- t1$p
    Cc <- t1$bicor

  } else if (method == "categorical") {  

      if (verbose) {message(method)}

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

    	  Cc[varname, lev] <- GKtau(xvec, yvec) # 

        }
      }
   } else if (method == "mi") {

      if (verbose) {message(method)}
    
      InstallMarginal("minet")

      Cc <- matrix(NA, nrow = ncol(x), ncol = ncol(y))
      rownames(Cc) <- colnames(x)
      colnames(Cc) <- colnames(y)

      for (i in 1:ncol(x)) {
        for (j in 1:ncol(y)) {

          Cc[i,j] <- minet::build.mim(cbind(x[,i], y[,j]), estimator = "spearman")[1, 2]

        }
      }
   }

  if (!all(is.na(Pc))) {

     if (verbose) {message("p adjustment")}

     rownames(Pc) <- xnames
     colnames(Pc) <- ynames

     rownames(Cc) <- xnames
     colnames(Cc) <- ynames

     # Corrected p-values
     qv <- array(NA, dim = dim(Pc))

     if (p.adj.method == "qvalue") {
       if (((prod(dim(Pc)) - sum(is.na(Pc))) >= 100)) {
         qv <- matrix.qvalue(Pc)
       } else {
         warning("Not enough p-values for qvalue calculation - q-value calculation skipped. Try BH method instead for multiple correction?")
       }
     } else {
       if (verbose) {message(paste("Multiple testing correction with", p.adj.method))}

       qv <- matrix(p.adjust(Pc, method = p.adj.method), nrow = nrow(Pc))
       dimnames(qv) <- dimnames(Pc)

     }
  }

  if (verbose) {message("OK")}

   # Filter
   if (!is.null(p.adj.threshold) || !is.null(cth)) {

     if (verbose) {message("Filtering")}

     # Replace NAs with extreme values for filtering purposes
     qv[is.na(qv)] <- 1
     Pc[is.na(qv)] <- 1
     Cc[is.na(Cc)] <- 0

     # Filter by qvalues and correlations
     inds1.q <- inds2.q <- inds1.c <- inds2.c <- NULL

     if (!is.null(p.adj.threshold)) {
       inds1.q <- apply(qv, 1, function(x) {sum(x < p.adj.threshold) >= n.signif}) 
       inds2.q <- apply(qv, 2, function(x) {sum(x < p.adj.threshold) >= n.signif}) 
     }

     if (!is.null(cth)) {
       inds1.c <- apply(abs(Cc), 1, function(x) { sum(x > cth) >= n.signif })
       inds2.c <- apply(abs(Cc), 2, function(x) { sum(x > cth) >= n.signif })
     }

     if (!is.null(p.adj.threshold) && !is.null(cth)) {

       inds1 <- inds1.q & inds1.c
       inds2 <- inds2.q & inds2.c

     } else if (is.null(p.adj.threshold) && !is.null(cth)) {
       inds1 <- inds1.c
       inds2 <- inds2.c
     } else if (!is.null(p.adj.threshold) && is.null(cth)) {
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


   res <- list(cor = Cc, pval = Pc, p.adj = qv)

   if (all(as.vector(x) == as.vector(y)) && filter.self.correlations){ 
     message("Ignore self-correlations in filtering")
     diag(res$cor) <- NA
     diag(res$pval) <- NA
     diag(res$p.adj) <- NA
   }

   if (mode == "matrix") {
     return(res)     
   } else if (mode == "table") {

     tab <- cmat2table(res)
     tab$X1 <- factor(tab$X1, levels = rownames(res$cor))
     tab$X2 <- factor(tab$X2, levels = colnames(res$cor))

     if (order) {

       message("Ordering factors")
       tab$X1 <- factor(as.character(tab$X1), levels = rownames(res$cor))
       tab$X2 <- factor(as.character(tab$X2), levels = colnames(res$cor))

     } 

     if (all(as.vector(x) == as.vector(y)) && filter.self.correlations) {
       # Remove self-correlations
       tab <- tab[!(tab$X1 == tab$X2),]
     }

     if ("p.adj" %in% colnames(tab)) {
       tab <- tab[order(tab$p.adj), ]
     } else if ("pvalue" %in% colnames(tab)) {
       tab <- tab[order(tab$pvalue), ]
     }

     return(tab)
   }
}


#' Description: Arrange correlation matrices from cross.correlate into a table format
#'              
#' Arguments:
#'   @param res Output from cross.correlate
#'   @param verbose verbose
#'
#' Returns:
#'   @return Correlation table
#'
#' @export
#'
#' @examples data(peerj32); cc <- cross.correlate(peerj32$microbes[1:20, 1:10], peerj32$lipids[1:20,1:10], mode = "matrix"); cmat <- cmat2table(cc)
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities


cmat2table <- function (res, verbose = FALSE) {

     ctab <- NULL

     if (!is.null(res$cor)) {
       ctab <- melt(res$cor)
       colnames(ctab) <- c("X1", "X2", "Correlation")
     }

     correlation <- NULL # circumwent warning on globabl vars

     if (!is.null(res$p.adj)) {

       if (verbose) {message("Arranging the table")}
       ctab <- cbind(ctab, melt(res$p.adj)$value)
       colnames(ctab) <- c("X1", "X2", "Correlation", "p.adj")
       ctab <- esort(ctab, ctab$p.adj, -abs(ctab$Correlation))
       colnames(ctab) <- c("X1", "X2", "Correlation", "p.adj")

     } else {
       message("No significant adjusted p-values")
       if (!is.null(ctab)) {
         ctab <- cbind(ctab, melt(res$pval)$value)
         ctab <- esort(ctab, -abs(ctab$Correlation))
         colnames(ctab) <- c("X1", "X2", "Correlation", "pvalue")
       }
     }

     ctab$X1 <- as.character(ctab$X1)
     ctab$X2 <- as.character(ctab$X2)

     ctab

}



#' Description: Stability analysis. Calculates average Pearson '
#'  correlation between samples in the input data and picks the lower '
#'  triangular matrix to avoid duplicating the correlations. Returns 
#'  correlations and stability estimate (average of the correlations). 
#'  Can also be used to calculate stability between two data sets. 
#'  Then provide two data sets as inputs.
#'
#' 
#' Arguments:
#'   @param dat1 data matrix phylotypes vs. samples (in log10 scale)
#'   @param dat2 Optional. Second data matrix phylotypes vs. samples. 
#'          Provide this to calculate stability between two (paired) 
#'          data sets.
#'   @param method Correlation method (see ?cor)
#'
#' Returns:
#'   @return List with correlations and astability estimate
#'
#' @export
#' @examples data(peerj32); s <- estimate.stability(t(peerj32$microbes)[, 1:5])
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

estimate.stability <- function (dat1, dat2 = NULL, method = "pearson") {

  if (is.null(dat2)) {
    # Within-matrix stability
    # NOTE: earlier this was calculated as 
    # the average of upper triangular correlation matrix
    # This is heavily biased since the values are dependent
    # Now replaced by calculating correlations against the
    # mean of the whole sample set
    #cors <- lower.triangle(cor(dat1))
    cors <- as.vector(cor(dat1, matrix(rowMeans(dat1)), method = method))
    names(cors) <- colnames(dat1)
    stab <- list(correlations = cors, stability = mean(cors))
  } else {
    # Between-matrices stability
    cors <- diag(cor(dat1, dat2, method = method))
    stab <- list(correlations = cors, stability = mean(cors))
  }

  stab

}



#' Description: Measure association between nominal (no order for levels) variables 
#' using Goodman and Kruskal tau. Code modified from the original source:
#' http://www.r-bloggers.com/measuring-associations-between-non-numeric-variables/
#' An important feature of this procedure is that it allows missing
#' values in either of the variables x or y, treating 'missing' as an
#' additional level.  In practice, this is sometimes very important since
#' missing values in one variable may be strongly associated with either
#' missing values in another variable or specific non-missing levels of
#' that variable. An important characteristic of Goodman and Kruskal's tau measure is
#' its asymmetry: because the variables x and y enter this expression
#' differently, the value of a(y,x) is not the same as the value of
#' a(x, y), in general.  This stands in marked contrast to either the
#' product-moment correlation coefficient or the Spearman rank
#' correlation coefficient, which are both symmetric, giving the same
#' association between x and y as that between y and x.  The fundamental
#' reason for the asymmetry of the general class of measures defined
#' above is that they quantify the extent to which the variable x is
#' useful in predicting y, which may be very different than the extent to
#' which the variable y is useful in predicting x.
#'
#' Arguments:
#'   @param x first variable
#'   @param y second variable
#'
#' Returns:
#'   @return Dependency measure
#'
#' @examples data(peerj32); tc <- GKtau(unlist(peerj32$microbes[,1]), unlist(peerj32$lipids[,1]))
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

GKtau <- function(x,y){

      #
      #  First, compute the IxJ contingency table between x and y
      #
      Nij = table(x,y,useNA="ifany")
      #
      #  Next, convert this table into a joint probability estimate
      #
      PIij = Nij/sum(Nij)
      #
      #  Compute the marginal probability estimates
      #
      PIiPlus = apply(PIij,MARGIN=1,sum)
      PIPlusj = apply(PIij,MARGIN=2,sum)
      #
      #  Compute the marginal variation of y
      #
      Vy = 1 - sum(PIPlusj^2)
      #
      #  Compute the expected conditional variation of y given x
      #
      InnerSum = apply(PIij^2,MARGIN=1,sum)
      VyBarx = 1 - sum(InnerSum/PIiPlus)
      #
      #  Compute and return Goodman and Kruskal's tau measure
      #
      tau = (Vy - VyBarx)/Vy
      tau
}

