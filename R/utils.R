# Copyright (C) 2011-2014 Leo Lahti and Jarkko Salojarvi Contact:
# <microbiome-admin@googlegroups.com>. All rights reserved.

# This file is a part of the microbiome R package
# http://microbiome.github.com/

# This program is open source software; you can redistribute it and/or
# modify it under the terms of the FreeBSD License (keep this notice):
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.



#' Description: Arrange correlation matrices from cross.correlate into 
#'         a table format
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
#' @examples data(peerj32); 
#'          cc <- cross.correlate(peerj32$microbes[1:20, 1:10], 
#'                               peerj32$lipids[1:20,1:10], 
#'                   mode = 'matrix'); 
#'                    cmat <- cmat2table(cc)
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities


cmat2table <- function(res, verbose = FALSE) {
    
    ctab <- NULL
    
    if (!is.null(res$cor)) {
        ctab <- melt(res$cor)
        colnames(ctab) <- c("X1", "X2", "Correlation")
    }
    
    correlation <- NULL  # circumwent warning on globabl vars
    
    if (!is.null(res$p.adj)) {
        
        if (verbose) {
            message("Arranging the table")
        }
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



#' Description: List color scales
#'
#' Arguments:
#'
#' Returns:
#'   @return list of color scales
#'
#' @export
#' @examples list.color.scales()
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

list.color.scales <- function() {
    ## Different colour scales
    list(`white/blue` = colorRampPalette(c("white", "darkblue"), 
         interpolate = "linear")(100), 
         `white/black` = colorRampPalette(c("white", "black"), 
         interpolate = "linear")(100), 
        `black/yellow/white` = colorRampPalette(c("black", "yellow", "white"), 
         bias = 0.5, interpolate = "linear")(100))
}




#' calculate.hclust
#' 
#' Description: Calculate hierarchical clustering for standard selections in 
#' profiling script
#'
#' Arguments:
#'   @param dat data matrix (use log10 with pearson!)
#'   @param method hierarchical clustering method (see ?hclust)
#'   @param metric clustering metrics (spearman / pearson / euclidean)
#'
#' Returns:
#'   @return hclust object for log10 and for absolute scale data
#'
#' @export
#' @examples data(peerj32); 
#'            dat <- peerj32$microbes;
#'           hc <- calculate.hclust(dat, 'complete', 'pearson') 
#' @references See citation('microbiome')
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

calculate.hclust <- function(dat, method = "complete", metric = "pearson") {
    
    if (metric == "euclidean") {
        hc <- hclust(dist(t(dat)), method = method)
    } else if (metric %in% c("spearman", "pearson")) {
        hc <- hclust(as.dist(1 - cor(dat, use = "complete.obs", 
                     method = metric)), 
            method = method)
    } else {
        stop("Provide proper metric for calculate.hclust!")
    }
    
    hc
    
}


#' Description: get probeset data matrix
#' 
#' Arguments:
#'   @param name name
#'   @param level taxonomic level
#'   @param phylogeny.info phylogeny.info
#'   @param oligo.matrix oligos vs. samples preprocessed data matrix; 
#'                    absolute scale
#'   @param log10 Logical. Logarithmize the data TRUE/FALSE
#'
#' Returns:
#'   @return probeset data matrix
#'
#' @export
#' @examples 
#'      phylogeny.info <- GetPhylogeny('HITChip', 'filtered')
#'      # ps <- get.probeset('Vibrio', 'L2', phylogeny.info, oligo.matrix)
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

get.probeset <- function(name, level, phylogeny.info, oligo.matrix, 
                         log10 = TRUE) {
    
    # Pick probes for this entity
    probes <- retrieve.probesets(phylogeny.info, level, name)
    
    sets <- vector(length = length(probes), mode = "list")
    names(sets) <- names(probes)
    
    for (nam in names(probes)) {
        
        # Pick expression for particular probes (absolute scale)
        p <- intersect(probes[[nam]], rownames(oligo.matrix))
        dat <- NULL
        if (length(p) > 0) {
            dat <- oligo.matrix[p, ]
            
            dat <- matrix(dat, nrow = length(probes[[nam]]))
            rownames(dat) <- probes[[nam]]
            colnames(dat) <- colnames(oligo.matrix)
            
            # Logarithmize probeset?
            if (log10) {
                dat <- log10(dat)
            }
        }
        sets[[nam]] <- dat
        
    }
    
    if (length(sets) == 1) {
        sets <- sets[[1]]
    }
    
    # Return
    sets
    
}





#' PhylotypeRatios
#'
#' Calculate phylotype ratios (eg. Bacteroides vs. Prevotella etc.) 
#'          for a given phylotypes vs. samples matrix
#'
#' @param dat phylotypes vs. samples data matrix in log10 scale
#'
#' @return phylotype pairs x samples matrix indicating the ratio 
#'                 (in log10 domain) between each unique pair
#' @export 
#' @examples data(peerj32); 
#'          dat <- peerj32$microbes; 
#'         ratios <- PhylotypeRatios(dat)
#' @references
#' See citation('microbiome')
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

PhylotypeRatios <- function(dat) {
    
    phylogroups <- rownames(dat)
    Nratios <- (length(phylogroups)^2 - length(phylogroups))/2
    Nsamples <- ncol(dat)
    ratios <- list()
    for (i in 1:(length(phylogroups) - 1)) {
        for (j in (i + 1):length(phylogroups)) {
            pt1 <- phylogroups[[i]]
            pt2 <- phylogroups[[j]]
            ratios[[paste(pt1, pt2, sep = "-")]] <- dat[pt1, ] - dat[pt2, ]
        }
    }
    ratios <- do.call(cbind, ratios)
    
    t(ratios)
}



#' matrix.padjust
#'
#' Calculate adjusted p-values for a matrix of pvalues 
#' which may contain missing values.
#' @param pvals p-value matrix
#' @param p.adjust.method p-value adjustment method: for options, see ?p.adjust
#' @return Adjusted p-value matrix
#' @export 
#' @references 
#'    JD Storey 2003. Ann. Statist. 31(6):2013-2035. 
#'    The positive false discovery rate: 
#'          a Bayesian interpretation and the q-value. 
#'    To cite the microbiome R package, see citation('microbiome')
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples qvals <- matrix.padjust(matrix(runif(1000), nrow = 100))
#' @keywords utilities

matrix.padjust <- function(pvals, p.adjust.method = "BH") {
    
    pvec <- as.vector(pvals)
    nai <- is.na(pvec)
    qvec <- rep(NA, length(pvec))
    qvec[!nai] <- p.adjust(pvec[!nai], method = p.adjust.method)
    qmat <- matrix(qvec, nrow = nrow(pvals))
    dimnames(qmat) <- dimnames(pvals)
    qmat
    
}

#' polish.phylogeny.info
#'
#' Ensure phylogeny.info is in correct format
#' 
#' @param phylogeny.info phylogeny.info data frame
#'
#' @return polished phylogeny.info
#' @export 
#' @references See citation('microbiome')
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples phylogeny.info <- GetPhylogeny('HITChip', 'filtered');
#'           phylogeny.info <- polish.phylogeny.info(phylogeny.info)
#' @keywords utilities

polish.phylogeny.info <- function(phylogeny.info) {
    
    colnames(phylogeny.info)[
          which(colnames(phylogeny.info) == "level.0")] <- "L0"
    colnames(phylogeny.info)[
          which(colnames(phylogeny.info) == "level.1")] <- "L1"
    colnames(phylogeny.info)[
          which(colnames(phylogeny.info) == "level.2")] <- "L2"
    
    colnames(phylogeny.info)[
          which(colnames(phylogeny.info) == "level 0")] <- "L0"
    colnames(phylogeny.info)[
          which(colnames(phylogeny.info) == "level 1")] <- "L1"
    colnames(phylogeny.info)[
          which(colnames(phylogeny.info) == "level 2")] <- "L2"
    
    phylogeny.info
    
}


#' Impute missing values from a Gaussian. 
#' 
#' @param X data matrix (features x samples)
#'
#' @return imputed data matrix
#' @export 
#' @references
#' See citation('microbiome')
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples data(peerj32)
#'          x <- peerj32$microbes
#'          xi <- impute(x) 
#' @keywords utilities

impute <- function(X) {
    
    for (i in 1:ncol(X)) {
        x <- X[, i]
        nas <- is.na(x)
        X[nas, i] <- rnorm(sum(nas), mean(x[!is.na(x)]), sd(x[!is.na(x)]))
    }
    
    X
    
}



#' Strip string i.e. remove spaces from the beginning and end
#' @param s string or character vector
#'
#' @return Stripped string
#' @export 
#' @references
#' See citation('microbiome')
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples Strip(' aa b c ') 
#' @keywords utilities
Strip <- function(s) {
    
    ss <- c()
    
    for (i in 1:length(s)) {
        
        si <- s[[i]]
        if (!is.na(si)) {
            # Strip string i.e. remove spaces from the beginning and end
            while (substr(si, 1, 1) == " ") {
                si <- substr(si, 2, nchar(si))
            }
            while (substr(si, nchar(si), nchar(si)) == " ") {
                si <- substr(si, 1, nchar(si) - 1)
            }
        }
        ss[[i]] <- si
    }
    
    ss
}



#' Description: Sort data frame dd by columns like: esort(dd, -z, b)
#'
#' Arguments:
#'   @param x data frame to sort
#'   @param sortvar sorted variable/s
#'   @param ... further parameters to pass
#'
#' Returns:
#'   @return sorted data frame
#'
#' @export
#' @examples data(peerj32)
#'          esort(peerj32$meta, -time, gender)
#'
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

esort <- function(x, sortvar, ...) {
    
    attach(x, warn.conflicts = FALSE)
    x <- x[with(x, order(sortvar, ...)), ]
    return(x)
    detach(x)
}




#' Description: 
#' Get lower triangle of a square matrix 
#' as a numeric vector such that
#' row-by-row, picking elements in the order
#' 2,1;3,1;3,2;4,1,...
#'        
#' Arguments:
#'   @param mat data matrix
#' Returns:
#'   @return lower triangle as vector 
#'
#' @export
#' @examples 
#'          mat <- rbind(c(1,2,3), c(4,5,6), c(7,8,9))
#'          vec <- lower.triangle(mat)
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
lower.triangle <- function(mat) {
    
    elements <- c()
    nr <- dim(mat)[[1]]
    nc <- dim(mat)[[2]]
    
    for (i in 2:nr) {
        for (j in 1:(i - 1)) {
            elements <- c(elements, mat[i, j])
        }
    }
    elements
} 
