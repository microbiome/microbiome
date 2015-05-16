#' Arrange correlation matrices from cross.correlate into a table format
#'              
#' @param res Output from cross.correlate
#' @param verbose verbose
#'
#' @return Correlation table
#'
#' @export
#'
#' @examples 
#'   data(peerj32)
#'   d1 <- peerj32$microbes[1:20, 1:10]
#'   d2 <- peerj32$lipids[1:20,1:10]
#'   cc <- cross.correlate(d1, d2, mode = 'matrix')
#'   cmat <- cmat2table(cc)
#'
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
        ctab <- ctab[order(ctab$p.adj), ]
        colnames(ctab) <- c("X1", "X2", "Correlation", "p.adj")
        
    } else {
        message("No significant adjusted p-values")
        if (!is.null(ctab)) {
            ctab <- cbind(ctab, melt(res$pval)$value)
            ctab <- ctab[order(-abs(ctab$Correlation)), ]
            colnames(ctab) <- c("X1", "X2", "Correlation", "pvalue")
        }
    }
    
    ctab$X1 <- as.character(ctab$X1)
    ctab$X2 <- as.character(ctab$X2)

    # Keep the original order of factor levels
    ctab$X1 <- factor(as.character(ctab$X1), levels = rownames(res$cor))
    ctab$X2 <- factor(as.character(ctab$X2), levels = colnames(res$cor))

    # Remove NAs
    ctab <- ctab[!is.na(ctab$Correlation),]

    # Order the table by p-value
    if ("p.adj" %in% colnames(ctab)) {
      ctab <- ctab[order(ctab$p.adj), ]
    } else if ("pvalue" %in% colnames(ctab)) {
      ctab <- ctab[order(ctab$pvalue), ]
    }

    ctab
    
}



#' List color scales
#'
#' 
#'
#' 
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
#' Calculate hierarchical clustering for standard selections in 
#' profiling script
#'
#' 
#'   @param dat data matrix (use log10 with pearson!)
#'   @param method hierarchical clustering method (see ?hclust)
#'   @param metric clustering metrics (spearman / pearson / euclidean)
#'
#' 
#'   @return hclust object for log10 and for absolute scale data
#'
#' @export
#' @examples 
#'   data(peerj32)
#'   dat <- peerj32$microbes
#'   hc <- calculate.hclust(dat, 'complete', 'pearson')
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

#' get probeset data matrix
#' 
#' 
#'   @param name name
#'   @param level taxonomic level
#'   @param phylogeny.info phylogeny.info
#'   @param probedata oligos vs. samples preprocessed data matrix; 
#'                    absolute scale
#'   @param log10 Logical. Logarithmize the data TRUE/FALSE
#'
#' 
#'   @return probeset data matrix
#'
#' @export
#' @examples 
#'   phylogeny.info <- GetPhylogeny('HITChip', 'filtered')
#'   data.dir <- system.file("extdata", package = "microbiome")
#'   probedata <- read.profiling("frpa", data.dir = data.dir)$oligo
#'   ps <- get.probeset('Akkermansia', 'L2', phylogeny.info, probedata)
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
get.probeset <- function(name, level, phylogeny.info, probedata, log10 = TRUE) {
    
    # Pick probes for this entity
    probes <- retrieve.probesets(phylogeny.info, level, name)
    
    sets <- vector(length = length(probes), mode = "list")
    names(sets) <- names(probes)
    
    for (nam in names(probes)) {
        
        # Pick expression for particular probes (absolute scale)
        p <- intersect(probes[[nam]], rownames(probedata))
        dat <- NULL
        if (length(p) > 0) {
            dat <- probedata[p, ]
            
            dat <- matrix(dat, nrow = length(probes[[nam]]))
            rownames(dat) <- probes[[nam]]
            colnames(dat) <- colnames(probedata)
            
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
#' @examples 
#'   data(peerj32)
#'   ratios <- PhylotypeRatios(peerj32$microbes)
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
#'   JD Storey 2003. Ann. Statist. 31(6):2013-2035. The positive false 
#'   discovery rate: a Bayesian interpretation and the q-value. 
#'
#'   To cite the microbiome R package, see citation('microbiome')
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




#' Impute missing values from a Gaussian. 
#' 
#' @param X data matrix (features x samples)
#'
#' @return imputed data matrix
#' @export 
#' @references
#' See citation('microbiome')
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples 
#'   data(peerj32)
#'   x <- peerj32$microbes
#'   xi <- impute(x) 
#' @keywords utilities

impute <- function(X) {

  # TODO Replace with standard R functions
    
    for (i in 1:ncol(X)) {
        x <- X[, i]
        nas <- is.na(x)
        X[nas, i] <- rnorm(sum(nas), mean(x[!is.na(x)]), sd(x[!is.na(x)]))
    }
    
    X
    
}


#' Get lower triangle of a square matrix 
#' as a numeric vector such that
#' row-by-row, picking elements in the order
#' 2,1;3,1;3,2;4,1,...
#'        
#'   @param mat data matrix
#'
#'   @return lower triangle as vector 
#'
#' @export
#' @examples 
#'   mat <- rbind(c(1,2,3), c(4,5,6), c(7,8,9))
#'   vec <- lower.triangle(mat)
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
lower.triangle <- function(mat) {
    
  # TODO is this easy replace with standard R functions ?

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
