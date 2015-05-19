#' Description: Check number of matching phylotypes for each probe
#' 
#' Arguments:
#'   @param taxonomy oligo - phylotype matching data.frame
#'   @param level phylotype level
#'
#' Returns:
#'   @return number of matching phylotypes for each probe
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
n.phylotypes.per.oligo <- function (taxonomy, level) {
  sapply(split(as.vector(taxonomy[, level]), as.vector(taxonomy[, "oligoID"])), function(x) length(unique(x)))
}


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
