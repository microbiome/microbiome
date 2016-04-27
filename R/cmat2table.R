#' @title Convert cross correlation results to a table
#' @description Arrange correlation matrices from cross.correlate into a table format.
#' @param res Output from cross.correlate
#' @param verbose verbose
#' @return Correlation table
#' @export
#' @examples 
#'   data(peerj32)
#'   d1 <- peerj32$microbes[1:20, 1:10]
#'   d2 <- peerj32$lipids[1:20,1:10]
#'   cc <- cross.correlate(d1, d2, mode = 'matrix')
#'   cmat <- cmat2table(cc)
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
cmat2table <- function(res, verbose = FALSE) {
    
    ctab <- ID <- NULL
    
    if (!is.null(res$cor)) {
  	ctab <- as.data.frame(res$cor)
  	ctab$ID <- rownames(res$cor)
  	ctab <- gather(ctab, ID)    
        colnames(ctab) <- c("X1", "X2", "Correlation")
  	ctab$Correlation <- as.numeric(as.character(ctab$Correlation))  	
    }
    
    correlation <- NULL  # circumwent warning on globabl vars
    
    if (!is.null(res$p.adj)) {
        
        if (verbose) {
            message("Arranging the table")
        }

	ctab2 <- as.data.frame(res$p.adj)
  	ctab2$ID <- rownames(res$p.adj)
  	ctab2 <- gather(ctab2, ID)    
        colnames(ctab2) <- c("X1", "X2", "p.adj")
  	ctab2$p.adj <- as.numeric(as.character(ctab2$p.adj))  	

        ctab <- cbind(ctab, ctab2$p.adj)
        colnames(ctab) <- c("X1", "X2", "Correlation", "p.adj")
        ctab <- ctab[order(ctab$p.adj), ]
        colnames(ctab) <- c("X1", "X2", "Correlation", "p.adj")
        
    } else {
        message("No significant adjusted p-values")
        if (!is.null(ctab)) {
	
	    ctab2 <- as.data.frame(res$pval)
  	    ctab2$ID <- rownames(res$pval)
  	    ctab2 <- gather(ctab2, ID)    
            colnames(ctab2) <- c("X1", "X2", "value")
  	    ctab2$value <- as.numeric(as.character(ctab2$value))
	
            ctab <- cbind(ctab, ctab2$value)
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





