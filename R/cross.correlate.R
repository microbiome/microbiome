
#' Description: Cross-correlate columns of the input matrices
#'              
#' Arguments:
#'   @param x matrix (samples x features if annotation matrix)
#'   @param y matrix (samples x features if cross-correlated with annotations)
#'   @param method association method ('pearson', 'spearman', or 'bicor' 
#'                  for continuous; categorical for discrete)
#'   @param p.adj.threshold q-value threshold to include features 
#'   @param cth correlation threshold to include features 
#'   @param order order the results
#'   @param n.signif mininum number of significant correlations for each 
#'                    element
#'   @param mode Specify output format ('table' or 'matrix')
#'   @param p.adj.method p-value multiple testing correction method. 
#'                    One of the methods in p.adjust 
#'             function ('BH' and others; see help(p.adjust)). 
#'             Default: 'fdr'
#'   @param verbose verbose
#'   @param filter.self.correlations Filter out correlations between 
#'                            identical items.
#'
#' Returns:
#'   @return List with cor, pval, pval.adjusted
#'
#' @examples data(peerj32); 
#'           cc <- cross.correlate(peerj32$microbes[1:20, 1:10], 
#'                           peerj32$lipids[1:20,1:10])
#' @export
#' @import WGCNA
#'
#' @details As the method=categorical (discrete) association measure
#'          for nominal (no order for levels) variables we use Goodman and
#'          Kruskal tau based on
#' r-bloggers.com/measuring-associations-between-non-numeric-variables/
#'          The 'bicor' method is from the WGCNA package.
#'
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

cross.correlate <- function(x, y = NULL, method = "pearson", 
                            p.adj.threshold = Inf, 
    cth = NULL, order = FALSE, n.signif = 0, mode = "table", 
    p.adj.method = "fdr", 
    verbose = FALSE, filter.self.correlations = FALSE) {
    
    if (is.null(y)) {
        message("Cross-correlating the data with itself")
        y <- x
        
        if (filter.self.correlations) {
            # Ignore self-correlations in filtering
            n.signif <- n.signif + 1
        }
    }
    
    x <- as.data.frame(x)  # numeric or discrete
    y <- y  # numeric
    
    if (is.null(colnames(y))) {
        colnames(y) <- paste("column-", 1:ncol(y), sep = "")
    }
    
    xnames <- colnames(x)
    ynames <- colnames(y)
    qv <- NULL    
    numeric.methods <- c("spearman", "pearson", "bicor")
    categorical.methods <- c("categorical")
    
    # Rows paired.
    if (method %in% numeric.methods) {
        inds <- sapply(x, is.numeric)
        if (any(!inds)) {
            warning("Considering only numeric annotations for \n       
                     pearson/spearman/bicor/mi")
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
    
    if (method %in% c("pearson", "spearman")) {
        
        for (j in 1:ncol(y)) {
            jc <- apply(x, 2, function(xi) {
                if (sum(!is.na(xi)) >= 8) {
                  res <- cor.test(xi, y[, j], method = method, 
                            use = "pairwise.complete.obs")
                  res <- c(res$estimate, res$p.value)
                } else {
                  warning(paste("Not enough observations; \n   
                          (",  
                    sum(!is.na(xi)), ") \n \n 
                           - skipping correlation estimation"))
                  res <- c(NA, NA)
                }
                res
            })
            
            Cc[, j] <- jc[1, ]
            Pc[, j] <- jc[2, ]
            
        }
       
    } else if (method == "bicor") {
        
        if (verbose) {
            message(method)
        }
        
        t1 <- bicorAndPvalue(x, y, use = "pairwise.complete.obs")
        Pc <- t1$p
        Cc <- t1$bicor
        
    } else if (method == "categorical") {
        
        if (verbose) {
            message(method)
        }
        
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
                
                # Number of data-annotation samples for calculating
                # the correlations
                n <- sum(keep)
                
                Cc[varname, lev] <- GKtau(xvec, yvec)  # 
                
            }
        }
    }

    
    if (!all(is.na(Pc))) {
        
        if (verbose) {
            message("p adjustment")
        }
        
        rownames(Pc) <- xnames
        colnames(Pc) <- ynames
        
        rownames(Cc) <- xnames
        colnames(Cc) <- ynames
        
        # Corrected p-values
        qv <- array(NA, dim = dim(Pc))            
        qv <- matrix(p.adjust(Pc, method = p.adj.method), nrow = nrow(Pc))
        dimnames(qv) <- dimnames(Pc)    
    }
    
    # Filter
    if (!is.null(p.adj.threshold) || !is.null(cth)) {
        
        # Replace NAs with extreme values for filtering purposes
        qv[is.na(qv)] <- 1
        Pc[is.na(qv)] <- 1
        Cc[is.na(Cc)] <- 0
        
        # Filter by adjusted pvalues and correlations
        inds1.q <- inds2.q <- inds1.c <- inds2.c <- NULL
        
        if (!is.null(p.adj.threshold)) {
            inds1.q <- apply(qv, 1, function(x) {
                sum(x < p.adj.threshold) >= n.signif
            })
            inds2.q <- apply(qv, 2, function(x) {
                sum(x < p.adj.threshold) >= n.signif
            })
        }
        
        if (!is.null(cth)) {
            inds1.c <- apply(abs(Cc), 1, function(x) {
                sum(x > cth) >= n.signif
            })
            inds2.c <- apply(abs(Cc), 2, function(x) {
                sum(x > cth) >= n.signif
            })
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

        if (sum(inds1) >= n.signif && sum(inds2) >= n.signif) {
            
            rnams <- rownames(Cc)[inds1]
            cnams <- colnames(Cc)[inds2]
            
            Cc <- matrix(Cc[inds1, inds2, drop = FALSE], nrow = sum(inds1))
            Pc <- matrix(Pc[inds1, inds2, drop = FALSE], nrow = sum(inds1))
            qv <- matrix(qv[inds1, inds2, drop = FALSE], nrow = sum(inds1))
            
            rownames(qv) <- rownames(Pc) <- rownames(Cc) <- rnams
            colnames(qv) <- colnames(Pc) <- colnames(Cc) <- cnams
            
            if (order && sum(inds1) >= 2 && sum(inds2) >= 2) {
                # Order in visually appealing order
                tmp <- Cc
                rownames(tmp) <- NULL
                colnames(tmp) <- NULL
                
                #h <- heatmap(tmp, xlab = NULL, 
                #   ylab = NULL, xaxt = "n", yaxt = "n")                
                #rind <- h$rowInd
		#cind <- h$colInd

	    	rind <- hclust(as.dist(1-cor(t(tmp), use = "pairwise.complete.obs")))$order
	    	cind <- hclust(as.dist(1-cor(tmp, use = "pairwise.complete.obs")))$order
                rnams <- rownames(Cc)[rind]
                cnams <- colnames(Cc)[cind]
                Cc <- Cc[rind, cind]
                Pc <- Pc[rind, cind]
                qv <- qv[rind, cind]
                
                rownames(qv) <- rownames(Pc) <- rownames(Cc) <- rnams
                colnames(qv) <- colnames(Pc) <- colnames(Cc) <- cnams
                
            }
            
        } else {
            message("No significant correlations with the given criteria\n")
            Cc <- Pc <- qv <- NULL
        }
    }
        
    res <- list(cor = Cc, pval = Pc, p.adj = qv)
    
    if (all(as.vector(x) == as.vector(y)) && filter.self.correlations) {
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
            tab <- tab[!(tab$X1 == tab$X2), ]
        }
        
        if ("p.adj" %in% colnames(tab)) {
            tab <- tab[order(tab$p.adj), ]
        } else if ("pvalue" %in% colnames(tab)) {
            tab <- tab[order(tab$pvalue), ]
        }
        
        return(tab)
    }
}


