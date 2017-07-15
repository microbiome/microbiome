#' @importFrom WGCNA bicorAndPvalue

#' @title Cross Correlation Wrapper
#' @description Cross-correlate columns of the input matrices.
#' @param x matrix (samples x features if annotation matrix)
#' @param y matrix (samples x features if cross-correlated with annotations)
#' @param method association method ('pearson', 'spearman', or 'bicor' 
#' for continuous; categorical for discrete)
#' @param p.adj.threshold q-value threshold to include features 
#' @param cth correlation threshold to include features 
#' @param order order the results
#' @param n.signif mininum number of significant correlations for each 
#' element
#' @param mode Specify output format ('table' or 'matrix')
#' @param p.adj.method p-value multiple testing correction method. 
#' One of the methods in p.adjust 
#' function ('BH' and others; see help(p.adjust)). Default: 'fdr'
#' @param verbose verbose
#' @param filter.self.correlations Filter out correlations between 
#' identical items.
#' @return List with cor, pval, pval.adjusted
#' @examples 
#' data(peerj32)
#' d1 <- peerj32$microbes[1:20, 1:10]
#' d2 <- peerj32$lipids[1:20,1:10]
#' cc <- associate(d1, d2, method='pearson')
#' @export
#' @details As the method=categorical (discrete) association measure
#' for nominal (no order for levels) variables we use Goodman and
#' Kruskal tau based on
#' r-bloggers.com/measuring-associations-between-non-numeric-variables/
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @aliases cross_correlate
#' @keywords utilities
associate <- function(x, y=NULL,
method="spearman", p.adj.threshold=Inf,
cth=NULL, order=FALSE, n.signif=0, mode="table",
p.adj.method="fdr",
verbose=FALSE, filter.self.correlations=FALSE) {
    
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
        colnames(y) <- paste("column-", 1:ncol(y), sep="")
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
    
    if (!is.vector(x)) {
        x <- suppressWarnings(as.matrix(x[, inds], ncol=length(inds)))
    } else {
        x <- as.matrix(x[inds], ncol=length(inds))
    }
    
    colnames(x) <- xnames
    
    Pc <- matrix(NA, ncol(x), ncol(y))
    Cc <- matrix(NA, ncol(x), ncol(y))
    rownames(Cc) <- colnames(x)
    colnames(Cc) <- colnames(y)
    rownames(Pc) <- colnames(x)
    colnames(Pc) <- colnames(y)
    
    if (verbose) {
        message(method)
    }
    
    if (method %in% c("pearson", "spearman")) {
        
        minobs <- 8
        
        for (j in 1:ncol(y)) {
            
            jc <- apply(x, 2, function(xi) {
                
                if (sum(!is.na(xi)) >= minobs) {

                    res <- suppressWarnings(
                        cor.test(xi, unlist(y[, j], use.names=FALSE), 
                    method=method, use="pairwise.complete.obs"))

                    res <- c(res$estimate, res$p.value)

                } else {
                
                    warning(paste("Not enough observations (",
                minobs, "required); \n   
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
        
        t1 <- suppressWarnings(
        bicorAndPvalue(x, y, use="pairwise.complete.obs"))
        Pc <- t1$p
        Cc <- t1$bicor
        
        
    } else if (method == "categorical") {
        
        if (verbose) {
            message(method)
        }
        
        Cc <- matrix(NA, nrow=ncol(x), ncol=ncol(y))
        rownames(Cc) <- colnames(x)
        colnames(Cc) <- colnames(y)
        
        for (varname in colnames(x)) {
            
            for (lev in colnames(y)) {
                
                xvec <- x[, varname]
                yvec <- y[, lev]
                keep <- rowSums(is.na(cbind(xvec, yvec))) == 0
                xvec <- xvec[keep]
                yvec <- yvec[keep]
                
                # Number of data-annotation samples for
                # calculating the correlations
                n <- sum(keep)
                Cc[varname, lev] <- GKtau(xvec, yvec) 
                
            }
        }
    }
    
    if (!all(is.na(Pc))) {
        
        rownames(Pc) <- xnames
        colnames(Pc) <- ynames
        
        rownames(Cc) <- xnames
        colnames(Cc) <- ynames
        
        # Corrected p-values
        qv <- array(NA, dim=dim(Pc))
        qv <- matrix(p.adjust(Pc, method=p.adj.method), nrow=nrow(Pc))
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
            
            Cc <- matrix(Cc[inds1, inds2, drop=FALSE], nrow=sum(inds1))
            Pc <- matrix(Pc[inds1, inds2, drop=FALSE], nrow=sum(inds1))
            qv <- matrix(qv[inds1, inds2, drop=FALSE], nrow=sum(inds1))
            
            rownames(qv) <- rownames(Pc) <- rownames(Cc) <- rnams
            colnames(qv) <- colnames(Pc) <- colnames(Cc) <- cnams
            
            if (order && sum(inds1) >= 2 && sum(inds2) >= 2) {
                
                # Order in visually appealing order
                tmp <- Cc
                rownames(tmp) <- NULL
                colnames(tmp) <- NULL
                
                rind <- hclust(as.dist(1 - cor(t(tmp),
            use="pairwise.complete.obs")))$order
                cind <- hclust(as.dist(1 - cor(tmp,
            use="pairwise.complete.obs")))$order
                
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
    
    res <- list(cor=Cc, pval=Pc, p.adj=qv)
    
    # message('Ignore self-correlations in filtering')
    
    if (nrow(x) == nrow(y) && ncol(x) == ncol(y) && filter.self.correlations) {
        diag(res$cor) <- diag(res$pval) <- diag(res$p.adj) <- NA
    }
    
    if (mode == "table") {
        res <- cmat2table(res)
    }
    
    res
    
}


#' @title Convert Correlation Matrix into a Table
#' @description Arrange correlation matrices from associate into a table format.
#' @param res Output from associate
#' @param verbose verbose
#' @return Correlation table
#' @examples 
#' data(peerj32)
#' d1 <- peerj32$microbes[1:20, 1:10]
#' d2 <- peerj32$lipids[1:20,1:10]
#' cc <- associate(d1, d2, mode='matrix', method='pearson')
#' cmat <- associate(d1, d2, mode='table', method='spearman')
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
cmat2table <- function(res, verbose=FALSE) {
    
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
    ctab$X1 <- factor(as.character(ctab$X1), levels=rownames(res$cor))
    ctab$X2 <- factor(as.character(ctab$X2), levels=colnames(res$cor))
    
    # Remove NAs
    ctab <- ctab[!is.na(ctab$Correlation), ]
    
    # Order the table by p-value
    if ("p.adj" %in% colnames(ctab)) {
        ctab <- ctab[order(ctab$p.adj), ]
    } else if ("pvalue" %in% colnames(ctab)) {
        ctab <- ctab[order(ctab$pvalue), ]
    }
    
    ctab
    
}


