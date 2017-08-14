#' @title project.data
#' @description Project high-dimensional data on two-dimensional plane
#'              by various methods
#' 
#' @param amat data matrix (samples x features)
#' @param type projection type 
#'           (options: PCA, MDS.classical, MDS.nonmetric, Sammon)
#'
#' @return projected data matrix
#'
#' @export
#' @importFrom MASS isoMDS
#' @importFrom MASS sammon
#' @importFrom mixOmics spca
#'
#' @examples 
#'   data(peerj32)
#'   xy <- project.data(peerj32$microbes[,1:3])
#'
#' @references 
#'    
#'    D. Hand and H. Mannila and P. Smyth: 
#'    Principles of Data Mining. MIT Press. Cambridge, MA, US (2001).
#'    
#'    To cite microbiome R package, see citation('microbiome') 
#'
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

project.data <- function(amat, type = "PCA") {
    
    if (type == "PCA") {
        if (nrow(amat) < ncol(amat)) {
            
            message("More samples than features, using sparse PCA")
            
            ## Spca example: we are selecting 50 variables on each of the PCs
            result <- spca(amat, ncomp = 2, center = TRUE, scale = TRUE, 
                           keepX = rep(50, 2))
            scores <- result$x
        } else {
            message("PCA")
            pca <- princomp(amat)  # Classical PCA
            scores <- pca$scores
        }
        tab <- data.frame(scores[, 1:2])
        rownames(tab) <- rownames(amat)
    } else if (type == "Sammon") {
        
        d <- as.dist(1 - cor(t(amat)))
        # This gave the clearest visualization.  
        # Tuning magic parameter could still
        # improve.  Try for instance magic = 0.05.
        fit <- sammon(d, k = 2)
        # Plot solution
        tab <- data.frame(list(Comp.1 = fit$points[, 1], 
                               Comp.2 = fit$points[, 2]))
        rownames(tab) <- rownames(amat)
    } else if (type == "MDS.classical") {
        d <- as.dist(1 - cor(t(amat)))
        fit <- cmdscale(d, eig = TRUE, k = 2)  # classical MDS
        tab <- data.frame(list(Comp.1 = fit$points[, 1], 
                               Comp.2 = fit$points[, 2]))
    } else if (type == "MDS.nonmetric") {
        d <- as.dist(1 - cor(t(amat)))
        fit <- isoMDS(d, k = 2)  # nonmetric MDS
        tab <- data.frame(list(Comp.1 = fit$points[, 1], 
                               Comp.2 = fit$points[, 2]))
    }
    
    # TODO Kernel-PCA kpc <- kpca(~., data=as.data.frame(x.train),
    # kernel='rbfdot', features = 2) Print the principal component
    # vectors pcv(kpc) Plot the data projection on the components
    # par(mfrow=c(2,2)) plot(rotated(kpc), col =
    # as.integer(as.factor(ann[rownames(x.train),'time'])), xlab='1st
    # Principal Component', ylab='2nd Principal Comp onent')
    # plot(rotated(kpc), col =
    # as.integer(as.factor(ann[rownames(x.train),'lipids.group'])),
    # xlab='1st Principal Component', ylab='2nd Principal Component')
    # embed remaining points emb <- predict(kpc, x.test)
    # plot(rotated(kpc), col =
    # as.integer(as.factor(ann[rownames(x.train),'lipids.group'])),
    # xlab='1st Principal Component', ylab='2nd Principal Component')
    # points(emb, col =
    # as.integer(as.factor(ann[rownames(x.train),'lipids.group'])))
    
    colnames(tab) <- c("Comp.1", "Comp.2")
    
    tab
}