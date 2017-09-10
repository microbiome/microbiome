#' @title Plot RDA
#' @description Visualize rda_bagged output.
#' @param x Output from rda_bagged
#' @param which.bac TBA
#' @param ptype Plot type. 'spider' or 'hull'
#' @param comp TBA
#' @param cex.bac Plot size.
#' @param plot.names Plot names
#' @param group.cols Group colors.
#' @param ... Other arguments to be passed
#' @return TBA
#' @examples
#' data(peerj32)
#' x <- t(peerj32$microbes)
#' y <- factor(peerj32$meta$time); names(y) <- rownames(peerj32$meta)
#' # Use more iterations in practice
#' res <- rda_bagged(x, y, bs.iter=5)
#' tmp <- plot_rda_bagged(res)
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
plot_rda_bagged <- function(x, which.bac=1:nrow(x$loadings),
    ptype="spider", comp=1:2, cex.bac=0.5, plot.names=TRUE,
    group.cols=as.numeric(unique(Y)), ...) {

    # TODO: can we speed up this function ?
    # For instance by switching from for loops to vectorization etc

    y <- cluster <- x.centroid <- y.centroid <- NULL
    
    bag <- x  #[which(!names(x) == 'variable')] #$bagged.rda
    Y <- x$variable
    
    scaled.loadings <- (bag$loadings/max(abs(bag$loadings)))[, comp]
    scaled.scores <- (bag$scores/max(abs(bag$scores)))[, comp]
    
    plot(rbind(scaled.scores, scaled.loadings), type="n",
        xlab=paste(names(bag$R2)[1], 
        " (", format(100 * bag$R2[1], digits=2), "%)", sep=""),
    ylab=paste(names(bag$R2)[2], 
        " (", format(100 * bag$R2[2], digits=2), "%)", sep=""))
    
    if (ptype == "spider") 
        s.class(scaled.scores, factor(Y), grid=FALSE, col=group.cols,
        cellipse=0.5, cpoint=0, add.plot=TRUE)
    
    # TODO: Same with ggplot (not ready)
    skip <- TRUE
    if (!skip) {
        # build ggplot dataframe with points (x,y) and
        # corresponding groups (cluster)
        gg <- as.data.frame(scaled.scores)
        names(gg) <- c("x", "y")
        gg$cluster <- factor(Y)
        # calculate group centroid locations
        centroids <- aggregate(cbind(x, y) ~ cluster, data=gg, mean)
        # merge centroid locations into ggplot dataframe
        gg <- merge(gg, centroids, by="cluster",
        suffixes=c("", ".centroid"))
        # generate star plot...
        ggplot(gg) +
        geom_point(aes(x=x, y=y, color=cluster), size=3) +
        geom_point(data=centroids, 
                aes(x=x, y=y, color=cluster), size=4) +
        geom_segment(aes(x=x.centroid, 
                y=y.centroid, xend=x, yend=y, color=cluster))
    }
    
    if (ptype == "hull") {
        
        ll <- split(rownames(scaled.scores), Y)
        hulls <- lapply(ll, function(ii) ii[chull(scaled.scores[ii, ])])
        for (i in 1:length(hulls)) {
            polygon(scaled.scores[hulls[[i]], ], border=group.cols[i])
        }
    }
    if (plot.names) {
        text(scaled.scores, rownames(scaled.scores), cex=0.5, ...)
    } else {
        points(scaled.scores, ...)
    }
    text(scaled.loadings[which.bac, ],
        rownames(scaled.loadings)[which.bac], cex=cex.bac)
}

