#' @title Association Heatmap
#' @description Visualizes n x m association table as heatmap. 
#' @param df Data frame. Each row corresponds to a pair of associated 
#' variables. The columns give variable names, association scores and 
#' significance estimates.
#' @param Xvar X axis variable column name. For instance 'X'.
#' @param Yvar Y axis variable column name. For instance 'Y'.
#' @param fill Column to be used for heatmap coloring. 
#' For instance 'association'.
#' @param star Column to be used for cell highlighting. For instance 'p.adj'.
#' @param p.adj.threshold Significance threshold for the stars.
#' @param association.threshold Include only elements that have absolute 
#' association higher than this value
#' @param step color interval
#' @param colours heatmap colours
#' @param limits colour scale limits
#' @param legend.text legend text
#' @param order.rows Order rows to enhance visualization
#' interpretability. If this is logical, then hclust is applied. If this
#' is a vector then the rows are ordered using this index.
#' @param order.cols Order columns to enhance visualization
#' interpretability. If this is logical, then hclust is applied. If this
#' is a vector then the rows are ordered using this index.
#' @param filter.significant Keep only the elements with at least one 
#' significant entry
#' @param star.size NULL Determine size of the highlight symbols
#' @param plot.values Show values as text
#' @return ggplot2 object
#' @examples
#' data(peerj32)
#' d1 <- peerj32$lipids[, 1:10]
#' d2 <- peerj32$microbes[, 1:10]
#' cc <- associate(d1, d2, method='pearson') 
#' p <- heat(cc, 'X1', 'X2', 'Correlation', star='p.adj')
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
heat <- function(df, Xvar = names(df)[[1]], Yvar = names(df)[[2]],
    fill = names(df)[[3]], star = NULL, p.adj.threshold=1,
    association.threshold=0, 
    step=0.2, colours=c("darkblue", "blue", "white", "red", "darkred"),
    limits=NULL, legend.text="", order.rows=TRUE, order.cols=TRUE,
    filter.significant=TRUE, 
    star.size=NULL, plot.values=FALSE) {

    if (is.null(limits)) {
        maxval <- max(abs(df[[fill]]))
        if (maxval <= 1) {
            limits <- c(-1, 1)
        } else {
            xmaxval <- ceiling(maxval)
            limits <- c(-maxval, maxval)
        }
    }

    # Truncate anything exceeding the limits
    df[df[[fill]] < limits[[1]], fill] <- limits[[1]]
    df[df[[fill]] > limits[[2]], fill] <- limits[[2]]

    if (nrow(df) == 0) {
        warning("Input data frame is empty.")
        return(NULL)
    }

    if (filter.significant & !is.null(star)) {
        keep.X <- as.character(unique(df[((df[[star]] < p.adj.threshold) &
        (abs(df[[fill]]) > association.threshold)), Xvar]))
        keep.Y <- as.character(unique(df[((df[[star]] < p.adj.threshold) &
        (abs(df[[fill]]) > 
            association.threshold)), Yvar]))
        df <- df[((df[[Xvar]] %in% keep.X) & (df[[Yvar]] %in% keep.Y)), ]
    }

    if (any(c("XXXX", "YYYY", "ffff") %in% names(df))) {
        stop("XXXX, YYYY, ffff are not allowed in df")
    }
    
    df[[Xvar]] <- factor(df[[Xvar]])
    df[[Yvar]] <- factor(df[[Yvar]])

    # TODO neatmap with options
    # xo <- neat(mat, method = "NMDS", distance = "euclidean",
    # first.row = "VDP.03231", first.col = "Paraprevotella") 
    
    if (is.logical(order.rows) || is.logical(order.cols)) {
        
        rnams <- unique(as.character(df[[Xvar]]))
        cnams <- unique(as.character(df[[Yvar]]))    

        # Rearrange into matrix
        # FIXME: could be better done with cast
        mat <- matrix(0, nrow=length(rnams), ncol=length(cnams))
        rownames(mat) <- rnams
        colnames(mat) <- cnams

        for (i in seq_len(nrow(df))) {
            mat[as.character(df[i, Xvar]),
                as.character(df[i, Yvar])] <- df[i, fill]
        }
        
        mat <- t(mat)
        cind <- seq_len(ncol(mat))
        rind <- seq_len(nrow(mat))
    
    } 

    if (is.logical(order.rows)) {
    
        if (order.rows) {
        
            if (nrow(mat) > 1 && ncol(mat) > 1) {
                rind <- hclust(as.dist(1 - cor(t(mat),
                    use="pairwise.complete.obs")))$order
            }

            if (nrow(mat) > 1 && ncol(mat) == 1) {
                rind <- order(mat[, 1])
            }

            order.rows <- rownames(mat)[rind]

        } else {

            order.rows <- rownames(mat)[rind]
    
        }
    
    }

    if (is.logical(order.cols)) {

        if (order.cols) {

            if (ncol(mat) > 1 && nrow(mat) > 1) {
    
                cind <- hclust(as.dist(1 - cor(mat,
                    use="pairwise.complete.obs")))$order
        
            } else {
    
                cind <- order(mat[1, ])
        
            }

            order.cols <- colnames(mat)[cind]
    
        } else {

            order.cols <- colnames(mat)[cind]

        }
        
    }

    df[[Yvar]] <- factor(df[[Yvar]], levels=order.rows)
    df[[Xvar]] <- factor(df[[Xvar]], levels=order.cols)
    
    XXXX <- YYYY <- ffff <- NULL
    df[["XXXX"]] <- df[[Xvar]]
    df[["YYYY"]] <- df[[Yvar]]
    df[["ffff"]] <- df[[fill]]
    
    p <- ggplot(df, aes(x=XXXX, y=YYYY, fill=ffff)) +
            geom_tile()
    
    p <- p + scale_fill_gradientn(legend.text,
            breaks=seq(from=min(limits), to=max(limits), 
            by=step), colours=colours, limits=limits) +
            labs(x = "", y = "") +
            theme(axis.text.x=element_text(angle=90))
    
    # Mark significant cells with stars
    if (!is.null(star)) {
        inds <- which((df[[star]] < p.adj.threshold) &
            (abs(df[[fill]]) > association.threshold))
        if (!is.null(star) & length(inds) > 0) {
            df.sub <- df[inds, ]
        
        if (is.null(star.size)) {
            star.size <- 1 # max(1, floor(text.size/2))
        }
        
        p <- p + geom_text(data=df.sub, aes(x=XXXX, y=YYYY, label="+"),
                           col="white", size=star.size)
        }
    }
    if (plot.values) {
        p <- p + geom_text(aes(label=round(ffff, 2)), size=3)
    }
    
    p
    
}


