#' @title Landscape Plot
#' @description Wrapper for visualizing sample similarity landscape
#' ie. sample density in various 2D projections.
#' @param x \code{\link{phyloseq-class}} object or a data matrix
#' (samples x features; eg. samples vs. OTUs). If the input x is a 2D matrix
#' then it is plotted as is.
#' @param method Ordination method, see phyloseq::plot_ordination; or "PCA",
#'        or "t-SNE" (from the \pkg{Rtsne} package)
#' @param distance Ordination distance, see phyloseq::plot_ordination; for
#'        method = "PCA", only euclidean distance is implemented now.
#' @param transformation Transformation applied on the input object x
#' @param col Variable name to highlight samples (points) with colors
#' @param main title text
#' @param x.ticks Number of ticks on the X axis
#' @param rounding Rounding for X axis tick values
#' @param add.points Plot the data points as well
#' @param point.opacity Transparency-level for points
#' @param adjust Kernel width adjustment
#' @param size point size
#' @param legend plot legend TRUE/FALSE
#' @param shading Add shading in the background.
#' @param shading.low Color for shading low density regions
#' @param shading.high Color for shading high density regions
#' @return A \code{\link{ggplot}} plot object.
#' @export
#' @details For consistent results, set random seet (set.seed) before
#' function call. Note that the distance and transformation arguments may
#' have a drastic effect on the outputs.
#' @examples
#'
#' data(dietswap)
#'
#' # PCoA
#' p <- plot_landscape(transform(dietswap, "compositional"),
#'    distance = "bray", method = "PCoA")
#'
#  # t-SNE
#' p <- plot_landscape(dietswap, method = "t-SNE", distance = "bray",
#'        transformation = "compositional")
#'
#' # PCA
#' p <- plot_landscape(dietswap, method = "PCA", transformation = "clr")
#'
#' @keywords utilities
plot_landscape <- function(x, method="PCoA",
                           distance="bray",
                           transformation = "identity",
                           col=NULL,
                           main=NULL,
                           x.ticks=10,
                           rounding=0,
                           add.points=TRUE,
                           adjust=1, size=1,
                           legend=FALSE,
                           shading=TRUE,
                           shading.low="#ebf4f5",
                           shading.high="#e9b7ce",
                           point.opacity=0.75) {

    if (is.matrix(x) || is.data.frame(x)) {

        if (ncol(x) == 2) {
            proj <- as.data.frame(x)
        } else if (ncol(x) > 2) {
            # Convert the matrix into phyloseq object
            x <- phyloseq(otu_table(t(x), taxa_are_rows = TRUE))
        }
    }

    if (is.phyloseq(x)) {

        x <- transform(x, transformation)

        if (method == "PCA") {

            # TODO: add option to calculate PCA with different
            # distances using the distance argument
            # (now assumes the default and cannot be altered)
            d <- t(abundances(x))
            pca <- princomp(d)
            proj <- pca$scores[, c(1,2)]
            rownames(proj) <- sample_names(x)

            # TODO add robust PCA
            #(pc.rob <- princomp(stackloss, covmat = MASS::cov.rob(stackloss)))

        } else if (method == "t-SNE") {

            dm <- vegdist(t(otu_table(x)), distance)

            ## Run TSNE
            tsne_out <- Rtsne(dm, dims = 2)
            proj <- tsne_out$Y
            rownames(proj) <- sample_names(x)

        } else {

            #quiet(proj <- get_ordination(x, method, distance))
            quiet(x.ord <- ordinate(x, method, distance))
            # Pick the projected data (first two columns + metadata)
            quiet(proj <- phyloseq::plot_ordination(x, x.ord, justDF=TRUE))
            # Rename the projection axes
        }

    }

    proj <- as.data.frame(proj)
    colnames(proj) <- paste0("PC", c(1,2))

    guide.title <- "color"
    if (is.null(col)) {
        proj$col <- as.factor(rep("black", nrow(proj)))
    } else if (length(col) == 1 && col %in% names(meta(x))) {
        proj$col <- meta(x)[, col]
        guide.title <- col
    } else {
        proj$col <- col
    }

    p <- densityplot(proj[, c(1,2)], main=NULL, x.ticks=10,
            rounding=0, add.points=add.points,
            adjust=1, size=size, col=proj$col,
            legend = TRUE,
            shading=shading,
            shading.low=shading.low,
            shading.high=shading.high,
            point.opacity=point.opacity) +
            guides(color = guide_legend(title = guide.title))

    p

}





#' @title Density Plot
#' @description Density visualization for data points overlaid on cross-plot.
#' @param x Data matrix to plot. The first two columns will be visualized as a
#'    cross-plot.
#' @param main title text
#' @param x.ticks Number of ticks on the X axis
#' @param rounding Rounding for X axis tick values
#' @param add.points Plot the data points as well
#' @param point.opacity Transparency-level for points
#' @param col Color of the data points. NAs are marked with darkgray.
#' @param adjust Kernel width adjustment
#' @param size point size
#' @param legend plot legend TRUE/FALSE
#' @param shading Shading
#' @param shading.low Color for shading low density regions
#' @param shading.high Color for shading high density regions
#' @return ggplot2 object
#' @examples
#' # p <- densityplot(cbind(rnorm(100), rnorm(100)))
#' @references See citation('microbiome')
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
densityplot <- function(x, main=NULL, x.ticks=10, rounding=0,
    add.points=TRUE, col="black",
    adjust=1, size=1,
    legend=FALSE, shading=TRUE,
    shading.low="white",
    shading.high="black",
    point.opacity=0.75) {

    df <- x
    if (!is.data.frame(df)) {
        df <- as.data.frame(as.matrix(df))
    }

    # Avoid warnings
    x <- y <- ..density.. <- color <- NULL

    # If colors are NA then mark them dark gray
    lev <- levels(col)
    if (!is.numeric(col) & any(is.na(col))) {
        col <- as.character(col)
        col[unname(which(is.na(col)))] <- "darkgray"
        col <- factor(col, levels = unique(c(lev, "darkgray")))
    }

    xvar <- colnames(df)[[1]]
    yvar <- colnames(df)[[2]]
    df[["x"]] <- df[, 1]
    df[["y"]] <- df[, 2]

    df[["color"]] <- col
    df[["size"]] <- size

    # Remove NAs
    df <- df[!(is.na(df[["x"]]) | is.na(df[["y"]])), ]

    # Determine bandwidth for density estimation
    bw <- adjust * c(bwi(df[["x"]]), bwi(df[["y"]]))
    if (any(bw == 0)) {
        warning("Zero bandwidths
                (possibly due to small number of observations).
                Using minimal bandwidth.")
        bw[bw == 0]=bw[bw == 0] + min(bw[!bw == 0])
    }

    # Construct the figure
    p <- ggplot(df)

    if (shading) {
        p <- p + stat_density2d(aes(x, y, fill=..density..),
                geom="raster", h=bw,
                contour=FALSE)
        p <- p + scale_fill_gradient(low=shading.low, high=shading.high)
    }

    if (add.points) {

        if (length(unique(df$color)) == 1 && length(unique(df$size)) == 1) {
            p <- p + geom_point(aes(x=x, y=y),
            col=unique(df$color), size=unique(df$size),
            alpha=point.opacity)
        } else if (length(unique(df$color)) == 1 &&
            length(unique(df$size)) > 1) {
            p <- p + geom_point(aes(x=x, y=y, size=size),
            col=unique(df$color),
            alpha=point.opacity)
        } else if (length(unique(df$color)) > 1 &&
            length(unique(df$size)) == 1) {
            p <- p + geom_point(aes(x=x, y=y, col=color),
            size=unique(df$size),
            alpha=point.opacity)
        } else {
            p <- p + geom_point(aes(x=x, y=y, col=color, size=size),
                                alpha=point.opacity)
        }
    }

    p <- p + xlab(xvar) + ylab(yvar)

    if (!legend) {
        p <- p + theme(legend.position="none")
    }

    p <- p + scale_x_continuous(breaks=round(seq(floor(min(df[["x"]])),
        ceiling(max(df[["x"]])), length=x.ticks), rounding))

    if (!is.null(main)) {
        p <- p + ggtitle(main)
    }

    p + theme(panel.grid = element_blank(),
              panel.background = element_rect(fill = "white", color="grey70"),
              legend.key=element_blank())

}



# Bandwidth
# As in MASS::bandwidth.nrd but rewritten. Internal.
bwi <- function (x) {
    r <- quantile(x, c(0.25, 0.75))
    4 * 1.06 * min(sd(x), (r[[2]] - r[[1]])/1.34) * length(x)^(-.2)
}
