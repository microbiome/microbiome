#' @title Visualize OTU Core
#' @description Core visualization (2D).
#' @param x A \code{\link{phyloseq}} object or a core matrix
#' @param prevalences a vector of prevalence percentages in [0,1]
#' @param detections a vector of intensities around the data range,
#' or a scalar indicating the number of intervals in the data range.
#' @param plot.type Plot type ('lineplot' or 'heatmap')
#' @param colours colours for the heatmap
#' @param min.prevalence If minimum prevalence is set, then filter out those
#' rows (taxa) and columns (detections) that never exceed this
#' prevalence. This helps to zoom in on the actual core region
#' of the heatmap. Only affects the plot.type='heatmap'.
#' @param taxa.order Ordering of the taxa: a vector of names.
#' @param horizontal Logical. Horizontal figure.
#' @return A list with three elements: the ggplot object and the data.
#' The data has a different form for the lineplot and heatmap.
#' Finally, the applied parameters are returned.
#' @examples 
#' data(dietswap)
#' p <- plot_core(transform(dietswap, "compositional"),
#'   prevalences=seq(0.1, 1, .1), detections=seq(0.01, 1, length = 10))
#' @export 
#' @references 
#' A Salonen et al. The adult intestinal core microbiota is determined by 
#' analysis depth and health status. Clinical Microbiology and Infection 
#' 18(S4):16 20, 2012. 
#' To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
plot_core <- function(x, prevalences=seq(.1, 1, 0.1), detections=20,
    plot.type="lineplot", colours=NULL, # gray(seq(0, 1, length=5)),
    min.prevalence=NULL, taxa.order=NULL, horizontal=FALSE) {
    
    if (length(detections) == 1) {
        detections <- 10^seq(log10(0.001), log10(max(abundances(x),
        na.rm=TRUE)), length=detections)
    }

    if (!is_compositional(x)) {
        warning("The plot_core function is typically used with compositional 
                data. The data is not compositional. Make sure that you
                intend to operate on non-compositional data.")
    }

    if (plot.type == "lineplot") {

        # Calculate the core matrix (prevalences x abundance thresholds)
        coremat <- core_matrix(x, prevalences, detections)
        res <- core_lineplot(coremat)

    } else if (plot.type == "heatmap") {

        # Here we use taxon x abundance thresholds table
        #  indicating prevalences
        res <- core_heatmap(
                abundances(x),
                    dets=detections,
                    cols=colours, 
                    min.prev=min.prevalence,
                    taxa.order=taxa.order)
    }

    p <- res$plot
    
    if (horizontal) {
        p <- p + coord_flip() + theme(axis.text.x=element_text(angle=90))
    }
    
    p
    
}


#' @title Core Matrix 
#' @description Creates the core matrix.
#' @param x \code{\link{phyloseq}} object or a taxa x samples abundance matrix
#' @param prevalences a vector of prevalence percentages in [0,1]
#' @param detections a vector of intensities around the data range
#' @return Estimated core microbiota
#' @examples
#' # Not exported
#' #data(peerj32)
#' #core <- core_matrix(peerj32$phyloseq)
#' @references 
#' A Salonen et al. The adult intestinal core microbiota is determined by 
#' analysis depth and health status. Clinical Microbiology and Infection 
#' 18(S4):16 20, 2012. 
#' To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
core_matrix <- function(x, prevalences=seq(0.1, 1, , 1), detections=NULL) {
    
    # Pick abundances
    data <- abundances(x)
    
    # Convert prevalences from percentages to sample counts
    p.seq <- 0.01 * prevalences * ncol(data)
    
    ## Intensity vector
    if (is.null(detections)) {
        detections <- seq(min(data), max(data), length=10)
    }
    i.seq <- detections
    
    coreMat <- matrix(NA, nrow=length(i.seq), ncol=length(p.seq),
        dimnames=list(i.seq, p.seq))
    
    n <- length(i.seq) * length(p.seq)
    cnt <- 0
    for (i in i.seq) {
        for (p in p.seq) {
            # Number of OTUs above a given prevalence
            coreMat[as.character(i), as.character(p)] <-
            sum(rowSums(data > i) >= p)                
        }
    }
    
    # # Convert Prevalences to percentages
    colnames(coreMat) <- as.numeric(colnames(coreMat))/ncol(data)
    rownames(coreMat) <- as.character(as.numeric(rownames(coreMat)))
    
    coreMat
    
}


#' @title Core Heatmap
#' @description Core heatmap.
#' @param x OTU matrix
#' @param dets A vector or a scalar indicating the number of intervals
#' in (0, log10(max(data))). The dets are calculated for relative
#' abundancies.
#' @param cols colours for the heatmap
#' @param min.prev If minimum prevalence is set, then filter out those
#' rows (taxa) and columns (dets) that never exceed this prevalence.
#' This helps to zoom in on the actual core region of the heatmap.
#' @param taxa.order Ordering of the taxa.
#' @return Used for its side effects
#' @references 
#' A Salonen et al. The adult intestinal core microbiota is determined by 
#' analysis depth and health status. Clinical Microbiology and Infection 
#' 18(S4):16 20, 2012. 
#' To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
core_heatmap <- function(x, dets, cols, min.prev, taxa.order)
{

    data <- x
    DetectionThreshold <- Taxa <- Prevalence <- NULL
    
    # Prevalences with varying dets
    prev <- lapply(dets, function(th) {
        prevalence(data, detection=th)
    })
    prev <- do.call("cbind", prev)
    colnames(prev) <- as.character(dets)

    # Exclude rows and cols that never exceed the given prevalence
    if (!is.null(min.prev)) {
        rinds <- rowMeans(prev > min.prev) > 0
	cinds <- colMeans(prev > min.prev) > 0
        prev <- prev[rinds, cinds, drop=FALSE]	
    }

    df <- as.data.frame(prev)
    if (nrow(df) == 0) {stop("Too few taxa fulfil the criteria on detection and prevalence. Apply less conservative limits.")}
    
    df$ID <- rownames(prev)

    df <- melt(df, "ID")
    names(df) <- c("Taxa", "DetectionThreshold", "Prevalence")
    df$DetectionThreshold <- as.numeric(as.character(df$DetectionThreshold))
    df$Prevalence <- as.numeric(as.character(df$Prevalence))
    df$DetectionThreshold <- factor(df$DetectionThreshold)

    if (is.null(taxa.order)) {
        o <- names(sort(rowSums(prev)))
    } else {
        o <- taxa.order
    }
    df$Taxa <- factor(df$Taxa, levels=o)

    p <- ggplot(df, aes(x=DetectionThreshold, y=Taxa, fill=Prevalence)) +
            geom_tile() +
            labs(y = "")

    if (is_compositional(x)) {

        lab <- paste0(100 *
             as.numeric(as.character(unique(df$DetectionThreshold))), "%")
	
        p <- p + scale_x_discrete(labels=lab)

        if (!is.null(cols)) {
            p <- p + scale_fill_gradientn("Prevalence",
                    breaks=seq(from=0, to=1, by=0.1),
                    labels=scales::percent,                
                    colours=cols,
                    limits=c(0, 1))
        }

    } else {

        if (!is.null(cols)) {

            p <- p + scale_fill_gradientn("Prevalence",
                breaks=seq(from=0, to=1, by=0.1),
                colours=cols,
                limits=c(0, 1))
        }

    }
    p <- p + labs(x = "Detection Threshold")
        
    return(list(plot=p, data=df))
    
}


core_lineplot <- function(x, 
    xlabel="Abundance", ylabel="Core size (N)") {

    Abundance <- Prevalence <- Count <- NULL

    df <- as.data.frame(x)
    df$ID <- rownames(x)
    df <- melt(df, "ID")    
    names(df) <- c("Abundance", "Prevalence", "Count")
    
    df$Abundance <- as.numeric(as.character(df$Abundance))
    df$Prevalence <- as.numeric(as.character(df$Prevalence))
    df$Count <- as.numeric(as.character(df$Count))
    
    p <- ggplot(df, aes(x=Abundance, y=Count,
        color=Prevalence, group=Prevalence))
    
    p <- p + geom_line()
    p <- p + geom_point()
    p <- p + scale_x_log10()
    p <- p + xlab(xlabel)
    p <- p + ylab(ylabel)
    
    list(plot=p, data=x)
}





