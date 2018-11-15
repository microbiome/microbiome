#' @title Taxonomic Composition Plot
#' @description Plot taxon abundance for samples.
#' @param x \code{\link{phyloseq-class}} object
#' @param sample.sort Order samples. Various criteria are available:
#' \itemize{
#' \item NULL or 'none': No sorting
#' \item A single character string: indicate the metadata field to be used
#'  for ordering
#' \item A character vector: sample IDs indicating the sample ordering.
#' \item 'neatmap' Order samples based on the neatmap approach.
#' See \code{\link{neatsort}}. By default, 'NMDS' method with 'bray'
#' distance is used. For other options, arrange the samples manually
#' with the function.
#' }
#' @param otu.sort Order taxa. Same options as for the sample.sort argument
#' but instead of metadata, taxonomic table is used.
#' Also possible to sort by 'abundance'.
#' @param x.label Specify how to label the x axis.
#' This should be one of the variables in sample_variables(x).
#' @param plot.type Plot type: 'barplot' or 'heatmap'
#' @param verbose verbose
#'  (but not in sample/taxon ordering).
#' The options are 'Z-OTU', 'Z-Sample', 'log10' and 'compositional'.
#' See the \code{\link{transform}} function.
#' @param mar Figure margins
#' @param average_by Average the samples by the average_by variable 
#' @param ... Arguments to be passed (for \code{\link{neatsort}} function)
#' @return A \code{\link{ggplot}} plot object.
#' @export
#' @examples
#' data(dietswap)
#' pseq <- subset_samples(dietswap, group == 'DI' & nationality == 'AFR' & sex == "female")
#' p <- plot_composition(pseq)
#' @keywords utilities
plot_composition <- function(x, sample.sort=NULL,
    otu.sort=NULL, x.label="sample", plot.type="barplot",
    verbose=FALSE, 
    mar=c(5, 12, 1, 1), average_by=NULL, ...) {
    
    # Avoid warnings
    Sample <- Abundance <- Taxon <- horiz <- value <- scales <- ID <- 
        meta <- OTU <- taxic <- otu.df <- taxmat <-  new.tax<- NULL
    if (!is.null(x@phy_tree)){
        x@phy_tree= NULL
    }
    taxic <- as.data.frame(x@tax_table) 
    otu.df <- as.data.frame(otu_table(x))
    taxic$OTU <- row.names(otu.df)
    taxmat <- as.matrix(taxic) # convert it into a matrix.
    new.tax <- tax_table(taxmat)  # convert into phyloseq compaitble file.
    tax_table(x) <- new.tax # incroporate into phyloseq Object
    
    xorig <- x

    # Pick the abundance matrix taxa x samples
    abu <- abundances(x)
    
    # Average the samples by group
    group <- NULL
    if (!is.null(average_by)) {
        dff <- as.data.frame(t(abu))
        dff$group <- sample_data(x)[[average_by]]
        if (is.numeric(dff$group)) {
            dff$group <- factor(dff$group, levels=sort(unique(dff$group)))
        }
        # Remove samples with no group info
        dff <- dff %>% filter(!is.na(group))
        dff$group <- droplevels(dff$group)
        # av <- ddply(dff, "group", colwise(mean))
        av <- aggregate(. ~ group, data = dff, mean)
        rownames(av) <- as.character(av$group)
        av$group <- NULL
        abu <- t(av)  # taxa x groups
    }
    
    # Sort samples
    if (is.null(sample.sort) || sample.sort == "none" || !is.null(average_by)) {
        # No sorting sample.sort <- sample_names(x)
        sample.sort <- colnames(abu)
    } else if (length(sample.sort) == 1 &&
        sample.sort %in% names(sample_data(x)) && 
        is.null(average_by)) {
        # Sort by metadata field
        sample.sort <- rownames(sample_data(x))[order(
        sample_data(x)[[sample.sort]])]
    } else if (all(sample.sort %in% sample_names(x)) & is.null(average_by)) {
        # Use predefined order
        sample.sort <- sample.sort
    } else if (length(sample.sort) == 1 && sample.sort == "neatmap") {
    
        sample.sort <- neatsort(x, method="NMDS", distance="bray",
        target="sites", first=NULL)

    } else if (is.vector(sample.sort) && length(sample.sort) > 1) {
        sample.sort <- sample_names(x)[sample.sort]        	
    } else if (!sample.sort %in% names(sample_data(x))) {
        warning(paste("The sample.sort argument", sample.sort,
        "is not included in sample_data(x). 
            Using original sample ordering."))
        sample.sort <- sample_names(x)
    }
    
    # Sort taxa
    if (is.null(otu.sort) || otu.sort == "none") {
        # No sorting
        otu.sort <- taxa(x)
    } else if (length(otu.sort) == 1 && otu.sort == "abundance") {
        otu.sort <- rev(names(sort(rowSums(abu))))
    } else if (length(otu.sort) == 1 && otu.sort %in% names(tax_table(x))) {
        # Sort by phylogenetic group
        otu.sort <- rownames(sample_data(x))[order(tax_table(x)[[otu.sort]])]
    } else if (all(otu.sort %in% taxa(x))) {
        # Use predefined order
        otu.sort <- otu.sort
    } else if (length(otu.sort) == 1 && otu.sort == "neatmap") {
        otu.sort <- neatsort(x, method="NMDS", distance="bray",
        target="species", first=NULL)
    }
    
    if (verbose) {
        message("Prepare data.frame.")
    }
    # Abundances as data.frame dfm <- psmelt(x)
    #dfm <- melt(abu)
    dfm <- psmelt(otu_table(abu, taxa_are_rows = TRUE))
    names(dfm) <- c("OTU", "Sample", "Abundance")

    dfm$Sample <- factor(dfm$Sample, levels=sample.sort)
    dfm$OTU <- factor(dfm$OTU, levels=otu.sort)
    
    # SampleIDs for plotting
    if (x.label %in% colnames(sample_data(x)) & is.null(average_by)) {
        
        meta <- sample_data(x)
        dfm$xlabel <- as.vector(unlist(meta[as.character(dfm$Sample), x.label]))
        
        # Sort the levels as in the original metadata
        if (is.factor(meta[, x.label])) {
            
            lev <- levels(meta[, x.label])
            
        } else {
            
            lev <- unique(as.character(unname(unlist(meta[, x.label]))))
            
        }
        
        dfm$xlabel <- factor(dfm$xlabel, levels=lev)
        
    } else {
        
        dfm$xlabel <- dfm$Sample
        
    }
    
    if (verbose) {
        message("Construct the plots")
    }
    
    if (plot.type == "barplot") {
        
        # Provide barplot
        dfm <- dfm %>% arrange(OTU)  # Show OTUs always in the same order

        p <- ggplot(dfm, aes(x=Sample, y=Abundance, fill=OTU))
        p <- p + geom_bar(position="stack", stat="identity")
        p <- p + scale_x_discrete(labels=dfm$xlabel, breaks=dfm$Sample)

        # Name appropriately
        p <- p + ylab("Abundance")
        
        # Rotate horizontal axis labels, and adjust
        p <- p + theme(axis.text.x=element_text(angle=90, vjust=0.5,
            hjust=0))
        p <- p + guides(fill=guide_legend(reverse=FALSE))

    } else if (plot.type == "heatmap") {
        
        if (verbose) {
            message("Constructing the heatmap.")
        }
        
        # Taxa x samples otu matrix
        otu <- abundances(x)
        
        # Remove NAs after the transform
        otu <- otu[rowMeans(is.na(otu)) < 1, colMeans(is.na(otu)) < 1]
        
        otu.sort <- otu.sort[otu.sort %in% rownames(otu)]
        sample.sort <- sample.sort[sample.sort %in% colnames(otu)]
        
        # Plot TODO: move it in here from netresponse and return the
        # #ggplot object as well
        # p <- plot_matrix(otu[otu.sort, sample.sort], type="twoway", mar=mar)
        tmp <- melt(otu[otu.sort, sample.sort])
        p <- heat(tmp, colnames(tmp)[[1]], colnames(tmp)[[2]],
            colnames(tmp)[[3]])
        
    } else if (plot.type == "lineplot") {

        dfm <- dfm %>% arrange(OTU)  # Show OTUs always in the same order
        p <- ggplot(dfm, aes(x=Sample, y=Abundance, color=OTU, group = OTU))
        p <- p + geom_point()
        p <- p + geom_line()    
        p <- p + scale_x_discrete(labels=dfm$xlabel, breaks=dfm$Sample)
        
        # Name appropriately
        if (!is.null(transform) && transform == "relative.abundance") {
            p <- p + ylab("Relative abundance (%)")
        } else {
            p <- p + ylab("Abundance")
        }
        
        # Rotate horizontal axis labels, and adjust
        p <- p + theme(axis.text.x=element_text(angle=90, vjust=0.5,
        hjust=0))
        p <- p + guides(fill=guide_legend(reverse=FALSE))

    }

    p
    
}

