#' @title plot_composition
#' @description Plot taxon abundance for samples
#' @param x \code{\link{phyloseq-class}} object or an OTU matrix (samples x phylotypes)
#' @param taxonomic.level Merge the OTUs (for phyloseq object) into a higher taxonomic level. This has to be one from \code{colnames(tax_table(x))}.
#' @param relative.abundance Logical. Show relative abundances or not.
#' @param sort.by Sort by sample data column. Or provide vector of sample IDs.
#' @param x.label Specify how to label the x axis. This should be one of the variables in 
#'        \code{sample_variables(x)}
#' @return A \code{\link{ggplot}} plot object.
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom phyloseq tax_glom
#' @export
#' @examples \dontrun{
#'   # Example data
#'     library(microbiome)
#'     pseq0 <- download_microbiome("dietswap")
#'     pseq <- subset_samples(pseq0, group == "DI" & nationality == "AFR")
#'     plot_composition(pseq, taxonomic.level = "Phylum")
#'           }
#' @keywords utilities
plot_composition <- function (x, taxonomic.level = NULL, relative.abundance = FALSE, sort.by = NULL, x.label = "sample") {

  # Avoid warnings
  Sample <- Abundance <- Taxon <- horiz <- value <- NULL

  # Merge the taxa at a higher taxonomic level
  if (!is.null(taxonomic.level)) {	         
    x <- tax_glom(x, taxonomic.level)
    # Fix the taxon names; tax_glom assigns wrong names
  }

  # Pick the OTU data
  if (class(x) == "phyloseq") {
    otu <- otu_table(x)@.Data
    # FIXME: Remove this when this has been fixed to phyloseq package - pending 5/2015
    if (!is.null(taxonomic.level)) {
      rownames(otu) <- as.vector(tax_table(x)[, taxonomic.level])
    }
  } else {
    otu <- t(x)
  }

  # Estimate relative abundances
  if (relative.abundance) {
    otu <- relative.abundance(otu)
  }

  # Define sample ordering based on the sort.by column
  meta <- sample_data(x)
  if (is.null(sort.by)) {
    sort.by <- rownames(meta)    
  } else if (length(sort.by) == 1) {
    sort.by <- rownames(meta)[order(meta[[sort.by]])]
  }

  # Prepare data.frame
  dfm <- melt(otu)
  colnames(dfm) <- c("Taxon", "Sample", "Abundance")
  dfm$Sample <- factor(as.character(dfm$Sample), levels = sort.by)

  # SampleIDs used in plotting
  if (x.label %in% colnames(meta)) {
    dfm$xlabel <- as.vector(unlist(meta[as.character(dfm$Sample), x.label]))
    # Sort the levels as in the original metadata
    if (is.factor(meta[, x.label])) {
      lev <- levels(meta[, x.label])
    } else {
      lev <- unique(as.character(meta[, x.label]))
    }
    dfm$xlabel <- factor(dfm$xlabel, levels = lev)
  } else {
    dfm$xlabel <- dfm$Sample
  }
  
  # Provide barplot of relative abundances
  p <- ggplot(dfm, aes(x = Sample, y = Abundance, fill = Taxon))
  p <- p + geom_bar(position = "stack", stat = "identity")
  p <- p + scale_x_discrete(labels = dfm$xlabel, breaks = dfm$Sample)

  # Rotate horizontal axis labels, and adjust
  p <- p + theme(axis.text.x=element_text(angle=-90, vjust=0.5, hjust=0))

  p

}

