#' @title Taxonomic Abundance Plot
#' @description Plot taxon abundance for samples.
#' @param x \code{\link{phyloseq-class}} object
#' @param level Taxonomic level to show.
#' @param tax Taxonomic group to use in sorting.
#' @return \code{\link{ggplot}} object
#' @examples
#' data(dietswap)
#' pseq <- subset_samples(dietswap, group == 'DI' & nationality == 'AFR' &
#'    sex == "female")
#' p <- plot_abundances(pseq)
#' @keywords utilities
#' @export
plot_abundances <- function(x, level = "Phylum", tax = "Bacteroidetes") {

    x <- aggregate_taxa(x, level = level)

    xc <- transform(x, "compositional")
    sample.sort <- rev(sample_names(x)[order(abundances(xc)[tax,])])
    
    p <- plot_composition(xc,
            sample.sort = sample.sort,
            otu.sort = "abundance") +
            scale_y_continuous(labels = scales::percent)

    p
    
}

