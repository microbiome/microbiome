#' @title Visualize Prevalence Distributions for Taxa
#' @description Create taxa prevalence plots at various taxonomic levels.
#' @details This helps to obtain first insights into how is the taxa
#' distribution in the data. It also gives an idea about the taxonomic
#' affiliation of rare and abundant taxa in the data.
#' This may be helpful for data filtering or other downstream analysis.
#' @param x \code{\link{phyloseq-class}} object, OTU data must be counts and
#'        not relative abundance or other transformed data.
#' @param level Phylum/Order/Class/Family
#' @param detection Detection threshold for presence (prevalance)
#' @return A \code{\link{ggplot}} plot object.
#' @export
#' @examples
#' data(atlas1006)
#' # Pick data subset just to speed up example
#' p0 <- subset_samples(atlas1006, DNA_extraction_method == "r")
#' p0 <- prune_taxa(taxa(p0)[grep("Bacteroides", taxa(p0))], p0)
#' # Detection threshold (0 by default; higher especially with HITChip)
#' p <- plot_taxa_prevalence(p0, 'Phylum', detection = 1)
#' print(p)
#' @keywords utilities
#' @author Sudarshan A. Shetty \email{sudarshanshetty9@@gmail.com}
plot_taxa_prevalence <- function(x, level, detection = 0) {

    abundance <- NULL
    prevalence <- NULL

    x <- check_phyloseq(x, fill_na_taxa = TRUE)

    if (level == "Phylum") {
        tax.abun <- apply(abundances(x), 1, mean)
    
        #tax.prev <- rowSums(abundances(x) != 0)/nsamples(x)
        tax.prev <- prevalence(x, detection = detection)
    
        Phylum <- as.vector(data.frame(tax_table(x))$Phylum)
        Phylum <- as.vector(Phylum)
        Phylum <- as.factor(Phylum)
        xdf <- data.frame(abundance=log(tax.abun), prevalence=tax.prev, 
            Phylum)
        plot.phylum <- ggplot(xdf,
            aes(x=abundance, y=prevalence, color=Phylum)) + 
            geom_point(shape=16, alpha=0.9) +
        xlab("Average count abundance (log scale)") + 
            ylab("Taxa prevalence") +
        theme_bw() +
        theme(plot.background=element_blank(), 
                panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(), 
                panel.background=element_blank(),
        axis.text.x=element_text(angle=90, 
                    vjust=0.5, size=6)) +
        geom_vline(xintercept=0, linetype="dashed", 
                    color="grey") +
        facet_wrap(~Phylum) +
        theme(strip.background=element_rect(fill="white"))
        return(plot.phylum)
        
    } else if (level == "Family") {
        tax.abun <- apply(abundances(x), 1, mean)
        tax.prev <- rowSums(abundances(x) != 0)/nsamples(x)
        Family <- as.vector(data.frame(tax_table(x))$Family)
        Family <- as.vector(Family)
        Family <- as.factor(Family)
        xdf <- data.frame(abundance=log(tax.abun), prevalence=tax.prev, 
            Family)
        plot.fam <- ggplot(xdf, aes(x=abundance, y=prevalence,
        color=Family)) + 
            geom_point(shape=16, alpha=0.9) +
        xlab("Average count abundance (log scale)") + 
            ylab("Taxa prevalence") +
        theme_bw() +
        theme(plot.background=element_blank(), 
            panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
            panel.background=element_blank(),
        axis.text.x=element_text(angle=90, 
                vjust=0.5, size=6)) +
        geom_vline(xintercept=0, linetype="dashed", 
            color="grey") +
        facet_wrap(~Family) +
        theme(strip.background=element_rect(fill="white"))
        return(plot.fam)
    } else if (level == "Order") {
        tax.abun <- apply(abundances(x), 1, mean)
        tax.prev <- rowSums(abundances(x) != 0)/nsamples(x)
        Order <- as.vector(data.frame(tax_table(x))$Order)
        Order <- as.vector(Order)
        Order <- as.factor(Order)
        # FIXME: Most of the time log10 is used for microbial abundances (not
        # log) - should change here for consistency?
        xdf <- data.frame(abundance=log(tax.abun), prevalence=tax.prev, 
            Order)
        plot.order <- ggplot(xdf,
        aes(x=abundance, y=prevalence, color=Order)) + 
            geom_point(shape=16, alpha=0.9) +
        xlab("Average count abundance (log scale)") + 
            ylab("Taxa prevalence") +
        theme_bw() +
        theme(plot.background=element_blank(), 
            panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
            panel.background=element_blank(),
        axis.text.x=element_text(angle=90, 
                vjust=0.5, size=6)) +
        geom_vline(xintercept=0, linetype="dashed", 
            color="grey") +
        facet_wrap(~Order) +
        theme(strip.background=element_rect(fill="white"))
        return(plot.order)
    } else if (level == "Class") {
        tax.abun <- apply(abundances(x), 1, mean)
        tax.prev <- rowSums(abundances(x) != 0)/nsamples(x)
        Class <- as.vector(data.frame(tax_table(x))$Class)
        Class <- as.vector(Class)
        Class <- as.factor(Class)
        xdf <- data.frame(abundance=log(tax.abun), prevalence=tax.prev, 
            Class)
        plot.class <- ggplot(xdf,
        aes(x=abundance, y=prevalence, color=Class)) + 
            geom_point(shape=16, alpha=0.9) +
        xlab("Average count abundance (log scale)") + 
            ylab("Taxa prevalence") +
        theme_bw() +
        theme(plot.background=element_blank(), 
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(), 
            panel.background=element_blank(),
            axis.text.x=element_text(angle=90, 
                vjust=0.5, size=6)) +
        geom_vline(xintercept=0, linetype="dashed", 
                color="grey") +
        facet_wrap(~Class) +
        theme(strip.background=element_rect(fill="white"))
        return(plot.class)
    }
    
}

