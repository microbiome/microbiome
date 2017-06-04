#' @title Visualize Prevalence Distributions for Taxa
#' @description Create taxa prevalence plots at various taxonomic levels.
#' @details This helps to obtain first insights into how is the taxa
#' distribution in the data. It also gives an idea about the taxonomic
#' affiliation of rare and abundant taxa in the data.
#' This may be helpful for data filtering or other downstream analysis.
#' @param x \code{\link{phyloseq-class}} object
#' @param level Phylum/Order/Class/Family
#' @return A \code{\link{ggplot}} plot object.
#' @export
#' @examples
#' data(DynamicsIBD)
#' # Pick data subset to speed up example
#' p0 <- subset_samples(DynamicsIBD, sex == 'male' & timepoint == 1)
#' # Check the names of the taxonomic level 
#' colnames(tax_table(p0)) 
#' # Change the "Rank" label to taxonomic levels
#' colnames(tax_table(p0)) <- c("Kingdom", "Phylum", "Class", "Order",
#'     "Family", "Genus", "Species")
#' p <- plot_taxa_prevalence(p0, 'Phylum')
#' print(p)
#' @keywords utilities
#' @author Sudarshan A. Shetty \email{sudarshanshetty9@@gmail.com}
plot_taxa_prevalence <- function(x, level) {
    abundance <- NULL
    prevalence <- NULL
    if (level == "Phylum") {
        tax.abun <- apply(otu_table(x), 1, mean)
        tax.prev <- rowSums(otu_table(x) != 0)/nsamples(x)
        Phylum <- as.vector(data.frame(tax_table(x))$Phylum)
        Phylum <- as.vector(Phylum)
        Phylum <- as.factor(Phylum)
        xdf <- data.frame(abundance=log(tax.abun), prevalence=tax.prev, 
            Phylum)
        plot.phylum <- ggplot(xdf,
            aes(x=abundance, y=prevalence, color=Phylum)) + 
            geom_point(shape=16, alpha=0.9) +
        xlab("Average relative abundance (log scale)") + 
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
        tax.abun <- apply(otu_table(x), 1, mean)
        tax.prev <- rowSums(otu_table(x) != 0)/nsamples(x)
        Family <- as.vector(data.frame(tax_table(x))$Family)
        Family <- as.vector(Family)
        Family <- as.factor(Family)
        xdf <- data.frame(abundance=log(tax.abun), prevalence=tax.prev, 
            Family)
        plot.fam <- ggplot(xdf, aes(x=abundance, y=prevalence,
        color=Family)) + 
            geom_point(shape=16, alpha=0.9) +
        xlab("Average relative abundance (log scale)") + 
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
        tax.abun <- apply(otu_table(x), 1, mean)
        tax.prev <- rowSums(otu_table(x) != 0)/nsamples(x)
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
        xlab("Average relative abundance (log scale)") + 
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
        tax.abun <- apply(otu_table(x), 1, mean)
        tax.prev <- rowSums(otu_table(x) != 0)/nsamples(x)
        Class <- as.vector(data.frame(tax_table(x))$Class)
        Class <- as.vector(Class)
        Class <- as.factor(Class)
        xdf <- data.frame(abundance=log(tax.abun), prevalence=tax.prev, 
            Class)
        plot.class <- ggplot(xdf,
        aes(x=abundance, y=prevalence, color=Class)) + 
            geom_point(shape=16, alpha=0.9) +
        xlab("Average relative abundance (log scale)") + 
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

