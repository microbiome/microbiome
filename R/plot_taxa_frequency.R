#' @title Plotting the Taxa Frequency
#' @description A simplified wrapper for creating taxa frequency plots at Phylum/Order/Class/Family level.
#' @details This plot is useful to have a first insights into how is the taxa distribution in your data.
#'          This also gives an idea about the taxonomic affiliation of rare and abundant taxa in your dataset. This may be helpful in some data filtering or other downstream analysis.
#' @param x \code{\link{phyloseq-class}} object
#' @param level Phylum/Order/Class/Family
#' @return A \code{\link{ggplot}} plot object.
#' @export
#' @examples \dontrun{
#'  library(microbiome)
#'  data(DynamicsIBD)
#'  p0 <- DynamicsIBD
#'  p0.f <- format_phyloseq(p0)
#'  plot <- plot_taxa_frequency(p0.f, 'Phylum')
#'     }
#' @keywords utilities
#' @author Sudarshan A. Shetty \email{sudarshanshetty9@@gmail.com}
plot_taxa_frequency <- function(x, level) {
    if (level == "Phylum") {
        tax.abun <- apply(otu_table(x), 1, mean)
        tax.freq <- rowSums(otu_table(x) != 0)/nsamples(x)
        Phylum <- as.vector(data.frame(tax_table(x))$Phylum)
        Phylum <- as.vector(Phylum)
        Phylum <- as.factor(Phylum)
        xdf <- data.frame(abundance = log(tax.abun), frequency = tax.freq, 
            Phylum)
        plot.phylum <- ggplot(xdf, aes(x = abundance, y = frequency, color = Phylum)) + 
            geom_point(shape = 16, alpha = 0.9) + xlab("Average relative abundance (log scale)") + 
            ylab("Taxa Frequency") + theme_bw() + theme(plot.background = element_blank(), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.text.x = element_text(angle = 90, 
                vjust = 0.5, size = 6)) + geom_vline(xintercept = 0, linetype = "dashed", 
            color = "grey") + facet_wrap(~Phylum) + theme(strip.background = element_rect(fill = "white"))
        return(plot.phylum)
        
    } else if (level == "Family") {
        tax.abun <- apply(otu_table(x), 1, mean)
        tax.freq <- rowSums(otu_table(x) != 0)/nsamples(x)
        Family <- as.vector(data.frame(tax_table(x))$Family)
        Family <- as.vector(Family)
        Family <- as.factor(Family)
        xdf <- data.frame(abundance = log(tax.abun), frequency = tax.freq, 
            Family)
        plot.fam <- ggplot(xdf, aes(x = abundance, y = frequency, color = Family)) + 
            geom_point(shape = 16, alpha = 0.9) + xlab("Average relative abundance (log scale)") + 
            ylab("Taxa Frequency") + theme_bw() + theme(plot.background = element_blank(), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.text.x = element_text(angle = 90, 
                vjust = 0.5, size = 6)) + geom_vline(xintercept = 0, linetype = "dashed", 
            color = "grey") + facet_wrap(~Family) + theme(strip.background = element_rect(fill = "white"))
        return(plot.fam)
    } else if (level == "Order") {
        tax.abun <- apply(otu_table(x), 1, mean)
        tax.freq <- rowSums(otu_table(x) != 0)/nsamples(x)
        Order <- as.vector(data.frame(tax_table(x))$Order)
        Order <- as.vector(Order)
        Order <- as.factor(Order)
        xdf <- data.frame(abundance = log(tax.abun), frequency = tax.freq, 
            Order)
        plot.order <- ggplot(xdf, aes(x = abundance, y = frequency, color = Order)) + 
            geom_point(shape = 16, alpha = 0.9) + xlab("Average relative abundance (log scale)") + 
            ylab("Taxa Frequency") + theme_bw() + theme(plot.background = element_blank(), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.text.x = element_text(angle = 90, 
                vjust = 0.5, size = 6)) + geom_vline(xintercept = 0, linetype = "dashed", 
            color = "grey") + facet_wrap(~Order) + theme(strip.background = element_rect(fill = "white"))
        return(plot.order)
    } else if (level == "Class") {
        tax.abun <- apply(otu_table(x), 1, mean)
        tax.freq <- rowSums(otu_table(x) != 0)/nsamples(x)
        Class <- as.vector(data.frame(tax_table(x))$Class)
        Class <- as.vector(Class)
        Class <- as.factor(Class)
        xdf <- data.frame(abundance = log(tax.abun), frequency = tax.freq, 
            Class)
        plot.class <- ggplot(xdf, aes(x = abundance, y = frequency, color = Class)) + 
            geom_point(shape = 16, alpha = 0.9) + xlab("Average relative abundance (log scale)") + 
            ylab("Taxa Frequency") + theme_bw() + theme(plot.background = element_blank(), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.text.x = element_text(angle = 90, 
                vjust = 0.5, size = 6)) + geom_vline(xintercept = 0, linetype = "dashed", 
            color = "grey") + facet_wrap(~Class) + theme(strip.background = element_rect(fill = "white"))
        return(plot.class)
    }
    
}
