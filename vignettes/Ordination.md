<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - ordination}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->
Ordination examples
-------------------

Some examples with HITChip data. See also [phyloseq ordination
tutorial](http://joey711.github.io/phyloseq/plot_ordination-examples.html).

Load example data:

    library(microbiome)
    library(phyloseq)
    library(ggplot2)
    data(dietswap)
    pseq <- dietswap

    # Convert signal to compositionals
    pseq.rel <- transform_phyloseq(pseq, "compositional")

    # Pick OTUs that are present with >1 percent compositional 
    # in >10 percent of the samples
    pseq2 <- filter_taxa(pseq.rel,
               function(x) sum(x > 1) > (0.1*nsamples(pseq.rel)), TRUE)

### Sample ordination

Project the samples with the given method and distance. See also
plot\_ordination from the phyloseq package.

    set.seed(423542)
    pseq.ord <- ordinate(pseq2, "NMDS", "bray")
    # Just pick the projected data (first two columns + metadata)
    proj <- plot_ordination(pseq2, pseq.ord, justDF = T)

Then visualize the projected data:

    # Highlighting nationality
    p <- microbiome::densityplot(proj[, 1:2], col = proj$nationality, legend = T)
    print(p)

    # Projection with sample names:
    ax1 <- names(proj)[[1]]
    ax2 <- names(proj)[[2]]
    p <- ggplot(aes_string(x = ax1, y = ax2, label = "sample"), data = proj) +
           geom_text(size = 2)
    print(p)

<img src="Ordination_files/figure-markdown_strict/ordination4-1.png" width="400px" /><img src="Ordination_files/figure-markdown_strict/ordination4-2.png" width="400px" />

Ordinate the taxa in NMDS plot with Bray-Curtis distances

    p <- plot_ordination(pseq2, pseq.ord, type = "taxa", color = "Phylum", title = "Taxa ordination")
    print(p)

![](Ordination_files/figure-markdown_strict/ordination-pca-ordination21-1.png)

Grouping by Phyla

    p + facet_wrap(~Phylum, 5)

![](Ordination_files/figure-markdown_strict/ordination-pca-ordination22-1.png)

### Multidimensional scaling (MDS / PCoA)

    plot_ordination(pseq, ordinate(pseq, "MDS"), color = "nationality") +
                    geom_point(size = 5)

![](Ordination_files/figure-markdown_strict/ordination-ordinate23-1.png)

### RDA

See a separate page on [RDA](RDA.md).

### Canonical correspondence analysis (CCA)

    # With samples
    p <- plot_ordination(pseq, ordinate(pseq, "CCA"),
           type = "samples", color = "nationality")
    p <- p + geom_point(size = 4)
    print(p)

    # With taxa:
    p <- plot_ordination(pseq, ordinate(pseq, "CCA"),
           type = "taxa", color = "Phylum")
    p <- p + geom_point(size = 4)
    print(p)

<img src="Ordination_files/figure-markdown_strict/ordination-ordinate24a-1.png" width="400px" /><img src="Ordination_files/figure-markdown_strict/ordination-ordinate24a-2.png" width="400px" />

### Split plot

    plot_ordination(pseq, ordinate(pseq, "CCA"),
                  type = "split", shape = "nationality", 
                      color = "Phylum", label = "nationality")

![](Ordination_files/figure-markdown_strict/ordination-ordinate25-1.png)

### Ordination biplot

    plot_ordination(pseq, ordinate(pseq, "CCA"), type = "biplot", color = "Phylum")

![](Ordination_files/figure-markdown_strict/ordination-ordinate26-1.png)
