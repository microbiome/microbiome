<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - ordination}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->
Ordination examples
-------------------

Full examples for standard ordination techniques applied to phyloseq
data, based on the [phyloseq ordination
tutorial](http://joey711.github.io/phyloseq/plot_ordination-examples.html).
For handy wrappers for some common ordination tasks in microbiome
analysis, see [landscaping examples](Landscaping.md)

Load example data:

    library(microbiome)
    library(phyloseq)
    library(ggplot2)
    data(dietswap)
    pseq <- dietswap

    # Convert to compositional data
    pseq.rel <- transform(pseq, "compositional")

    # Pick core taxa with with >10 percent prevalence of the samples
    # at a >1 percent relative abundance detection limit
    pseq.core <- core(pseq.rel, detection = 1, prevalence = 10)

### Sample ordination

Project the samples with the given method and dissimilarity measure.

    # Ordinate the data
    set.seed(423542)
    proj <- get_ordination(pseq.core, "NMDS", "bray")

Ordinate the taxa in NMDS plot with Bray-Curtis distances

    # "quiet" is here used to suppress intermediate output
    quiet(p <- plot_ordination(pseq.core, ordinate(pseq.core, "NMDS", "bray"), type = "taxa", color = "Phylum", title = "Taxa ordination"))
    print(p)

![](Ordination_files/figure-markdown_strict/ordination-pca-ordination21-1.png)

Grouping by Phyla

    p + facet_wrap(~Phylum, 5)

![](Ordination_files/figure-markdown_strict/ordination-pca-ordination22-1.png)

### Multidimensional scaling (MDS / PCoA)

    plot_ordination(pseq, ordinate(pseq, "MDS"), color = "nationality") +
                    geom_point(size = 5)

<img src="Ordination_files/figure-markdown_strict/ordination-ordinate23-1.png" width="200px" />

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

### RDA

See a separate page on [RDA](RDA.md).
