<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - composition}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->

Also see [phyloseq barplot
examples](http://joey711.github.io/phyloseq/plot_bar-examples.html).

Read example data from a [diet swap
study](http://dx.doi.org/10.1038/ncomms7342):

    # Example data
    library(microbiome)
    library(dplyr)
    data(dietswap)

    # Just use prevalent taxa to speed up examples
    # (not absolute counts used in this example)
    pseq <- core(dietswap, detection = 10^3, prevalence = 95/100)

    # Pick sample subset
    library(phyloseq)
    pseq2 <- subset_samples(pseq, group == "DI" & nationality == "AFR" & timepoint.within.group == 1)

### Barplots for composition

Same with compositional (relative) abundances; for each sample (left),
or averafged by group (right).

    # Try another theme
    # from https://github.com/hrbrmstr/hrbrthemes
    library(hrbrthemes)
    library(gcookbook)
    library(tidyverse)

    # Limit the analysis on core taxa and specific sample group
    p <- plot_composition(pseq2,
                  taxonomic.level = "OTU",
                          sample.sort = "nationality",
                          x.label = "nationality",
                          transform = "compositional") +
         guides(fill = guide_legend(ncol = 1)) +
         scale_y_percent() +
         labs(x = "Samples", y = "Relative abundance (%)",
                                       title = "Relative abundance data",
                                       subtitle = "Subtitle",
                                       caption = "Caption text.") + 
         theme_ipsum(grid="Y")
    print(p)  

    # Averaged by group
    p <- plot_composition(pseq2,
                          average_by = "bmi_group", transform = "compositional")
    print(p)

<img src="Composition_files/figure-markdown_strict/composition-example4b-1.png" width="400px" /><img src="Composition_files/figure-markdown_strict/composition-example4b-2.png" width="400px" />

### Composition heatmaps

Heatmap for CLR-transformed abundances, with samples and OTUs sorted
with the neatmap method:

    tmp <- plot_composition(pseq2, plot.type = "heatmap", transform = "clr",
                          sample.sort = "neatmap", otu.sort = "neatmap",
                          mar = c(6, 13, 1, 1))

    ## Ordering columns

    ## Ordering rows

### Plot taxa prevalence

We use the Dynamics IBD data set from [Halfvarson J., et al. Nature
Microbiology, 2017](http://www.nature.com/articles/nmicrobiol20174) as
downloaded from [Qiita ID
1629](https://qiita.ucsd.edu/study/description/1629). This function
allows you to have an overview of OTU prevalences alongwith their
taxonomic affiliations. This will aid in checking if you filter OTUs
based on prevalence, then what taxonomic affliations will be lost.

    data(DynamicsIBD)

    # Use sample and taxa subset to speed up example
    p0 <- subset_samples(DynamicsIBD, sex == "male" & timepoint == 1)

    # Check the names of the taxonomic level 
    colnames(tax_table(p0)) 

    ## [1] "Rank1" "Rank2" "Rank3" "Rank4" "Rank5" "Rank6" "Rank7"

    # Change the "Rank" label to taxonomic levels
    colnames(tax_table(p0)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

    p0 <- core(p0, detection = 10, prevalence = 0)

    # For the available taxonomic levels
    plot_taxa_prevalence(p0, "Phylum")

![](Composition_files/figure-markdown_strict/plot_prev-1.png)
