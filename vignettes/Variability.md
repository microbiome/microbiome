<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - variability}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->
### Inter-individual homogeneity (within group of samples)

Assess 'inter-individual stability', or homogeneity, as in [Salonen et
al. ISME J
2014](http://www.nature.com/ismej/journal/v8/n11/full/ismej201463a.html).
This is defined as the average correlation between the samples and their
mean for a given samples vs phylotypes matrix. For illustration,
calculate inter-individual homogeneity separately for Placebo and LGG
groups. Note that this homogeneity measure is affected by sample size.

Load example data

    library(microbiome)
    data("dietswap")
    x <- dietswap

    # Add time field (two time points needed within each group for the 
    # intraindividual method)
    sample_data(x)$time <- sample_data(x)$timepoint.within.group

Heterogeneity across subjects within a group

    res <- estimate_homogeneity(x, "interindividual")

Visualize

    library(ggplot2)
    theme_set(theme_bw(20))
    p <- ggplot(res$data, aes(x = group, y = correlation))
    p <- p + geom_boxplot()
    p <- p + ggtitle(paste("Inter-individual homogeneity (p=", round(res$p.value, 6), ")", sep = ""))
    p <- p + ylab("Correlation")
    print(p)

![](Variability_files/figure-markdown_strict/homogeneity-example2d-1.png)

### Intra-individual stability

Homogeneity within subjects over time (also called intra-individual
stability in [Salonen et al. ISME J
2014](http://www.nature.com/ismej/journal/v8/n11/full/ismej201463a.html)).
Defined as the average correlation between two time points within
subjects within each group. For illustration, check intra-individual
stability (homogeneity) separately for Placebo and LGG groups.

    res <- estimate_homogeneity(x, "intraindividual")

Visualize

    library(ggplot2)
    theme_set(theme_bw(20))
    p <- ggplot(res$data, aes(x = group, y = correlation))
    p <- p + geom_boxplot()
    p <- p + ggtitle(paste("Intra-individual homogeneity (p=", round(res$p.value, 6), ")"))
    p <- p + ylab("Correlation")
    print(p)

![](Variability_files/figure-markdown_strict/homogeneity-intra-1.png)

### Time series

    data("atlas1006")
    pseq <- atlas1006
    pseq <- subset_samples(pseq, DNA_extraction_method == "r")
    pseq <- transform_phyloseq(pseq, "compositional")
    p <- plot_timeseries(pseq, "Dialister", subject = "831", tipping.point = 0.5)
    print(p)

![](Variability_files/figure-markdown_strict/homogeneity-timeseries-1.png)

Pick samples at the baseline time points only:

    data("atlas1006")
    pseq0 <- pick_baseline(atlas1006)

Further visualization tools
---------------------------

Draw regression curve with smoothed error bars based on the
[Visually-Weighted
Regression](http://www.fight-entropy.com/2012/07/visually-weighted-regression.html)
by Solomon M. Hsiang. The sorvi implementation extends [Felix
Schonbrodt's original
code](http://www.nicebread.de/visually-weighted-watercolor-plots-new-variants-please-vote/).

    data(atlas1006)
    p <- plot_regression(diversity ~ age, sample_data(atlas1006))
    print(p)

![](Variability_files/figure-markdown_strict/variability-regression-1.png)
