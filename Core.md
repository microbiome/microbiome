<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - core}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->
Core microbiota analysis
------------------------

See also related functions for the analysis of rare and variable taxa
(noncore\_members; noncore\_abundance; rare\_members; rare\_abundance;
low\_abundance).

Load example data:

    # Load data
    library(microbiome)
    data(peerj32)

    # Rename the data
    pseq <- peerj32$phyloseq

    # Calculate compositional version of the data
    # (relative abundances)
    pseq.rel <- microbiome::transform(pseq, "compositional")

### Prevalence of taxonomic groups

Relative population frequencies; at 1% compositional abundance
threshold:

    head(prevalence(pseq.rel, detection = 1, sort = TRUE))

    ##  Yersinia et rel.  Xanthomonadaceae  Wissella et rel.            Vibrio 
    ##                 0                 0                 0                 0 
    ## Weissella et rel.       Veillonella 
    ##                 0                 0

Absolute population frequencies (sample count):

    head(prevalence(pseq.rel, detection = 1, sort = TRUE, count = TRUE))

    ##  Yersinia et rel.  Xanthomonadaceae  Wissella et rel.            Vibrio 
    ##                 0                 0                 0                 0 
    ## Weissella et rel.       Veillonella 
    ##                 0                 0

### Core microbiota analysis

If you only need the names of the core taxa, do as follows. This returns
the taxa that exceed the given prevalence and detection thresholds.

    core.taxa.standard <- core_members(pseq.rel, detection = 0, prevalence = 50/100)

A full phyloseq object of the core microbiota is obtained as follows:

    pseq.core <- core(pseq.rel, detection = 0, prevalence = .5)

Retrieving the associated taxa names from the phyloseq object:

    core.taxa <- taxa(pseq.core)

### Core abundance and diversity

Total core abundance in each sample (sum of abundances of the core
members):

    core.abundance <- sample_sums(core(pseq.rel, detection = .01, prevalence = .95))

Core visualization
------------------

### Core line plots

Determine core microbiota across various abundance/prevalence thresholds
with the blanket analysis [(Salonen et al. CMI,
2012)](http://onlinelibrary.wiley.com/doi/10.1111/j.1469-0691.2012.03855.x/abstract)
based on various signal and prevalences.

    # With compositional (relative) abundances
    det <- c(0, 0.1, 0.5, 2, 5, 20)/100
    prevalences <- seq(.05, 1, .05)
    plot_core(pseq.rel, prevalences = prevalences, detections = det, plot.type = "lineplot") + xlab("Relative Abundance (%)")

<img src="Core_files/figure-markdown_strict/core2-1.png" width="400px" />

### Core heatmaps

This visualization method has been used for instance in [Intestinal
microbiome landscaping: Insight in community assemblage and implications
for microbial modulation
strategies](https://academic.oup.com/femsre/article/doi/10.1093/femsre/fuw045/2979411/Intestinal-microbiome-landscaping-insight-in#58802539).
Shetty et al. *FEMS Microbiology Reviews* fuw045, 2017.

Note that you can order the taxa on the heatmap with the order.taxa
argument.

    # Core with compositionals:
    prevalences <- seq(.05, 1, .05)
    detections <- 10^seq(log10(1e-3), log10(.2), length = 10)

    # Also define gray color palette
    gray <- gray(seq(0,1,length=5))
    p <- plot_core(pseq.rel, plot.type = "heatmap", colours = gray,
        prevalences = prevalences, detections = detections) +
        xlab("Detection Threshold (Relative Abundance (%))")
    print(p)    


    # Same with the viridis color palette
    # color-blind friendly and uniform
    # options: viridis, magma, plasma, inferno
    # https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html
    # Also discrete=TRUE versions available
    library(viridis)

    ## Loading required package: viridisLite

    print(p + scale_fill_viridis())

    ## Scale for 'fill' is already present. Adding another scale for 'fill',
    ## which will replace the existing scale.

    # Core with absolute counts and horizontal view:
    # and minimum population prevalence (given as percentage)
    detections <- 10^seq(log10(1), log10(max(abundances(pseq))/10), length = 10)

    library(RColorBrewer)
    p <- plot_core(pseq, plot.type = "heatmap", 
                 prevalences = prevalences,
                 detections = detections,
             colours = rev(brewer.pal(5, "Spectral")),
             min.prevalence = .2, horizontal = TRUE)
    print(p)

<img src="Core_files/figure-markdown_strict/core-example3-1.png" width="200px" /><img src="Core_files/figure-markdown_strict/core-example3-2.png" width="200px" /><img src="Core_files/figure-markdown_strict/core-example3-3.png" width="200px" />
