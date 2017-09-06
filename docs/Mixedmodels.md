<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - comparisons}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->
Mixed models for univariate comparisons
---------------------------------------

Load example data:

    # Load libraries
    library(microbiome)
    library(ggplot2)
    library(dplyr)

    # Probiotics intervention example data 
    data(peerj32) # Source: https://peerj.com/articles/32/
    pseq <- peerj32$phyloseq # Rename the example data

Abundance boxplot

    p <- boxplot_abundance(pseq, x = "time", y = "Akkermansia", line = "subject") +
        scale_y_log10()
    print(p)

<img src="Mixedmodels_files/figure-markdown_strict/boxplot2-1.png" width="300px" />

### Linear model comparison with random effect subject term

Test individual taxonomic group

    # Get sample metadata
    dfs <- meta(pseq)

    # Add abundance as the signal to model
    dfs$signal <- abundances(pseq)["Akkermansia", rownames(dfs)]

    # Paired comparison
    # with fixed group effect and random subject effect
    library(lme4)

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'Matrix'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    out <- lmer(signal ~ group + (1|subject), data = dfs)
    out0 <- lmer(signal ~ (1|subject), data = dfs)
    comp <- anova(out0, out)

    ## refitting model(s) with ML (instead of REML)

    pv <- comp[["Pr(>Chisq)"]][[2]]
    print(pv)

    ## [1] 0.4556962
