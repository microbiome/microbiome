---
title: "Comparisons of microbiome community composition"
bibliography: 
- bibliography.bib
- references.bib
output: 
  prettydoc::html_pretty:
    theme: cayman
    highlight: github
---
<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - comparisons}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->


## Mixed models for univariate comparisons


Load example data:


```r
# Load libraries
library(microbiome)
library(ggplot2)
library(dplyr)

# Probiotics intervention example data 
data(peerj32) # Source: https://peerj.com/articles/32/
pseq <- peerj32$phyloseq # Rename the example data
```


Abundance boxplot


```r
p <- boxplot_abundance(pseq, x = "time", y = "Akkermansia", line = "subject", color = "gender") + scale_y_log10()
print(p)
```

<img src="figure/boxplot2-1.png" title="plot of chunk boxplot2" alt="plot of chunk boxplot2" width="300px" />


### Linear model comparison with random effect subject term

Test individual taxonomic group


```r
# Get sample metadata
dfs <- meta(pseq)

# Add abundance as the signal to model
dfs$signal <- abundances(pseq)["Akkermansia", rownames(dfs)]

# Paired comparison
# with fixed group effect and random subject effect
library(lme4)
```

```
## Loading required package: Matrix
```

```r
out <- lmer(signal ~ group + (1|subject), data = dfs)
out0 <- lmer(signal ~ (1|subject), data = dfs)
comp <- anova(out0, out)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
pv <- comp[["Pr(>Chisq)"]][[2]]
print(pv)
```

```
## [1] 0.4556962
```


```
## [1] "maptree"
## [1] "Biobase"
## [1] "dynamicTreeCut"
## [1] "tidyr"
## [1] "jsonlite"
## [1] "splines"
## [1] "foreach"
## [1] "moments"
## [1] "Formula"
## [1] "assertthat"
## [1] "highr"
## [1] "stats4"
## [1] "phyloseq"
## [1] "latticeExtra"
## [1] "grDevices"
## [1] "tensorA"
## [1] "robustbase"
## [1] "impute"
## [1] "RSQLite"
## [1] "backports"
## [1] "lattice"
## [1] "base"
## [1] "digest"
## [1] "RColorBrewer"
## [1] "XVector"
## [1] "checkmate"
## [1] "minqa"
## [1] "colorspace"
## [1] "htmltools"
## [1] "preprocessCore"
## [1] "Matrix"
## [1] "plyr"
## [1] "microbiome"
## [1] "zlibbioc"
## [1] "tgp"
## [1] "GO.db"
## [1] "scales"
## [1] "lme4"
## [1] "htmlTable"
## [1] "tibble"
## [1] "mgcv"
## [1] "datasets"
## [1] "IRanges"
## [1] "ggplot2"
## [1] "fastcluster"
## [1] "nnet"
## [1] "BiocGenerics"
## [1] "lazyeval"
## [1] "survival"
## [1] "magrittr"
## [1] "memoise"
## [1] "evaluate"
## [1] "methods"
## [1] "doParallel"
## [1] "nlme"
## [1] "MASS"
## [1] "foreign"
## [1] "compositions"
## [1] "utils"
## [1] "vegan"
## [1] "tools"
## [1] "data.table"
## [1] "energy"
## [1] "matrixStats"
## [1] "stringr"
## [1] "S4Vectors"
## [1] "munsell"
## [1] "cluster"
## [1] "AnnotationDbi"
## [1] "stats"
## [1] "Biostrings"
## [1] "ade4"
## [1] "nloptr"
## [1] "rhdf5"
## [1] "grid"
## [1] "iterators"
## [1] "biomformat"
## [1] "graphics"
## [1] "htmlwidgets"
## [1] "WGCNA"
## [1] "igraph"
## [1] "base64enc"
## [1] "boot"
## [1] "gtable"
## [1] "codetools"
## [1] "multtest"
## [1] "DBI"
## [1] "R6"
## [1] "reshape2"
## [1] "bayesm"
## [1] "gridExtra"
## [1] "knitr"
## [1] "dplyr"
## [1] "Hmisc"
## [1] "permute"
## [1] "ape"
## [1] "stringi"
## [1] "parallel"
## [1] "Rcpp"
## [1] "rpart"
## [1] "acepack"
## [1] "DEoptimR"
```



