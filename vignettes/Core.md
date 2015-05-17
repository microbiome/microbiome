### Prevalence of taxonomic groups


```r
# Example data
library(microbiome)
pseq <- download_microbiome("peerj32")$physeq
```

```
## Downloading data set from Lahti et al. PeerJ, 2013: https://peerj.com/articles/32/
```

```r
# List prevalence measure for each group with a given detection threshold.
# Also sort the taxa by prevalence.
head(prevalence(pseq, detection.threshold = 10, sort = FALSE))
```

```
##             Actinomycetaceae                   Aerococcus 
##                   0.13636364                   0.25000000 
##                    Aeromonas                  Akkermansia 
##                   0.31818182                   1.00000000 
## Alcaligenes faecalis et rel.           Allistipes et rel. 
##                   0.04545455                   1.00000000
```

```r
# Just list the names of taxa that are present over abundance threshold 2
# in over 20 percent of the samples:
prevalent.taxa <- prevalent_taxa(pseq, detection.threshold = 50, prevalence.threshold = 0.2)
```


### Core microbiota

Determine core microbiota with the [blanket
analysis](http://onlinelibrary.wiley.com/doi/10.1111/j.1469-0691.2012.03855.x/abstract)
based on various signal and prevalence thresholds.
 

```r
core <- core_matrix(pseq, prevalence.intervals = seq(10, 100, 10), intensity.intervals = c(0, 10^(0:4)))
```

```
## Error in core_matrix(pseq, prevalence.intervals = seq(10, 100, 10), intensity.intervals = c(0, : unused argument (intensity.intervals = c(0, 10^(0:4)))
```

Two alternative ways to visualize the core microbiota:


```r
# Core 2D line plots
p <- plot_core(pseq, prevalence.intervals = seq(10, 100, 10), detection.thresholds = c(0, 10^(0:4)), plot.type = "lineplot")
```

![plot of chunk core-example2](figure/core-example2-1.png) 

```r
# Core heatmap
p <- plot_core(pseq, prevalence.intervals = seq(10, 100, 10), detection.thresholds = c(0, 2^(0:14)), plot.type = "heatmap")
```

![plot of chunk core-example2](figure/core-example2-2.png) 
