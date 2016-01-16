### Prevalence of taxonomic groups


```r
# Load example data
library(microbiome)
pseq <- download_microbiome("peerj32")$physeq
```

```
## Downloading data set from Lahti et al. PeerJ, 2013: https://peerj.com/articles/32/
```

List prevalence measure for each group with a given detection threshold:


```r
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

List the taxa that are present over detection threshold in given
fraction of the samples:


```r
prevalent.taxa <- prevalent_taxa(pseq, detection.threshold = 50, prevalence.threshold = 0.2)
```


### Core microbiota

Determine core microbiota with the [blanket
analysis](http://onlinelibrary.wiley.com/doi/10.1111/j.1469-0691.2012.03855.x/abstract)
based on various signal and prevalence thresholds.
 

```r
core <- core_matrix(pseq, prevalence.intervals = seq(10, 100, 10), detection.thresholds = c(0, 10^(0:4)))
```

### Core 2D line plots


```r
p <- plot_core(pseq, prevalence.intervals = seq(10, 100, 10), detection.thresholds = c(0, 10^(0:4)), plot.type = "lineplot")
```

![plot of chunk core-example2](figure/core-example2-1.png)

### Core heatmaps


```r
p <- plot_core(pseq, plot.type = "heatmap", palette = "bw")
```

![plot of chunk core-example3](figure/core-example3-1.png)


```r
p <- plot_core(pseq, plot.type = "heatmap", palette = "spectral")
```

![plot of chunk core-example4](figure/core-example4-1.png)

