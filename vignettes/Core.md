### Prevalence of taxonomic groups



```r
# Load example data
library(microbiome)
data("peerj32")
pseq <- peerj32$phyloseq

# Calculate relative abundances
pseq.rel <- transform_phyloseq(pseq, "relative.abundance", "OTU")
```


List prevalence for each group at 1 percent relative abundance abundance threshold:


```r
head(prevalence(pseq.rel, detection.threshold = 1, sort = FALSE))
```

```
## sample-1 sample-2 sample-3 sample-4 sample-5 sample-6 
## 20.76923 19.23077 20.00000 20.00000 20.00000 20.76923
```


List the taxa that are present at the given detection threshold (1% relative abundance) at a given prevalence (80%) (fraction of the samples):


```r
prevalent.taxa <- prevalent_taxa(pseq.rel, detection.threshold = 1, prevalence.threshold = 80)
```


### Core microbiota

Determine core microbiota with the [blanket
analysis](http://onlinelibrary.wiley.com/doi/10.1111/j.1469-0691.2012.03855.x/abstract)
based on various signal and prevalence thresholds. See also the the
bootstrap_microbes function.
 

```r
det <- c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20)
prev <- seq(10, 100, 10)
core <- core_matrix(pseq.rel, prevalence.intervals = prev, detection.thresholds = det)
```

### Core 2D line plots


```r
# Core lineplot with absolute read counts
det <- c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 1e4)
res <- plot_core(pseq, prevalence.intervals = prev, detection.thresholds = det, plot.type = "lineplot", plot = FALSE)
```

```
## Error in plot_core(pseq, prevalence.intervals = prev, detection.thresholds = det, : unused argument (plot = FALSE)
```

```r
res$plot + xlab("Abundance (OTU read count)")
```

```
## Error in res$plot + xlab("Abundance (OTU read count)"): non-numeric argument to binary operator
```

```r
# Core lineplot with relative abundances
det <- c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20)
res <- plot_core(pseq.rel, prevalence.intervals = prev, detection.thresholds = det, plot.type = "lineplot", plot = FALSE)
```

```
## Error in plot_core(pseq.rel, prevalence.intervals = prev, detection.thresholds = det, : unused argument (plot = FALSE)
```

```r
res$plot + xlab("Relative Abundance (%)")
```

```
## Error in res$plot + xlab("Relative Abundance (%)"): non-numeric argument to binary operator
```

```r
# Retrieve the core count data matrix
coremat <- res$data
print(coremat)
```

```
## NULL
```


### Core heatmaps




```r
# Core with relative abundances:
prevalence.intervals <- seq(5, 100, 5)
detection.thresholds <- 10^seq(log10(1e-3), log10(20), length = 20)		 
res <- plot_core(pseq.rel, plot.type = "heatmap", palette = "bw", prevalence.intervals = prevalence.intervals, detection.thresholds = detection.thresholds, plot = FALSE) 
```

```
## Error in plot_core(pseq.rel, plot.type = "heatmap", palette = "bw", prevalence.intervals = prevalence.intervals, : unused argument (plot = FALSE)
```

```r
print(res$plot + xlab("Detection Threshold (Relative Abundance (%))"))
```

```
## Error in print(res$plot + xlab("Detection Threshold (Relative Abundance (%))")): error in evaluating the argument 'x' in selecting a method for function 'print': Error in res$plot + xlab("Detection Threshold (Relative Abundance (%))") : 
##   non-numeric argument to binary operator
```

```r
# Core with absolute counts:
prevalence.intervals = seq(5, 100, 5)
detection.thresholds <- 10^seq(log10(1), log10(max(otu_table(pseq))/10), length = 20)		 
res <- plot_core(pseq, plot.type = "heatmap", palette = "bw", prevalence.intervals = prevalence.intervals,
       		       detection.thresholds = detection.thresholds, min.prevalence = NULL)$plot
```


Zoom in on the core region by filtering out rows and columns not passing min prevalence (given as percentages):


```r
res <- plot_core(pseq, plot.type = "heatmap", palette = "bw", prevalence.intervals = prevalence.intervals,
		detection.thresholds = detection.thresholds, min.prevalence = 10, plot = TRUE)
```

```
## Error in plot_core(pseq, plot.type = "heatmap", palette = "bw", prevalence.intervals = prevalence.intervals, : unused argument (plot = TRUE)
```

```r
res <- plot_core(pseq, plot.type = "heatmap", palette = "spectral", prevalence.intervals = prevalence.intervals,
		detection.thresholds = detection.thresholds, min.prevalence = 0, plot = TRUE)
```

```
## Error in plot_core(pseq, plot.type = "heatmap", palette = "spectral", : unused argument (plot = TRUE)
```



Retrieve the core prevalence data matrix


```r
prevalences <- res$data
kable(head(prevalences))
```



|Taxa                         | DetectionThreshold| Prevalence|
|:----------------------------|------------------:|----------:|
|Actinomycetaceae             |                  1|   61.36364|
|Aerococcus                   |                  1|   61.36364|
|Aeromonas                    |                  1|   65.90909|
|Akkermansia                  |                  1|  100.00000|
|Alcaligenes faecalis et rel. |                  1|   13.63636|
|Allistipes et rel.           |                  1|  100.00000|

