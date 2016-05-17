### Prevalence of taxonomic groups

Load example data:


```r
library(microbiome)
data("peerj32")
pseq <- peerj32$phyloseq

# Calculate relative abundances
pseq.rel <- transform_phyloseq(pseq, "relative.abundance", "OTU")
```

Taxa prevalences at 1 percent relative abundance abundance threshold:


```r
head(prevalence(pseq.rel, detection.threshold = 1, sort = FALSE))
```

```
##             Actinomycetaceae                   Aerococcus 
##                      0.00000                      0.00000 
##                    Aeromonas                  Akkermansia 
##                      0.00000                     52.27273 
## Alcaligenes faecalis et rel.           Allistipes et rel. 
##                      0.00000                     34.09091
```


### Core microbiota

List the core taxa above a given detection and prevalence thresholds:


```r
core.taxa <- core(pseq.rel, detection.threshold = 1, prevalence.threshold = 95)
```

Total core abundance:


```r
core.abundance <- core_abundance(pseq.rel, detection.threshold = 1, prevalence.threshold = 95)
```

```
## Error in rownames(meta2): object 'meta2' not found
```


### Core 2D line plots

Determine core microbiota across various abundance/prevalence
thresholds with the [blanket
analysis](http://onlinelibrary.wiley.com/doi/10.1111/j.1469-0691.2012.03855.x/abstract) based on various signal and prevalence thresholds. See also the the
bootstrap_microbes function.


```r
# With absolute read counts
det <- c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 1e4)
res <- plot_core(pseq, prevalence.intervals = prev, detection.thresholds = det, plot.type = "lineplot")
res$plot + xlab("Abundance (OTU read count)")

# With relative abundances
det <- c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20)
res <- plot_core(pseq.rel, prevalence.intervals = prev, detection.thresholds = det, plot.type = "lineplot")
res$plot + xlab("Relative Abundance (%)")

# Retrieve the core count data matrix
coremat <- res$data
```

<img src="figure/core-example2-1.png" title="plot of chunk core-example2" alt="plot of chunk core-example2" width="430px" /><img src="figure/core-example2-2.png" title="plot of chunk core-example2" alt="plot of chunk core-example2" width="430px" />


### Core heatmaps


```r
# Core with relative abundances:
prevalence.intervals <- seq(5, 100, 5)
detection.thresholds <- 10^seq(log10(1e-3), log10(20), length = 20)
# Also define gray color palette
gray <- gray(seq(0,1,length=5))
res <- plot_core(pseq.rel, plot.type = "heatmap", colours = gray,
    prevalence.intervals = prevalence.intervals, detection.thresholds = detection.thresholds) 
print(res$plot + xlab("Detection Threshold (Relative Abundance (%))"))

res <- plot_core(pseq.rel, plot.type = "heatmap", colours = gray, 
    prevalence.intervals = prevalence.intervals, detection.thresholds = detection.thresholds) 

# Core with absolute counts:
prevalence.intervals = seq(5, 100, 5)
detection.thresholds <- 10^seq(log10(1), log10(max(otu_table(pseq))/10), length = 20)		 
plot_core(pseq, plot.type = "heatmap", colours = gray,
       		 prevalence.intervals = prevalence.intervals,
       		 detection.thresholds = detection.thresholds,
		 min.prevalence = NULL)$plot
```

<img src="figure/core-example3-1.png" title="plot of chunk core-example3" alt="plot of chunk core-example3" width="430px" /><img src="figure/core-example3-2.png" title="plot of chunk core-example3" alt="plot of chunk core-example3" width="430px" />


Zoom in on the core region by filtering out rows and columns not passing min prevalence (given as percentages):


```r
res <- plot_core(pseq, plot.type = "heatmap", colours = gray,
                 prevalence.intervals = prevalence.intervals,
		 detection.thresholds = detection.thresholds,
		 min.prevalence = 10)
print(res$plot)		 


library(RColorBrewer)
res <- plot_core(pseq, plot.type = "heatmap",
                 #colours = colorRampPalette(rev(brewer.pal(11, "Spectral")))(5),
		 colours = rev(brewer.pal(5, "Spectral")),
                 prevalence.intervals = prevalence.intervals,
		 detection.thresholds = detection.thresholds,
		 min.prevalence = 0)		 
print(res$plot)		 
```

<img src="figure/core-example3bb-1.png" title="plot of chunk core-example3bb" alt="plot of chunk core-example3bb" width="430px" /><img src="figure/core-example3bb-2.png" title="plot of chunk core-example3bb" alt="plot of chunk core-example3bb" width="430px" />


Retrieve the core prevalence data matrix


```r
prevalences <- res$data
kable(head(prevalences), digits = 2)
```



|Taxa                         | DetectionThreshold| Prevalence|
|:----------------------------|------------------:|----------:|
|Actinomycetaceae             |                  1|      61.36|
|Aerococcus                   |                  1|      61.36|
|Aeromonas                    |                  1|      65.91|
|Akkermansia                  |                  1|     100.00|
|Alcaligenes faecalis et rel. |                  1|      13.64|
|Allistipes et rel.           |                  1|     100.00|

