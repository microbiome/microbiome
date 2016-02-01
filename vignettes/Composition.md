## Microbiota composition


### Barplots for composition

Also see [phyloseq barplot examples](http://joey711.github.io/phyloseq/plot_bar-examples.html) and [HITChip Barplots](Barplots.md)


Read example data from a [diet swap study](http://dx.doi.org/10.1038/ncomms7342):


```r
# Example data
library(microbiome)
data(dietswap)
pseq <- dietswap
```

Show OTU absolute abundance in each sample. Plot absolute taxon
abundances (Note: on HITChip data the Phylum level is only
approximate):


```r
# Pick sample subset
library(phyloseq)
pseq2 <- subset_samples(pseq, group == "DI" & nationality == "AFR")
plot_composition(pseq2, taxonomic.level = "Phylum")
```

```
## Aggregating the taxa.
```

```
## Check data transformations.
```

```
## Prepare data.frame.
```

```
## Construct the plots
```

```
## $plot
```

![plot of chunk composition-example1b](figure/composition-example1b-1.png)

```
## 
## $data
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 22 taxa and 32 samples ]
## sample_data() Sample Data:       [ 32 samples by 8 sample variables ]
## tax_table()   Taxonomy Table:    [ 22 taxa by 2 taxonomic ranks ]
```

Arrange by sample variable and use custom X axis labels. Africans have more Prevotella as expected. Absolute counts:


```r
# Pick the top OTUs only
top <- names(sort(taxa_sums(pseq), TRUE)[1:10])
pseq3 <- prune_taxa(top, pseq)
pseq3 <- subset_samples(pseq3, group == "DI" & timepoint.within.group == 1)

res <- plot_composition(pseq3, sample.sort = "nationality", x.label = "nationality")
```

```
## Check data transformations.
```

```
## Prepare data.frame.
```

```
## Construct the plots
```

```r
p <- res$plot
p <- p + guides(fill = guide_legend(ncol = 1))
p <- p + theme(legend.position = "bottom")
print(p)
```

![plot of chunk composition-example4](figure/composition-example4-1.png)


Same with relative abundances:


```r
res <- plot_composition(pseq3, sample.sort = "nationality", x.label = "nationality", transformation = "relative.abundance")
```

```
## Check data transformations.
```

```
## Prepare data.frame.
```

```
## Construct the plots
```

```r
p <- res$plot
p <- p + guides(fill = guide_legend(ncol = 1))
p <- p + theme(legend.position = "bottom")
p <- p + ylab("Relative Abundance (%)")
print(p)
```

![plot of chunk composition-example4b](figure/composition-example4b-1.png)



### Heatmaps for composition


Plain heatmap


```r
res <- plot_composition(pseq3, plot.type = "heatmap")
```

```
## Check data transformations.
```

```
## Prepare data.frame.
```

```
## Construct the plots
```

```
## Constructing the heatmap.
```

![plot of chunk composition-example5](figure/composition-example5-1.png)


Heatmap with Z-transformed OTUs


```r
res <- plot_composition(pseq3, plot.type = "heatmap", transformation = "Z-OTU")
```

```
## Check data transformations.
```

```
## Prepare data.frame.
```

```
## Construct the plots
```

```
## Constructing the heatmap.
```

![plot of chunk composition-example6](figure/composition-example6-1.png)


Same, but samples and OTUs sorted with the neatmap method


```r
res <- plot_composition(pseq3, plot.type = "heatmap", transformation = "Z-OTU", sample.sort = "neatmap", otu.sort = "neatmap")
```

```
## Check data transformations.
```

```
## Prepare data.frame.
```

```
## Construct the plots
```

```
## Constructing the heatmap.
```

![plot of chunk composition-example7](figure/composition-example7-1.png)


Same, but samples and OTUs sorted manually


```r
sample.sort <- order_neatmap(pseq3, method = "NMDS", distance = "bray", target = "sites", first = NULL) 
otu.sort <- order_neatmap(pseq3, method = "NMDS", distance = "bray", target = "species", first = NULL)
res <- plot_composition(pseq3, plot.type = "heatmap", transformation = "Z-OTU", sample.sort = sample.sort, otu.sort = otu.sort)
```

```
## Warning in if (otu.sort == "neatmap") {: the condition has length > 1 and
## only the first element will be used
```

```
## Check data transformations.
```

```
## Prepare data.frame.
```

```
## Construct the plots
```

```
## Constructing the heatmap.
```

![plot of chunk composition-example8](figure/composition-example8-1.png)


