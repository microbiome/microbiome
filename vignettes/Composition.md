## Microbiota composition


### Barplot visualizations

Also see [phyloseq barplot examples](http://joey711.github.io/phyloseq/plot_bar-examples.html) and [HITChip Barplots](Barplots.md)


Load example data:


```r
# Example data
library(microbiome)
pseq0 <- download_microbiome("dietswap")
```

```
## Downloading data set from O'Keefe et al. Nat. Comm. 6:6342, 2015 from Data Dryad: http://datadryad.org/resource/doi:10.5061/dryad.1mn1n
```

```r
# Pick sample subset
pseq <- subset_samples(pseq0, group == "DI" & nationality == "AFR")
```

Show OTU absolute abundance in each sample. Plot absolute taxon
abundances (Note: on HITChip data the Phylum level is only
approximate):


```r
plot_composition(pseq, taxonomic.level = "Phylum")
```

![plot of chunk composition-example1b](figure/composition-example1b-1.png) 

Same with relative abundances:


```r
p <- plot_composition(pseq, taxonomic.level = "Phylum", relative.abundance = TRUE)
p <- p + guides(fill = guide_legend(nrow = 12, byrow = TRUE))
p <- p + theme(legend.position = "bottom")
print(p)
```

![plot of chunk composition-example3](figure/composition-example3-1.png) 


Arrange by sample variable and use custom X axis labels. Africans have more Prevotella as expected:


```r
# Subset taxa and samples
pseq <- subset_samples(pseq0, group == "DI" & timepoint.within.group == 1)
# Pick the top OTUs only
pseq <- prune_taxa(names(sort(taxa_sums(pseq), TRUE)[1:5]), pseq)
p <- plot_composition(pseq, relative.abundance = TRUE, sort.by = "nationality", x.label = "nationality")
p <- p + guides(fill = guide_legend(ncol = 1))
p <- p + theme(legend.position = "bottom")
print(p)
```

![plot of chunk composition-example4](figure/composition-example4-1.png) 

### Coloured Barplots

The following example visualizes samples, colored by Phylum
percentages (in this example data the Phylum is approximated by 16S
sequence similarity, not exactly Phylum):


```r
pseq <- subset_samples(pseq0, group == "DI")
p <- plot_bar(pseq, x = "timepoint.within.group", fill = "Phylum", facet_grid = ~nationality)
```

```
## Error in .Method(..., na.last = na.last, decreasing = decreasing): argument 1 is not a vector
```

```r
print(p)
```

![plot of chunk composition-example5](figure/composition-example5-1.png) 

