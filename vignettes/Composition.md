## Microbiota composition


### Barplot visualizations

Also see [phyloseq barplot examples](http://joey711.github.io/phyloseq/plot_bar-examples.html) and [HITChip Barplots](Barplots.md)

Loading example data:


```r
# Example data
library(microbiome)
pseq0 <- download_microbiome("dietswap")
```

```
## Downloading data set from O'Keefe et al. Nat. Comm. 6:6342, 2015 from Data Dryad: http://datadryad.org/resource/doi:10.5061/dryad.1mn1n
```

```
## Error in curl::curl_fetch_memory(url, handle = handle): Timeout was reached
```

```r
# Pick sample subset
pseq <- subset_samples(pseq0, group == "DI" & nationality == "AFR")
```

```
## Error in sample_data(physeq): error in evaluating the argument 'object' in selecting a method for function 'sample_data': Error: object 'pseq0' not found
```

Show OTU absolute abundance in each sample. Plot absolute taxon
abundances (Note: on HITChip data the Phylum level is only
approximate):


```r
plot_composition(pseq, taxonomic.level = "Phylum")
```

```
## Error in match.fun(FUN): argument "FUN" is missing, with no default
```

Composition with relative abundances:


```r
p <- plot_composition(pseq, taxonomic.level = "Phylum", relative.abundance = TRUE)
```

```
## Error in match.fun(FUN): argument "FUN" is missing, with no default
```

```r
p <- p + guides(fill = guide_legend(nrow = 7, byrow = TRUE))
p <- p + theme(legend.position = "bottom")
print(p)
```

```
## Error in eval(expr, envir, enclos): object 'Correlation' not found
```

![plot of chunk composition-example3](figure/composition-example3-1.png)

Arrange by sample variable and use custom X axis labels. Africans have more Prevotella as expected:


```r
# Subset taxa and samples
pseq <- subset_samples(pseq0, group == "DI" & timepoint.within.group == 1)
```

```
## Error in sample_data(physeq): error in evaluating the argument 'object' in selecting a method for function 'sample_data': Error: object 'pseq0' not found
```

```r
# Pick the top OTUs only
pseq <- prune_taxa(names(sort(taxa_sums(pseq), TRUE)[1:5]), pseq)
p <- plot_composition(pseq, relative.abundance = TRUE, sort.by = "nationality", x.label = "nationality")
```

```
## Error in order(meta[[sort.by]]): argument 1 is not a vector
```

```r
p <- p + guides(fill = guide_legend(ncol = 1))
p <- p + theme(legend.position = "bottom")
print(p)
```

```
## Error in eval(expr, envir, enclos): object 'Correlation' not found
```

![plot of chunk composition-example4](figure/composition-example4-1.png)

