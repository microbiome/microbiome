### Coloured Barplots

The following example visualizes samples, colored by Phylum
percentages. For further examples, see [phyloseq
package](http://joey711.github.io/phyloseq/plot_bar-examples.html).


```r
library(microbiome)
library(ggplot2)

# Microbiota profiling data. Read as: bacteria x samples matrix
# From https://peerj.com/articles/32/
x <- download_microbiome("peerj32")$physeq

# Visualize samples by Phyla (note: in HITChip data it is only approximately Phylum level)
plot_bar(x, x = "sample", fill = "Phylum")
```

```
## Error in .Method(..., na.last = na.last, decreasing = decreasing): argument 1 is not a vector
```


[Phyloseq](http://joey711.github.io/phyloseq/plot_bar-examples.html) example, filling by Phylum:


```r
p <- plot_bar(x, fill = "Phylum")
```

```
## Error in .Method(..., na.last = na.last, decreasing = decreasing): argument 1 is not a vector
```

```r
print(p)
```

![plot of chunk taxbar](figure/taxbar-1.png) 


Top OTU plot


```r
library(microbiome)
data.dietswap <- download_microbiome("dietswap")
TopNOTUs <- names(sort(taxa_sums(x), TRUE)[1:3])
tops <- prune_taxa(TopNOTUs, x)
plot_bar(tops, "group", fill = "gender", facet_grid = ~Phylum)
```

```
## Error in .Method(..., na.last = na.last, decreasing = decreasing): argument 1 is not a vector
```

```r
#plot_bar(ent10, "Genus", fill = "Genus", facet_grid = SeqTech ~ Enterotype)
```


