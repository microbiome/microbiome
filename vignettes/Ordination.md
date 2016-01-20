## Ordination examples

Here some examples with HITChip data. See also [phyloseq ordination tutorial](joey711.github.io/phyloseq/plot_ordination-examples.html).

Load example data:


```r
library(microbiome)
library(ggplot2)
pseq <- download_microbiome("atlas1006")

# Convert signal to relative abundances
pseq.rel <- transform_sample_counts(pseq, function (x) {x/sum(x)})
```

```
## Error in eval(expr, envir, enclos): could not find function "transform_sample_counts"
```

```r
# Pick OTUs that are present with >1 percent relative abundance 
# in >10 percent of the samples
pseq2 <- filter_taxa(pseq.rel, function(x) sum(x > .01) > (0.1*nsamples(pseq.rel)), TRUE)
```

```
## Error in eval(expr, envir, enclos): could not find function "filter_taxa"
```


### Principal component analysis (PCA)

Project data on the first two PCA axes:


```r
set.seed(423542)
x <- t(otu_table(pseq2)@.Data)
```

```
## Error in t(otu_table(pseq2)@.Data): could not find function "otu_table"
```

```r
proj <- project.data(log10(x), type = "PCA") # MDS.classical etc.
```

```
## Error in if (nrow(amat) < ncol(amat)) {: argument is of length zero
```


Visualize and highlight:


```r
# Highlight gender
p <- densityplot(proj, col = sample_data(pseq2)$gender, legend = T)
```

```
## Error in as.vector(x, mode): cannot coerce type 'closure' to vector of type 'any'
```

```r
print(p)
```

![plot of chunk ordination4](figure/ordination4-1.png)

```r
# Highlight low/high Prevotella subjects
prevotella.abundance  <- log10(x[, "Prevotella melaninogenica et rel."]) 
```

```
## Error in x[, "Prevotella melaninogenica et rel."]: incorrect number of dimensions
```

```r
p <- densityplot(proj, col = prevotella.abundance, legend = T)
```

```
## Error in as.vector(x, mode): cannot coerce type 'closure' to vector of type 'any'
```

```r
print(p)
```

![plot of chunk ordination4](figure/ordination4-2.png)

Projection with sample names:


```r
ggplot(aes(x = Comp.1, y = Comp.2, label = rownames(proj)), data = proj) + geom_text(size = 2)
```

```
## Error: ggplot2 doesn't know how to deal with data of class function
```


### Unifrac with PCoA

Not implemented with HITChip but popular with sequencing data. See [phyloseq tutorial](http://joey711.github.io/phyloseq/plot_ordination-examples.html). 


### NMDS with Bray-Curtis distances


```r
pseq.ord <- ordinate(pseq2, "NMDS", "bray")
```

```
## Error in eval(expr, envir, enclos): could not find function "ordinate"
```

Ordinate the taxa in NMDS plot with Bray-Curtis distances


```r
p <- plot_ordination(pseq2, pseq.ord, type = "taxa", color = "Genus", title = "taxa")
```

```
## Error in eval(expr, envir, enclos): could not find function "plot_ordination"
```

```r
print(p)
```

![plot of chunk pca-ordination21](figure/pca-ordination21-1.png)

Grouping the plots by Phylum


```r
p + facet_wrap(~Phylum, 5)
```

```
## Error in layout_base(data, vars, drop = drop): At least one layer must contain all variables used for facetting
```

![plot of chunk pca-ordination22](figure/pca-ordination22-1.png)


### Multidimensional scaling (MDS)


```r
plot_ordination(pseq, ordinate(pseq, "MDS"), color = "DNA_extraction_method") + geom_point(size = 5)
```

```
## Error in eval(expr, envir, enclos): could not find function "plot_ordination"
```


### Canonical correspondence analysis (CCA)

With samples:


```r
p <- plot_ordination(pseq, ordinate(pseq, "CCA"), type = "samples", color = "gender")
```

```
## Error in eval(expr, envir, enclos): could not find function "plot_ordination"
```

```r
p + geom_point(size = 5)
```

```
## Error: geom_point requires the following missing aesthetics: y
```

![plot of chunk ordinate24a](figure/ordinate24a-1.png)

```r
#p <- p + geom_polygon(aes(fill = gender))
```

With taxa:


```r
p <- plot_ordination(pseq, ordinate(pseq, "CCA"), type = "taxa", color = "Phylum")
```

```
## Error in eval(expr, envir, enclos): could not find function "plot_ordination"
```

```r
p <- p + geom_point(size = 4)
print(p)
```

```
## Error: geom_point requires the following missing aesthetics: y
```

![plot of chunk ordinate24b](figure/ordinate24b-1.png)


### Split plot


```r
plot_ordination(pseq, ordinate(pseq, "CCA"), type = "split", shape = "gender", 
    color = "Phylum", label = "gender")
```

```
## Error in eval(expr, envir, enclos): could not find function "plot_ordination"
```


### Ordination biplot


```r
plot_ordination(pseq, ordinate(pseq, "CCA"), type = "biplot", color = "Phylum")
```

```
## Error in eval(expr, envir, enclos): could not find function "plot_ordination"
```






