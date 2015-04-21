# HITChip data analysis in phyloseq format

The [phyloseq](https://github.com/joey711/phyloseq) is an external high-quality R package that provides additional tools for microbiome data analysis. These examples show how to convert HITChip data into phyloseq format and perform some standard analyses. For more info on phyloseq tools, see [phyloseq demo page](http://joey711.github.io/phyloseq-demo/).


## Example data

Loading example data (L2 data and metadata; you can replace these with your own data). Make sure that the L2 datamatrix


```r
library(microbiome)
data(peerj32)
data <- peerj32$microbes # Samples x L2 groups; L2 data matrix; ABSOLUTE scale, not log10
meta <- peerj32$meta # Samples x features metadata
```

## Converting HITChip L2 data into [phyloseq](https://github.com/joey711/phyloseq) format


```r
library("phyloseq")
physeq <- hitchip2physeq(data, meta)
```

```
## Reading /home/antagomir/R/x86_64-pc-linux-gnu-library/3.1/microbiome/extdata/phylogeny.full.tab
```

## Barplots


```r
p <- plot_bar(physeq, fill = "Phylum")
print(p)
```

![plot of chunk taxbar](figure/taxbar-1.png) 


## Heatmaps

[Heatmaps](http://joey711.github.io/phyloseq/plot_heatmap-examples) and [Neatmaps](http://www.biomedcentral.com/1471-2105/11/45)


```r
plot_heatmap(physeq, taxa.label = "Phylum")
```

![plot of chunk heatmap](figure/heatmap-1.png) 

```r
plot_heatmap(physeq, "NMDS", "bray", "gender", "Phylum", low="#000033", high="#FF3300", na.value="white")
```

![plot of chunk heatmap](figure/heatmap-2.png) 

```r
#plot_heatmap(physeq, "NMDS", "bray", "gender", "Phylum", trans = log_trans(10))
#plot_heatmap(physeq, "NMDS", "bray", "gender", "Phylum", trans = identity_trans())
#plot_heatmap(physeq, "NMDS", "bray", "gender", "Phylum", trans = boxcox_trans(0.15))
```

## Richness


```r
library(ggplot2)
plot_richness(physeq, x = "gender", color = "group") + geom_boxplot()
```

![plot of chunk richness](figure/richness-1.png) 

## Top OTU plot


```r
TopNOTUs <- names(sort(taxa_sums(physeq), TRUE)[1:3])
tops <- prune_taxa(TopNOTUs, physeq)
plot_bar(tops, "group", fill = "gender", facet_grid = ~Genus)
```

![plot of chunk topotu](figure/topotu-1.png) 

## Ordination


```r
plot_ordination(physeq, ordinate(physeq, "MDS"), color = "group") + geom_point(size = 5)
```

![plot of chunk ordinate](figure/ordinate-1.png) 

```r
nmds <- ordinate(physeq, "NMDS", "bray")
```

```
## Square root transformation
## Wisconsin double standardization
## Run 0 stress 0.1747778 
## Run 1 stress 0.1754901 
## Run 2 stress 0.1754901 
## Run 3 stress 0.1747798 
## ... procrustes: rmse 0.0005217542  max resid 0.002997529 
## *** Solution reached
```

```r
require("ggplot2")
p2 <- plot_ordination(physeq, ordinate(physeq, "CCA"), type = "samples", color = "gender")
p2 + geom_point(size = 5) + geom_polygon(aes(fill = gender))
```

![plot of chunk ordinate](figure/ordinate-2.png) 

```r
plot_ordination(physeq, ordinate(physeq, "CCA"), type = "taxa", color = "Phylum") + 
    geom_point(size = 4)
```

![plot of chunk ordinate](figure/ordinate-3.png) 

```r
plot_ordination(physeq, ordinate(physeq, "CCA"), type = "split", color = "gender")
```

![plot of chunk ordinate](figure/ordinate-4.png) 

```r
plot_ordination(physeq, ordinate(physeq, "CCA"), type = "biplot", shape = "Phylum")
```

![plot of chunk ordinate](figure/ordinate-5.png) 

```r
plot_ordination(physeq, ordinate(physeq, "CCA"), type = "split", color = "gender", 
    shape = "Phylum", label = "gender")
```

![plot of chunk ordinate](figure/ordinate-6.png) 

## Filtering and pruning


```r
f1 <- filterfun_sample(topp(0.1))
wh1 <- genefilter_sample(physeq, f1, A = round(0.5 * nsamples(physeq)))
ex2 <- prune_taxa(wh1, physeq)
r <- transform_sample_counts(physeq, function(x) x/sum(x))
f <- filter_taxa(r, function(x) var(x) > 1e-05, TRUE)
```

## Networks

[Networks](http://joey711.github.io/phyloseq/plot_network-examples)


```r
plot_net(physeq, maxdist = 0.45, point_label = "group")
```

![plot of chunk networks](figure/networks-1.png) 

```r
ig <- make_network(physeq, max.dist = 0.45)
plot_network(ig, physeq, color = "gender", shape = "group", line_weight = 0.4, label = NULL)
```

![plot of chunk networks](figure/networks-2.png) 



For more, see [phyloseq demo page](http://joey711.github.io/phyloseq-demo/phyloseq-demo.html)


```r
#ntaxa(physeq)
#rank_names(physeq)
#sample_names(physeq)
#taxa_names(physeq)
#sample_variables(physeq)
#merge_phyloseq
#get_variable(physeq, sample_variables(physeq)[5])[1:10]
#get_taxa_unique(physeq, "Phylum")
#sample_sums(physeq)
#taxa_sums(physeq)
#get_taxa(physeq, sample_names(physeq)[5])
#get_sample(physeq, taxa_names(physeq)[5])
#plot_bar(ent10, "Genus", fill = "Genus", facet_grid = SeqTech ~ Enterotype)
#GP <- prune_taxa(taxa_sums(physeq) > 0, physeq)
#sample_data(GP)$human = ...
#eso <- rarefy_even_depth(physeq)
# UniFrac(physeq) # This requires tree that we do not have
#GP.chl = subset_taxa(GlobalPatterns, Phylum == "Chlamydiae")
```
