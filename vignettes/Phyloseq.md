<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{Project Template}
  %\usepackage[utf8]{inputenc}
-->
Converting HITChip data in phyloseq format
------------------------------------------

The [phyloseq](https://github.com/joey711/phyloseq) is an external high-quality R package that provides additional tools for microbiome data analysis. These examples show how to convert HITChip data into phyloseq format and perform some standard analyses. For more info on phyloseq tools, see [phyloseq demo page](http://joey711.github.io/phyloseq-demo/).

Loading example data (L2 data and metadata; you can replace these with your own data). Make sure that the L2 datamatrix

``` r
library(microbiome)
```

    ## Loading required package: ade4
    ## Loading required package: e1071
    ## Loading required package: phyloseq
    ## Loading required package: reshape2
    ## Loading required package: vegan
    ## Loading required package: permute
    ## Loading required package: lattice
    ## This is vegan 2.2-1
    ## 
    ## Attaching package: 'vegan'
    ## 
    ## The following object is masked from 'package:ade4':
    ## 
    ##     cca
    ## 
    ## Loading required package: DBI
    ## Loading required package: AnnotationDbi
    ## Loading required package: stats4
    ## Loading required package: BiocGenerics
    ## Loading required package: parallel
    ## 
    ## Attaching package: 'BiocGenerics'
    ## 
    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB
    ## 
    ## The following object is masked from 'package:stats':
    ## 
    ##     xtabs
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, as.vector, cbind,
    ##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
    ##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unlist, unsplit
    ## 
    ## Loading required package: Biobase
    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.
    ## 
    ## 
    ## Attaching package: 'Biobase'
    ## 
    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     sampleNames
    ## 
    ## Loading required package: GenomeInfoDb
    ## Loading required package: S4Vectors
    ## Loading required package: IRanges
    ## 
    ## Attaching package: 'IRanges'
    ## 
    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     distance
    ## 
    ## The following object is masked from 'package:ade4':
    ## 
    ##     score
    ## 
    ## 
    ## Attaching package: 'AnnotationDbi'
    ## 
    ## The following object is masked from 'package:GenomeInfoDb':
    ## 
    ##     species
    ## 
    ## 
    ## 
    ## microbiome R package (microbiome.github.com)
    ##           
    ## 
    ## 
    ##  Copyright (C) 2011-2015
    ##           Leo Lahti and Jarkko Salojarvi 
    ## 
    ##         
    ##           <microbiome-admin@googlegroups.com>
    ## 
    ## 
    ## Attaching package: 'microbiome'
    ## 
    ## The following object is masked from 'package:lattice':
    ## 
    ##     densityplot
    ## 
    ## The following object is masked from 'package:e1071':
    ## 
    ##     impute

``` r
data(peerj32)
data <- peerj32$microbes # Samples x L2 groups; L2 data matrix; ABSOLUTE scale, not log10
meta <- peerj32$meta # Samples x features metadata
```

Convert the HITChip L2 data into [phyloseq](https://github.com/joey711/phyloseq) format

``` r
library("phyloseq")
physeq <- hitchip2physeq(data, meta)
```

    ## Reading /home/antagomir/R/x86_64-pc-linux-gnu-library/3.1/microbiome/extdata/phylogeny.full.tab

Some phyloseq examples

``` r
plot_bar(physeq, fill = "Phylum")
```

![](Phyloseq_files/figure-markdown_github/taxbar-1.png)

[Heatmaps](http://joey711.github.io/phyloseq/plot_heatmap-examples) and [Neatmaps](http://www.biomedcentral.com/1471-2105/11/45)

``` r
plot_heatmap(physeq, taxa.label = "Phylum")
```

![](Phyloseq_files/figure-markdown_github/heatmap-1.png)

``` r
plot_heatmap(physeq, "NMDS", "bray", "gender", "Phylum", low="#000033", high="#FF3300", na.value="white")
```

![](Phyloseq_files/figure-markdown_github/heatmap-2.png)

``` r
#plot_heatmap(physeq, "NMDS", "bray", "gender", "Phylum", trans = log_trans(10))
#plot_heatmap(physeq, "NMDS", "bray", "gender", "Phylum", trans = identity_trans())
#plot_heatmap(physeq, "NMDS", "bray", "gender", "Phylum", trans = boxcox_trans(0.15))
```

``` r
library(ggplot2)
plot_richness(physeq, x = "gender", color = "group") + geom_boxplot()
```

![](Phyloseq_files/figure-markdown_github/richness-1.png)

Top OTU plot:

``` r
TopNOTUs <- names(sort(taxa_sums(physeq), TRUE)[1:3])
tops <- prune_taxa(TopNOTUs, physeq)
plot_bar(tops, "group", fill = "gender", facet_grid = ~Genus)
```

![](Phyloseq_files/figure-markdown_github/topotu-1.png)

Ordination

``` r
plot_ordination(physeq, ordinate(physeq, "MDS"), color = "group") + geom_point(size = 5)
```

![](Phyloseq_files/figure-markdown_github/ordinate-1.png)

``` r
nmds <- ordinate(physeq, "NMDS", "bray")
```

    ## Square root transformation
    ## Wisconsin double standardization
    ## Run 0 stress 0.1747778 
    ## Run 1 stress 0.17549 
    ## Run 2 stress 0.1747778 
    ## ... New best solution
    ## ... procrustes: rmse 7.119887e-05  max resid 0.0004098436 
    ## *** Solution reached

``` r
require("ggplot2")
p2 <- plot_ordination(physeq, ordinate(physeq, "CCA"), type = "samples", color = "gender")
p2 + geom_point(size = 5) + geom_polygon(aes(fill = gender))
```

![](Phyloseq_files/figure-markdown_github/ordinate-2.png)

``` r
plot_ordination(physeq, ordinate(physeq, "CCA"), type = "taxa", color = "Phylum") + 
    geom_point(size = 4)
```

![](Phyloseq_files/figure-markdown_github/ordinate-3.png)

``` r
plot_ordination(physeq, ordinate(physeq, "CCA"), type = "split", color = "gender")
```

![](Phyloseq_files/figure-markdown_github/ordinate-4.png)

``` r
plot_ordination(physeq, ordinate(physeq, "CCA"), type = "biplot", shape = "Phylum")
```

    ## Warning: The shape palette can deal with a maximum of 6 discrete values
    ## because more than 6 becomes difficult to discriminate; you have
    ## 23. Consider specifying shapes manually. if you must have them.

    ## Warning: The shape palette can deal with a maximum of 6 discrete values
    ## because more than 6 becomes difficult to discriminate; you have
    ## 23. Consider specifying shapes manually. if you must have them.

![](Phyloseq_files/figure-markdown_github/ordinate-5.png)

``` r
plot_ordination(physeq, ordinate(physeq, "CCA"), type = "split", color = "gender", 
    shape = "Phylum", label = "gender")
```

    ## Warning: The shape palette can deal with a maximum of 6 discrete values
    ## because more than 6 becomes difficult to discriminate; you have
    ## 23. Consider specifying shapes manually. if you must have them.

    ## Warning: The shape palette can deal with a maximum of 6 discrete values
    ## because more than 6 becomes difficult to discriminate; you have
    ## 23. Consider specifying shapes manually. if you must have them.

    ## Warning: The shape palette can deal with a maximum of 6 discrete values
    ## because more than 6 becomes difficult to discriminate; you have
    ## 23. Consider specifying shapes manually. if you must have them.

![](Phyloseq_files/figure-markdown_github/ordinate-6.png)

Filtering and pruning

``` r
f1 <- filterfun_sample(topp(0.1))
wh1 <- genefilter_sample(physeq, f1, A = round(0.5 * nsamples(physeq)))
ex2 <- prune_taxa(wh1, physeq)
r <- transform_sample_counts(physeq, function(x) x/sum(x))
f <- filter_taxa(r, function(x) var(x) > 1e-05, TRUE)
```

[Networks](http://joey711.github.io/phyloseq/plot_network-examples)

``` r
plot_net(physeq, maxdist = 0.45, point_label = "group")
```

![](Phyloseq_files/figure-markdown_github/networks-1.png)

``` r
ig <- make_network(physeq, max.dist = 0.45)
plot_network(ig, physeq, color = "gender", shape = "group", line_weight = 0.4, label = NULL)
```

    ## Warning: attributes are not identical across measure variables; they will
    ## be dropped

![](Phyloseq_files/figure-markdown_github/networks-2.png)

For more, see [phyloseq demo page](http://joey711.github.io/phyloseq-demo/phyloseq-demo.html)

``` r
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
