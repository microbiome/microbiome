## Intestinal microbiota diversity in 1006 western adults

Let us investigate an example data set from [Lahti et al. Nat. Comm. 5:4344, 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html). This contains microbiota profiling of 130 genus-like taxa across 1006 normal western adults from [Data Dryad](http://doi.org/10.5061/dryad.pk75d).


### Download HITChip Atlas data

[Load the HITChip Atlas microbiome profiling data in R](Data.md)


```r
# Download the required R packages and then the HITChip Atlas data set
library("rdryad")
library("microbiome")
```

```
## Loading required package: phyloseq
## Creating a generic function for 'nchar' from package 'base' in package 'S4Vectors'
## Loading required package: RPA
## Loading required package: parallel
## Loading required package: affy
## Loading required package: BiocGenerics
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
## 
## RPA Copyright (C) 2008-2014 Leo Lahti.
## This program comes with ABSOLUTELY NO WARRANTY.
## This is free software, and you are welcome to redistribute it under the FreeBSD open source license.
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
```

```r
pseq <- download_microbiome("atlas1006")
```

```
## Downloading data set from Lahti et al. Nat. Comm. 5:4344, 2014 from Data Dryad: http://doi.org/10.5061/dryad.pk75d
```


### Diversity 

### Estimating microbial diversity with different diversity measures


```r
library(phyloseq)
div <- estimate_diversity(pseq, measures = c("Observed", "Shannon", "Simpson"))
kable(head(div))
```



|         | Observed|  Shannon|   Simpson|
|:--------|--------:|--------:|---------:|
|Sample.1 |      130| 3.189726| 0.9230387|
|Sample.2 |      130| 3.396135| 0.9397719|
|Sample.3 |      130| 2.866104| 0.8850959|
|Sample.4 |      130| 3.058653| 0.9066459|
|Sample.5 |      130| 3.076850| 0.9184565|
|Sample.6 |      130| 2.945709| 0.8966565|


### Diversity vs. obesity


```r
p <- plot_diversity(pseq, x = "bmi_group", measures = c("Observed", "Shannon", "Simpson"), det.th = 250)
print(p)
```

![plot of chunk div-example2](figure/div-example2-1.png) 


### Diversity vs. age


```r
# Pick the subset of RBB-preprocessed samples from time point 0
pseq <- subset_samples(pseq, time == 0 & DNA_extraction_method == "r")

# Visualize
library(sorvi)
p <- sorvi::regression_plot(diversity~age, sample_data(pseq))
print(p)
```

![plot of chunk atlas-example3](figure/atlas-example3-1.png) 


## Further resources

For further examples, see [microbiome tutorial](https://github.com/microbiome/microbiome/blob/master/vignettes/vignette.md)
