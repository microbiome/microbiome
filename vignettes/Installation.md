### Installing R

If you do not already have R/RStudio installed, do the
following. 

  1. Install [R](http://www.r-project.org/) 
  1. Consider installing [RStudio](http://rstudio.org); GUI for R
  1. If you use Windows, install also [RTools](http://cran.r-project.org/bin/windows/Rtools/) (version corresponding to your R version)


### Install dependencies

Open R and install dependencies from the Tools panel, or run the
following commands:


```r
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("ade4", "fastcluster", "df2json", "rjson", "gplots", "devtools", 
  	   "ggplot2","MASS","minet","mixOmics", "plyr","qvalue","reshape2",
	   "RPA","svDialogs","vegan","WGCNA", "rpart"))
```

If some of these installations fail, ensure from the RStudio tools
panel that you have access to CRAN and Bioconductor repositories. If
you cannot install some packages, some functionality in microbiome may
not work.


### Install/update the microbiome package

To install microbiome package and recommended dependencies, run in R:


```r
library(devtools) # Load the devtools package
install_github("microbiome/microbiome") # Install the package
install_github("ropensci/rdryad") # Install proposed packages
```

### Loading the package

Once the package has been installed, you can load it in R with:


```r
library(microbiome)  
```

```
## Loading required package: e1071
## Loading required package: vegan
## Loading required package: permute
## Loading required package: lattice
## This is vegan 2.2-1
## Loading required package: reshape
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
## Loading required package: GenomeInfoDb
## Loading required package: S4Vectors
## 
## Attaching package: 'S4Vectors'
## 
## The following object is masked from 'package:reshape':
## 
##     rename
## 
## Loading required package: IRanges
## 
## Attaching package: 'IRanges'
## 
## The following object is masked from 'package:reshape':
## 
##     expand
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
```



## Install HITChipDB package

The HITChipDB package contains additional routines to fetch and
preprocess HITChip (or MIT/PITChip) data from the MySQL database. Note
that this package is not needed by most users and the data is
protected by password/IP combinations. Ask details from
admins. Install the package in R with:


```r
library(devtools) # Load the devtools package
install_github("microbiome/HITChipDB") # Install the package
# Also install RMySQL, multicore and tcltk !
source("http://www.bioconductor.org/biocLite.R")
biocLite("RMySQL") # multicore, tcltk?
# Test installation by loading the microbiome package in R
library("HITChipDB")
```

