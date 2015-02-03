---
title: "Project Template"
author: "Your Name"
date: "2015-02-03"
output:
  html_document:
    toc: true
    number_sections: true
    theme: united
    highlight: pygments
---
<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{Project Template}
  %\usepackage[utf8]{inputenc}
-->


Example template
===========


```r
library(devtools)
```

```
## 
## Attaching package: 'devtools'
## 
## The following object is masked from 'package:permute':
## 
##     check
```

```r
install_github("microbiome/microbiome")
```

```
## Downloading github repo microbiome/microbiome@master
## Installing microbiome
## '/usr/lib/R/bin/R' --vanilla CMD INSTALL  \
##   '/tmp/RtmpGM2MfC/devtools3bfbf407a3/microbiome-microbiome-d3b0a8e'  \
##   --library='/home/antagomir/R/x86_64-pc-linux-gnu-library/3.1'  \
##   --install-tests 
## 
## Reloading installed microbiome
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

```r
plot(c(1,2,3))
```

![plot of chunk install2](figure/install2-1.png) 

