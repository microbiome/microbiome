### Abundance histograms

Different sample sets have different population distributions in
microbial abundance. It is also important to consider whether to use
absolute or logarithmic abundances!

Start by loading the [HITChip Atlas data set](Data.Rmd)



```r
# Load microbiome package
library(microbiome)  

# Load tools
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
## 
## The following object is masked from 'package:graph':
## 
##     union
## 
## The following object is masked from 'package:AnnotationDbi':
## 
##     select
## 
## The following objects are masked from 'package:BiocGenerics':
## 
##     intersect, setdiff, union
## 
## The following object is masked from 'package:reshape':
## 
##     rename
## 
## The following object is masked from 'package:stats':
## 
##     filter
## 
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
# 1. List all samples (all time points and DNA extraction methods)
all.samples <- meta$SampleID

# 2. List samples at time point 0 that have specific DNA extraction method 
rbb.samples <- filter(meta, Time == "0" & DNA_extraction_method == "r")$SampleID

# Visualize
#tax <- "Prevotella.melaninogenica.et.rel."
tax <- "Bifidobacterium"
d <- data[all.samples, tax]
par(mfrow = c(1, 2))
plot(density(d), main = paste(tax, "(All samples)"), xlab = "Abundance (Absolute HITChip signal)")
plot(density(log10(d)), main = paste(tax, "(All samples)"), xlab = "Abundance (Log10 HITChip signal)")
```

![plot of chunk hist](figure/hist-1.png) 

```r
d <- data[rbb.samples, tax]
par(mfrow = c(1, 2))
plot(density(d), main = paste(tax, "(RBB samples)"), xlab = "Abundance (Absolute HITChip signal)")
plot(density(log10(d)), main = paste(tax, "(RBB samples)"), xlab = "Abundance (Log10 HITChip signal)")
```

![plot of chunk hist](figure/hist-2.png) 
