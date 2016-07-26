---
title: "HITChip Atlas Overview"
author: "Leo Lahti, Jarkko Salojarvi, Anne Salonen, Willem M de Vos"
date: "2016-07-26"
bibliography: 
- bibliography.bib
- references.bib
output: 
  md_document:
    variant: markdown_github
---
<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial}
  %\usepackage[utf8]{inputenc}
-->




HITChip Atlas Overview
===========

We use the [phyloseq](http://joey711.github.io/phyloseq/import-data)
class, a standard representation format in R for taxonomic
profiling. This package provides extensions and convenient wrappers
for many standard tasks encountered in microbiome studies. 


## Intestinal microbiota diversity 

This data set has microbiota profiling of 130 genus-like taxa across over 10,000 samples.

Download the data in R:


```r
# Download the required R packages and then the HITChip Atlas data set
library("microbiome")
#data("atlas1006")
#pseq <- atlas1006
# Load the data (not yet public)
atlas = readRDS("/home/lei/proj/hitchip-atlas-annotation-scripts/Atlas.RData") # atlas
```


### Data overview

 * ``10763`` samples
 * ``5007`` unique subjects


```r
p = ggplot(sample_data(atlas), aes(x = sample_type)) + geom_bar() + coord_flip() + ggtitle("Sample type")
print(p)

p = ggplot(sample_data(atlas), aes(x = DNA_extraction_method)) + geom_bar() + coord_flip() + ggtitle("DNA extraction method")
print(p)

p = ggplot(sample_data(atlas), aes(x = gender)) + geom_bar() + coord_flip() + ggtitle("Sex")
print(p)

p = ggplot(sample_data(atlas), aes(x = health_status)) + geom_bar() + coord_flip() + ggtitle("Health status")
print(p)

p = ggplot(sample_data(atlas), aes(x = health_info)) + geom_bar() + coord_flip() + ggtitle("Health info")
print(p)

p = ggplot(sample_data(atlas), aes(x = antibio)) + geom_bar() + coord_flip() + ggtitle("Antibiotics use")
print(p)

p = ggplot(sample_data(atlas), aes(x = probiotics)) + geom_bar() + coord_flip() + ggtitle("Probiotics")
print(p)

p = ggplot(sample_data(atlas), aes(x = medication)) + geom_bar() + coord_flip() + ggtitle("Medication")
print(p)

p = ggplot(sample_data(atlas), aes(x = nationality)) + geom_bar() + coord_flip() + ggtitle("Nationality")
print(p)

p = ggplot(sample_data(atlas), aes(x = age)) + geom_histogram() + ggtitle("Age") + xlab("Age")
print(p)

p = ggplot(sample_data(atlas), aes(x = bmi_group)) + geom_bar() + ggtitle("BMI group") + xlab("BMI")
print(p)
```

<img src="figure/hatlas-sampletype-1.png" title="plot of chunk hatlas-sampletype" alt="plot of chunk hatlas-sampletype" width="200px" /><img src="figure/hatlas-sampletype-2.png" title="plot of chunk hatlas-sampletype" alt="plot of chunk hatlas-sampletype" width="200px" /><img src="figure/hatlas-sampletype-3.png" title="plot of chunk hatlas-sampletype" alt="plot of chunk hatlas-sampletype" width="200px" /><img src="figure/hatlas-sampletype-4.png" title="plot of chunk hatlas-sampletype" alt="plot of chunk hatlas-sampletype" width="200px" /><img src="figure/hatlas-sampletype-5.png" title="plot of chunk hatlas-sampletype" alt="plot of chunk hatlas-sampletype" width="200px" /><img src="figure/hatlas-sampletype-6.png" title="plot of chunk hatlas-sampletype" alt="plot of chunk hatlas-sampletype" width="200px" /><img src="figure/hatlas-sampletype-7.png" title="plot of chunk hatlas-sampletype" alt="plot of chunk hatlas-sampletype" width="200px" /><img src="figure/hatlas-sampletype-8.png" title="plot of chunk hatlas-sampletype" alt="plot of chunk hatlas-sampletype" width="200px" /><img src="figure/hatlas-sampletype-9.png" title="plot of chunk hatlas-sampletype" alt="plot of chunk hatlas-sampletype" width="200px" /><img src="figure/hatlas-sampletype-10.png" title="plot of chunk hatlas-sampletype" alt="plot of chunk hatlas-sampletype" width="200px" /><img src="figure/hatlas-sampletype-11.png" title="plot of chunk hatlas-sampletype" alt="plot of chunk hatlas-sampletype" width="200px" />


### Diversity vs. age


```r
# Pick the subset of RBB-preprocessed samples from time point 0
pseq <- subset_samples(atlas, time == 0 & DNA_extraction_method == "rbb")

# Visualize
p <- plot_regression(diversity~age, sample_data(pseq))
print(p)
```

![plot of chunk hatlas-example3](figure/hatlas-example3-1.png)



### Licensing and Citations

This material can be freely used, modified and distributed under the
[Two-clause FreeBSD
license](http://en.wikipedia.org/wiki/BSD\_licenses). Kindly cite as
'Leo Lahti, Jarkko Salojarvi, Anne Salonen and Willem M de
Vos. HITChip Atlas. URL: http://microbiome.github.com'.


### References





