---
title: "Format Phyloseq"
author: "Sudarshan A. Shetty"
date: "2017-03-05"
output: 
  rmarkdown::html_vignette
---
<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - Format Phyloseq}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->


## Formatting the Phyloseq Object

Load [example data](Data.md):  
For this example we will use data from [Halfvarson J., et al. Nature Microbiology, 2017](http://www.nature.com/articles/nmicrobiol20174). It was downloaded from [Qitta](https://qiita.ucsd.edu/study/description/1629).  



```r
library(microbiome)
data(DynamicsIBD)
p0 <- DynamicsIBD
```


### Check the taxonomy 

We will check the taxonomy information stored in the phyloseq object.  


```r
library(knitr)
kable(head(tax_table(p0)))
```



|        |Rank1       |Rank2         |Rank3         |Rank4              |Rank5               |Rank6               |Rank7 |
|:-------|:-----------|:-------------|:-------------|:------------------|:-------------------|:-------------------|:-----|
|577110  |k__Bacteria |p__Firmicutes |c__Clostridia |o__Clostridiales   |f__                 |g__                 |s__   |
|181342  |k__Bacteria |p__Firmicutes |c__Clostridia |o__Clostridiales   |f__Clostridiaceae   |g__02d06            |s__   |
|581609  |k__Bacteria |p__Firmicutes |c__Clostridia |o__Clostridiales   |f__Ruminococcaceae  |g__                 |s__   |
|4341234 |k__Bacteria |p__Firmicutes |c__Clostridia |o__Clostridiales   |f__Peptococcaceae   |g__Desulfotomaculum |s__   |
|181348  |k__Bacteria |p__Firmicutes |c__Clostridia |o__Clostridiales   |f__Lachnospiraceae  |g__Coprococcus      |s__   |
|4467992 |k__Bacteria |p__Firmicutes |c__Bacilli    |o__Lactobacillales |f__Streptococcaceae |g__Streptococcus    |s__   |

It can be observed that the not all the OTUs are classified until the lowest taxonomic level (here, species level). This is especially the case with high throughput sequencing data sets. In doing OTU level testing for differential abundance, you may need information regading the specific otu number or taxonomy of the otu. This can help in easily tracing back the sequence and also make the plots with best taxonomic classification possible. Additionally, the names of taxonomic ranks are corrected using this function.


```r
p0.f <- format_phyloseq(p0)

#Check the taxonomy again with the formatted phyloseq object.
kable(head(tax_table(p0.f)))
```



|        |Domain   |Phylum     |Class      |Order           |Family                  |Genus                     |Species                     |
|:-------|:--------|:----------|:----------|:---------------|:-----------------------|:-------------------------|:---------------------------|
|577110  |Bacteria |Firmicutes |Clostridia |Clostridiales   |o__Clostridiales_577110 |o__Clostridiales_577110   |o__Clostridiales_577110     |
|181342  |Bacteria |Firmicutes |Clostridia |Clostridiales   |Clostridiaceae          |02d06                     |f__02d06_181342             |
|581609  |Bacteria |Firmicutes |Clostridia |Clostridiales   |Ruminococcaceae         |f__Ruminococcaceae_581609 |f__Ruminococcaceae_581609   |
|4341234 |Bacteria |Firmicutes |Clostridia |Clostridiales   |Peptococcaceae          |Desulfotomaculum          |f__Desulfotomaculum_4341234 |
|181348  |Bacteria |Firmicutes |Clostridia |Clostridiales   |Lachnospiraceae         |Coprococcus               |f__Coprococcus_181348       |
|4467992 |Bacteria |Firmicutes |Bacilli    |Lactobacillales |Streptococcaceae        |Streptococcus             |f__Streptococcus_4467992    |

Check [Preprocessing](Preprocessing.md) for additional utilitiy functions.
