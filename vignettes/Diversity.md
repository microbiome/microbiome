---
title: "Diversity"
bibliography: 
- bibliography.bib
- references.bib
output: 
  prettydoc::html_pretty:
    theme: cayman
    highlight: github
---
<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - diversity}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->


## Global Ecosystem State Variables 

Load example data:


```r
library(microbiome)
data(atlas1006)
pseq <- atlas1006
```


### Global indicators


A comprehensive list of global indicators of the ecosystem state can be obtained as follows. This includes various measures of richness, evenness, diversity, dominance, and rarity with default parameters. See the individual functions for more options regarding parameter tuning.


```r
global.inds <- global(pseq, index = "all")
head(kable(global.inds))
```

```
## [1] "|            | richness_0%| richness_20%| richness_50%| richness_80%| diversities_inverse_simpson| diversities_gini_simpson| diversities_shannon| diversities_fisher| diversities_observed| diversities_inverse_gini| diversities_coverage| evenness_camargo| evenness_pielou| evenness_simpson| evenness_evar| evenness_bulla| dominance_DBP| dominance_DMN| dominance_absolute| dominance_relative| dominance_simpson| dominance_core_abundance| rarity_log_modulo_skewness| rarity_low_abundance| rarity_rare_abundance|"
## [2] "|:-----------|-----------:|------------:|------------:|------------:|---------------------------:|------------------------:|-------------------:|------------------:|--------------------:|------------------------:|--------------------:|----------------:|---------------:|----------------:|-------------:|--------------:|-------------:|-------------:|------------------:|------------------:|-----------------:|------------------------:|--------------------------:|--------------------:|---------------------:|"
## [3] "|Sample-1    |         109|           99|           60|           25|                   12.993537|                0.9230387|            3.189726|           12.29785|                  109|                 1.178316|                    5|        0.1513312|       0.6553063|        0.0999503|     0.0665828|      0.2902021|     0.1758679|     0.3377921|              84316|          0.1758679|         0.0001472|                0.9597792|                   2.060324|            0.0246043|             0.0164759|"
## [4] "|Sample-2    |         110|           99|           63|           34|                   16.603545|                0.9397719|            3.396135|           11.93702|                  110|                 1.221472|                    7|        0.1813160|       0.6977115|        0.1277196|     0.0599369|      0.3293462|     0.1716273|     0.2820096|             109940|          0.1716273|         0.0001472|                0.9015118|                   2.058091|            0.0199587|             0.0141982|"
## [5] "|Sample-3    |         109|           99|           59|           13|                    8.702908|                0.8850959|            2.866104|           12.38015|                  109|                 1.135699|                    4|        0.1194850|       0.5888204|        0.0669454|     0.0738372|      0.2381299|     0.2793253|     0.3985027|             125664|          0.2793253|         0.0001472|                0.9391221|                   2.056895|            0.0393057|             0.0192827|"
## [6] "|Sample-4    |         111|          100|           67|           25|                   10.711903|                0.9066459|            3.058653|           11.85667|                  111|                 1.162582|                    4|        0.1398459|       0.6283784|        0.0823993|     0.0644078|      0.2896032|     0.1957585|     0.3813535|             134094|          0.1957585|         0.0001472|                0.9509151|                   2.059333|            0.0249986|             0.0141300|"
```


### Alpha diversity

This returns a table with selected diversity indicators. See a separate page on [Beta diversity](Betadiversity.html).


```r
tab <- diversities(pseq, index = "all")
head(kable(tab))
```

```
## [1] "|            | inverse_simpson| gini_simpson|  shannon|   fisher| observed| inverse_gini| coverage|"
## [2] "|:-----------|---------------:|------------:|--------:|--------:|--------:|------------:|--------:|"
## [3] "|Sample-1    |       12.993537|    0.9230387| 3.189726| 12.29785|      109|     1.178316|        5|"
## [4] "|Sample-2    |       16.603545|    0.9397719| 3.396135| 11.93702|      110|     1.221472|        7|"
## [5] "|Sample-3    |        8.702908|    0.8850959| 2.866104| 12.38015|      109|     1.135699|        4|"
## [6] "|Sample-4    |       10.711903|    0.9066459| 3.058653| 11.85667|      111|     1.162582|        4|"
```


### Richness

This returns observed richness with given detection threshold(s).


```r
tab <- richness(pseq)
head(kable(tab))
```

```
## [1] "|            |  0%| 20%| 50%| 80%|"
## [2] "|:-----------|---:|---:|---:|---:|"
## [3] "|Sample-1    | 109|  99|  60|  25|"
## [4] "|Sample-2    | 110|  99|  63|  34|"
## [5] "|Sample-3    | 109|  99|  59|  13|"
## [6] "|Sample-4    | 111| 100|  67|  25|"
```


### Dominance 

The dominance index refers to the abundance of the most abundant species. Various dominance indices are available (see the function help for a list of options).


```r
# Absolute abundances for the single most abundant taxa in each sample
do <- dominance(pseq, index = "all")
kable(head(do))
```



|         |       DBP|       DMN| absolute|  relative|   simpson| core_abundance|
|:--------|---------:|---------:|--------:|---------:|---------:|--------------:|
|Sample-1 | 0.1758679| 0.3377921|    84316| 0.1758679| 0.0001472|      0.9597792|
|Sample-2 | 0.1716273| 0.2820096|   109940| 0.1716273| 0.0001472|      0.9015118|
|Sample-3 | 0.2793253| 0.3985027|   125664| 0.2793253| 0.0001472|      0.9391221|
|Sample-4 | 0.1957585| 0.3813535|   134094| 0.1957585| 0.0001472|      0.9509151|
|Sample-5 | 0.1685621| 0.3332058|   127719| 0.1685621| 0.0001472|      0.9438707|
|Sample-6 | 0.2271913| 0.3799754|   113490| 0.2271913| 0.0001472|      0.9557649|



### Rarity and low abundance

The rarity indices quantify the concentration of rare or low abundance taxa. Various rarity indices are available (see the function help for a list of options).


```r
ra <- rarity(pseq, index = "all")
kable(head(ra))
```



|         | log_modulo_skewness| low_abundance| rare_abundance|
|:--------|-------------------:|-------------:|--------------:|
|Sample-1 |            2.060324|     0.0246043|      0.0164759|
|Sample-2 |            2.058091|     0.0199587|      0.0141982|
|Sample-3 |            2.056895|     0.0393057|      0.0192827|
|Sample-4 |            2.059333|     0.0249986|      0.0141300|
|Sample-5 |            2.059235|     0.0237021|      0.0166293|
|Sample-6 |            2.060673|     0.0384898|      0.0151221|



### Coverage

The coverage index gives the number of groups needed to have a given proportion of the ecosystem occupied (by default 0.5 ie 50%).


```r
do <- coverage(pseq, threshold = 0.5)
```


### Core abundance

The core_abundance function refers to the relative proportion of the core species. Rare abundance provides the complement.


```r
ra <- rare_abundance(pseq, detection = .1/100, prevalence = 50/100)
co <- core_abundance(pseq, detection = .1/100, prevalence = 50/100)
```



### Gini index

Gini index is a common measure for inequality in economical income, but can also be used as a community diversity measure.


```r
gi <- inequality(pseq)
```


### Evenness

Various evenness measures are also available.


```r
ev <- evenness(pseq, "all")
kable(head(ev))
```



|         |   camargo|    pielou|   simpson|      evar|     bulla|
|:--------|---------:|---------:|---------:|---------:|---------:|
|Sample-1 | 0.1513312| 0.6553063| 0.0999503| 0.0665828| 0.2902021|
|Sample-2 | 0.1813160| 0.6977115| 0.1277196| 0.0599369| 0.3293462|
|Sample-3 | 0.1194850| 0.5888204| 0.0669454| 0.0738372| 0.2381299|
|Sample-4 | 0.1398459| 0.6283784| 0.0823993| 0.0644078| 0.2896032|
|Sample-5 | 0.1332354| 0.6321168| 0.0943338| 0.0611757| 0.2662343|
|Sample-6 | 0.1270281| 0.6051747| 0.0744343| 0.0701696| 0.2649482|




### Visualization

Show indicators:


```r
library(ggplot2)
theme_set(theme_bw(20)) # Set bw color scheme
p <- ggplot(global.inds, aes(x = diversities_shannon)) + geom_histogram()
print(p)
```

![plot of chunk div-example2](figure/div-example2-1.png)

