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
data(dietswap)
pseq <- dietswap
```


### Global indicators


A comprehensive list of global indicators of the ecosystem state can be obtained as follows. This includes various measures of richness, evenness, diversity, dominance, and rarity with default parameters. See the individual functions for more options regarding parameter tuning.


```r
global.inds <- global(pseq, index = "all")
head(kable(global.inds))
```

```
## [1] "|           | richness_0%| richness_20%| richness_50%| richness_80%| diversities_inverse_simpson| diversities_gini_simpson| diversities_shannon| diversities_fisher| diversities_observed| diversities_inverse_gini| diversities_coverage| evenness_camargo| evenness_pielou| evenness_simpson| evenness_evar| evenness_bulla| dominance_DBP| dominance_DMN| dominance_absolute| dominance_relative| dominance_simpson| dominance_core_abundance| rarity_log_modulo_skewness| rarity_low_abundance| rarity_rare_abundance|"
## [2] "|:----------|-----------:|------------:|------------:|------------:|---------------------------:|------------------------:|-------------------:|------------------:|--------------------:|------------------------:|--------------------:|----------------:|---------------:|----------------:|-------------:|--------------:|-------------:|-------------:|------------------:|------------------:|-----------------:|------------------------:|--------------------------:|--------------------:|---------------------:|"
## [3] "|Sample-1   |         112|          104|           61|           27|                    7.562984|                0.8677771|            2.942723|           12.16148|                  112|                 1.159830|                    4|        0.1378045|       0.6045614|        0.0581768|     0.0736465|      0.2925991|     0.3279166|     0.4296966|             175035|          0.3279166|         0.0012093|                0.9274756|                   2.059086|            0.0291825|             0.0160722|"
## [4] "|Sample-2   |         118|          110|           68|           39|                    8.105283|                0.8766237|            2.824184|           11.11824|                  118|                 1.131139|                    3|        0.1159349|       0.5802083|        0.0623483|     0.0722394|      0.2506029|     0.2428268|     0.4655585|             323085|          0.2428268|         0.0012093|                0.9328351|                   2.058747|            0.0302304|             0.0167326|"
## [5] "|Sample-3   |         113|          104|           71|           37|                    4.292701|                0.7670464|            2.409584|           10.80073|                  113|                 1.101253|                    2|        0.0919433|       0.4950318|        0.0330208|     0.0608332|      0.2233591|     0.4593873|     0.5602856|             837328|          0.4593873|         0.0012093|                0.9513098|                   2.056009|            0.0341229|             0.0125840|"
## [6] "|Sample-4   |         114|          106|           73|           30|                    7.937365|                0.8740136|            2.994672|           11.62450|                  114|                 1.167401|                    4|        0.1433967|       0.6152338|        0.0610567|     0.0692447|      0.2829995|     0.3229230|     0.3956421|             269963|          0.3229230|         0.0012093|                0.8617545|                   2.054981|            0.0349690|             0.0140909|"
```


### Alpha diversity

This returns a table with selected diversity indicators. See a separate page on [Beta diversity](Betadiversity.html).


```r
tab <- diversities(pseq, index = "all")
head(kable(tab))
```

```
## [1] "|           | inverse_simpson| gini_simpson|  shannon|   fisher| observed| inverse_gini| coverage|"
## [2] "|:----------|---------------:|------------:|--------:|--------:|--------:|------------:|--------:|"
## [3] "|Sample-1   |        7.562984|    0.8677771| 2.942723| 12.16148|      112|     1.159830|        4|"
## [4] "|Sample-2   |        8.105283|    0.8766237| 2.824184| 11.11824|      118|     1.131139|        3|"
## [5] "|Sample-3   |        4.292701|    0.7670464| 2.409584| 10.80073|      113|     1.101253|        2|"
## [6] "|Sample-4   |        7.937365|    0.8740136| 2.994672| 11.62450|      114|     1.167401|        4|"
```


### Richness

This returns observed richness with given detection threshold(s).


```r
tab <- richness(pseq)
head(kable(tab))
```

```
## [1] "|           |  0%| 20%| 50%| 80%|" "|:----------|---:|---:|---:|---:|"
## [3] "|Sample-1   | 112| 104|  61|  27|" "|Sample-2   | 118| 110|  68|  39|"
## [5] "|Sample-3   | 113| 104|  71|  37|" "|Sample-4   | 114| 106|  73|  30|"
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
|Sample-1 | 0.3279166| 0.4296966|   175035| 0.3279166| 0.0012093|      0.9274756|
|Sample-2 | 0.2428268| 0.4655585|   323085| 0.2428268| 0.0012093|      0.9328351|
|Sample-3 | 0.4593873| 0.5602856|   837328| 0.4593873| 0.0012093|      0.9513098|
|Sample-4 | 0.3229230| 0.3956421|   269963| 0.3229230| 0.0012093|      0.8617545|
|Sample-5 | 0.5448817| 0.6315785|   596658| 0.5448817| 0.0012093|      0.9533809|
|Sample-6 | 0.5692406| 0.6428424|   709407| 0.5692406| 0.0012093|      0.9255292|



### Rarity and low abundance

The rarity indices quantify the concentration of rare or low abundance taxa. Various rarity indices are available (see the function help for a list of options).


```r
ra <- rarity(pseq, index = "all")
kable(head(ra))
```



|         | log_modulo_skewness| low_abundance| rare_abundance|
|:--------|-------------------:|-------------:|--------------:|
|Sample-1 |            2.059086|     0.0291825|      0.0160722|
|Sample-2 |            2.058747|     0.0302304|      0.0167326|
|Sample-3 |            2.056009|     0.0341229|      0.0125840|
|Sample-4 |            2.054981|     0.0349690|      0.0140909|
|Sample-5 |            2.060330|     0.0443278|      0.0156928|
|Sample-6 |            2.060149|     0.0378701|      0.0193888|



### Coverage

The coverage index gives the number of groups needed to have a given proportion of the ecosystem occupied (by default 0.5 ie 50%).


```r
do <- coverage(pseq, threshold = 0.5)
```


### Core abundance

The core_abundance function refers to the relative proportion of the core species. Rare abundance provides the complement (1-x; see rare_abundance).


```r
co <- core_abundance(pseq, detection = .1/100, prevalence = 50/100)
```


### Gini index

Gini index is a common measure for inequality in economical income. The inverse gini index (1/x) can also be used as a community diversity measure.


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
|Sample-1 | 0.1378045| 0.6045614| 0.0581768| 0.0736465| 0.2925991|
|Sample-2 | 0.1159349| 0.5802083| 0.0623483| 0.0722394| 0.2506029|
|Sample-3 | 0.0919433| 0.4950318| 0.0330208| 0.0608332| 0.2233591|
|Sample-4 | 0.1433967| 0.6152338| 0.0610567| 0.0692447| 0.2829995|
|Sample-5 | 0.0790664| 0.4331198| 0.0244025| 0.0709762| 0.2066856|
|Sample-6 | 0.0811396| 0.4259505| 0.0227324| 0.0685637| 0.2084814|


