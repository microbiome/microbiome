---
title: "Introduction to the microbiome R package"
author: "Leo Lahti, Sudarshan Shetty, et al."
bibliography: 
- bibliography.bib
date: "2017-06-03"
output:
  rmarkdown::pdf_document:
    toc: true
  rmarkdown::md_document:
    toc: true    
  rmarkdown::html_vignette:
    fig_caption: yes
vignette: >
  %\VignetteIndexEntry{microbiome R package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




## Introduction

The [microbiome R package](http://microbiome.github.io/microbiome) facilitates exploration and analysis of microbiome profiling data, in particular 16S taxonomic profiling. Here we provide a brief overview to the package functionality using example data sets from published studies [@lahti14natcomm, @Lahti13provasI, @OKeefe15]. A comprehensive tutorial is available [on-line](http://microbiome.github.io/microbiome). 

The microbiome R package has been utilized in recent publications [@salonen14, @Faust15, @Shetty2017], and is released under the [Two-clause FreeBSD license](http://en.wikipedia.org/wiki/BSD\_licenses). The package extends the independent [phyloseq](http://joey711.github.io/phyloseq/import-data) class structure and [tools](http://joey711.github.io/phyloseq/).

Input from the user community is very welcome via the [issue tracker](https://github.com/microbiome/microbiome/issues) or [pull requests](http://microbiome.github.io/microbiome/Contributing.html). See the [Github site](https://github.com/microbiome/microbiome) for source code and other details. 

Kindly cite the work as follows: "Leo Lahti [et al.](https://github.com/microbiome/microbiome/graphs/contributors) (2017). Tools for microbiome analysis in R. Microbiome package version 0.99.1. URL: [http://microbiome.github.com/microbiome](http://microbiome.github.com/microbiome). 


## Installation

To install microbiome package (the latest development version) in R, use


```r
library(devtools) # Load the devtools package
install_github("microbiome/microbiome") # Install the package
```

Then load the package in R


```r
library(microbiome)  
```


## Data

A [phyloseq](http://joey711.github.io/phyloseq) data object typically contains an OTU table (taxa abundances), sample metadata, taxonomy table (mapping between OTUs and higher-level taxonomic classifications), and a phylogenetic tree (relations between the taxa). 


### Example data

HITChip Atlas data set [Lahti et al. Nat. Comm. 5:4344, 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html)) contains 130 genus-level taxonomic groups across 1006 western adults. Get example data with


```r
data(atlas1006) 
```

A two-week diet swap study between western (USA) and traditional (rural Africa) diets, reported in [O'Keefe et al. Nat. Comm. 6:6342, 2015](http://dx.doi.org/10.1038/ncomms7342)


```r
data(dietswap) 
```

A gut microbiome profiling data set, which includes parallel profiling of intestinal microbiota versus blood metabolites from [Lahti et al. PeerJ 1:e32, 2013](https://peerj.com/articles/32/) to characterize associations between human intestinal microbiota and blood serum lipids. 


```r
data(peerj32) # Data from https://peerj.com/articles/32/
```


### Data import

The [tutorial](http://microbiome.github.io/microbiome/Data.html) describes how to import data from standard formats (Mother, BIOM, CSV).


### Data manipulation

A phyloseq object can be subsetted, filtered, aggregated, transformed, and otherwise manipulated with the [phyloseq](http://joey711.github.io/phyloseq/) and [microbiome](http://microbiome.github.io/microbiome/Preprocessing.html) tools. The microbiome package provides a wrapper for many standard transformations such  as Z, centered log-ratio, hellinger, log10, and other transformations. 

To convert absolute counts to compositional (relative) abundances, for instance, use


```r
pseq.comp <- transform(dietswap, "compositional")
```



## Diversity and other ecosystem indices 

Standard ecosystem state variables include richness, evenness, alpha diversity, dominance, and rarity. The function `global` calls these indicators with default parameters. For further options and [beta diversity](http://microbiome.github.io/microbiome/Betadiversity.html) quantification, see the [tutorial](http://microbiome.github.io/microbiome/Diversity.html). 



```r
g <- global(dietswap, index = c("shannon", "gini", "DMN"))
print(kable(head(g, 10), format = "html"))
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> diversities_shannon </th>
   <th style="text-align:right;"> dominance_gini </th>
   <th style="text-align:right;"> dominance_dmn </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Sample-1 </td>
   <td style="text-align:right;"> 2.942723 </td>
   <td style="text-align:right;"> 0.8621955 </td>
   <td style="text-align:right;"> 0.3279166 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sample-2 </td>
   <td style="text-align:right;"> 2.824184 </td>
   <td style="text-align:right;"> 0.8840651 </td>
   <td style="text-align:right;"> 0.2428268 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sample-3 </td>
   <td style="text-align:right;"> 2.409584 </td>
   <td style="text-align:right;"> 0.9080567 </td>
   <td style="text-align:right;"> 0.4593873 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sample-4 </td>
   <td style="text-align:right;"> 2.994672 </td>
   <td style="text-align:right;"> 0.8566033 </td>
   <td style="text-align:right;"> 0.3229230 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sample-5 </td>
   <td style="text-align:right;"> 2.108225 </td>
   <td style="text-align:right;"> 0.9209336 </td>
   <td style="text-align:right;"> 0.5448817 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sample-6 </td>
   <td style="text-align:right;"> 2.073329 </td>
   <td style="text-align:right;"> 0.9188604 </td>
   <td style="text-align:right;"> 0.5692406 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sample-7 </td>
   <td style="text-align:right;"> 3.012372 </td>
   <td style="text-align:right;"> 0.8722852 </td>
   <td style="text-align:right;"> 0.1726627 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sample-8 </td>
   <td style="text-align:right;"> 2.170369 </td>
   <td style="text-align:right;"> 0.9174670 </td>
   <td style="text-align:right;"> 0.5372549 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sample-9 </td>
   <td style="text-align:right;"> 3.337931 </td>
   <td style="text-align:right;"> 0.8318029 </td>
   <td style="text-align:right;"> 0.1097180 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sample-10 </td>
   <td style="text-align:right;"> 2.829861 </td>
   <td style="text-align:right;"> 0.8840996 </td>
   <td style="text-align:right;"> 0.2682443 </td>
  </tr>
</tbody>
</table>



## Core microbiota analysis


**Prevalence**

Relative population frequencies; at 1% compositional abundance threshold:


```r
head(prevalence(pseq.rel, detection = 1, sort = TRUE))
```

```
##  Yersinia et rel.  Xanthomonadaceae  Wissella et rel. Weissella et rel. 
##                 0                 0                 0                 0 
##            Vibrio       Veillonella 
##                 0                 0
```


This returns names of the core taxa 


```r
core.taxa.standard <- core_members(pseq.rel, detection = 0, prevalence = 50/100)
```

A full phyloseq object of the core microbiota is obtained as follows:


```r
pseq.core <- core(pseq.rel, detection = 0, prevalence = .5)
```

Core papers: [@Jalanka-Tuovinen11, @Salonen12cmi]

See also core abundance.

See also related functions for the analysis of rare and variable taxa (noncore_members; noncore_abundance; rare_members; rare_abundance; low_abundance).

This visualization method has been used for instance in [Intestinal microbiome landscaping: Insight in community assemblage and implications for microbial modulation strategies](https://academic.oup.com/femsre/article/doi/10.1093/femsre/fuw045/2979411/Intestinal-microbiome-landscaping-insight-in#58802539). Shetty et al. _FEMS Microbiology Reviews_ fuw045, 2017.

Note that you can order the taxa on the heatmap with the order.taxa argument.


```r
# Core with compositionals:
prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-3), log10(.2), length = 10)

# Also define gray color palette
gray <- gray(seq(0,1,length=5))
p <- plot_core(pseq.rel, plot.type = "heatmap", colours = gray,
    prevalences = prevalences, detections = detections) +
    xlab("Detection Threshold (Relative Abundance (%))")
print(p)    
```

<img src="figure/core-example3-1.png" title="plot of chunk core-example3" alt="plot of chunk core-example3" width="200px" />


## Microbiome Landscaping

See also composition.

[Microbiome Landscaping](https://academic.oup.com/femsre/article/doi/10.1093/femsre/fuw045/2979411/Intestinal-microbiome-landscaping-insight-in#58802539) refers to the analysis and illustration of population frequencies. Typically, these are wrappers based on standard ordination methods (for more examples, see [ordination examples](http://microbiome.github.io/microbiome/Ordination.html))


### Two-dimensional microbiome landscape

Load example data:


```r
library(phyloseq)
library(ggplot2)

data(dietswap)
pseq <- dietswap

# Convert to compositional data
pseq.rel <- microbiome::transform(pseq, "compositional")

# Pick core taxa
pseq.core <- core(pseq.rel, detection = 5/100, prevalence = 50/100)
pseq.core <- subset_samples(pseq.core, sex == "Female" &
	                               bmi_group == "overweight")
```


Visualize the microbiome landscape (sample similarities on two-dimensional projection). For direct access to the ordination coordinates, use the following:


```r
# Project the samples with the given method and dissimilarity measure. 
# Ordinate the data; note that some ordinations are sensitive to random seed
# "quiet" is used to suppress intermediate outputs
set.seed(423542)
quiet(proj <- get_ordination(pseq.core, "NMDS", "bray"))

# Same with a generic data.frame
# (note that random seed will affect the exact ordination)
p <- plot_landscape(proj[, 1:2], col = proj$nationality, legend = T)
print(p)

# Visualize sample names:
ax1 <- names(proj)[[1]]
ax2 <- names(proj)[[2]]
p <- ggplot(aes_string(x = ax1, y = ax2, label = "sample"), data = proj) +
       geom_text(size = 2)
print(p)
```

<img src="figure/landscape4-1.png" title="plot of chunk landscape4" alt="plot of chunk landscape4" width="400px" /><img src="figure/landscape4-2.png" title="plot of chunk landscape4" alt="plot of chunk landscape4" width="400px" />




### Cross-correlating data sets

Cross-correlate columns of two data sets from related to microbiome and blood serum lipids associations ([PeerJ 1:e32](https://peerj.com/articles/32/)).

The function returns correlations, raw p-values, and fdr estimates (not strictly proper as the comparisons are not independent). Here robust biweight midcorrelation ('bicor') from the [WGCNA package](http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/). Keep only those elements that have at least only one significant correlation (n.signif):


```r
# Load example data 
otu <- peerj32$microbes 
lipids <- peerj32$lipids 

# Define data sets to cross-correlate
x <- log10(otu) # OTU Log10 (44 samples x 130 genera)
y <- as.matrix(lipids) # Lipids (44 samples x 389 lipids)

# Cross correlate data sets
correlations <- associate(x, y, method = "bicor", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)

# Or, alternatively, the same output is also available in a handy table format
correlation.table <- associate(x, y, method = "bicor", mode = "table", p.adj.threshold = 0.05, n.signif = 1)

kable(head(correlation.table))
```



|    |X1                               |X2         | Correlation|     p.adj|
|:---|:--------------------------------|:----------|-----------:|---------:|
|833 |Ruminococcus gnavus et rel.      |TG(54:5).2 |   0.7207818| 0.0017385|
|547 |Ruminococcus gnavus et rel.      |TG(52:5)   |   0.6996301| 0.0031929|
|141 |Eubacterium cylindroides et rel. |PC(40:3)   |  -0.6771286| 0.0038006|
|144 |Helicobacter                     |PC(40:3)   |  -0.6838424| 0.0038006|
|437 |Ruminococcus gnavus et rel.      |TG(50:4)   |   0.6852226| 0.0038006|
|525 |Ruminococcus gnavus et rel.      |TG(52:4).1 |   0.6716223| 0.0038006|

### Association heatmaps

Rearrange the data and plot the heatmap and mark significant correlations with stars to reproduce microbiota-lipidome heatmap from [Lahti et al. PeerJ (2013)](https://peerj.com/articles/32/) (the ordering of rows and columns may be different): 


```r
heat(correlation.table, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05) 
```

![plot of chunk heatmap-example-stars3](figure/heatmap-example-stars3-1.png)



### Bagged RDA

Bagged RDA provides added robustness in the analysis compared to the standard RDA. Fit bagged (bootstrap aggregated) RDA on a phyloseq object.

For PCoA, RDA, and other popular techniques, see [tutorial](http://microbiome.github.io/microbiome/Ordination.html).


```r
pseq <- peerj32$phyloseq # phyloseq data

# Core taxa to speed up examples
pseq <- core(pseq, detection = 10^2, prevalence = 95/100)
pseq.trans <- microbiome::transform(pseq, "hell") # Hellinger transform

# In any real study, use bs.iter = 100 or higher
# to achieve meaningful benefits from the bagged version.
# In this example we use bs.iter = 2 just to speed up the
# example code for educational purposes
res <- rda_bagged(pseq.trans, "group", bs.iter=2)
```

Visualizing bagged RDA:


```r
plot_rda_bagged(res)
```

![plot of chunk rda6](figure/rda6-1.png)



## Regression plots

Regression curve with smoothed error bars based on the [Visually-Weighted Regression](http://www.fight-entropy.com/2012/07/visually-weighted-regression.html) by Solomon M. Hsiang. The sorvi implementation extends [Felix Schonbrodt's original code](http://www.nicebread.de/visually-weighted-watercolor-plots-new-variants-please-vote/). See also [potential analysis](Potential.html).


```r
pseq <- atlas1006
p <- plot_regression(diversity ~ age, meta(pseq))
print(p)
```

![plot of chunk variability-regression](figure/variability-regression-1.png)


## Microbiome stability analysis

It has been reported that certain microbial groups exhibit bi-stable abundance distributions with an instable intermediate abundance range. [Lahti et al. 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html) used correlation analysis to quantify how the distance from the intermediate abundance region (50% quantile) is associated with the observed shifts between consecutive time points. 


```r
intermediate.stability <- intermediate_stability(pseq, output = "scores")
```


## Bimodality analysis

Plot the subject abundance variation lineplot (**Variation lineplot**) and **Bimodality hotplot** for a given taxon as in [Lahti et al. 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html). The bi-stable Dialister has bimodal population distribution and reduced temporal stability within subjects at intermediate abundances. 76 subjects in the HITChip Atlas data (1006) set have also short time series available for temporal stability analysis:


```r
# Use relative abundances
pseq <- transform(atlas1006, "compositional")

# For cross-sectional analysis, include only the zero time point:
pseq0 <- subset_samples(pseq, time == 0)
```


### Bimodality indicators

Bimodality of the abundance distribution provides an indirect indicator of bistability. [Multimodality test with bootstrap](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html). Check the [bimodality page](Bimodality.html) for more examples on bimodality indicators. Bimodality is better estimated from log10 abundances.


```r
bimodality <- bimodality(
    transform(pseq0, "log10"), 
    method = "potential_analysis", bs.iter = 20)
```


## Variation lineplot and bimodality hotplot

Pick subset of the [HITChip Atlas data set](http://doi.org/10.5061/dryad.pk75d) and plot the subject abundance variation lineplot (**Variation tip plot**) and **Bimodality hotplot** for a given taxon as in [Lahti et al. 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html). The bi-stable Dialister has bimodal population distribution and reduced temporal stability within subjects at intermediate abundances.


```r
# Bimodality hotplot:
# Consider the baseline time point for each subject
p <- hotplot(pseq0, "Dialister", tipping.point = 0.004)
print(p)

# Set tipping point manually in this example
pv <- plot_tipping(pseq, "Dialister", tipping.point = 0.004)
print(pv)
```

<img src="figure/stability-variationplot2-1.png" title="plot of chunk stability-variationplot2" alt="plot of chunk stability-variationplot2" width="430px" /><img src="figure/stability-variationplot2-2.png" title="plot of chunk stability-variationplot2" alt="plot of chunk stability-variationplot2" width="430px" />


## Other tools

 * [Heatmaps for microbiome analysis](Composition.html) 
 * [Networks](Comparisons.html)
 * [Group-wise community comparisons](Comparisons.html)


### Acknowledgements

Thanks to all [contributors](https://github.com/microbiome/microbiome/graphs/contributors). Financial support has been provided by Academy of Finland (grants 256950 and 295741), [University of Turku](http://www.utu.fi/en/Pages/home.aspx), Department of Mathematics and Statistics, [VIB lab for Bioinformatics and (eco-)systems biology](http://www.vib.be/en/research/scientists/Pages/Jeroen-Raes-Lab.aspx), VIB/KULeuven, Belgium, [Molecular Ecology group](http://www.mib.wur.nl/UK/), Laboratory of Microbiology, Wageningen University, Netherlands, and [Department of Veterinary Bioscience](http://www.vetmed.helsinki.fi/apalva/index.htm), University of Helsinki, Finland. This work relies on the independent [phyloseq](https://github.com/joey711/phyloseq) package and data structures for R-based microbiome analysis developed by Paul McMurdie and Susan Holmes. This work also utilizes a number of independent R extensions, including ade4 [@ade4], dplyr [@dplyr], ggplot2 [@ggplot2], MASS [@MASS], moments [@moments], phyloseq [@phyloseq], RColorBrewer [@RColorBrewer], scales [@scales], stats [@stats], tidyr [@tidyr], and vegan [@vegan].



# References



[1] D. Chessel, A. Dufour and J. Thioulouse. "The ade4 package-I-
One-table methods". In: _R News_ 4 (2004), pp. 5-10.
[1] S. Dray and A. Dufour. "The ade4 package: implementing the
duality diagram for ecologists". In: _Journal of Statistical
Software_ 22.4 (2007), pp. 1-20.
[1] S. Dray, A. Dufour and D. Chessel. "The ade4 package-II:
Two-table and K-table methods." In: _R News_ 7.2 (2007), pp.
47-52.
[1] J. Jalanka-Tuovinen, A. Salonen, J. N. ä, et al. "Intestinal
microbiota in healthy adults: Temporal analysis reveals individual
and common core and relation to intestinal symptoms". In: _PLoS
One_ 6.7 (2011), p. e23035.
[1] L. Komsta and F. Novomestky. _moments: Moments, cumulants,
skewness, kurtosis and related tests_. R package version 0.14.
2015. <URL: https://CRAN.R-project.org/package=moments>.
[1] L. Lahti, A. Salonen, R. A. Kekkonen, et al. "Associations
between the human intestinal microbiota, Lactobacillus rhamnosus
GG and serum lipids indicated by integrated analysis of
high-throughput profiling data". In: _PeerJ_ 1 (2013), p. e32.
<URL: http://dx.doi.org/10.7717/peerj.32>.
[1] P. J. McMurdie and S. Holmes. "phyloseq: An R package for
reproducible interactive analysis and graphics of microbiome
census data". In: _PLoS ONE_ 8.4 (2013), p. e61217. <URL:
http://dx.plos.org/10.1371/journal.pone.0061217>.
[1] E. Neuwirth. _RColorBrewer: ColorBrewer Palettes_. R package
version 1.1-2. 2014. <URL:
https://CRAN.R-project.org/package=RColorBrewer>.
[1] J. Oksanen, F. G. Blanchet, M. Friendly, et al. _vegan:
Community Ecology Package_. R package version 2.4-3. 2017. <URL:
https://CRAN.R-project.org/package=vegan>.
[1] R Core Team. _R: A Language and Environment for Statistical
Computing_. R Foundation for Statistical Computing. Vienna,
Austria, 2017. <URL: https://www.R-project.org/>.
[1] A. Salonen, J. Salojärvi, L. Lahti, et al. "The adult
intestinal core microbiota is determined by analysis depth and
health status". In: _Clinical Microbiology and Infection_
18.Suppl. 4 (2012), p. 16–20. <URL:
http://onlinelibrary.wiley.com/doi/10.1111/j.1469-0691.2012.03855.x/abstract>.
[1] W. N. Venables and B. D. Ripley. _Modern Applied Statistics
with S_. Fourth. ISBN 0-387-95457-0. New York: Springer, 2002.
<URL: http://www.stats.ox.ac.uk/pub/MASS4>.
[1] H. Wickham. _ggplot2: Elegant Graphics for Data Analysis_.
Springer-Verlag New York, 2009. ISBN: 978-0-387-98140-6. <URL:
http://ggplot2.org>.
[1] H. Wickham. _scales: Scale Functions for Visualization_. R
package version 0.4.1. 2016. <URL:
https://CRAN.R-project.org/package=scales>.
[1] H. Wickham. _tidyr: Easily Tidy Data with 'spread()' and
'gather()' Functions_. R package version 0.6.3. 2017. <URL:
https://CRAN.R-project.org/package=tidyr>.
[1] H. Wickham and R. Francois. _dplyr: A Grammar of Data
Manipulation_. R package version 0.5.0. 2016. <URL:
https://CRAN.R-project.org/package=dplyr>.


