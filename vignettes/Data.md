---
title: "Data"
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
  %\VignetteIndexEntry{microbiome tutorial - data}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->


## Importing microbiome data in R

To import your own data in R phyloseq format, you can use the [import functions](http://joey711.github.io/phyloseq/import-data) (from the independent phyloseq R package) for standard formats (mothur, qiime etc.).

Alternatively, you can read your data in R (read.table, read.csv or other standard functions). You can then pick the OTU table (taxa x samples matrix x), metadata (samples x variables data.frame y), and taxonomic table (taxa x phylogenetic levels matrix z) in R. Make sure that these objects do not contain other fields. They should also have row and column names (samples, taxa, phylogenetic levels), and the rows and columns should match between these objects. Then you can create your own phyloseq object pseq as follows:


```r
library(phyloseq)

# x: taxa x samples OTU matrix;
# put taxa and samples as row and col names, respectively
pseq <- phyloseq(otu_table(x, taxa_are_rows = TRUE))

# Add taxonomy
# y: taxa x levels taxonomy matrix
# Create tax_table
TAX <- tax_table(y[rownames(x), ])
# Now the tax_table rownames are newly created and different from the
# original rownames in y
# and should be renamed so that they match with the otu_table above
rownames(TAX) <- rownames(y)
# Merge the taxonomy with the phyloseq object
pseq <- merge_phyloseq(pseq, TAX)
 
# Add metadata
# z: samples x variables data.frame
pseq <- merge_phyloseq(pseq, z[colnames(x),])
```


See also examples on [filtering, subsetting and other data processing](Preprocessing.html) for phyloseq objects.



## Microbiome example data sets

### Intestinal microbiota profiling of 1006 Western adults

[The HITChip Atlas](Atlas.html) data set is available via the microbiome R package in phyloseq format, and via [Data Dryad](http://doi.org/10.5061/dryad.pk75d) in tabular format. This data set from [Lahti et al. Nat. Comm. 5:4344, 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html) comes with 130 genus-like taxonomic groups across 1006 western adults with no reported health complications. Some subjects have also short time series. Load the data in R with:


```r
library(microbiome)
data(atlas1006) 
print(atlas1006)
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 130 taxa and 1172 samples ]
## sample_data() Sample Data:       [ 1172 samples by 10 sample variables ]
## tax_table()   Taxonomy Table:    [ 130 taxa by 2 taxonomic ranks ]
```


### Diet swap between Rural and Western populations

A two-week diet swap study between western (USA) and traditional
(rural Africa) diets, reported in [O'Keefe et al. Nat. Comm. 6:6342,
2015](http://dx.doi.org/10.1038/ncomms7342). The data is also
available for download from [Data
Dryad](http://dx.doi.org/10.5061/dryad.1mn1n). Load in R with:


```r
data(dietswap)
print(dietswap)
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
## sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
## tax_table()   Taxonomy Table:    [ 130 taxa by 2 taxonomic ranks ]
```


### Dynamics of the human gut microbiome in inflammatory bowel disease

Data set from [Halfvarson et al. Nature Microbiology 2, 2017](http://www.nature.com/articles/nmicrobiol20174) characterizes longitudinal fluctuations of human intestinal microbiota in Inflammatory bowel disease (IBD) subjects and compare it with those of healthy individuals. The authors use Illumina HiSeq 2000 for the V4 region with 515F/806RBC. This data was downloaded from [Qiita](https://qiita.ucsd.edu/study/description/1629) and the Qiita study ID is 1629.
Load the data in R with:


```r
data(DynamicsIBD)
print(DynamicsIBD)
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 10996 taxa and 683 samples ]
## sample_data() Sample Data:       [ 683 samples by 32 sample variables ]
## tax_table()   Taxonomy Table:    [ 10996 taxa by 7 taxonomic ranks ]
```

### Intestinal microbiota versus blood metabolites

Data set from [Lahti et al. PeerJ 1:e32,
2013](https://peerj.com/articles/32/) characterizes associations
between human intestinal microbiota and blood serum lipids. Note that
this data set contains an additional data matrix of lipid
species. Load the data in R with:


```r
data(peerj32)
print(names(peerj32))
```

```
## [1] "lipids"   "microbes" "meta"     "phyloseq"
```
