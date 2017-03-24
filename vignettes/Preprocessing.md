---
title: "Preprocessing of taxonomic profiling data"
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
  %\VignetteIndexEntry{microbiome tutorial - Preprocessing}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->


## Processing taxonomic profiling data

Instructions to manipulate microbiome data sets using tools from the [phyloseq package](http://joey711.github.io/phyloseq/) and some extensions from the [microbiome package](https://github.com/microbiome/microbiome), including subsetting, aggregating and filtering.


Load example data:


```r
library(phyloseq)
library(microbiome)

data(atlas1006)   # Load the data
pseq <- atlas1006 # Rename the data
```


### Retrieving data elements from a phyloseq object  

A phyloseq object contains OTU table (taxa abundances), sample
metadata, taxonomy table (mapping between OTUs and higher-level
taxonomic classifications), and phylogenetic tree (relations between
the taxa). Some of these are optional.


Pick metadata:


```r
meta <- sample_data(pseq)
```

Taxonomy table:


```r
taxonomy <- tax_table(pseq)
```


Abundances for taxonomic groups ('OTU table') as a TaxaxSamples matrix:


```r
# Absolute abundances
otu.absolute <- abundances(pseq)

# Relative abundances
otu.relative <- abundances(pseq, "compositional")
```


Melt phyloseq data for easier plotting:


```r
df <- psmelt(pseq)
kable(head(df))
```



|       |OTU                               |Sample     | Abundance| age|gender |nationality   |DNA_extraction_method |project | diversity|bmi_group   |subject | time|sample     |Phylum        |Genus                             |
|:------|:---------------------------------|:----------|---------:|---:|:------|:-------------|:---------------------|:-------|---------:|:-----------|:-------|----:|:----------|:-------------|:---------------------------------|
|113110 |Prevotella melaninogenica et rel. |Sample-448 |    944002|  54|female |CentralEurope |o                     |18      |      5.98|lean        |448     |    0|Sample-448 |Bacteroidetes |Prevotella melaninogenica et rel. |
|113015 |Prevotella melaninogenica et rel. |Sample-360 |    902034|  45|female |CentralEurope |o                     |13      |      5.49|severeobese |360     |    0|Sample-360 |Bacteroidetes |Prevotella melaninogenica et rel. |
|112747 |Prevotella melaninogenica et rel. |Sample-190 |    862870|  34|female |CentralEurope |r                     |7       |      6.06|lean        |190     |    0|Sample-190 |Bacteroidetes |Prevotella melaninogenica et rel. |
|113109 |Prevotella melaninogenica et rel. |Sample-743 |    852350|  52|male   |US            |NA                    |19      |      5.21|obese       |743     |    0|Sample-743 |Bacteroidetes |Prevotella melaninogenica et rel. |
|112944 |Prevotella melaninogenica et rel. |Sample-366 |    851147|  52|female |CentralEurope |o                     |15      |      5.63|obese       |366     |    0|Sample-366 |Bacteroidetes |Prevotella melaninogenica et rel. |
|113639 |Prevotella melaninogenica et rel. |Sample-375 |    844482|  45|female |CentralEurope |o                     |16      |      5.64|severeobese |375     |    0|Sample-375 |Bacteroidetes |Prevotella melaninogenica et rel. |



### Sample operations

Sample names and variables


```r
head(sample_names(pseq))
```

```
## [1] "Sample-1" "Sample-2" "Sample-3" "Sample-4" "Sample-5" "Sample-6"
```

Total OTU abundance in each sample


```r
head(sample_sums(pseq))
```

```
## Sample-1 Sample-2 Sample-3 Sample-4 Sample-5 Sample-6 
##   479428   640574   449884   684997   757697   499535
```

Abundance of a given species in each sample


```r
head(get_sample(pseq, "Akkermansia"))
```

```
## Sample-1 Sample-2 Sample-3 Sample-4 Sample-5 Sample-6 
##     1319     2299    29980     3824     2133      864
```

Filter samples


```r
f1 <- filterfun_sample(topp(0.1))
taxa <- genefilter_sample(pseq, f1, A = round(0.5 * nsamples(pseq)))
```

Select samples by specific metadata fields


```r
pseq.subset <- subset_samples(pseq, nationality == "US")
```


Pick samples at the baseline time points only:


```r
pseq0 <- baseline(pseq)
```



### Data transformations

The microbiome package provides a wrapper for standard sample/OTU transforms. For arbitrary transforms, use the transform_sample_counts function in the phyloseq package.

Log10 transform (log(1+x) if the data contains zeroes). Also "Z",
"clr", and "hellinger" are available as common transforms.


```r
pseq.log <- microbiome::transform(pseq, "log10")
```

Relative abundances (the input data needs to be in absolute scale, not logarithmic!):


```r
pseq1 <- microbiome::transform(pseq, "compositional", "OTU")
pseq2 <- phyloseq::transform_sample_counts(pseq, function(x) x/sum(x))
```


### Variable operations

Sample variable names


```r
sample_variables(pseq)
```

```
##  [1] "age"                   "gender"               
##  [3] "nationality"           "DNA_extraction_method"
##  [5] "project"               "diversity"            
##  [7] "bmi_group"             "subject"              
##  [9] "time"                  "sample"
```

Pick values for a given variable


```r
head(get_variable(pseq, sample_variables(pseq)[1]))
```

```
## [1] 28 24 52 22 25 42
```

Assign new fields to metadata


```r
# Calculate diversity for samples
div <- microbiome::diversity(pseq, measures = "Shannon")$Shannon

# Assign this to sample metadata
sample_data(pseq)$diversity <- div
```

### Taxa operations


Number of taxa


```r
n <- ntaxa(pseq)
```


Most abundant taxa


```r
topx <- top_taxa(pseq, n = 10)
```


Names


```r
ranks <- rank_names(pseq)  # Taxonomic levels
taxa  <- taxa(pseq)        # Taxa names at the analysed level
```


Subset taxa:


```r
pseq.bac <- subset_taxa(pseq, Phylum == "Bacteroidetes")
```


Prune (select) taxa:


```r
# List of Genera in the Bacteroideted Phylum
taxa <- map_levels(NULL, "Phylum", "Genus", pseq)$Bacteroidetes

# With given taxon names
ex2 <- prune_taxa(taxa, pseq)

# Taxa with positive sum across samples
ex3 <- prune_taxa(taxa_sums(pseq) > 0, pseq)
```


Filter by user-specified function values (here variance):


```r
f <- filter_taxa(pseq, function(x) var(x) > 1e-05, TRUE)
```


List unique phylum-level groups: 


```r
head(get_taxa_unique(pseq, "Phylum"))
```

```
## [1] "Actinobacteria"         "Bacilli"               
## [3] "Proteobacteria"         "Verrucomicrobia"       
## [5] "Bacteroidetes"          "Clostridium cluster XV"
```

Pick the taxa abundances for a given sample:


```r
samplename <- sample_names(pseq)[[1]]

# Two ways to pick abundances for a particular taxon
tax.abundances <- get_taxa(pseq, samplename)
tax.abundances2 <- abundances(pseq)[, samplename]
```


### Merging operations

Aggregate taxa to higher taxonomic levels. This is particularly useful if the phylogenetic tree is missing. When it is available, see [merge_samples, merge_taxa and tax_glom](http://joey711.github.io/phyloseq/merge.html))


```r
pseq2 <- aggregate_taxa(pseq, "Phylum") 
```


Merging phyloseq objects


```r
merge_phyloseq(pseqA, pseqB)
```


### Rarification


```r
pseq.rarified <- rarefy_even_depth(pseq)
```


### Taxonomy 

Convert between taxonomic levels (here from Genus (Akkermansia) to
Phylum (Verrucomicrobia):


```r
m <- map_levels("Akkermansia", "Genus", "Phylum", tax_table(pseq))
print(m)
```

```
## [1] "Verrucomicrobia"
```


### Metadata

Visualize frequencies of given factor (sex) levels within the
indicated groups (group):


```r
res <- plot_frequencies(sample_data(pseq), "bmi_group", "gender")
print(res$plot)
```

![plot of chunk phylogeny-example3](figure/phylogeny-example3-1.png)

```r
# Retrieving the actual data values:
kable(head(res$data), digits = 2)
```



|Groups      |Factor |   n|   pct|
|:-----------|:------|---:|-----:|
|underweight |female |  21| 91.30|
|underweight |male   |   2|  8.70|
|lean        |female | 304| 61.66|
|lean        |male   | 189| 38.34|
|overweight  |female | 102| 50.00|
|overweight  |male   | 102| 50.00|
