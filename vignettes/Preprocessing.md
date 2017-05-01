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
pseq <- core(subset_samples(atlas1006, nationality == "EasternEurope"), detection = 10^2, prevalence = 50/100) # Rename the data and pick subset for faster examples
```


### Retrieving data elements from a phyloseq object  

A phyloseq object contains OTU table (taxa abundances), sample
metadata, taxonomy table (mapping between OTUs and higher-level
taxonomic classifications), and phylogenetic tree (relations between
the taxa). Some of these are optional.


Pick metadata as data.frame:


```r
meta <- meta(pseq)
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


Melting phyloseq data for easier plotting:


```r
df <- psmelt(pseq)
kable(head(df))
```



|     |OTU                              |Sample     | Abundance| age|gender |nationality   |DNA_extraction_method |project | diversity|bmi_group |subject | time|sample     |Phylum                 |Genus                            |
|:----|:--------------------------------|:----------|---------:|---:|:------|:-------------|:---------------------|:-------|---------:|:---------|:-------|----:|:----------|:----------------------|:--------------------------------|
|597  |Escherichia coli et rel.         |Sample-910 |    179473|  56|male   |EasternEurope |NA                    |27      |      5.51|NA        |910     |    0|Sample-910 |Proteobacteria         |Escherichia coli et rel.         |
|1211 |Subdoligranulum variable at rel. |Sample-911 |    162402|  45|male   |EasternEurope |NA                    |27      |      5.62|NA        |911     |    0|Sample-911 |Clostridium cluster IV |Subdoligranulum variable at rel. |
|1215 |Subdoligranulum variable at rel. |Sample-919 |    144757|  64|male   |EasternEurope |NA                    |28      |      5.47|NA        |919     |    0|Sample-919 |Clostridium cluster IV |Subdoligranulum variable at rel. |
|1201 |Subdoligranulum variable at rel. |Sample-908 |    123448|  53|male   |EasternEurope |NA                    |27      |      5.72|NA        |908     |    0|Sample-908 |Clostridium cluster IV |Subdoligranulum variable at rel. |
|223  |Bifidobacterium                  |Sample-917 |    109982|  43|male   |EasternEurope |NA                    |28      |      5.80|NA        |917     |    0|Sample-917 |Actinobacteria         |Bifidobacterium                  |
|1209 |Subdoligranulum variable at rel. |Sample-909 |     97965|  64|female |EasternEurope |NA                    |27      |      5.66|NA        |909     |    0|Sample-909 |Clostridium cluster IV |Subdoligranulum variable at rel. |



### Sample operations

Sample names and variables


```r
head(sample_names(pseq))
```

```
## [1] "Sample-312" "Sample-907" "Sample-908" "Sample-909" "Sample-910"
## [6] "Sample-911"
```

Total OTU abundance in each sample


```r
s <- sample_sums(pseq)
```

Abundance of a given species in each sample


```r
head(abundances(pseq)["Akkermansia",])
```

```
## Sample-312 Sample-907 Sample-908 Sample-909 Sample-910 Sample-911 
##       3649       7446       1461       2633       1052       2023
```


Filtering samples


```r
f1 <- filterfun_sample(topp(0.1))
taxa <- genefilter_sample(pseq, f1, A = round(0.5 * nsamples(pseq)))
```


Select a subset by metadata fields:


```r
pseq.subset <- subset_samples(pseq, nationality == "AFR")
```


Select a subset by providing sample names: 


```r
# Check sample names for African Females in this phyloseq object
s <- rownames(subset(meta(pseq), nationality == "AFR" & sex == "Female"))

# Pick the phyloseq subset with these sample names
pseq.subset2 <- prune_samples(s, pseq)
```


Pick samples at the baseline time points only:


```r
pseq0 <- baseline(pseq)
```



### Data transformations

The microbiome package provides a wrapper for standard sample/OTU transforms. For arbitrary transforms, use the transform_sample_counts function in the phyloseq package. Log10 transform is log(1+x) if the data contains zeroes. Also "Z",
"clr", "hellinger", and "shift" are available as common
transformations. Relative abundances (note that the input data needs
to be in absolute scale, not logarithmic!):


```r
pseq.compositional <- microbiome::transform(pseq, "compositional")
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
## [1] 36 40 53 64 56 45
```

Assign new fields to metadata


```r
# Calculate diversity for samples
div <- global(pseq, index = "shannon")

# Assign the estimated diversity to sample metadata
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


Subset taxa


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
## [1] "Verrucomicrobia"          "Proteobacteria"          
## [3] "Bacteroidetes"            "Clostridium cluster XIVa"
## [5] "Clostridium cluster IV"   "Clostridium cluster XI"
```

Pick the taxa abundances for a given sample:


```r
samplename <- sample_names(pseq)[[1]]

# Pick abundances for a particular taxon
tax.abundances <- abundances(pseq)[, samplename]
```


### Merging operations

Aggregate taxa to higher taxonomic levels. This is particularly useful if the phylogenetic tree is missing. When it is available, see [merge_samples, merge_taxa and tax_glom](http://joey711.github.io/phyloseq/merge.html)).

Put the less abundant taxa together in the "Other" category:


```r
pseq2 <- aggregate_taxa(pseq, "Phylum", top = 5) 
```


Merging phyloseq objects with the phyloseq package:


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



|Groups |Factor |  n|    pct|
|:------|:------|--:|------:|
|lean   |female |  1| 100.00|
|NA     |female |  6|  42.86|
|NA     |male   |  8|  57.14|


## Formatting the Phyloseq Object 

For this example, we use example data from [Halfvarson J., et al. Nature Microbiology, 2017](http://www.nature.com/articles/nmicrobiol20174) as downloaded from [Qiita](https://qiita.ucsd.edu/study/description/1629).  



```r
library(microbiome)
data(DynamicsIBD)
p0 <- DynamicsIBD

# Check rank names
rank_names(p0) 
```

```
## [1] "Rank1" "Rank2" "Rank3" "Rank4" "Rank5" "Rank6" "Rank7"
```

```r
# Rename rank names
# If you have used the default parsing function for greengenes taxonomy
# for creating the phyloseq object then this step can be omitted
colnames(tax_table(p0)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species") 

# Check the taxonomy information stored in the phyloseq object.  
library(knitr)
kable(head(tax_table(p0)))
```



|        |Domain      |Phylum        |Class         |Order              |Family              |Genus               |Species |
|:-------|:-----------|:-------------|:-------------|:------------------|:-------------------|:-------------------|:-------|
|577110  |k__Bacteria |p__Firmicutes |c__Clostridia |o__Clostridiales   |f__                 |g__                 |s__     |
|181342  |k__Bacteria |p__Firmicutes |c__Clostridia |o__Clostridiales   |f__Clostridiaceae   |g__02d06            |s__     |
|581609  |k__Bacteria |p__Firmicutes |c__Clostridia |o__Clostridiales   |f__Ruminococcaceae  |g__                 |s__     |
|4341234 |k__Bacteria |p__Firmicutes |c__Clostridia |o__Clostridiales   |f__Peptococcaceae   |g__Desulfotomaculum |s__     |
|181348  |k__Bacteria |p__Firmicutes |c__Clostridia |o__Clostridiales   |f__Lachnospiraceae  |g__Coprococcus      |s__     |
|4467992 |k__Bacteria |p__Firmicutes |c__Bacilli    |o__Lactobacillales |f__Streptococcaceae |g__Streptococcus    |s__     |

Note that the not all the OTUs are classified down to the lowest taxonomic level (here, species level). For OTU level testing of differential abundance, you may need information on the specific otu number or taxonomy of the OTU. This can help in easily tracing back the sequence and also make the plots with most fine-resolution taxonomic classification possible. To analyze other taxonomic levels, you can use the phyloseq::tax_glom function (which is relatively slow).


```r
# using the p0 object, agglomerate at genus level
# commented out but shown:
# p0 <- tax_glom(p0, "Genus")
p1 <- p0

# Improve the classification
pf <- format_phyloseq(p1)

# Check the taxonomy with the formatted phyloseq object
kable(head(tax_table(pf)))
```



|        |Domain   |Phylum     |Class      |Order           |Family                  |Genus                     |Species                     |
|:-------|:--------|:----------|:----------|:---------------|:-----------------------|:-------------------------|:---------------------------|
|577110  |Bacteria |Firmicutes |Clostridia |Clostridiales   |o__Clostridiales_577110 |o__Clostridiales_577110   |o__Clostridiales_577110     |
|181342  |Bacteria |Firmicutes |Clostridia |Clostridiales   |Clostridiaceae          |02d06                     |f__02d06_181342             |
|581609  |Bacteria |Firmicutes |Clostridia |Clostridiales   |Ruminococcaceae         |f__Ruminococcaceae_581609 |f__Ruminococcaceae_581609   |
|4341234 |Bacteria |Firmicutes |Clostridia |Clostridiales   |Peptococcaceae          |Desulfotomaculum          |f__Desulfotomaculum_4341234 |
|181348  |Bacteria |Firmicutes |Clostridia |Clostridiales   |Lachnospiraceae         |Coprococcus               |f__Coprococcus_181348       |
|4467992 |Bacteria |Firmicutes |Bacilli    |Lactobacillales |Streptococcaceae        |Streptococcus             |f__Streptococcus_4467992    |

