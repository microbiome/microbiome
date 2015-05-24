## Filtering and pruning

The external high-quality [phyloseq package](http://joey711.github.io/phyloseq/) provides a complete set of tools for data preprocessing and filtering purposes. Download example data from [O'Keefe et al. Nat. Comm. 6:6342, 2015](http://dx.doi.org/10.1038/ncomms7342) from [Data Dryad](http://dx.doi.org/10.5061/dryad.1mn1n). This is a two-week diet swap study between western (USA) and traditional (rural Africa) diets including microbiota profiling:


```r
library(microbiome)
pseq <- download_microbiome("dietswap")
```


### Estimating relative abundancies

Estimate relative abundance of the taxa in each sample. Note: the
input data set needs to be in absolute scale (not logarithmic).


```r
dat <- otu_table(pseq)@.Data
rel <- relative.abundance(dat, det.th = 0)
```

Calculating relative abundances for phyloseq objects:


```r
pseq <- transform_sample_counts(pseq, function(x) x/sum(x))
```


### Sample operations

Transform sample counts. Here, calculate relative abundances:


```r
r <- transform_sample_counts(pseq, function(x) x/sum(x))
```

Sample names and variables


```r
head(sample_names(pseq))
```

```
## [1] "Sample-1" "Sample-2" "Sample-3" "Sample-4" "Sample-5" "Sample-6"
```

Sample sums


```r
head(sample_sums(pseq))
```

```
## Sample-1 Sample-2 Sample-3 Sample-4 Sample-5 Sample-6 
##        1        1        1        1        1        1
```

Abundance for species ‘i’ in each sample


```r
head(get_sample(pseq, taxa_names(pseq)[1]))
```

```
##     Sample-1     Sample-2     Sample-3     Sample-4     Sample-5 
## 2.060778e-05 5.035640e-05 1.152133e-05 5.023935e-05 1.461157e-05 
##     Sample-6 
## 1.604835e-05
```





```r
# Pick samples by specific metadata fields
pseq.subset2 <- subset_samples(pseq.subset, nationality == "US")
```


### Variable operations

Sample variable names


```r
sample_variables(pseq)
```

```
## [1] "subject"                "gender"                
## [3] "nationality"            "group"                 
## [5] "sample"                 "timepoint"             
## [7] "timepoint.within.group" "bmi_group"
```

Pick variable values for a given variable


```r
head(get_variable(pseq, sample_variables(pseq)[3]))
```

```
## [1] AAM AFR AFR AFR AFR AFR
## Levels: AAM AFR
```

```r
# Assign fields to sample metadata
# sample_data(GP)$human <- ..
```

### Taxa operations


Filter samples


```r
f1 <- filterfun_sample(topp(0.1))
taxa <- genefilter_sample(pseq, f1, A = round(0.5 * nsamples(pseq)))
```

Prune taxa


```r
taxa <- levelmap(NULL, "Phylum", "Genus", tax_table(pseq))$Bacteroidetes

# With given taxon names
ex2 <- prune_taxa(taxa, pseq)

# Taxa with positive sum across samples
ex3 <- prune_taxa(taxa_sums(pseq) > 0, pseq)
```

Subset taxa


```r
pseq <- subset_taxa(pseq, Phylum == "Bacteroidetes")
```


Filter by user-specified function values (here variance)


```r
f <- filter_taxa(r, function(x) var(x) > 1e-05, TRUE)
```

Number of taxa


```r
ntaxa(pseq)
```

```
## [1] 16
```


Names


```r
rank_names(pseq)
```

```
## [1] "Phylum" "Genus"
```

```r
taxa_names(pseq)
```

```
##  [1] "Allistipes et rel."                
##  [2] "Bacteroides fragilis et rel."      
##  [3] "Bacteroides intestinalis et rel."  
##  [4] "Bacteroides ovatus et rel."        
##  [5] "Bacteroides plebeius et rel."      
##  [6] "Bacteroides splachnicus et rel."   
##  [7] "Bacteroides stercoris et rel."     
##  [8] "Bacteroides uniformis et rel."     
##  [9] "Bacteroides vulgatus et rel."      
## [10] "Parabacteroides distasonis et rel."
## [11] "Prevotella melaninogenica et rel." 
## [12] "Prevotella oralis et rel."         
## [13] "Prevotella ruminicola et rel."     
## [14] "Prevotella tannerae et rel."       
## [15] "Tannerella et rel."                
## [16] "Uncultured Bacteroidetes"
```


Pick taxa by rank


```r
# List unique phylum-level groups in the data
head(get_taxa_unique(pseq, "Phylum"))
```

```
## [1] "Bacteroidetes"
```

Pick taxa by sample name:


```r
get_taxa(pseq, sample_names(pseq)[[1]])
```

```
##                 Allistipes et rel.       Bacteroides fragilis et rel. 
##                       0.0397580272                       0.0523156587 
##   Bacteroides intestinalis et rel.         Bacteroides ovatus et rel. 
##                       0.0014200634                       0.0504197430 
##       Bacteroides plebeius et rel.    Bacteroides splachnicus et rel. 
##                       0.0076304988                       0.0053111868 
##      Bacteroides stercoris et rel.      Bacteroides uniformis et rel. 
##                       0.0025235163                       0.0366855946 
##       Bacteroides vulgatus et rel. Parabacteroides distasonis et rel. 
##                       0.3279166097                       0.0220540711 
##  Prevotella melaninogenica et rel.          Prevotella oralis et rel. 
##                       0.0029281781                       0.0025965802 
##      Prevotella ruminicola et rel.        Prevotella tannerae et rel. 
##                       0.0001311404                       0.0074262944 
##                 Tannerella et rel.           Uncultured Bacteroidetes 
##                       0.0064596022                       0.0002116981
```


Taxa sums


```r
head(taxa_sums(pseq))
```

```
##               Allistipes et rel.     Bacteroides fragilis et rel. 
##                        4.8226896                        3.4909888 
## Bacteroides intestinalis et rel.       Bacteroides ovatus et rel. 
##                        0.3137315                        2.0882334 
##     Bacteroides plebeius et rel.  Bacteroides splachnicus et rel. 
##                        0.8411087                        1.0741668
```


### Merging operations for phyloseq objects


```r
# merge_phyloseq()
# merge_taxa()
# tax_glom()
```


### Rarification


```r
pseq.rarified <- rarefy_even_depth(pseq)
```

```
## Error in validObject(.Object): invalid class "otu_table" object: 
##  OTU abundance data must have non-zero dimensions.
```

