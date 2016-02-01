## Standard data processing operations

The high-quality [phyloseq package](http://joey711.github.io/phyloseq/) provides a complete set of tools for data subsetting, aggregating and filtering.

Download example data from [O'Keefe et al. Nat. Comm. 6:6342, 2015](http://dx.doi.org/10.1038/ncomms7342) from [Data Dryad](http://dx.doi.org/10.5061/dryad.1mn1n). This is a two-week diet swap study between western (USA) and traditional (rural Africa) diets including microbiota profiling:


```r
library(microbiome)
pseq <- download_microbiome("dietswap")
```

```
## Error in file(file, "rt"): cannot open the connection
```


### Sample operations

Sample names and variables


```r
head(sample_names(pseq))
```

```
## Error in sample_names(pseq): error in evaluating the argument 'physeq' in selecting a method for function 'sample_names': Error: object 'pseq' not found
```

Sample sums


```r
head(sample_sums(pseq))
```

```
## Error in otu_table(x): error in evaluating the argument 'object' in selecting a method for function 'otu_table': Error: object 'pseq' not found
```

Abundance of a given species in each sample


```r
head(get_sample(pseq, taxa_names(pseq)[1]))
```

```
## Error in get_sample(pseq, taxa_names(pseq)[1]): error in evaluating the argument 'physeq' in selecting a method for function 'get_sample': Error: object 'pseq' not found
```

Filter samples


```r
f1 <- filterfun_sample(topp(0.1))
taxa <- genefilter_sample(pseq, f1, A = round(0.5 * nsamples(pseq)))
```

```
## Error in genefilter_sample(pseq, f1, A = round(0.5 * nsamples(pseq))): error in evaluating the argument 'X' in selecting a method for function 'genefilter_sample': Error: object 'pseq' not found
```

Select samples by specific metadata fields:


```r
pseq.subset <- subset_samples(pseq, nationality == "AFR")
```

```
## Error in sample_data(physeq): error in evaluating the argument 'object' in selecting a method for function 'sample_data': Error: object 'pseq' not found
```


### Data transformations

The microbiome package provides a wrapper for standard sample/OTU transformations. For arbitrary transformations, use the transform_sample_counts function in the phyloseq package.

Log10 transformation (log(1+x) if the data contains zeroes)


```r
pseq.log <- transform_phyloseq(pseq, "log10")
```

```
## Error in otu_table(x): error in evaluating the argument 'object' in selecting a method for function 'otu_table': Error: object 'pseq' not found
```

Z transformation:


```r
# Z-transform OTUs
pseq.zotu <- transform_phyloseq(pseq, "Z", "OTU")
```

```
## Error in otu_table(x): error in evaluating the argument 'object' in selecting a method for function 'otu_table': Error: object 'pseq' not found
```

```r
# Z-transform samples
pseq.zsample <- transform_phyloseq(pseq, "Z", "sample")
```

```
## Error in otu_table(x): error in evaluating the argument 'object' in selecting a method for function 'otu_table': Error: object 'pseq' not found
```

Relative abundances (the input data needs to be in absolute scale, not logarithmic!):


```r
pseq1 <- transform_phyloseq(pseq, "relative.abundance", "OTU")
```

```
## Error in otu_table(x): error in evaluating the argument 'object' in selecting a method for function 'otu_table': Error: object 'pseq' not found
```

```r
pseq2 <- transform_sample_counts(pseq, function(x) x/sum(x))
```

```
## Error in taxa_are_rows(physeq): error in evaluating the argument 'physeq' in selecting a method for function 'taxa_are_rows': Error: object 'pseq' not found
```

Relative abundance for plain abundance matrix:


```r
dat <- otu_table(pseq)@.Data
```

```
## Error in otu_table(pseq): error in evaluating the argument 'object' in selecting a method for function 'otu_table': Error: object 'pseq' not found
```

```r
rel <- relative.abundance(dat, det.th = 0)
```



### Variable operations

Sample variable names


```r
sample_variables(pseq)
```

```
## Error in colnames(sample_data(physeq, errorIfNULL)): error in evaluating the argument 'x' in selecting a method for function 'colnames': Error in sample_data(physeq, errorIfNULL) : 
##   error in evaluating the argument 'object' in selecting a method for function 'sample_data': Error: object 'pseq' not found
```

Pick variable values for a given variable


```r
head(get_variable(pseq, sample_variables(pseq)[1]))
```

```
## Error in sample_data(physeq, FALSE): error in evaluating the argument 'object' in selecting a method for function 'sample_data': Error: object 'pseq' not found
```

```r
# Assign fields to sample metadata
# sample_data(GP)$human <- ..
```

### Taxa operations


Number of taxa


```r
ntaxa(pseq)
```

```
## Error in ntaxa(pseq): error in evaluating the argument 'physeq' in selecting a method for function 'ntaxa': Error: object 'pseq' not found
```


Names


```r
ranks <- rank_names(pseq)
```

```
## Error in colnames(tax_table(physeq, errorIfNULL)): error in evaluating the argument 'x' in selecting a method for function 'colnames': Error in tax_table(physeq, errorIfNULL) : 
##   error in evaluating the argument 'object' in selecting a method for function 'tax_table': Error: object 'pseq' not found
```

```r
taxa <- taxa_names(pseq)
```

```
## Error in taxa_names(pseq): error in evaluating the argument 'physeq' in selecting a method for function 'taxa_names': Error: object 'pseq' not found
```


Prune taxa:


```r
taxa <- levelmap(NULL, "Phylum", "Genus", tax_table(pseq))$Bacteroidetes
```

```
## Error in tax_table(pseq): error in evaluating the argument 'object' in selecting a method for function 'tax_table': Error: object 'pseq' not found
```

```r
# With given taxon names
ex2 <- prune_taxa(taxa, pseq)
```

```
## Error in prune_taxa(taxa, pseq): error in evaluating the argument 'taxa' in selecting a method for function 'prune_taxa': Error: object 'taxa' not found
```

```r
# Taxa with positive sum across samples
ex3 <- prune_taxa(taxa_sums(pseq) > 0, pseq)
```

```
## Error in prune_taxa(taxa_sums(pseq) > 0, pseq): error in evaluating the argument 'taxa' in selecting a method for function 'prune_taxa': Error in otu_table(x) : 
##   error in evaluating the argument 'object' in selecting a method for function 'otu_table': Error: object 'pseq' not found
```


Subset taxa:


```r
pseq <- subset_taxa(pseq, Phylum == "Bacteroidetes")
```

```
## Error in tax_table(physeq): error in evaluating the argument 'object' in selecting a method for function 'tax_table': Error: object 'pseq' not found
```


Filter by user-specified function values (here variance):


```r
f <- filter_taxa(pseq, function(x) var(x) > 1e-05, TRUE)
```

```
## Error in inherits(x, get.component.classes()): object 'pseq' not found
```


List unique phylum-level groups: 


```r
head(get_taxa_unique(pseq, "Phylum"))
```

```
## Error in unique(as(tax_table(physeq, errorIfNULL)[, taxonomic.rank], "character")): error in evaluating the argument 'x' in selecting a method for function 'unique': Error in tax_table(physeq, errorIfNULL) : 
##   error in evaluating the argument 'object' in selecting a method for function 'tax_table': Error: object 'pseq' not found
```

Pick detected taxa by sample name:


```r
samplename <- sample_names(pseq)[[1]]
```

```
## Error in sample_names(pseq): error in evaluating the argument 'physeq' in selecting a method for function 'sample_names': Error: object 'pseq' not found
```

```r
tax.abundances <- get_taxa(pseq, samplename)
```

```
## Error in get_taxa(pseq, samplename): error in evaluating the argument 'physeq' in selecting a method for function 'get_taxa': Error: object 'pseq' not found
```


Taxa sums


```r
head(taxa_sums(pseq))
```

```
## Error in otu_table(x): error in evaluating the argument 'object' in selecting a method for function 'otu_table': Error: object 'pseq' not found
```


### Merging operations

Aggregate OTUs to higher taxonomic levels, use (on HITChip we use L1/L2 instead of Phylum/Genus). See also [merge_samples and and merge_taxa](http://joey711.github.io/phyloseq/merge.html):


```r
pseq <- read_hitchip(data.directory, method = "frpa")$pseq
```

```
## Loading pre-calculated RPA preprocessing parameters
```

```r
pseq.L2 <- aggregate_taxa(pseq, level = "L2")
```

Merging phyloseq objects


```r
merge_phyloseq(pseq1, pseq2)
```


### Rarification


```r
pseq.rarified <- rarefy_even_depth(pseq)
```


