## Tools focusing on the HITChip phylogenetic microarray analysis

  * [Extracting data from HITChip database](https://github.com/microbiome/HITChipDB/blob/master/vignettes/vignette.md)
  * [Cross hybridization (phylogenetic microarrays)](Crosshyb.md)
  * [Probe level studies (phylogenetic microarrays)](Probelevel.md)


### Importing HITChip data to phyloseq format

Define the data folder. 


```r
# Define example data path (replace here data.directory with your own path)
library(microbiome)
data.directory <- system.file("extdata", package = "microbiome")
print(data.directory)
```

```
## [1] "/home/lei/R/x86_64-pc-linux-gnu-library/3.2/microbiome/extdata"
```

With HITChip,
[fRPA](http://www.computer.org/csdl/trans/tb/2011/01/ttb2011010217-abs.html)
is the recommended preprocessing method. You can add new metadata
fields in the template metadata file in your HITChip data folder and
exporting it again to tab-separated .tab format. Some standard,
self-explanatory field names include 'sample', 'time', 'subject',
'group', 'gender', 'diet', 'age'. You can leave these out or include
further fields. Import HITChip phylotype-level data in
[phyloseq](https://github.com/joey711/phyloseq) format (note: the
precalculated matrices are calculated with detection.threshold = 0):


```r
pseq <- read_hitchip(data.directory, method = "frpa")$pseq
```

```
## Error in eval(expr, envir, enclos): could not find function "read_hitchip"
```

Get higher taxonomic levels, use (on HITChip we use L1/L2 instead of Phylum/Genus):


```r
pseq.L2 <- aggregate_taxa(pseq, level = "L2")
```

```
## Error in tax_glom(pseq, level): Bad taxrank argument. Must be among the values of rank_names(physeq)
```

```r
pseq.L1 <- aggregate_taxa(pseq, level = "L1")
```

```
## Error in tax_glom(pseq, level): Bad taxrank argument. Must be among the values of rank_names(physeq)
```

Importing HITChip probe-level data and taxonomy from HITChip
output directory (these are not available in the phyloseq object):


```r
probedata <- read_hitchip(data.directory, method = "frpa")$probedata
```

```
## Error in eval(expr, envir, enclos): could not find function "read_hitchip"
```

```r
taxonomy.full <- read_hitchip(data.directory, method = "frpa")$taxonomy.full
```

```
## Error in eval(expr, envir, enclos): could not find function "read_hitchip"
```

Convert your own data into phyloseq format as follows:


```r
# We need to choose the HITChip data level to be used in the analyses
# In this example use HITChip L2 data (note: this is in absolute scale)
res <- read_hitchip(method = "frpa", data.dir = data.directory)

# Species-level data matrix
otu <- otu_table(res$pseq)@.Data 

# Corresponding sample metadata
meta <- res$meta

# Taxonomy
taxonomy <- GetPhylogeny("HITChip", "filtered")
taxonomy <- unique(as.data.frame(taxonomy[, c("L1", "L2", "species")]))
rownames(taxonomy) <- as.vector(taxonomy[, "species"])

# Merging data matrices into phyloseq format:
pseq <- hitchip2physeq(t(otu), meta, taxonomy)
```


