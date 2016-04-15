## Example data sets for microbiome analysis

Use [phyloseq](http://joey711.github.io/phyloseq/import-data) tools to import standard microbiome data formats (mothur, qiime etc.) into an R phyloseq object. Here we provide some published example data sets in R. For further microbiome data sets, check [this](http://joey711.github.io/phyloseq/download-microbio.me.html). For data preprocessing (filtering, subsetting etc.), see the [preprocessing tutorial](Preprocessing.md). 


### HITChip Atlas data 

Data from [Lahti et al. Nat. Comm. 5:4344, 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html) contains large-scale profiling of 130 genus-like taxa across 1006 normal western adults. Some subjects have also short time series. This data set is available in [Data Dryad](http://doi.org/10.5061/dryad.pk75d) and can be downloaded with the download_microbiome function. However [the HITChip Atlas](Atlas.md) data set is also readily available via the microbiome R package in phyloseq format:


```r
library(microbiome)
data("atlas1006") 
pseq <- atlas1006
```

### Diet swap data set

Data from [O'Keefe et al. Nat. Comm. 6:6342, 2015](http://dx.doi.org/10.1038/ncomms7342) from [Data Dryad](http://dx.doi.org/10.5061/dryad.1mn1n). This is a two-week diet swap study between western (USA) and traditional (rural Africa) diets including microbiota profiling. 


```r
data(dietswap) 
```

### Intestinal microbiota and blood serum lipid metabolites

Data from [Lahti et al. PeerJ 1:e32, 2013](https://peerj.com/articles/32/) characterizes associations between human intestinal microbiota and blood serum lipids. This data set is not readily provided in phyloseq format since it also contains additional data matrix of lipid species. Loading the data in R:


```r
data(peerj32)
```


## HITChip to phyloseq 

### Importing HITChip data

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
## Warning in file(file, "rt"): cannot open file '/home/lei/R/x86_64-pc-
## linux-gnu-library/3.2/microbiome/extdata/oligoprofile.tab': No such file or
## directory
```

```
## Error in file(file, "rt"): cannot open the connection
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
## Warning in file(file, "rt"): cannot open file '/home/lei/R/x86_64-pc-
## linux-gnu-library/3.2/microbiome/extdata/oligoprofile.tab': No such file or
## directory
```

```
## Error in file(file, "rt"): cannot open the connection
```

```r
taxonomy.full <- read_hitchip(data.directory, method = "frpa")$taxonomy.full
```

```
## Warning in file(file, "rt"): cannot open file '/home/lei/R/x86_64-pc-
## linux-gnu-library/3.2/microbiome/extdata/oligoprofile.tab': No such file or
## directory
```

```
## Error in file(file, "rt"): cannot open the connection
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


### Picking data from phyloseq  

Assuming your data 'pseq' is in the phyloseq format, many standard tools can directly operate on that data. If you need to pick specific data sets separately, you can mimic these examples.


Sample metadata:


```r
library(phyloseq)
meta <- sample_data(pseq)
```

Taxonomy table


```r
tax.table <- tax_table(pseq)
```

Pick taxa abundance data matrix. In this example the OTU level corresponds to genus-like groups (the function name otu_table is somewhat misleading):


```r
otu <- otu_table(pseq)@.Data
```

Aggregate the abundance matrix to higher-level taxa on HITChip:


```r
pseq2 <- aggregate_taxa(pseq, "Phylum") # Aggregate into phyloseq object
dat <- otu_table(pseq2)@.Data # Pick aggregated abundance table
```

Melted data for plotting is readily obtained with the phyloseq psmelt function:


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
