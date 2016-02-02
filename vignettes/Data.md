## Example data sets

This page shows how to import microbiome profiling data into phyloseq format and to load some published example data sets in R. For further microbiome data sets in phyloseq format, check [this](http://joey711.github.io/phyloseq/download-microbio.me.html). For data preprocessing (filtering, subsetting etc.), see the [preprocessing tutorial](Preprocessing.md).


### HITChip Atlas data 

Data from [Lahti et al. Nat. Comm. 5:4344, 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html) contains large-scale profiling of 130 genus-like taxa across 1006 normal western adults. Some subjects have also short time series. This data set is available in [Data Dryad](http://doi.org/10.5061/dryad.pk75d). [Download the HITChip Atlas in R phyloseq format](Atlas.md):


```r
library(microbiome)
atlas1006 <- download_microbiome("atlas1006")
```

```
## Downloading data set from Lahti et al. Nat. Comm. 5:4344, 2014 from Data Dryad: http://doi.org/10.5061/dryad.pk75d
```

```
## Warning in file(file, "rt"): cannot open file '/home/antagomir/R/x86_64-pc-
## linux-gnu-library/3.2/microbiome/extdata/taxonomy.full.tab': No such file
## or directory
```

```
## Error in file(file, "rt"): cannot open the connection
```

Also available as a ready-made example data set in microbiome package:


```r
data(atlas1006)
```

### Diet swap data set

Data from [O'Keefe et al. Nat. Comm. 6:6342, 2015](http://dx.doi.org/10.1038/ncomms7342) from [Data Dryad](http://dx.doi.org/10.5061/dryad.1mn1n). This is a two-week diet swap study between western (USA) and traditional (rural Africa) diets including microbiota profiling.


```r
library(microbiome)
dietswap <- download_microbiome("dietswap")
```

```
## Downloading data set from O'Keefe et al. Nat. Comm. 6:6342, 2015 from Data Dryad: http://datadryad.org/resource/doi:10.5061/dryad.1mn1n
```

```
## Warning in file(file, "rt"): cannot open file '/home/antagomir/R/x86_64-pc-
## linux-gnu-library/3.2/microbiome/extdata/taxonomy.full.tab': No such file
## or directory
```

```
## Error in file(file, "rt"): cannot open the connection
```

Also available as a ready-made example data set in microbiome package:


```r
data(dietswap)
```

### Intestinal microbiota and blood serum lipid metabolites

Data from [Lahti et al. PeerJ 1:e32, 2013](https://peerj.com/articles/32/) characterizes associations between human intestinal microbiota and blood serum lipids. This data set is not readily provided in phyloseq format since it also contains additional data matrix of lipid species. Loading the data in R:


```r
library(microbiome)
data.peerj32 <- download_microbiome("peerj32")
```

```
## Downloading data set from Lahti et al. PeerJ, 2013: https://peerj.com/articles/32/
```

```
## Warning in file(file, "rt"): cannot open file '/home/antagomir/R/x86_64-pc-
## linux-gnu-library/3.2/microbiome/extdata/taxonomy.full.tab': No such file
## or directory
```

```
## Error in file(file, "rt"): cannot open the connection
```


Also available as a ready-made example data set in microbiome package:


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
## [1] "/home/antagomir/R/x86_64-pc-linux-gnu-library/3.2/microbiome/extdata"
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
## Warning in file(file, "rt"): cannot open file '/home/antagomir/R/x86_64-pc-
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
## Error in inherits(x, get.component.classes()): object 'pseq' not found
```

```r
pseq.L1 <- aggregate_taxa(pseq, level = "L1")
```

```
## Error in inherits(x, get.component.classes()): object 'pseq' not found
```

Importing HITChip probe-level data and taxonomy from HITChip
output directory (these are not available in the phyloseq object):


```r
probedata <- read_hitchip(data.directory, method = "frpa")$probedata
```

```
## Warning in file(file, "rt"): cannot open file '/home/antagomir/R/x86_64-pc-
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
## Warning in file(file, "rt"): cannot open file '/home/antagomir/R/x86_64-pc-
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

```
## Error in sample_data(pseq): error in evaluating the argument 'object' in selecting a method for function 'sample_data': Error: object 'pseq' not found
```

Taxonomy table


```r
tax.table <- tax_table(pseq)
```

```
## Error in tax_table(pseq): error in evaluating the argument 'object' in selecting a method for function 'tax_table': Error: object 'pseq' not found
```

Pick taxa abundance data matrix. In this example the OTU level corresponds to genus-like groups (the function name otu_table is somewhat misleading):


```r
otu <- otu_table(pseq)@.Data
```

```
## Error in otu_table(pseq): error in evaluating the argument 'object' in selecting a method for function 'otu_table': Error: object 'pseq' not found
```

Aggregate the abundance matrix to higher-level taxa on HITChip:


```r
pseq2 <- aggregate_taxa(pseq, "Phylum") # Aggregate into phyloseq object
```

```
## Error in inherits(x, get.component.classes()): object 'pseq' not found
```

```r
dat <- otu_table(pseq2)@.Data # Pick aggregated abundance table
```

```
## Error in otu_table(pseq2): error in evaluating the argument 'object' in selecting a method for function 'otu_table': Error: object 'pseq2' not found
```

Melted data for plotting is readily obtained with the phyloseq psmelt function:


```r
df <- psmelt(pseq)
```

```
## Error in inherits(physeq, "phyloseq"): object 'pseq' not found
```

```r
kable(head(df))
```



|                                         |
|:----------------------------------------|
|function (x, df1, df2, ncp, log = FALSE) |
|{                                        |
|if (missing(ncp))                        |
|.Call(C_df, x, df1, df2, log)            |
|else .Call(C_dnf, x, df1, df2, ncp, log) |
|}                                        |
