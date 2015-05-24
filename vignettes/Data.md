## Example data sets

This page shows how to import HITChip data in R, how to convert HITChip data into phyloseq format, and how to load some published example data sets for microbiome analyses in R. For further microbiome data sets in phyloseq format, check [this](http://joey711.github.io/phyloseq/download-microbio.me.html).

For examples on preprocessing the data (filtering, subsetting etc.), see the [preprocessing tutorial](Preprocessing.md).


### HITChip Atlas data 

Data from [Lahti et al. Nat. Comm. 5:4344, 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html) contains large-scale profiling of 130 genus-like taxa across 1006 normal western adults. Some subjects have also short time series. This data set is available in [Data Dryad](http://doi.org/10.5061/dryad.pk75d). [Downloading the HITChip Atlas in R phyloseq format](Atlas.md):


```r
library(microbiome)
data.atlas <- download_microbiome("atlas1006")
```

```
## Downloading data set from Lahti et al. Nat. Comm. 5:4344, 2014 from Data Dryad: http://doi.org/10.5061/dryad.pk75d
```


### Diet swap data set

Data from [O'Keefe et al. Nat. Comm. 6:6342, 2015](http://dx.doi.org/10.1038/ncomms7342) from [Data Dryad](http://dx.doi.org/10.5061/dryad.1mn1n). This is a two-week diet swap study between western (USA) and traditional (rural Africa) diets including microbiota profiling.


```r
library(microbiome)
data.dietswap <- download_microbiome("dietswap")
```

```
## Downloading data set from O'Keefe et al. Nat. Comm. 6:6342, 2015 from Data Dryad: http://datadryad.org/resource/doi:10.5061/dryad.1mn1n
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


### Importing HITChip data

Importing HITChip data from data folder. With HITChip,
[fRPA](http://www.computer.org/csdl/trans/tb/2011/01/ttb2011010217-abs.html)
is the recommended preprocessing method. You can provide sample
metadata by adding new fields in the template metadata file your
HITChip data folder and exporting it again to tab-separated .tab
format. Some standard, self-explanatory field names include 'sample',
'time', 'subject', 'group', 'gender', 'diet', 'age'. You can leave
these out or include further fields. See this [example
file](https://raw.github.com/microbiome/microbiome/master/inst/extdata/meta.tab).



```r
# Define example data path (replace here data.directory with your own path)
library(microbiome)
data.directory <- system.file("extdata", package = "microbiome")
print(data.directory)
```

```
## [1] "/home/lei/R/x86_64-unknown-linux-gnu-library/3.2/microbiome/extdata"
```

Read HITChip data from the specified folder and probe summarization
method (returns the precalculated data matrices for all levels, sample
metadata and taxonomy; without detection thresholding):


```r
# Read precalculated HITChip data and check available data entries
chipdata <- read.profiling(method = "frpa", data.dir = data.directory)
print(names(chipdata))
```

```
## [1] "probedata"     "L1"            "L2"            "species"      
## [5] "taxonomy"      "taxonomy.full" "meta"
```

```r
# Pick specific data field (for instance oligo level or L2-level data)
probedata <- chipdata[["probedata"]]
dat <- chipdata[["L2"]]

# Check the output
kable(head(dat))
```



|                             |   Sample.1|    Sample.2|    Sample.3|   Sample.4|   Sample.5|   Sample.6|   Sample.7|   Sample.8|    Sample.9|   Sample.10|  Sample.11|  Sample.12| Sample.13| Sample.14| Sample.15| Sample.16|  Sample.17| Sample.18|  Sample.19|  Sample.20|
|:----------------------------|----------:|-----------:|-----------:|----------:|----------:|----------:|----------:|----------:|-----------:|-----------:|----------:|----------:|---------:|---------:|---------:|---------:|----------:|---------:|----------:|----------:|
|Actinomycetaceae             |   84.76163|    78.29504|   134.13189|  109.30349|  102.27557|   91.90667|  123.30067|   76.80188|    92.14881|   100.80771|  105.48939|   75.44185|  80.54118|  79.29750| 385.23378|  93.70611|  145.87980|  84.28841|   74.86562|   99.49987|
|Aerococcus                   |   39.46335|    37.15725|    57.04694|   49.92221|   50.12097|   45.20561|   39.18982|   40.00816|    48.15633|    48.88911|   47.43598|   40.01769|  38.47069|  51.66101|  57.82417|  59.01372|   74.22681|  40.46727|   38.66126|   49.77266|
|Aeromonas                    |   50.05587|    39.16444|    65.74415|   61.35530|   54.83985|   44.42839|   46.46214|   37.13851|    48.41742|    49.84928|   49.77854|   38.50172|  54.63650|  36.81774|  59.30448|  57.19805|   76.79569|  41.01024|   39.13462|   57.59464|
|Akkermansia                  | 3241.64899| 16118.25264|  3545.40458| 2695.79856| 1373.24405| 3092.14005| 6116.98477| 5479.95341|  2785.46574|  2980.10595| 1997.84190| 2480.54241| 723.33085| 894.08274| 681.38990| 559.76203| 1457.17800| 850.12411| 2778.85817| 2479.55738|
|Alcaligenes faecalis et rel. |  193.64745|   274.46270|   242.94563|  185.10122|  188.41686|  184.45777|  147.32411|  144.44846|   257.87644|   179.00981|  496.61390|  153.35142| 149.67784| 149.41375| 177.97378| 165.68556|  294.35603| 162.99679|  152.82539|  181.06110|
|Allistipes et rel.           | 5275.98593|  2885.60708| 26293.81603| 2176.93497| 3004.86107|  679.42221| 1355.59886|  609.17888| 14413.70966| 14237.53229| 7305.26704| 4402.88386| 579.25640| 665.34332| 740.46958| 973.90064| 4646.47323| 906.95959| 4966.18710| 5052.53186|


## HITChip to phyloseq format


The [phyloseq](https://github.com/joey711/phyloseq) R package provides
many additional tools for microbiome analyses. See [phyloseq demo
page](http://joey711.github.io/phyloseq-demo/).

Import HITChip phylotype-level data in
[phyloseq](https://github.com/joey711/phyloseq) format (note: the
precalculated matrices are calculated with detection.threshold = 0):


```r
pseq <- read_hitchip(data.directory, method = "frpa", detection.threshold = 10^1.8)$pseq
```

```
## Reading Chip data from /home/lei/R/x86_64-unknown-linux-gnu-library/3.2/microbiome/extdata
## Loading pre-calculated RPA preprocessing parameters
## Reading Chip data from /home/lei/R/x86_64-unknown-linux-gnu-library/3.2/microbiome/extdata
```

Get higher taxonomic levels, use (on HITChip we use L1/L2 instead of Phylum/Genus):


```r
pseq.L2 <- aggregate_taxa(pseq, level = "L2")
pseq.L1 <- aggregate_taxa(pseq, level = "L1")
```

Importing HITChip probe-level data and taxonomy from HITChip
output directory (these are not available in the phyloseq object):


```r
probedata <- read_hitchip(data.directory, method = "frpa", detection.threshold = 10^1.8)$probedata
```

```
## Reading Chip data from /home/lei/R/x86_64-unknown-linux-gnu-library/3.2/microbiome/extdata
## Loading pre-calculated RPA preprocessing parameters
## Reading Chip data from /home/lei/R/x86_64-unknown-linux-gnu-library/3.2/microbiome/extdata
```

```r
taxonomy.full <- read_hitchip(data.directory, method = "frpa", detection.threshold = 10^1.8)$taxonomy.full
```

```
## Reading Chip data from /home/lei/R/x86_64-unknown-linux-gnu-library/3.2/microbiome/extdata
## Loading pre-calculated RPA preprocessing parameters
## Reading Chip data from /home/lei/R/x86_64-unknown-linux-gnu-library/3.2/microbiome/extdata
```


Convert your own data matrices into phyloseq format as follows:


```r
# We need to choose the HITChip data level to be used in the analyses
# In this example use HITChip L2 data (note: this is in absolute scale)
otu <- read.profiling(method = "frpa", data.dir = data.directory)$L2
meta <- read.profiling(method = "frpa", data.dir = data.directory)$meta
taxonomy <- GetPhylogeny("HITChip", "filtered")
taxonomy <- unique(as.data.frame(taxonomy[, c("L1", "L2")]))
rownames(taxonomy) <- as.vector(taxonomy[, "L2"])

# Convert to phyloseq
pseq <- hitchip2physeq(t(otu), meta, taxonomy, detection.limit = 10^1.8)
```


### Picking data from phyloseq  

Assuming your data is in the phyloseq format, many standard tools can directly operate on that data. If you need to pick specific data sets separately, you can mimic these examples:


```r
library(microbiome)

# Get some HITChip data in phyloseq format
pseq <- download_microbiome("atlas1006")
```

```
## Downloading data set from Lahti et al. Nat. Comm. 5:4344, 2014 from Data Dryad: http://doi.org/10.5061/dryad.pk75d
```

```r
# Sample metadata
meta <- sample_data(pseq)

# Taxonomy table
tax.table <- tax_table(pseq)

# OTU data (the zero point has been moved to the detection threshold;
# typically signal 1.8 at HITChip log10 scale). In this example
# the OTU level corresponds to genus-like groups
otu <- otu_table(pseq)@.Data

# Higher-level taxa on HITChip
pseq2 <- aggregate_taxa(pseq, "Phylum")
dat <- otu_table(pseq2)@.Data

# Alternatively, you can use phyloseq functions for this:
level <- "Phylum"
tg <- tax_glom(pseq, level) # Agglomerate taxa
x <- tg@otu_table # Pick the agglomerated data
# On HITChip, we are missing the taxonomic tree and 
# need to use the following to provide correct names for the
# agglomerated taxa:
rownames(x) <- as.character(as.data.frame(tax_table(tg))[[level]]) 
```

### Subsetting phyloseq data


```r
bacteroidetes <- levelmap(NULL, "Phylum", "Genus", tax_table(pseq))$Bacteroidetes

# Keep only the given taxa 
pseq.subset <- prune_taxa(bacteroidetes, pseq)

# Pick samples by specific metadata fields
pseq.subset2 <- subset_samples(pseq.subset, nationality == "US")
```


### Estimating relative abundancies

Estimate relative abundance of the taxa in each sample. Note: the
input data set needs to be in absolute scale (not logarithmic).


```r
pseq <- download_microbiome("dietswap")
```

```
## Downloading data set from O'Keefe et al. Nat. Comm. 6:6342, 2015 from Data Dryad: http://datadryad.org/resource/doi:10.5061/dryad.1mn1n
```

```r
dat <- otu_table(pseq)@.Data
rel <- relative.abundance(dat, det.th = 0)
```

Calculating relative abundances for phyloseq objects:


```r
pseq <- transform_sample_counts(pseq, function(x) x/sum(x))
```



<!--
### Long-term follow-up time series (David et al. 2014)

The data set from [David et al. Genome Biology 2014, 15:R89](http://genomebiology.com/2014/15/7/R89):


```r
library(microbiome)
data.david2014 <- download_microbiome("david2014")
```
-->

