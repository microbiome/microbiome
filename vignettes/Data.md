## Example data sets

Example data sets for microbiome analyses.


### HITChip Atlas data 


The data from [Lahti et al. Nat. Comm. 5:4344, 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html) contains large-scale profiling of 130 genus-like taxa across 1006 normal western subjects. Some subjects have also short time series. This data set is available in [Data Dryad](http://doi.org/10.5061/dryad.pk75d). [Download the HITChip Atlas microbiome profiling data in R](Atlas.md):


```r
library(microbiome)
library(rdryad)
data.atlas <- download_microbiome("atlas1006")
```

```
## Downloading data set from Lahti et al. Nat. Comm. from Data Dryad: http://doi.org/10.5061/dryad.pk75d
```


### Diet swap data set

An example data set from [O'Keefe et al. Nat. Comm. 6:6342, 2015](http://dx.doi.org/10.1038/ncomms7342) from [Data Dryad](http://dx.doi.org/10.5061/dryad.1mn1n). This is a follow-up of microbial profiles during a 2-week diet swap of western (USA) and traditional (rural Africa) diets.


```r
library(microbiome)
library(rdryad)
data.dietswap <- download_microbiome("dietswap")
```

```
## Downloading data set from O'Keefe et al. Nat. Comm. 6:6342, 2015 from Data Dryad: http://datadryad.org/resource/doi:10.5061/dryad.1mn1n
```

```
## Warning in harmonize_fields(meta): bmi information is a factor; renaming as
## bmi_group
```


### Intestinal microbiota and blood serum lipid metabolites

An example data set from [Lahti et al. PeerJ 1:e32, 2013](https://peerj.com/articles/32/) characterizes associations between human intestinal microbiota and blood serum lipids. Loading the data in R:


```r
library(microbiome)
data.peerj32 <- download_microbiome("peerj32")
```

```
## Downloading data set from Lahti et al. PeerJ, 2013: https://peerj.com/articles/32/
```


### Loadomg example data

Load simulated example data of the human gut microbiota. With HITChip,
[fRPA](http://www.computer.org/csdl/trans/tb/2011/01/ttb2011010217-abs.html)
is the recommended preprocessing method.


```r
library(microbiome)
# Define data path (you can replace data.directory with your own path)
data.directory <- system.file("extdata", package = "microbiome")
print(data.directory)
```

```
## [1] "/home/lei/R/x86_64-unknown-linux-gnu-library/3.2/microbiome/extdata"
```

```r
# Read HITChip data matrix (genus-level (L2) log10 values)
level <- "L1"
method <- "frpa"
l1.data <- read.profiling(level = level, 
	     		       method = method, 
              		       data.dir = data.directory, 
	      	       	       log10 = TRUE)  

# Read HITChip probe level data (absolute values - no log10)
oligo.data <- read.profiling(level = "oligo", 
                             data.dir = data.directory, 
			     log10 = FALSE)  

# Probe-taxon mapping table
phylogeny.info <- read.profiling(level = "phylogeny.full", 
                           	 data.dir = data.directory)

# Phylogeny that is used to summarize the probes to phylotype/genus/phylum levels
phylogeny.info.filtered <- read.profiling(level = "phylogeny.filtered", 
                           	 data.dir = data.directory)
```


### Reading metadata

An easy way to provide sample metadata is to create a tab-separated metadata file. You can create the file in Excel and export it to tab-separated csv format. The standard (and self-explanatory) field names include 'sampleID', 'time', 'subjectID', 'group', 'gender', 'diet', 'age'. You can leave these out or include further fields. See this [example file](https://raw.github.com/microbiome/microbiome/master/inst/extdata/metadata.xls). Read the metadata with:


```r
# Read simulated example metadata
library(gdata)
metadata.file <- paste(data.directory, "/metadata.xls", sep = "")
metadata <- read.xls(metadata.file, as.is = TRUE)
rownames(metadata) <- metadata$sampleID
```


### Estimating relative abundancies

Estimate relative abundance of the taxa in each sample. Note: the
input data set needs to be in absolute scale (not logarithmic).


```r
rel <- relative.abundance(oligo.data, det.th = min(na.omit(oligo.data)))
```

