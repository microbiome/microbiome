
### PeerJ example data set

An example data set from Lahti et al. [PeerJ 1:e32, 2013](https://peerj.com/articles/32/) concerns associations between human intestinal microbiota and blood serum lipids. Load the data in R:


```r
library(microbiome)
data(peerj32)
names(peerj32)
```

```
## [1] "lipids"   "microbes" "meta"
```

### Load example data

Load simulated example data of the human gut microbiota. With HITChip,
[fRPA](http://www.computer.org/csdl/trans/tb/2011/01/ttb2011010217-abs.html)
is the recommended preprocessing method.


```r
# Define data path (you can replace data.directory with your own path)
data.directory <- system.file("extdata", package = "microbiome")
print(data.directory)
```

```
## [1] "/home/antagomir/R/x86_64-pc-linux-gnu-library/3.1/microbiome/extdata"
```

```r
# Read HITChip data matrix (genus-level (L2) log10 values)
level <- "L2"
method <- "frpa"
genus.data <- read.profiling(level = level, 
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


### Read metadata

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
rel <- relative.abundance(oligo.data, det.th = NULL)
```

```
## Warning in relative.abundance(oligo.data, det.th = NULL): Applying
## detection threshold at 0.8 quantile: 232.026771597465
```

