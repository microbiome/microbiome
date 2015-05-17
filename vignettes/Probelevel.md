### Probe summarization

Summarize (preprocessed) oligo-level data into phylotype level; examples with simulated data; see [read.profiling](reading) to use your own data.


```r
library(microbiome)
data.directory <- system.file("extdata", package = "microbiome")

# Read oligo-level data (here: simulated example data)
oligo.matrix.nolog.simulated <- read.profiling(level = "oligo", 
                                   data.dir = data.directory, log10 = FALSE)
```

```
## Error in read.profiling(level = "oligo", data.dir = data.directory, log10 = FALSE): unused arguments (level = "oligo", log10 = FALSE)
```

```r
# Read phylogeny map
# NOTE: use phylogeny.filtered for species/L1/L2 summarization
# Load phylogeny.info from output directory
phylogeny.info <- GetPhylogeny("HITChip", "filtered")

# Summarize oligos into higher level phylotypes
dat <- summarize.probesets(phylogeny.info, 
                 log10(oligo.matrix.nolog.simulated), 
                 "frpa", "species", verbose = TRUE, 
                  phylotype.rm.list("HITChip"))
```

```
## Summarizing through species...
```

```
## Error in phylogeny.info$oligoID: $ operator is invalid for atomic vectors
```


### Retrieve probe-level data

Get oligos for each probeset:


```r
sets <- retrieve.probesets(phylogeny.info, level = "species", name = NULL)
```

Get probeset data matrix/matrices:


```r
set <- get.probeset("Actinomyces naeslundii", "species", 
             phylogeny.info, oligo.matrix.nolog.simulated, log10 = TRUE)
```

```
## Error in rownames(probedata): object 'oligo.matrix.nolog.simulated' not found
```
