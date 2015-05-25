### Probe summarization

Summarize (preprocessed) oligo-level data into phylotype level; examples with simulated data; see [read.profiling](reading) to use your own data. We use and recommend the [Robust Probabilistic Averaging (RPA)](https://github.com/antagomir/RPA/wiki) for probe summarization.



```r
library(microbiome)
data.directory <- system.file("extdata", package = "microbiome")

# Read oligo-level data (here: simulated example data)
probedata <- read.profiling(data.directory, method = "frpa")$probedata
```

```
## Reading Chip data from /home/lei/R/x86_64-unknown-linux-gnu-library/3.2/microbiome/extdata
```

```r
# Read phylogeny map
# NOTE: use phylogeny.filtered for species/L1/L2 summarization
# Load taxonomy from output directory
taxonomy <- GetPhylogeny("HITChip", "filtered")

# Summarize oligos into higher level phylotypes
dat <- summarize_probedata(
                 probedata = probedata,
		 taxonomy = taxonomy, 
                 method = "rpa",
		 level = "species")
```


### Retrieve probe-level data

Get probes for each probeset:


```r
sets <- retrieve.probesets(phylogeny.info, level = "species", name = NULL)
```


Get probeset data matrix/matrices:


```r
set <- get.probeset("Actinomyces naeslundii", "species",
       		    taxonomy, probedata, log10 = TRUE)
```





