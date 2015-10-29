### Probe summarization

Summarize (preprocessed) oligo-level data into phylotype level; examples with simulated data; see [read_hitchip](reading) to use your own data. We use and recommend the [Robust Probabilistic Averaging (RPA)](https://github.com/antagomir/RPA/wiki) for probe summarization.



```r
library(microbiome)
data.directory <- system.file("extdata", package = "microbiome")

# Read oligo-level data (here: simulated example data)
probedata <- read_hitchip(data.directory, method = "frpa")$probedata
```

```
## Loading pre-calculated RPA preprocessing parameters
```

```r
# Read phylogeny map
# NOTE: use phylogeny.filtered for species/L1/L2 summarization
# Load taxonomy from output directory
taxonomy <- GetPhylogeny("HITChip", "filtered")
```

```
## Error in eval(expr, envir, enclos): could not find function "GetPhylogeny"
```

```r
# Summarize oligos into higher level phylotypes
dat <- summarize_probedata(
                 probedata = probedata,
		 taxonomy = taxonomy, 
                 method = "rpa",
		 level = "species")
```

```
## Error in summarize_probedata(probedata = probedata, taxonomy = taxonomy, : object 'taxonomy' not found
```


### Retrieve probe-level data

Get probes for each probeset:


```r
sets <- retrieve.probesets(taxonomy, level = "species", name = NULL)
```

```
## Error in as.data.frame(tax.table): object 'taxonomy' not found
```

```r
head(sets)
```

```
## Error in head(sets): error in evaluating the argument 'x' in selecting a method for function 'head': Error: object 'sets' not found
```


Get probeset data matrix/matrices:


```r
set <- get.probeset("Actinomyces naeslundii", "species",
       		     taxonomy, probedata, log10 = TRUE)
```

```
## Error in as.data.frame(taxonomy): object 'taxonomy' not found
```

```r
kable(set[, 1:5])
```

```
## Error in is.data.frame(x): object 'set' not found
```





