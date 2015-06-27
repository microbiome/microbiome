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
sets <- retrieve.probesets(taxonomy, level = "species", name = NULL)
head(sets)
```

```
## $`Achromobacter denitrificans`
## [1] "HIT 1497" "HIT 1498" "HIT 1500"
## 
## $`Acidaminococcus fermentans`
## [1] "HIT 1490" "HIT 1491" "HIT 1492" "HIT 1494"
## 
## $`Acinetobacter calcoaceticus`
## [1] "HIT 4156" "HIT 4158"
## 
## $`Actinomyces naeslundii`
## [1] "HIT 1589" "HIT 1590"
## 
## $`Aerococcus viridans`
## [1] "HIT 1072" "HIT 6678" "HIT 955"  "HIT 957"  "HIT 959"  "HIT 960" 
## 
## $`Aeromonas trota`
## [1] "HIT 1633" "HIT 1634" "HIT 1635" "HIT 1636" "HIT 1638"
```


Get probeset data matrix/matrices:


```r
set <- get.probeset("Actinomyces naeslundii", "species",
       		     taxonomy, probedata, log10 = TRUE)
kable(set[, 1:5])
```



|         | Sample.1| Sample.2| Sample.3| Sample.4| Sample.5|
|:--------|--------:|--------:|--------:|--------:|--------:|
|HIT 1589 | 1.606418| 1.569691| 1.805999| 1.666957| 1.670350|
|HIT 1590 | 1.599793| 1.594346| 1.890435| 1.913062| 1.732066|





