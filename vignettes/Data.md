## Description of data formats

This page shows how to define data matrices within R. To read data from files, see [this page](reading).


### Sample metadata

If you provide sample metadata in the standard format described below, you will have easier access to many [analysis tools](analysis). Sample metadata can be defined provided as data.frame, as in the simulated example:


```r
# Load the package
library(microbiome)

# Define file path for the simulated data
# (replace this with your own)
data.directory <- system.file("extdata", package = "microbiome")
metadata.example.file <- paste(data.directory, "/metadata.tab", sep = "")

# Read simulated metadata
metadata.simulated <- read.csv(metadata.example.file, sep = "\t", as.is = TRUE)

# Check the first metadata entries
print(head(metadata.simulated)) 
```

```
##   sampleID time   age   bmi   subjectID   group gender      diet
## 1 Sample.1    1 85.84 32.29 subjectID-1 group-3      M Beverages
## 2 Sample.2    2 94.63 23.11 subjectID-2 group-1      M Beverages
## 3 Sample.3    3 99.95 36.19 subjectID-3 group-1      F   Carrots
## 4 Sample.4    4 20.31 31.64 subjectID-4 group-3      M   Carrots
## 5 Sample.5    1 23.24 28.47 subjectID-1 group-2      F  Apricots
## 6 Sample.6    2 58.36 28.22 subjectID-2 group-1      F  Apricots
```


The metadata needs to contain the following fields: sampleID, subjectID, group, and time. Additional fields, such as age, gender, diet, etc. can be included. To create your own metadata, define the data.frame. Here we show an example with the simulated data:



```r
# Define data directory where the files can be found
data.directory <- system.file("extdata", package = "microbiome")  

# Read logarithmized L2-nmf data (phylotypes vs. samples):
genus.matrix.log10.simulated <- read.profiling(level = "L2", method = "frpa", 
			             data.dir = data.directory, log10 = TRUE)  

# Determine sample size
N <- ncol(genus.matrix.log10.simulated)

# Create some random sample information for the metadata
metadata.simulated <- data.frame(list(
              sampleID = colnames(genus.matrix.log10.simulated),
	      subjectID = paste("subjectID", rep(1:4, 5)),
	      group = sample(paste("group", rep(1:4, 5))),
	      time = paste("group", rep(1:4, 5)),
              age = runif(N, 0, 100),
              gender = sample(c("M", "F"), N, replace = TRUE),
              diet = sample(c("Apricots", "Beverages", "Carrots"), 
	      N, replace = TRUE)))
```


We can also add functions to read metadata directly from a text file. Kindly provide an example data set, and we can work it out..

### Generating simulated data

The simulated data for the package is generated into a given output directory as follows:


```r
library(microbiome)
output.dir <- GenerateSimulatedData(output.dir = ".")
```


The simulated data files are available as part of this package. Check the simulated data directory as follows:


```r
print(system.file("extdata", package = "microbiome"))
```

```
<<<<<<< HEAD
## [1] "/home/lei/R/x86_64-pc-linux-gnu-library/3.1/microbiome/extdata"
=======
## [1] "/home/antagomir/R/x86_64-pc-linux-gnu-library/3.1/microbiome/extdata"
>>>>>>> 19a59384b85504870f957ca3ff3e9be41803e8fc
```
