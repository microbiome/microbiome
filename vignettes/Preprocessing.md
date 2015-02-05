# Miscellaneous preprocessing functions

### Retrieve and preprocess data from the MySQL database

[Preprocessing array data ('profiling script')](profiling)  
[Preprocessing HITChip atlas collection](atlas)  

### Probe-level operations 

[Probe-level preprocessing](probelevel)

### Robust Probabilistic Averaging

[RPA](RPA.Rmd)

### Split data into training and test samples


```r
library(microbiome)
data.directory <- system.file("extdata", package = "microbiome")
metadata.example.file <- paste(data.directory, "/metadata.tab", sep = "")
metadata.simulated <- read.csv(metadata.example.file, sep = "\t", as.is = TRUE)
splitted.data <- pick.training.samples(metadata.simulated, 
	      	 			training.fraction = 0.80, rseed = 1463)
```

```
## Error in eval(expr, envir, enclos): could not find function "pick.training.samples"
```

