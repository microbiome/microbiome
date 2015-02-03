### Matrix visualization

Matrix heatmaps:


```r
library(microbiome)

data.directory <- system.file("extdata", package = "microbiome")

genus.matrix.log10.simulated <- read.profiling(level = "L2", method = "frpa", 
			              data.dir = data.directory, log10 = TRUE)  
```

```
## Reading /home/antagomir/R/x86_64-pc-linux-gnu-library/3.1/microbiome/extdata/L2-frpa.tab
## Logarithmizing the data
```

```r
tmp <- netresponse::plot_matrix(genus.matrix.log10.simulated, type = "twoway")
```

![plot of chunk matvisu-example](figure/matvisu-example-1.png) 


Project high-dimensional data on two-dimensional plane by various methods 
including PCA, MDS, Sammons mapping etc. (for visualization purposes):

* project.data

