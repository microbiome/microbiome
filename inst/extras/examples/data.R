# Load the package
library(microbiome)  

# Define data path (here we retrieve data from R package itself)
data.directory <- system.file("extdata", package = "microbiome")

# Read HITChip data matrix (log10 values)
level <- "L2"
method <- "frpa"
hitchip.data <- read.profiling(level = level, 
	     		       method = method, 
              		       data.dir = data.directory, 
	      	       	       log10 = TRUE)  

# Oligo level data (absolute values - no log10)
oligo.data <- read.profiling(level = "oligo", 
                             data.dir = data.directory, 
			     log10 = FALSE)  

# Oligo-phylogeny mapping table
phylogeny.info <- read.profiling(level = "phylogeny.full", 
                           	 data.dir = data.directory)

# --------------------------------------------------------------


