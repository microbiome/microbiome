## Install HITChipDB package 

The HITChipDB package contains additional routines to fetch and
preprocess HITChip (or MIT/PITChip) data from the MySQL database. Note
that this package is **not needed by most users** and the data is
protected by password/IP combinations. Ask details from
admins. Install the package in R with:



```r
# Install additional dependencies
source("http://www.bioconductor.org/biocLite.R")
biocLite("DBI")
biocLite("RPA")
biocLite("svDialogs")

library(devtools) # Load the devtools package
install_github("microbiome/HITChipDB") # Install the package
# Also install RMySQL, multicore and tcltk !
source("http://www.bioconductor.org/biocLite.R")
biocLite("RMySQL") # multicore, tcltk?
# Test installation by loading the microbiome package in R
library("HITChipDB")
```

