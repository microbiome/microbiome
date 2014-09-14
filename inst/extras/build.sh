/usr/bin/R CMD BATCH document.R
/usr/bin/R CMD build ../../
/usr/bin/R CMD check microbiome_0.99.34.tar.gz
/usr/bin/R CMD INSTALL microbiome_0.99.34.tar.gz
<<<<<<< HEAD

/usr/bin/R CMD BiocCheck microbiome_0.99.34.tar.gz
=======
>>>>>>> 413f7b12461595839da811dd5bd1ea945ce9d764
