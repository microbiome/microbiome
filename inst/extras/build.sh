/usr/bin/R CMD BATCH document.R
/usr/bin/R CMD build ../../
/usr/bin/R CMD check microbiome_0.99.34.tar.gz
/usr/bin/R CMD INSTALL microbiome_0.99.34.tar.gz
