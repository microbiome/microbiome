# https://support.rstudio.com/hc/en-us/articles/200626518-Customizing-Package-Bubuiild-Options
#/usr/bin/R CMD BATCH document.R
/usr/bin/R CMD build ../../ --resave-data #--no-examples  --no-build-vignettes 
/usr/bin/R CMD check microbiome_1.9.16.tar.gz #--no-build-vignettes --no-examples
/usr/bin/R CMD BiocCheck microbiome_1.9.16.tar.gz
/usr/bin/R CMD INSTALL microbiome_1.9.16.tar.gz 
