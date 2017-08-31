# https://support.rstudio.com/hc/en-us/articles/200526518-Customizing-Package-Build-Options
#/usr/bin/R CMD BATCH document.R
/usr/bin/R CMD build ../../ #--no-examples  --no-build-vignettes 
/usr/bin/R CMD check microbiome_0.99.63.tar.gz #--no-build-vignettes --no-examples
/usr/bin/R CMD BiocCheck microbiome_0.99.63.tar.gz
/usr/bin/R CMD INSTALL microbiome_0.99.63.tar.gz 

