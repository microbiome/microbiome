# https://support.rstudio.com/hc/en-us/articles/200626518-Customizing-Package-Build-Options
#/usr/bin/R CMD BATCH document.R
/usr/bin/R CMD build ../../ #--no-examples  --no-build-vignettes 
/usr/bin/R CMD check microbiome_1.1.2.tar.gz #--no-build-vignettes --no-examples
/usr/bin/R CMD BiocCheck microbiome_1.1.2.tar.gz
/usr/bin/R CMD INSTALL microbiome_1.1.2.tar.gz 
