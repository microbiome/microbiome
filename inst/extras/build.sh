# https://support.rstudio.com/hc/en-us/articles/200486518-Customizing-Package-Build-Options
/usr/bin/R CMD BATCH document.R
/usr/bin/R CMD build ../../ --no-build-vignettes --no-examples
/usr/bin/R CMD check microbiome_0.99.73.tar.gz --no-build-vignettes --no-examples
/usr/bin/R CMD INSTALL microbiome_0.99.73.tar.gz 
#/usr/bin/R CMD BiocCheck microbiome_0.99.32.tar.gz 
