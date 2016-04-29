# https://support.rstudio.com/hc/en-us/articles/200486518-Customizing-Package-Build-Options
~/bin/R-3.2.3/bin/R CMD BATCH document.R
~/bin/R-3.2.3/bin/R CMD build ../../ --no-build-vignettes --no-examples
~/bin/R-3.2.3/bin/R CMD check microbiome_0.99.68.tar.gz --no-build-vignettes --no-examples
~/bin/R-3.2.3/bin/R CMD INSTALL microbiome_0.99.68.tar.gz 
#/usr/bin/R CMD BiocCheck microbiome_0.99.32.tar.gz 
