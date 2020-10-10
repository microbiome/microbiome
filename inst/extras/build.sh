# https://support.rstudio.com/hc/en-us/articles/200626518-Customizing-Package-Build-Options
#~/bin/R-3.5.0/bin/R CMD BATCH document.R

~/bin/R-4.0.0/bin/R CMD build ../../ --resave-data #--no-examples  --no-build-vignettes 
~/bin/R-4.0.0/bin/R CMD check microbiome_2.1.0.tar.gz #--no-build-vignettes --no-examples
~/bin/R-4.0.0/bin/R CMD BiocCheck microbiome_2.1.0.tar.gz
~/bin/R-4.0.0/bin/R CMD INSTALL microbiome_2.1.0.tar.gz 
