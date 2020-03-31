# https://support.rstudio.com/hc/en-us/articles/200626518-Customizing-Package-Build-Options
#~/bin/R-3.5.0/bin/R CMD BATCH document.R

~/bin/R-patched/bin/R CMD build ../../ --resave-data #--no-examples  --no-build-vignettes 
~/bin/R-patched/bin/R CMD check microbiome_1.9.95.tar.gz #--no-build-vignettes --no-examples
~/bin/R-patched/bin/R CMD BiocCheck microbiome_1.9.95.tar.gz
~/bin/R-patched/bin/R CMD INSTALL microbiome_1.9.95.tar.gz 
