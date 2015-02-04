<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{Project Template}
  %\usepackage[utf8]{inputenc}
-->


Minimal example
===============

To test this example, do the following:

1.  Copy the
    [Template.Rmd](https://raw.githubusercontent.com/microbiome/microbiome/master/vignettes/Template.Rmd)
    (this file) on your computer. NOTE: change the file ending from .txt
    to .Rmd !!
2.  Start [RStudio](http://www.rstudio.com/)
3.  Open the
    [Template.Rmd](https://raw.githubusercontent.com/microbiome/microbiome/master/vignettes/Template.Rmd)
    file in RStudio
4.  Convert the Rmd file with the 'knit HTML' button
5.  Start expanding this file to make your own report
6.  [Examples to use unique microbiota profiling data
    set](https://github.com/microbiome/microbiome/blob/master/vignettes/Atlas.Rmd)
7.  Adapt further examples from [microbiome
    tutorial](https://github.com/microbiome/microbiome/blob/master/vignettes/vignette.md)

### Update the microbiome package

    library(devtools)
    install_github("microbiome/microbiome")

    ## Downloading github repo microbiome/microbiome@master
    ## Installing microbiome
    ## '/usr/lib/R/bin/R' --vanilla CMD INSTALL  \
    ##   '/tmp/Rtmpq61jQN/devtools5f121d293985/microbiome-microbiome-c5feb51'  \
    ##   --library='/home/antagomir/R/x86_64-pc-linux-gnu-library/3.1'  \
    ##   --install-tests 
    ## 
    ## Reloading installed microbiome
    ## 
    ## microbiome R package (microbiome.github.com)
    ##           
    ## 
    ## 
    ##  Copyright (C) 2011-2015
    ##           Leo Lahti and Jarkko Salojarvi 
    ## 
    ##         
    ##           <microbiome-admin@googlegroups.com>
    ## 
    ## 
    ## Attaching package: 'microbiome'
    ## 
    ## The following object is masked from 'package:lattice':
    ## 
    ##     densityplot
    ## 
    ## The following object is masked from 'package:e1071':
    ## 
    ##     impute

### Load the microbiome package

    library(microbiome)
