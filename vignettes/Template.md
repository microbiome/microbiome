<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{Project Template}
  %\usepackage[utf8]{inputenc}
-->


Minimal example
===============

To test this example, do the following:

1.  Copy the
    [Template.Rmd](https://github.com/microbiome/microbiome/tree/master/vignettes/Template.Rmd)
    (this file) on your computer
2.  Start [RStudio](http://www.rstudio.com/)
3.  open the Template.Rmd file in RStudio
4.  Convert the Rmd file with the knit button
5.  Start expanding this file to make your own report. You can check
    examples from the [example
    workflow](https://github.com/microbiome/microbiome/blob/master/vignettes/Atlas.Rmd)
    and [microbiome
    tutorial](https://github.com/microbiome/microbiome/blob/master/vignettes/vignette.md)

<!-- -->

    library(devtools)
    install_github("microbiome/microbiome")

    ## Downloading github repo microbiome/microbiome@master
    ## Installing microbiome
    ## '/usr/lib/R/bin/R' --vanilla CMD INSTALL  \
    ##   '/tmp/Rtmpq61jQN/devtools5f124c8f0655/microbiome-microbiome-69dfa32'  \
    ##   --library='/home/antagomir/R/x86_64-pc-linux-gnu-library/3.1'  \
    ##   --install-tests

    plot(c(1,2,3))

![](Template_files/figure-markdown_strict/install2-1.png)
