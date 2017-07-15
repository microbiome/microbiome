<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - stability}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->
Experimental functions
----------------------

### Time series for individual subjects

    source(system.file("extdata/plot_longitudinal.R", package = "microbiome"))
    p <- plot_longitudinal(pseq, "Dialister", subject = "831", tipping.point = 0.5)
    print(p)
