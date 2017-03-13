---
title: "Output"
author: "Leo Lahti"
date: "2017-03-05"
bibliography: 
- bibliography.bib
- references.bib
output: 
  rmarkdown::html_vignette
---
<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - output}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->


### Writing diversity table into file


```r
output.dir <- "./"
write.table(div.table, file = "DiversityTable.tab", sep = "\t")
```

### Save clustering image to a file

Save in PDF:


```r
pdf("myplot.pdf", width = 7, height = 7 * length(hc$order)/20)
plot(hc, hang=-1, main = "Hierarchical clustering")
dev.off()
```

Save in TIFF:


```r
tiff("myplot.tif", width = 480, height = 480 * length(hc$order)/20)
plot(hc, hang=-1, main = "Hierarchical clustering")
dev.off()
```

To save in Microsoft EMF format, try the following. If you find a
way to tune figure width for emf files kindly let the admins know.


```r
plot(hc, hang=-1, main = "Hierarchical clustering")
savePlot("myplot.emf", type = "emf")
dev.off()
```

