---
title: "ANCOM test for differential abundance"
bibliography: 
- bibliography.bib
- references.bib
output: 
  prettydoc::html_pretty:
    theme: cayman
    highlight: github
---
<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - comparisons}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->



## ANCOM

[ANCOM](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4450248/) has been claimed to outperform zero-inflated Gaussians and other recently popular models of differential abundance in microbiome studies. An R package is [available](https://www.niehs.nih.gov/research/resources/software/biostatistics/ancom/index.cfm) but not from standard repositories. Meanwhile, the ANCOM implementation is here modified from (Weiss, Xu, Peddada, Amir, Bittinger, Gonzalez, Lozupone, Zaneveld, Vázquez-Baeza, Birmingham, Hyde, and Knight, 2017). The R code was obtained from the first author and included with permission in the microbiome package. For reference to the original ANCOM method by (Mandal, Treuren, White, ø, Knight, and Peddada, 2015).

Load example data:


```r
# Load the example data
library(microbiome)
data(dietswap)
pseq <- dietswap
```

Calculate adjusted p-values on OTU abundance for the nationality:


```r
padj <- ancom(pseq, "nationality")
print(names(which(padj < 0.05)))
```


For validation purposes, compare ANCOM and Negative binomial adjusted p-values. Note that the ANCOM p-values in this example just take values 0 (significant) and 1 (non-significant). Note that this is a toy example, demonstrating that the significance estimates are correlated between the two tests. 


```r
library(MASS)
pvs <- c()
for (tax in taxa(pseq)) { 
  # Pick the signal (abundance) for this tax
  sample_data(pseq)$signal <- get_sample(pseq, tax)

  # Negative binomial test with group and gender included
  if (length(unique(sample_data(pseq)$signal)) > 1) {
    res <- glm.nb(signal ~ nationality, data = meta(pseq))

    # Show the results
    pvs[[tax]] <- anova(res)["nationality", "Pr(>Chi)"]
  } else {
    pvs[[tax]] <- 1
  }
}
pvs <- p.adjust(pvs)
pvs <- pvs[names(padj)]
boxplot(log10(pvs)  ~ padj, ylab = "Adj P (neg. binomial log10)", xlab = "Adj P (ANCOM)", main = "ANCOM vs. Negative binomial")
```

### References



