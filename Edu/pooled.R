---
title: "Pooled tests"
author: "`r Sys.Date()`"
bibliography: 
- bibliography.bib
output:
  rmdformats::readthedown:
    self_contained: true
    thumbnails: true
    lightbox: true
    gallery: true
    use_bookdown: false
    highlight: haddock
---
<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - atlas}
  %\usepackage[utf8]{inputenc}
-->


Example data:


```{r pooled1, warning=FALSE, message=FALSE}
library(microbiome)
data(dietswap)
d <- dietswap

# Pick microbial abundances for a given taxonomic group
taxa <- "Dialister"

# Construct a data.frame with the selected
# taxonomic group and grouping
df <- data.frame(Abundance = abundances(d)[taxa,],
                 Group = meta(d)$nationality,
		 Log10_Abundance = log10(1 + abundances(d)[taxa,])
		 )
```



Motivation: non-gaussianity


Counts

```{r pooled2, warning=FALSE, message=FALSE}
print(abundances(d)[1:5,1:3])
```


Sparsity

```{r pooled3, warning=FALSE, message=FALSE}
print(hist(log10(1 + abundances(d)), 100))

# - huge tail of rare taxa:
# median abundance per taxa
medians <- apply(abundances(d),1,median)/1e3
medians <- data.frame(Rank = rank(-medians),
                      Median = medians,
		      Taxon = names(medians))
p <- ggplot(medians, aes(x = Rank, y = Median)) +
        geom_point() +
	labs(y = "Median Abundance (reads)")
print(p)
```

Overdispersion

```{r pooled_overdispersion, warning=FALSE, message=FALSE}
means <- apply(abundances(d),1,mean)
variances <- apply(abundances(d),1,var)

# Calculate mean and variance over samples for each taxon
library(reshape2)
library(dplyr)
df <- melt(abundances(d))
names(df) <- c("Taxon", "Sample", "Reads")
df <- df %>% group_by(Taxon) %>%
             summarise(mean = mean(Reads),
	               variance = var(Reads))

# Illustrate overdispersion
library(scales)
p <- ggplot(df, aes(x = mean, y = variance)) +
       geom_point() +
       geom_abline(aes(intercept = 0, slope = 1)) +
       scale_x_log10(labels = scales::scientific) +
       scale_y_log10(labels = scales::scientific) +
       labs(title = "Overdispersion (variance > mean)")
print(p)
```


Running DESeq2. Running the DESeq2 analysis. Convert phyloseq object to deseq2 format.

```{r pooled_deseq, warning=FALSE, message=FALSE}
library(DESeq2)
ds2 <- phyloseq_to_deseq2(d, ~ group + nationality)
# Run DESeq2 analysis (all taxa at once!)
dds <- DESeq(ds2)
# Investigate results
res <- results(dds)
deseq.results <- as.data.frame(res)

df <- deseq.results
df$taxon <- rownames(df)
df <- df %>% arrange(log2FoldChange, padj)

# Print the results; fitered and sorted by pvalue and effectsize
library(knitr)
df <- df %>% filter(pvalue < 0.05 & log2FoldChange > 1) %>%
             arrange(pvalue, log2FoldChange)
print(kable(df, digits = 5))
```


Now identify top taxa based on Wilcoxon test.

```{r pooled_topw, warning=FALSE, message=FALSE}
all.taxa <- taxa(d)
pvalue.wilcoxon <- c()
foldchange <- c()
for (taxa in all.taxa) {
  # Create a new data frame for each taxonomic group
  df <- data.frame(Abundance = abundances(d)[taxa,],
                   Log10_Abundance = log10(1 + abundances(d)[taxa,]),
                   Group = meta(d)$nationality)
  # Calculate pvalue and effect size (difference beween log means)		 
  pvalue.wilcoxon[[taxa]] <- wilcox.test(Abundance ~ Group, data = df)$p.value
  foldchange[[taxa]] <- coef(lm(Log10_Abundance ~ Group, data = df))[[2]]
}
# Correct p-values for multiple testing
pvalue.wilcoxon.adjusted <- p.adjust(pvalue.wilcoxon)
```

Also the lowest p-values are well correlated (higher not so).

```{r pooled_pcomp, warning=FALSE, message=FALSE}
par(mfrow = c(1,2))
plot(deseq.results$padj, pvalue.wilcoxon.adjusted,
  xlab = "DESeq2 adjusted p-value",
  ylab = "Wilcoxon adjusted p-value",
  main = "P-value comparison")
abline(v = 0.05, h = 0.05, lty = 2)

plot(deseq.results$log2FoldChange, foldchange, 
  xlab = "DESeq2",
  ylab = "Linear model",
  main = "Effect size comparison")  )
abline(0,1)
```

