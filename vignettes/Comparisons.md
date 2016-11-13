---
title: "Comparisons"
author: "Leo Lahti"
date: "2016-11-13"
bibliography: 
- bibliography.bib
- references.bib
output: 
  rmarkdown::html_vignette
---
<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - comparisons}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->


## Group-wise comparisons

### Boxplots


```r
# Load libraries
library(microbiome)
library(ggplot2)
library(dplyr)

# Probiotics intervention example data from
# https://peerj.com/articles/32/
data("peerj32")

# Abundance boxplot
p <- boxplot_abundance(peerj32$phyloseq, x = "time", y = "Akkermansia", line = "subject", color = "gender")
print(p)
```

![plot of chunk boxplot-example](figure/boxplot-example-1.png)







### Negative binomial example

[Read more](http://www.ats.ucla.edu/stat/r/dae/nbreg.htm)


```r
library(MASS)
taxa <- taxa_names(x)[1:2]
x <- atlas1006
df <- as(sample_data(x), "data.frame")
for (tax in taxa) {
  df$signal <- get_sample(x, tax)
  res <- glm.nb(signal ~ bmi_group + gender, data = df)
  print(coef(summary(res)))
}
```


### Comparisons for individual taxa with random effect subject term


```r
# Get taxa x samples abundance matrix
x <- peerj32$phyloseq

# Get the data
mydata <- get_taxa(x)
tax <- "Dialister"
dfs <- sample_data(x)
dfs$signal <- mydata[tax, rownames(dfs)]
dfs$group <- dfs[[group]]

# Paired comparison
library(lme4)
out <- lmer(signal ~ group + (1|subject), data = dfs)
out0 <- lmer(signal ~ (1|subject), data = dfs)
comp <- anova(out0, out)
pv <- comp[["Pr(>Chisq)"]][[2]]
```


## Linear models with limma

Identify most significantly different taxa between males and females.

For further details, see [limma
homepage](http://bioinf.wehi.edu.au/limma/) and [limma User's
guide](http://www.lcg.unam.mx/~lcollado/R/resources/limma-usersguide.pdf). For
discussion on why limma is preferred over t-test, see [this
article](http://www.plosone.org/article/info:doi/10.1371/journal.pone.0012336).


```r
# Get example data
library(microbiome)
data("peerj32")
pseq <- peerj32$phyloseq
otu <- taxa_abundances(transform_phyloseq(pseq, "log10"))
meta <- sample_data(pseq)
grouping.variable <- "gender" 

# Compare the two groups with limma
library(limma)

# Prepare the design matrix which states the groups for each sample
# in the otu
design <- cbind(intercept = 1, Grp2vs1 = meta[[grouping.variable]])
rownames(design) <- rownames(meta)
design <- design[colnames(otu), ]

# NOTE: results and p-values are given for all groupings in the design matrix
# Now focus on the second grouping ie. pairwise comparison
coef.index <- 2
     
# Fit the limma model
fit <- lmFit(otu, design)
fit <- eBayes(fit)

# Limma P-values
pvalues.limma = fit$p.value[, 2]

# Limma effect sizes
efs.limma <-  fit$coefficients[, "Grp2vs1"]

# Summarise 
kable(topTable(fit, coef = coef.index, p.value=0.1), digits = 2)
```



|                               | logFC| AveExpr|     t| P.Value| adj.P.Val|     B|
|:------------------------------|-----:|-------:|-----:|-------:|---------:|-----:|
|Uncultured Clostridiales II    | -0.41|    1.37| -3.72|       0|      0.06| -0.24|
|Eubacterium siraeum et rel.    | -0.34|    1.67| -3.52|       0|      0.06| -0.77|
|Clostridium nexile et rel.     |  0.18|    2.84|  3.41|       0|      0.06| -1.04|
|Sutterella wadsworthia et rel. | -0.33|    1.50| -3.13|       0|      0.10| -1.74|

**Q-Q plot for limma**


```r
qqt(fit$t[, coef.index], df = fit$df.residual + fit$df.prior)
abline(0,1)
```

![plot of chunk limma-qq](figure/limma-qq-1.png)

**Volcano plot for limma**


```r
volcanoplot(fit, coef = coef.index, highlight = coef.index)
```

![plot of chunk limma-volcano](figure/limma-volcano-1.png)


### PERMANOVA

PERMANOVA can be also used to assess community-level differences
between groups. Here let us evaluate whether nationality has a
significant effect on gut microbiota.


```r
# Example data
data("dietswap")
x <- dietswap
group <- "nationality"

# Use relative abundances for simpler visualizations
x <- transform_phyloseq(x, "relative.abundance")
otu <- get_sample(x)
meta <- as(sample_data(x), "data.frame")
meta$group <- meta[[group]]

# PERMANOVA: samples x species as input
library(vegan)
permanova <- adonis(t(otu) ~ group, data=meta, permutations=99, method = "bray")
pv <- as.data.frame(permanova$aov.tab)["group", "Pr(>F)"]

# P-value
print(pv)
```

```
## [1] 0.01
```

```r
# Note the assumption of similar
# multivariate spread among the groups
# ie. analogous to variance homogeneity
# Here the groups have signif. different spreads and
# permanova result may be explained by that.
dist <- vegdist(t(otu))
anova(betadisper(dist,meta$group))
```

```
## Analysis of Variance Table
## 
## Response: Distances
##            Df  Sum Sq  Mean Sq F value    Pr(>F)    
## Groups      1 0.26114 0.261144  18.345 2.754e-05 ***
## Residuals 220 3.13168 0.014235                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# Coefs for the top taxa separating the groups
coef <- coefficients(permanova)["group1",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")
```

![plot of chunk comparisons-permanova](figure/comparisons-permanova-1.png)

