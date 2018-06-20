---
title: "Linear models"
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


# Linear models

This section provides a hands-on introduction to linear and
generalized linear models, including covariates.


Load example data.

```{r lm1, warning=FALSE, message=FALSE}
library(microbiome)
data(dietswap)
d <- dietswap

# Pick microbial abundances for a given taxonomic group
taxa <- "Dialister"
```


Construct a data.frame with the selected taxonomic group and grouping.

```{r lm2, warning=FALSE, message=FALSE}
df <- data.frame(Abundance = abundances(d)[taxa,],
                 Group = meta(d)$nationality,
		 Log10_Abundance = log10(1 + df$Abundance)		 
		 )

# Compare the groups with a linear model.
# Use Log10 abundances 
res <- glm(Log10_Abundance ~ Group, data = df, family = "gaussian")
```


Investigate model coefficients

```{r lm_coefs, warning=FALSE, message=FALSE}
print(summary(res)$coefficients)

# The intercept equals to the mean in the first group
# The group term equals to the difference between group means
print(mean(subset(df, Group == "AAM")$Log10_Abundance))
print(mean(subset(df, Group == "AFR")$Log10_Abundance) - mean(subset(df, Group == "AAM")$Log10_Abundance))
```


Significance with t-test assuming equal variances.

```{r lm_signif, warning=FALSE, message=FALSE}
print(t.test(Log10_Abundance ~ Group, data = df, var.equal=TRUE)$p.value)
```


Now, using the linear model allows incorporation of additional variables,
for instance potential confounders.

```{r lm_glm, warning=FALSE, message=FALSE}
df$sex <- meta(d)$sex
res <- glm(Log10_Abundance ~ Group + sex, data = df, family = "gaussian")
```


Generalized linear models

The GLM consists of three elements:
    1. A probability distribution (from exponential family)
    1. A linear predictor η = Xβ .
    1. A link function g such that E(Y) = μ = g−1(η).

We use Poisson with (its natural) log-link.

```{r lm_glm2, warning=FALSE, message=FALSE}
res <- glm(Abundance ~ Group, data = df, family = "poisson")
print(summary(res))
```

