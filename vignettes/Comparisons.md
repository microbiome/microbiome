<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - comparisons}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->
Group-wise comparisons
----------------------

### Boxplots

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

![](figure/boxplot-example-1.png)

### Negative binomial example

[Read more](http://www.ats.ucla.edu/stat/r/dae/nbreg.htm)

    library(MASS)
    taxa <- taxa_names(x)[1:2]
    x <- atlas1006
    df <- as(sample_data(x), "data.frame")
    for (tax in taxa) {
      df$signal <- get_sample(x, tax)
      res <- glm.nb(signal ~ bmi_group + gender, data = df)
      print(coef(summary(res)))
    }

### Comparisons for individual taxa with random effect subject term

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

Linear models with limma
------------------------

Identify most significantly different taxa between males and females.

For further details, see [limma
homepage](http://bioinf.wehi.edu.au/limma/) and [limma User's
guide](http://www.lcg.unam.mx/~lcollado/R/resources/limma-usersguide.pdf).
For discussion on why limma is preferred over t-test, see [this
article](http://www.plosone.org/article/info:doi/10.1371/journal.pone.0012336).

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

<table>
<thead>
<tr class="header">
<th></th>
<th align="right">logFC</th>
<th align="right">AveExpr</th>
<th align="right">t</th>
<th align="right">P.Value</th>
<th align="right">adj.P.Val</th>
<th align="right">B</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Uncultured Clostridiales II</td>
<td align="right">-0.41</td>
<td align="right">1.37</td>
<td align="right">-3.72</td>
<td align="right">0</td>
<td align="right">0.06</td>
<td align="right">-0.24</td>
</tr>
<tr class="even">
<td>Eubacterium siraeum et rel.</td>
<td align="right">-0.34</td>
<td align="right">1.67</td>
<td align="right">-3.52</td>
<td align="right">0</td>
<td align="right">0.06</td>
<td align="right">-0.77</td>
</tr>
<tr class="odd">
<td>Clostridium nexile et rel.</td>
<td align="right">0.18</td>
<td align="right">2.84</td>
<td align="right">3.41</td>
<td align="right">0</td>
<td align="right">0.06</td>
<td align="right">-1.04</td>
</tr>
<tr class="even">
<td>Sutterella wadsworthia et rel.</td>
<td align="right">-0.33</td>
<td align="right">1.50</td>
<td align="right">-3.13</td>
<td align="right">0</td>
<td align="right">0.10</td>
<td align="right">-1.74</td>
</tr>
</tbody>
</table>

**Q-Q plot for limma**

    qqt(fit$t[, coef.index], df = fit$df.residual + fit$df.prior)
    abline(0,1)

![](figure/limma-qq-1.png)

**Volcano plot for limma**

    volcanoplot(fit, coef = coef.index, highlight = coef.index)

![](figure/limma-volcano-1.png)

### PERMANOVA

PERMANOVA can be also used to assess community-level differences between
groups. Here let us evaluate whether nationality has a significant
effect on gut microbiota.

    # Example data
    data("dietswap")
    x <- dietswap
    group <- "nationality"

    # Use relative abundances for simpler visualizations
    x <- transform_phyloseq(x, "compositional")
    otu <- get_sample(x)
    meta <- as(sample_data(x), "data.frame")
    meta$group <- meta[[group]]

    # PERMANOVA: samples x species as input
    library(vegan)
    permanova <- adonis(t(otu) ~ group, data=meta, permutations=99, method = "bray")
    pv <- as.data.frame(permanova$aov.tab)["group", "Pr(>F)"]

    # P-value
    print(pv)

    ## [1] 0.01

    # Note the assumption of similar
    # multivariate spread among the groups
    # ie. analogous to variance homogeneity
    # Here the groups have signif. different spreads and
    # permanova result may be explained by that.
    dist <- vegdist(t(otu))
    anova(betadisper(dist,meta$group))

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##            Df  Sum Sq  Mean Sq F value    Pr(>F)    
    ## Groups      1 0.26114 0.261144  18.345 2.754e-05 ***
    ## Residuals 220 3.13168 0.014235                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    # Coefs for the top taxa separating the groups
    coef <- coefficients(permanova)["group1",]
    top.coef <- coef[rev(order(abs(coef)))[1:20]]
    par(mar = c(3, 14, 2, 1))
    barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")

![](Comparisons_files/figure-markdown_strict/comparisons-permanova-1.png)
