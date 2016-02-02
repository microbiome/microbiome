### TODO

 * Negative binomial for sample comparisons

  library(MASS)
  df <- sample_data(x)
  df$signal <- get_sample(x, "Dialister")
  summary(res <- glm.nb(daysabs ~ group + sex, data = df))

 * gower distance for sample similarities

 * SPIEC-EASE / SparCC network algorithms from the huge package -> Healthy vs other networks ?

 * Homogeneity vs. sample size dependence. Normalize.

 * Bagged and ordinary RDA for phyloseq objects with ggplot visualizations

 * Heatmap: indicate groups with a color bar

 * PERMANOVA shows the communities are significantly different across groups
   adonis(t(dat.full) ~ group, data=meta.full, permutations=99)

  * Licensing issues with phyloseq