### TODO

 * Homogeneity vs. sample size dependence. Normalize.

Mean squared error.

((x - mean(x))^2)/(x^2 * mean(x))

n <- 1000
x <- rnorm(n)
xm <- mean(x)
y <- rnorm(n)
ym <- mean(y)
# Correlation is normalized difference
cc <- sum((x - xm) * (y - ym)) / (sqrt(sum((x-xm)^2)) * sqrt(sum((y-ym)^2)))
cc
cor(x, y)


 * Bagged and ordinary RDA for phyloseq objects with ggplot visualizations

 * Heatmap: indicate groups with a color bar

 * PERMANOVA shows the communities are significantly different across groups
   adonis(t(dat.full) ~ group, data=meta.full, permutations=99)

  * Licensing issues with phyloseq