### TODO

 * Homogeneity vs. sample size dependence. Normalize.

Mean squared error.

x^2/n

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
sum((x - xm) * (y - xm)) / (sqrt(sum((x-xm)^2)) * sqrt(sum((y-xm)^2)))

  * Heatmap: indicate groups with a color bar

  * Licensing issues with phyloseq


