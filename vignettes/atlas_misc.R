library(plotly)
library(microbiome)
set.seed(100)
m <- sample_data(atlas)
x <- t(taxa_abundances(atlas))
pca = princomp(x)$scores[, 1:2]
d = data.frame(cbind(m, x))
plot_ly(d, x = Comp.1, y = Comp.2, mode = "markers", color = nationality, size = age)

