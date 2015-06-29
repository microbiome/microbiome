library(microbiome)
data(peerj32); 
x <- as.matrix(peerj32$lipids)
y <- as.matrix(peerj32$microbes[, 1])
cc <- cross.correlate(x, y); 



p <- correlation.heatmap(cc, "X1", "X2", "Correlation")
