library(microbiome)

# Example data
data(dietswap)
# Rename
pseq <- dietswap

# CLR transformation to remove compositionality biases
# Add small constant on the data to avoid problems from zero counts
pseq <- transform(pseq, "clr")

# Samples x taxa abundance matrix
otu <- t(abundances(pseq))

# PCA
proj <- princomp(otu)$scores[, 1:2]

# Pick the metadata
df <- meta(pseq)

# Add projection
df <- cbind(df, proj[rownames(df),])

# Visualize
library(ggplot2)
theme_set(theme_bw(20))
p <- ggplot(df, aes(x = Comp.1, y = Comp.2)) + geom_point(aes(color = group))

# Add trajectory for given subject
dfs <- subset(df, subject == "zaq")
p2 <- p + geom_path(data = dfs)

# Print figure
print(p2)


