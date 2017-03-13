load("~/data/Atlas/20160604/AtlasFull.RData")

library(plotly)
library(microbiome)
set.seed(100)
#m <- sample_data(atlas)
m = atlas.metadata
x = t(atlas$L2$rpa)
#x <- t(taxa_abundances(atlas))
pca = princomp(log10(x))$scores[, 1:2]
d = data.frame(cbind(m, x, pca))
d$nationality = as.character(d$nationality)
d$nationality[d$nationality %in% c("FIN", "SWE", "DNK")] = "NorthEurope"
d$nationality[d$nationality %in% c("BEL", "NLD", "DEU")] = "CentralEurope"
d$nationality[d$nationality %in% c("ITA", "ESP", "CC")] = "SouthEurope"
d$nationality[d$nationality %in% c("GBR", "IRL")] = "UK/IE"
d$nationality[d$nationality %in% c("IND", "PAK")] = "Asia"
d$nationality[d$nationality %in% c("ZAF")] = "AFR"
d$nationality[is.na(d$nationality)] = "Unknown"
d$nationality = factor(d$nationality)
#d = d[sample(1:nrow(d), 100), ]
d = subset(d, !is.na(age) & !is.na(nationality))
#plot_ly(d, x = Comp.1, y = Comp.2, mode = "markers", color = nationality, size = age)
#plot_ly(d, x = Comp.1, y = Comp.2, mode = "markers", color = log10(d$Prevotella_melaninogenica_et_rel), size = age)
#plot_ly(d, x = Comp.1, y = Comp.2, mode = "markers", color = log10(Bacteroides_fragilis_et_rel), size = age)
plot_ly(d, x = Comp.1, y = Comp.2, mode = "markers", color = age, size = age)

library(ggplot2)
theme_set(theme_bw(30))
ggplot(d, aes(x = Comp.1, y = Comp.2, color = DNA.extr.meth, size = max(age, na.rm = TRUE) - age)) + geom_point()
ggplot(d, aes(x = Comp.1, y = Comp.2, color = nationality, size = bmi)) + geom_point()
ggplot(d, aes(x = Comp.1, y = Comp.2, color = health.status, size = bmi)) + geom_point()
ggplot(d, aes(x = Comp.1, y = Comp.2, color = sample.type, size = age)) + geom_point()
ggplot(d, aes(x = Comp.1, y = Comp.2, color = Prevotella.melaninogenica.et.rel, size = Prevotella.melaninogenica.et.rel)) + geom_point()
ggplot(d, aes(x = Comp.1, y = Comp.2, color = Bacteroides.fragilis.et.rel, size = Bacteroides.fragilis.et.rel)) + geom_point()
ggplot(d, aes(x = Comp.1, y = Comp.2, color = Dialister, size = Dialister)) + geom_point()


