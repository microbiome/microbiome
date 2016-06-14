library(plotly)
library(microbiome)
set.seed(100)
m <- sample_data(atlas)
x <- t(taxa_abundances(atlas))
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
plot_ly(d, x = Comp.1, y = Comp.2, mode = "markers", color = nationality, size = age)
plot_ly(d, x = Comp.1, y = Comp.2, mode = "markers", color = log10(d$Prevotella_melaninogenica_et_rel), size = age)
plot_ly(d, x = Comp.1, y = Comp.2, mode = "markers", color = log10(Bacteroides_fragilis_et_rel), size = age)
plot_ly(d, x = Comp.1, y = Comp.2, mode = "markers", color = age, size = age)


library(ggplot2)
theme_set(theme_bw(30))
ggplot(d, aes(x = Comp.1, y = Comp.2, color = DNA_extraction_method, size = age)) + geom_point()
ggplot(d, aes(x = Comp.1, y = Comp.2, color = nationality, size = bmi_group)) + geom_point()
ggplot(d, aes(x = Comp.1, y = Comp.2, color = health_status, size = bmi_group)) + geom_point()
ggplot(d, aes(x = Comp.1, y = Comp.2, color = Prevotella_melaninogenica_et_rel, size = Prevotella_melaninogenica_et_rel)) + geom_point()
ggplot(d, aes(x = Comp.1, y = Comp.2, color = Bacteroides_fragilis_et_rel, size = Bacteroides_fragilis_et_rel)) + geom_point()
ggplot(d, aes(x = Comp.1, y = Comp.2, color = Dialister, size = Dialister)) + geom_point()

