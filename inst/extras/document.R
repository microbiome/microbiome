library(devtools)
document("../../")

library(knitr)
knit(input = "../../README.Rmd", output = "../../README.md")

#setwd("../../")
#library(pkgdown)
#build_site()


