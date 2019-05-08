library(microbiome)
library(phyloseq)
library(netresponse)
library(MASS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(devtools)
library(ggplot2)

#library(devtools)
#load_all()
library(microbiome)

library(rmarkdown)
library(knitr)
library(knitcitations)

# ---------------------------------------------

library(knitr)
#fs <- "Density.Rmd"
#fs <- "Experimental.Rmd"
#fs <- c("DMM.Rmd", "index.Rmd")
#fs <- "Diversity.Rmd"
#fs <- "Heatmap.Rmd"
#fs <- "Profiling.Rmd"
#fs <- "RDA.Rmd"
#fs <- "SQL.Rmd"
#fs <- "index.Rmd"
#fs <- "rstanarm.Rmd"
#fs <- "Landscaping.Rmd"
#fs <- sample(list.files(pattern = ".Rmd$"), 20)
fs <- sample(list.files(pattern = ".Rmd$")) # Random order

knitr::opts_chunk$set(fig.path = "figure/", dev="CairoPNG")

times <- c()
namespaces0 <-  loadedNamespaces()

#for (myfile in setdiff(fs, "Themes.Rmd")) {
for (myfile in fs) {

    rmarkdown::render(myfile)

}

# Time per index.page
# par(mar = c(3, 10, 1, 1)); barplot(sort(times), horiz = T, las = 1)

#system("git add *.md")
#system("git add figure/*")
#system("git add *.Rmd")
system("git add *.html")
system("git commit -a -m'homepage update'")
system("git push origin gh-pages")
#system("rm -rf cache")
#system("rm -rf figure")
