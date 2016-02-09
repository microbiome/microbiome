library(microbiome)
library(phyloseq)
library(netresponse)
library(MASS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sorvi)
library(limma)
library(devtools)
library(GGally)
library(ggnet)
library(network)
library(sna)
library(ggplot2)
library(intergraph) # ggnet2 works also with igraph with this

#library(rmarkdown)
#rmarkdown::render("vignette.Rmd")
#rmarkdown::render("Template.Rmd")
#rmarkdown::render("Atlas.Rmd", "all")

#render("vignette.Rmd", "html_document")
#rmarkdown::render("vignette.Rmd", "all")
library(knitr)
library(knitcitations)
#knit("vignette.Rmd")
#library(knitr)
#knit("Atlas.Rmd")

# ---------------------------------------------

library(knitr)
#fs <- c("Barplots.Rmd")
#fs <- c("Boxplots.Rmd")
#fs <- c("Clustering.Rmd")
#fs <- c("Comparisons.Rmd")
#fs <- c("Core.Rmd")
#fs <- c("Crosshyb.Rmd")
#fs <- "Density.Rmd"
#fs <- "Diversity.Rmd"
#fs <- "Heatmap.Rmd"
#fs <- "Profiling.Rmd"
#fs <- "RDA.Rmd"
#fs <- "SQL.Rmd"
#fs <- c("LatentClassAnalysis.Rmd", "NetResponse.Rmd")
#fs <- c("Phyloseq.Rmd")
#fs <- c("RPA.Rmd")
#fs <- c("wurcomputer.Rmd")
#fs <- c("Heatmap.Rmd")
#fs <- c("RPAtest.Rmd")
#fs <- c("Stability.Rmd")
#fs <- c("Installation.Rmd")
#fs <- c("ROC.Rmd")
#fs <- c("limma.Rmd")
fs <- c("Ordination.Rmd")
#fs <- c("Phylogeny.Rmd")
#fs <- "vignette.Rmd"
#fs <- sample(list.files(pattern = ".Rmd$"))
#fs <- sample(list.files(pattern = ".Rmd$"))
#for (f in setdiff(fs, c("misc.Rmd"))) { 
#for (f in setdiff(fs, c("Atlas.Rmd", "vignette.Rmd"))) {
# Motionchart as the last one
for (f in setdiff(fs, "misc.Rmd")) {
    print(f)
    knit(f) 
    #rmarkdown::render(f, "md_document")
}

system("git add *.md")
system("git add figure/*")
system("git add *.Rmd")
system("git commit -a -m'markdown update'")
system("git push")


