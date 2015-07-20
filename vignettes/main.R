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
#fs <- "vignette.Rmd"
#fs <- sample(list.files(pattern = ".Rmd$"))
#fs <- "Profiling.Rmd"
#fs <- "SQL.Rmd"
#fs <- "Heatmap.Rmd"
#fs <- c("LatentClassAnalysis.Rmd", "NetResponse.Rmd")
#fs <- c("RDA.Rmd")
#fs <- c("Phyloseq.Rmd")
#fs <- c("Barplots.Rmd")
#fs <- c("Boxplots.Rmd")
#fs <- c("RPA.Rmd")
#fs <- c("Clustering.Rmd")
#fs <- c("wurcomputer.Rmd")
#fs <- c("Heatmap.Rmd")
#fs <- c("RDA.Rmd")
#fs <- c("Comparisons.Rmd")
#fs <- c("RPAtest.Rmd")
#fs <- c("Stability.Rmd")
#fs <- c("Core.Rmd")
#fs <- "Diversity.Rmd"
#fs <- "RDA.Rmd"
#fs <- "Density.Rmd"
#fs <- c("Crosshyb.Rmd")
#fs <- c("Installation.Rmd")
#fs <- c("ROC.Rmd")
#fs <- c("limma.Rmd")
#fs <- c("Phylogeny.Rmd")
#for (f in setdiff(fs, "Installation.Rmd")) { 
fs <- list.files(pattern = ".Rmd$")
for (f in setdiff(fs, c("misc.Rmd"))) { 
#for (f in setdiff(fs, c("Atlas.Rmd", "vignette.Rmd"))) { 
#for (f in setdiff(fs, c("misc.Rmd", "Networks.Rmd"))) {
    print(f)
    knit(f) 
    #rmarkdown::render(f, "md_document")
}

system("git add *.md")
system("git add figure/*")
system("git add *.Rmd")
system("git commit -a -m'markdown update'")
system("git push")


