library(rmarkdown)
rmarkdown::render("vignette.Rmd")
rmarkdown::render("Atlas.Rmd")
rmarkdown::render("Template.Rmd")

#render("vignette.Rmd", "html_document")
#rmarkdown::render("vignette.Rmd", "all")
#library(knitr)
#knit("vignette.Rmd")
#library(knitr)
#knit("Atlas.Rmd")

# ---------------------------------------------

library(knitr)
fs <- list.files(pattern = ".Rmd$")
#fs <- "vignette.Rmd"
#fs <- sample(list.files(pattern = ".Rmd$"))
#fs <- "Profiling.Rmd"
#fs <- c("Metrics.Rmd", "Heatmap.Rmd")
#fs <- "SQL.Rmd"
#fs <- "Heatmap.Rmd"
#fs <- c("RDA.Rmd", "LatentClassAnalysis.Rmd", "NetResponse.Rmd")
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
for (f in setdiff(fs, "Installation.Rmd")) { 
    print(f)
    knit(f) 
}

system("git add *.md")
system("git add figure/*")
system("git add *.Rmd")
system("git commit -a -m'markdown update'")
system("git push")


