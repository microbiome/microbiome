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

#library(devtools)
#load_all()
library(microbiome)

#library(rmarkdown)
#rmarkdown::render("index.Rmd")
#rmarkdown::render("Template.Rmd")
#rmarkdown::render("Atlas.Rmd", "all")

#render("index.Rmd", "html_document")
#rmarkdown::render("index.Rmd", "all")
library(knitr)
library(knitcitations)
#knit("index.Rmd")
#library(knitr)
#knit("Atlas.Rmd")

# ---------------------------------------------

library(knitr)
#fs <- "Density.Rmd"
#fs <- "Diversity.Rmd"
#fs <- "Heatmap.Rmd"
#fs <- "Profiling.Rmd"
#fs <- "RDA.Rmd"
#fs <- "SQL.Rmd"
#fs <- "index.Rmd"
fs <- sample(list.files(pattern = ".Rmd$"))
knitr::opts_chunk$set(fig.path = "figure/", dev="CairoPNG")
times <- c()
namespaces0 <-  loadedNamespaces()

for (myfile in fs) {
    print(myfile)
    # times[[myfile]] <- system.time(knit(myfile))[["elapsed"]]
    # rmarkdown::render(myfile, "md_document")
    # rmarkdown::render(myfile, "all")
    rmarkdown::render(myfile, "html_document")    

    # Must do to clean up some space
    #for (i in 1:10) {
    #  pkgs0 <- setdiff(loadedNamespaces(), namespaces0)
    #  for (pkg in pkgs0) {
    #    print(c(pkg, length(pkgs0)))
    #    try(a <- unloadNamespace(pkg))
    #  }
    #}
}

# Time per index.page
# par(mar = c(3, 10, 1, 1)); barplot(sort(times), horiz = T, las = 1)

#system("git add *.md")
#system("git add figure/*")
#system("git add *.Rmd")
system("git add *.html")
system("git commit -a -m'homepage update'")
system("git push origin gh-pages")
system("rm -rf cache")
system("rm -rf figure")
