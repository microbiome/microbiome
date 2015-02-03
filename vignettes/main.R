library(rmarkdown)
render("vignette.Rmd", "html_document")

library(knitr)
knit("vignette.Rmd")

library(knitr)
knit("Atlas.Rmd")

