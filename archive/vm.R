#fs <- list.files("~/HITChipDB/R/", full.names = TRUE)
#for (f in fs)  { source(f) }
#fs <- list.files("~/microbiome/R/", full.names = TRUE)
#for (f in fs)  { source(f) }

# -------------------------------------------
library(devtools)
install_github("microbiome/microbiome")
install_github("microbiome/HITChipDB")

res <- run.profiling.script(dbuser, dbpwd, dbname, verbose = FALSE, host = host, port = port, summarization.methods = summarization.methods, which.projects = which.projects)

projs <- list.mysql.projects(dbuser, dbpwd, dbname, host = host, port = port)

res <- run.profiling.script(dbuser, dbpwd, dbname, verbose = FALSE, host = host, port = port, summarization.methods = summarization.methods, which.projects = "TURKU PET STUDY")

# -------------------------------



