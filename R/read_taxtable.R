read_taxtable <- function (taxonomy.file, sep = ",") {

    s.tax <- read.csv(taxonomy.file, row.names=1, check.names=FALSE, sep = sep)
    s.taxmat <- as.matrix(s.tax)

    tax_table(s.taxmat)

}

