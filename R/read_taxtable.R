read_taxtable <- function (taxonomy.file) {

    s.tax <- read.csv(taxonomy.file, row.names=1, check.names=FALSE)
    s.taxmat <- as.matrix(s.tax)
    s.tax_table <- tax_table(s.taxmat)
    s.tax_table
}

