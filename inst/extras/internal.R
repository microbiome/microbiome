
get.file.method <- function( f ) {
    
    method <- NULL
    
    if (f == "rpa") {
        method <- "rpa"
    }
    if (f == "frpa") {
        method <- "frpa"
    }
    if (!length(grep("sum", f)) == 0) {
        method <- "sum"
    }
    if (!length(grep("ave", f)) == 0) {
        method <- "ave"
    }
    if (!length(grep("nmf", f)) == 0) {
        method <- "nmf"
    }
    
    method
    
}

get.file.level <- function(f) {
    
    level <- NULL
    
    if (!length(grep("oligo", f)) == 0) {
        level <- "oligo"
    }
    if (!length(grep("species", f)) == 0) {
        level <- "species"
    }
    if (!length(grep("L0", f)) == 0) {
        level <- "L0"
    }
    if (!length(grep("L1", f)) == 0) {
        level <- "L1"
    }
    if (!length(grep("L2", f)) == 0) {
        level <- "L2"
    }
    if (!length(grep("phylogeny", f)) == 0) {
        level <- "phylogeny.info"
    }
    
    level
    
}

