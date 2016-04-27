tolong <- function (x, id.vars = NULL) {

    n <- ncol(x)
    d <- NULL
    
    if (is.null(colnames(x))) {
      colnames(x) <- as.character(1:ncol(x))
    }
    if (is.null(rownames(x))) {
      rownames(x) <- as.character(1:nrow(x))
    }
    
    cn <- setdiff(colnames(x), id.vars)
    
    for (k in rownames(x)) {
      if (is.null(id.vars)) {
        d[[k]] <- cbind(rep(k, n), cn, unlist(x[k, ]))
      } else {
        d[[k]] <- cbind(rep(id.vars, n-1), cn, unlist(x[k, cn]))
      }      
    }
    
    d <- do.call("rbind", d)
    rownames(d) <- NULL
    if (is.null(id.vars)) {id.vars <- "VarX"}
    colnames(d) <- c(id.vars, "VarY", "value")
    
    data.frame(d)

}    