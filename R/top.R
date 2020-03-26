#' @title Identify Top Entries
#' @description Identify top entries in a vector or given field in data frame.
#' @param x data.frame, matrix, or vector 
#' @param field Field or column to check for a data.frame or matrix
#' @param n Number of top entries to show
#' @param output Output format: vector or data.frame
#' @param round Optional rounding
#' @param na.rm Logical. Remove NA before calculating the statistics.
#' @param include.rank Include ranking if the output is data.frame. Logical.
#' @return Vector of top counts, named by the corresponding entries
#' @export
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("bibliographica")
#' @examples
#'   data(dietswap)
#'   p <- top(meta(dietswap), "group", 10)
#' @keywords utilities
top <- function (x, field = NULL, n = NULL, output = "vector", round = NULL, na.rm = FALSE, include.rank = FALSE) {

    if (is.factor(x)) {
        x <- as.character(x)
    }

    if (is.vector(x)) {
        if (na.rm) {
        inds <- which(x == "NA")
        if (length(inds) > 0) {
            x[inds] <- NA
            warning(paste("Interpreting NA string as missing value NA. 
            Removing", length(inds), "entries"))
        }
        x <- x[!is.na(x)]
    }
    
    s <- rev(sort(table(x)))
    N <- length(x)
    
    } else if (is.data.frame(x) || is.matrix(x)) {
        if (is.null(field)) {
        return(NULL)
    }

    x <- x[, field]
    if (na.rm) {
        inds <- which(x == "NA")
        if (length(inds) > 0) {
            x[inds] <- NA
            warning(
            paste("Interpreting NA string as missing value NA. Removing",
            length(inds), "entries"))
        }
        x <- x[!is.na(x)]
    }
    N <- length(x)
    s <- rev(sort(table(x)))
    }

    if (!is.null(n)) {
        s <- s[seq_len(min(n, length(s)))]
    }

    if (output == "data.frame") {
        s <- data_frame(name = names(s),
                    n = unname(s),
                    fraction = 100*unname(s)/N)
        if (is.null(field)) {field <- "Field"}
        names(s) <- c(field, "Entries (N)", "Fraction (%)")
        if (!is.null(round)) {
            s[,3] = round(s[,3], round)
        }
        if (include.rank) {
            s <- cbind(Rank = seq_len(nrow(s)), s)
        } 
    }

    s

}

