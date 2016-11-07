#' @title Chunk Reorder
#' @description Chunk re-order a vector so that specified newstart is first. Different than relevel.
#' @keywords internal
#' @details Borrowed from \pkg{phyloseq} package as needed here and not exported there.
#' @examples 
#' # Typical use-case
#' # chunk_reorder(1:10, 5)
#' # # Default is to not modify the vector
#' # chunk_reorder(1:10)
#' # # Another example not starting at 1
#' # chunk_reorder(10:25, 22)
#' # # Should silently ignore the second element of `newstart`
#' # chunk_reorder(10:25, c(22, 11))
#' # # Should be able to handle `newstart` being the first argument already
#' # # without duplicating the first element at the end of `x`
#' # chunk_reorder(10:25, 10)
#' # all(chunk_reorder(10:25, 10) == 10:25)
#' # # This is also the default
#' # all(chunk_reorder(10:25) == 10:25)
#' # # An example with characters
#' # chunk_reorder(LETTERS, "G") 
#' # chunk_reorder(LETTERS, "B") 
#' # chunk_reorder(LETTERS, "Z") 
#' # # What about when `newstart` is not in `x`? Return x as-is, throw warning.
#' # chunk_reorder(LETTERS, "g") 
chunk_reorder <- function(x, newstart = x[[1]]){

  pivot <- match(newstart[1], x, nomatch = NA)

  # If pivot is NA, then warn and return x as-is
  if(is.na(pivot)){
    warning("The `newstart` argument was not in `x`. Returning `x` without reordering.")
    newx <- x
  } else {
    newx <- c(tail(x, {length(x) - pivot + 1}), head(x, pivot - 1L))
  }
  
  newx

}
