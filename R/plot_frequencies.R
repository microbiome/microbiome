#' @title Plot Frequencies
#' @description Plot relative frequencies within each Group for the
#' levels of the given factor.
#' @param x \code{\link{data.frame}} 
#' @param Groups Name of the grouping variable
#' @param Factor Name of the frequency variable
#' @return \code{\link{ggplot}} plot object.
#' @details For table with the indicated frequencies, see the returned
#' phyloseq object.
#' @export
#' @examples
#' data(dietswap)
#' p <- plot_frequencies(meta(dietswap), 'group', 'sex')
#' @keywords utilities
plot_frequencies <- function(x, Groups, Factor) {
    
    # FIXME: do not capitalize arguments Rename Groups -> groups
    # Factor -> factor
    
    pct <- NULL
    
    x <- data.frame(x)
    x$Groups <- x[[Groups]]
    x$Factor <- x[[Factor]]
    x <- x %>% group_by(Groups, Factor) %>%
        summarise(n=n()) %>%
        mutate(pct=100 * n/sum(n))
    
    # Provide barplot of relative abundances
    p <- ggplot(x, aes(x=Groups, y=pct, fill=Factor)) +
        geom_bar(position="stack", stat="identity")
    
    # Rotate horizontal axis labels, and adjust
    p <- p + theme(axis.text.x =
        element_text(angle=-90, vjust=0.5, hjust=0)) +
        labs(x = "", y = "Proportion (%)") 
    
    p
    
}
