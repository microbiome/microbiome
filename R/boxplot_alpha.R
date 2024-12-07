#' @title Alpha Boxplot
#' @description Plot alpha index.
#' @param x \code{\link{phyloseq-class}} object
#' @param x_var Metadata variable to map to the horizontal axis.
#' @param index Alpha index to plot. See function \code{alpha}. 
#' @param zeroes Include zero counts in diversity estimation. Default is TRUE
#' @param violin Use violin version of the boxplot
#' @param na.rm Remove NAs
#' @param show.points Include data points in the figure
#' @param element.alpha Alpha value for plot elements. Controls the 
#'                      transparency of plots elements.
#' @param element.width Width value for plot elements. Controls the 
#'                      transparency of plots elements.
#' @param fill.colors Specify a list of colors passed on to ggplot2 
#'                    \code{scale_fill_manual}
#' @param outlier.fill If using boxplot and and points together how to deal with 
#'                     outliers. See ggplot2 outlier.fill argument in 
#'                     geom_ elements.  
#' @details A simple wrapper to visualize alpha diversity index.
#' @return A \code{\link{ggplot}} plot object
#' @export
#' @examples
#' data("dietswap")
#' p <- boxplot_alpha(dietswap, x_var = "sex", index="observed", violin=FALSE, 
#'                    na.rm=FALSE, show.points=TRUE, zeroes=TRUE, 
#'                    element.alpha=0.5, element.width=0.2, 
#'                    fill.colors= c("steelblue", "firebrick"),
#'                    outlier.fill="white")
#' p
#' 
#' @keywords utilities
boxplot_alpha <- function(x,
                          x_var = NULL,
                          index=NULL, 
                          violin=FALSE, 
                          na.rm=FALSE, 
                          show.points=TRUE,
                          zeroes=TRUE,
                          element.alpha=0.5,
                          element.width=0.2,
                          fill.colors= NA,
                          outlier.fill="grey50"){
    
    if(length(index) >1){
        stop("Please provide a single alpha index, e.g. index='shannon'")
    }
    
    d <- suppressMessages(alpha(x, index=index,zeroes=zeroes) )
    meta.df <-  cbind(meta(x),d) %>% 
        tibble::rownames_to_column(".sampleid")
    
    index.name <- colnames(meta.df)[ncol(meta.df)]
    
    if(is.null(x_var)){
        # Visualize example data with a boxplot
        p <- ggplot(meta.df, aes_string(".sampleid", index.name)) +
            geom_point() + 
            theme(axis.text.x = element_text(angle = 90)) +
            theme_bw()
        return(p)
    }
    
    if(!any(phyloseq::sample_variables(x) %in% x_var)){
      stop("'x_var' not available in `sample_data`")  
    }
    
    if(!is.na(fill.colors)[1]){
        .check.colors(meta.df, x_var, fill.colors)
    }
    
    p <- ggplot(meta.df, aes_string(x_var, index.name, fill=x_var)) 
    

    if (show.points) {
        p <- p + geom_point(size=2,
                            position=position_jitter(width=element.width), 
                            alpha=element.alpha,
                            shape=21)
    }
    
    # Box or Violin plot ?
    if (show.points & !violin) {
        p <- p + geom_boxplot(width=element.width,
                              alpha=element.alpha,
                              outlier.colour = outlier.fill,
                              outlier.fill = outlier.fill,
                              na.rm = na.rm)
    } else {
        p <- p + geom_violin(width=element.width,
                             alpha=element.alpha,
                             na.rm = na.rm)
    }  
    
    if(!is.na(fill.colors)[1]){
        p <- p + ggplot2::scale_fill_manual(values = fill.colors) +
            theme_bw()
    }
    p <- p + theme_bw()
    
    return(p)
}
    

.check.colors <- function(df, x_var=NULL, fill.colors=NULL){
    
    pl.vars = unique(df[,x_var])
    
    if(length(pl.vars) != length(fill.colors))
        stop("No. of fill.colors not equal to number of unique 'x_var'")
}



