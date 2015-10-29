core_lineplot <- function (x, title = "Common core",  
                   xlabel = "Abundance", 
                   ylabel = "Core size (number of taxa)", 
		   prevalence.intervals = NULL, 
		   detection.thresholds = NULL) {

    if (class(x) == "phyloseq") {
      x <- core_matrix(x, prevalence.intervals, detection.thresholds)
    }

    Abundance <- Prevalence <- Count <- NULL
    
    df <- melt(x)
    names(df) <- c("Abundance", "Prevalence", "Count")

    theme_set(theme_bw(20))
    p <- ggplot(df, aes(x = Abundance,
      	 	    	y = Count,
			color = Prevalence, 
                        group = Prevalence))
    
    p <- p + geom_line()
    p <- p + geom_point()
    p <- p + scale_x_log10()
    p <- p + xlab(xlabel)
    p <- p + ylab(ylabel)
    p <- p + ggtitle(title)
        
    return(p)
}


