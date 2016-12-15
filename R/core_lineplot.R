core_lineplot <- function (x, title = "Common core",  
                   xlabel = "Abundance", 
                   ylabel = "Core size (N)") {

    Abundance <- Prevalence <- Count <- NULL
    
    df <- as.data.frame(x)
    df$ID <- rownames(x)
    df <- gather(df, "ID")
    names(df) <- c("Abundance", "Prevalence", "Count")
    df$Abundance <- as.numeric(as.character(df$Abundance))
    df$Prevalence <- as.numeric(as.character(df$Prevalence))
    df$Count <- as.numeric(as.character(df$Count))    

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
        
    list(plot = p, data = x)
}


