### Cross-correlation example

    dat1 <- peerj32$lipids # Lipids (44 samples x 389 lipids)
    dat2 <- peerj32$microbes # Microbiota (44 samples x 130 bacteria)
    meta <- peerj32$meta

    correlations <- cross.correlate(dat1, dat2, 
                            method = "bicor", 
                mode = "matrix", 
                            n.signif = 1, 
                p.adj.threshold = 0.05, 
                            p.adj.method = "BH")

    ## Warning in as.vector(x) == as.vector(y): longer object length is not a
    ## multiple of shorter object length

    correlation.table <- cmat2table(correlations)
    head(correlation.table)

    ##             X1                               X2 Correlation       p.adj
    ## 278   PC(40:3) Eubacterium cylindroides et rel.  -0.6745063 0.008697453
    ## 356   PC(40:3)                     Helicobacter  -0.6807777 0.008697453
    ## 530 TG(54:5).2      Ruminococcus gnavus et rel.   0.6789932 0.008697453
    ## 525   TG(52:5)      Ruminococcus gnavus et rel.   0.6612294 0.013011631
    ## 521   TG(50:4)      Ruminococcus gnavus et rel.   0.6556532 0.013771802
    ## 474  PC(40:3e)                    Moraxellaceae  -0.6506371 0.014691099
