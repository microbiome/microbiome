# Just apply variance stabilizing transformation
# varianceStabilizingTransformation(ds2, blind = TRUE, fitType = "parametric")
# ds2 <- estimateSizeFactors(ds2)
# ds2 <- estimateDispersions(ds2)
# abund <- getVarianceStabilizedData(ds2)

#dfm <- as.data.frame(otu)
  #dfm$ID <- rownames(otu)
  #dfm <- gather(dfm, ID)    
  #colnames(dfm) <- c("Taxon", "Sample", "Abundance")
  #dfm$Abundance <- as.numeric(as.character(dfm$Abundance))  
  #dfm$Sample <- factor(as.character(dfm$Sample), levels = sort.samples)
  #dfm$ID <- factor(as.character(dfm$ID), levels = sort.taxa)  
