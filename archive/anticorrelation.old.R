anticorrelation.old <- function(x, type = "interindividual", group_by = "group", method = "spearman") {

    # Pick metadata		     
    meta <- sample_data(x)

    if (!"sample" %in% names(meta)) {
      warning("Using the metadata rownames as the sample ID")
      meta$sample = rownames(sample_data(x))
    }

    # OTU data Log10
    otu <- log10(t(abundances(x)@.Data))

    # Ensure compatiblity
    if (!nrow(otu) == nrow(meta)) {
      otu <- t(otu)
    }

    if (!all(rownames(otu) == rownames(meta))) {
      stop("OTU table and metadata do not match.")
    }

    correlation <- NULL
 		         
    # Split the data by group
    group <- NULL
    if (!group_by %in% names(meta)) {
      stop(paste("The group_by variable", group_by, "is not included in sample_data(x)."))
    }
    datasets <- split(as.data.frame(otu), meta[[group_by]], drop = TRUE)

    if (type == "interindividual") {

      tmp <- setdiff(c("sample", group_by), names(meta))
      if (length(tmp) > 0) {
        stop(paste("The following variables needed by betadiversity function type=interindividual are 
	            missing from sample metadata:", paste(tmp, collapse = ",")))
      }
  
      # Within-matrix stability NOTE: earlier this was calculated as
      # the average of upper triangular correlation matrix This is
      # heavily biased since the values are dependent Now replaced
      # by calculating correlations against the mean of the whole
      # sample set cors <- lower.triangle(cor(dat1))
      dfs <- NULL      
      for (ds in names(datasets)) {
        dat1 <- datasets[[ds]]
        cors <- as.vector(cor(t(dat1), matrix(colMeans(dat1)), method = method, use = "pairwise.complete.obs"))
        dfs <- rbind(dfs, data.frame(group = rep(ds, length(cors)),
                   	             sample = rownames(dat1),
  		   		     anticorrelation = 1 - cors))
      }
  
      pval <- anova(lm(anticorrelation ~ group, data = dfs))[["Pr(>F)"]][[1]]
      stats <- dfs %>% group_by(group) %>%
                       summarize(mean = mean(anticorrelation),
		       sd = sd(anticorrelation))
      homogeneity <- list(data = dfs, statistics = stats, p.value = pval)

    } else if (type == "intraindividual") {

      tmp <- setdiff(c("time", "subject", "sample", group_by), names(meta))
      if (length(tmp) > 0) {
        stop(paste("The following variables needed by betadiversity function 
               type=intraindividual are missing from sample metadata:", paste(tmp, collapse = ",")))
      }

      if (!all(sapply(split(meta$time, meta[[group_by]]), function (x) {length(unique(x))}) == 2)) {
        stop("Two time points needed for each group for the intraindividual type. 
              Some groups are having a different number of time points.")
      }

      homogeneity <- list()
      dfs <- NULL
      for (ds in names(datasets)) {
      
        # Pick the data and metadata for this group
        xsub <- datasets[[ds]]        
	msub <- meta[rownames(xsub),]

        # Use interindividual functionality to assess correlations
	# within subjects. Subjects are used as groups
        datasets2 <- split(xsub, droplevels(msub[["subject"]]))
	cors <- c()
        for (subj in names(datasets2)) {
          dats <- datasets2[[subj]]
          cors[[subj]] <- cor(unlist(dats[1,]), unlist(dats[2,]),
	                      method = method, use = "pairwise.complete.obs")
        }

        dfs <- rbind(dfs, data.frame(group = rep(ds, length(cors)),
                   	             subject = names(cors),
  		   		     anticorrelation = 1 - cors))

      }

      # Between time point correlations within subjects
      # and the mean over those correlations
      pval <- anova(lm(anticorrelation ~ group, data = dfs))[["Pr(>F)"]][[1]]      
      stats <- dfs %>%
                 group_by(group) %>%
                 summarize(mean = mean(anticorrelation, na.rm = T),
		           sd = sd(anticorrelation, na.rm = T))
			   
      anticor <- list(data = dfs, statistics = stats, p.value = pval)

    }
    
    anticor
    
}





