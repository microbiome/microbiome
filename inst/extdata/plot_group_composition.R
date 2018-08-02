
#library(microbiome)
#library(reshape2)
# Example data
#data(atlas1006)
#x <- atlas1006
# Grouping variable
#group <- "bmi_group"
# Selected OTUs (core taxa in this example)
#otus <- core_members(transform(x, "compositional"), detection = 1/100, prevalence = 50/100)
#x <- prune_taxa(otus, x)
# Visualize
#library(microbiome)
#source("~/Rpackages/microbiome/inst/extdata/plot_group_composition.R")
#p <- plot_group_composition(x, group, limits = c(-1,1), step = 0.5, order.cols = FALSE, type = "heatmap")
#p <- plot_group_composition(x, group, limits = c(-1,1), step = 0.5, order.cols = FALSE, type = "boxplot") 
#print(p)


plot_group_composition <- function (x, group, limits = NULL, step = 1, star = NULL, order.cols = TRUE, order.rows = TRUE, type = "heatmap") {

  # Grouping		       
  if (length(group) == 1) {
    groups <- meta(x)[, group]
  } else{
    groups <- group
  }
  if (is.null(names(groups))) {
    names(groups) <- colnames(A)
  }

  # Log10p abundances 
  # A <- log10(1 + 1e4 * abundances(x))
  # Log10 relative abundances
  # A <- log10(abundances(transform(x, "compositional")))
  # CLR transformed abundances
  A <- abundances(transform(x, "clr"))

  ## Scale and center to higlight differences of the log abundances
  d <- as.data.frame(scale(t(A)));

  # No scaling
  #d <- as.data.frame(t(A));
  d$group <- groups[colnames(A)]
  
  if (type == "heatmap") {

    # Mean z score per group
    a <- aggregate(. ~ group, d, mean)
    # Organize
    d <- melt(a, id = "group")
    names(d) <- c("group", "Genus", "Abundance")
    d$group <- factor(d$group)
    d$Genus <- factor(d$Genus)

    p <- heat(d, Xvar = "group", Yvar = "Genus", fill = "Abundance",
      star = star, step = step, limits = limits,
      order.cols = order.cols,
      order.rows = order.rows) +
      labs(title = "Z-transformed CLR scores") + 
      theme(axis.text.x = element_text(angle = 90))

  } else if (type == "boxplot") {

    # Then pick all abundances for these 
    dfs2 <- NULL
    for (dm in na.omit(unique(groups))) {
      a <- A[,which(groups == dm)]
      # a <- a + min(a[a>0])/2 # Add small constant for visualization purposes
      df <- melt(a)
      df$group <- dm
      dfs2 <- rbind(dfs2, df)
    }
    colnames(dfs2) <- c("Genus", "Sample", "Abundance", "Group")
    dfs2$Genus <- factor(dfs2$Genus)
    dfs2$Group   <- factor(dfs2$Group, levels = levels(groups))

    theme_set(theme_bw(20))
    #dfs2$Genus <- factor(dfs2$Genus, levels(figure.rf.discriminators$data$Genus))
    p <- ggplot(dfs2, aes(x = Genus, y = Abundance, fill = Group)) +
      geom_boxplot() +
      scale_y_continuous(label = scales::percent, trans = "log10", breaks = c(.1/100, 1/100, 10/100, 1)) +
      labs(y = "Relative abundance (%)", x = "") +
      # scale_fill_manual(values = dmmcolors[Group[as.character(dfs2$Sample)]]) +
      coord_flip()

  }

  p
  
}


