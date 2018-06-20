library(ggplot2)
library(rvg)
library(ggiraph)
library(devtools)
load_all()
x <- microbiome::transform(atlas1006, "compositional")

mytheme_main <- theme( panel.background = element_blank(), 
  panel.grid.major = element_line(colour = "#dddddd"), 
  axis.ticks = element_line(colour = "#dddddd") )

mytheme_map <- theme(
  panel.background = element_blank(), axis.title.x = element_blank(),
  axis.text = element_blank(), axis.line.x = element_blank(),
  axis.line.y = element_blank(), axis.title.y = element_blank(),
  axis.ticks.x = element_blank(), axis.ticks.y = element_blank() )

df <- as(sample_data(x), "data.frame")
df$Dialister <- get_sample(x, "Dialister")
df$Prevotella <- get_sample(x, "Prevotella melaninogenica et rel.")
df$sample <- row.names(df)

# geom_point_interactive example
gg_point_1 <- ggplot(df, aes(x = Prevotella, y = Dialister, 
        color = age, tooltip = sample) ) + 
    geom_point_interactive(size=3)

# htmlwidget call
ggiraph(code = {print(gg_point_1 + mytheme_main)}, width = 7, height = 6)
