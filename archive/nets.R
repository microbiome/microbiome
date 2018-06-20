```{r networks3, warning=FALSE, message=FALSE}
theme_set(theme_bw(20))
p <- plot_net(pseq, maxdist = 0.2,
              shape = "group", color = "nationality",
	      distance = "bray", laymeth = "auto") +
     scale_colour_brewer(palette = "Accent")
print(p)		 
```
