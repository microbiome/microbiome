
### Cross-correlation example


```r
dat1 <- peerj32$lipids # Lipids (44 samples x 389 lipids)
dat2 <- peerj32$microbes # Microbiota (44 samples x 130 bacteria)
meta <- peerj32$meta

correlations <- cross.correlate(dat1, dat2, 
                        method = "bicor", 
			mode = "matrix", 
                        n.signif = 1, 
			p.adj.threshold = 0.05, 
                        p.adj.method = "BH")
```

```
## Warning in as.vector(x) == as.vector(y): longer object length is not a
## multiple of shorter object length
```

```r
correlation.table <- cmat2table(correlations)
head(correlation.table)
```

```
##              X1                               X2 Correlation       p.adj
## 1100 TG(54:5).2      Ruminococcus gnavus et rel.   0.7207818 0.001738478
## 1087   TG(52:5)      Ruminococcus gnavus et rel.   0.6996301 0.003192887
## 479    PC(40:3) Eubacterium cylindroides et rel.  -0.6771286 0.003800575
## 656    PC(40:3)                     Helicobacter  -0.6838424 0.003800575
## 1082   TG(50:4)      Ruminococcus gnavus et rel.   0.6852226 0.003800575
## 1086 TG(52:4).1      Ruminococcus gnavus et rel.   0.6716223 0.003800575
```


