Sorghum miniheat maps
===================================


RPKM Functions
==================




```r
genebody_pval = function(file_name)
{
  
  x = read.table(file_name, col.names=c("pos","freq","mtype"))
  x$mtype <- tolower(x$mtype)
  x$mtype <- factor(x$mtype, levels = c(1,2,3,4))

  genebody = x[which(x$pos >= 0 & x$pos <= 1000),]
  d = kruskal.test(freq ~ mtype, data = genebody) 

  c <- posthoc.kruskal.nemenyi.test(freq ~ mtype, data=genebody, dist="Tukey")
  y <- data.frame(c$p.value)
  colnames(y) = c(1,2,3)
  ttype  <- c(2,3,4,2,3,4,2,3,4)

  z <- cbind.data.frame(melt(y),ttype)
  z$pval <- cut(z$value, breaks = c(-1, 0.0001, 0.001, 0.01, 1))
  z$ttype <- factor(z$ttype, levels = c(4,3,2))
  ggplot(z, aes(variable, ttype, fill = pval)) + geom_tile(colour = "black") + scale_fill_brewer(palette = "Blues",direction=-1, limits=c("(-1,0.0001]","(0.0001,0.01]","(0.001,0.01]")) + theme_classic() + theme(text=element_text(size=20))    +  
    theme(text=element_text(size=12), panel.margin = unit(0, "cm"), 
          plot.margin = unit(c(0, 0, 0.01, 0), "cm"),axis.ticks = element_line(size = 0.1), 
          axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + 
    xlab(NULL) + ylab(NULL) 
}
```


```r
upstream_pval = function(file_name)
{
  
  x = read.table(file_name, col.names=c("pos","freq","mtype"))
  x$mtype <- tolower(x$mtype)
  x$mtype <- factor(x$mtype, levels = c(1,2,3,4))

  up = x[which(x$pos < 0),]
  d = kruskal.test(freq ~ mtype, data = up)

  c <- posthoc.kruskal.nemenyi.test(freq ~ mtype, data=up, dist="Tukey")
  y <- data.frame(c$p.value)
  colnames(y) = c(1,2,3)
  ttype  <- c(2,3,4,2,3,4,2,3,4)

  z <- cbind.data.frame(melt(y),ttype)
  z$pval <- cut(z$value, breaks = c(-1, 0.0001, 0.001, 0.01, 1))
  z$ttype <- factor(z$ttype, levels = c(4,3,2))
  ggplot(z, aes(variable, ttype, fill = pval)) + geom_tile(colour = "black") + scale_fill_brewer(palette = "Blues",direction=-1, limits=c("(-1,0.0001]","(0.0001,0.01]","(0.001,0.01]")) + theme_classic() + theme(text=element_text(size=20))    +  
    theme(text=element_text(size=12), panel.margin = unit(0, "cm"), 
          plot.margin = unit(c(0, 0, 0.01, 0), "cm"),axis.ticks = element_line(size = 0.1), 
          axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + 
    xlab(NULL) + ylab(NULL) 
}
```


```r
downstream_pval  = function(file_name)
{
  
  x = read.table(file_name, col.names=c("pos","freq","mtype"))
  x$mtype <- tolower(x$mtype)
  x$mtype <- factor(x$mtype, levels = c(1,2,3,4))

  down = x[which(x$pos > 1000),]
  d = kruskal.test(freq ~ mtype, data = down) 

  c <- posthoc.kruskal.nemenyi.test(freq ~ mtype, data=down, dist="Tukey")
  y <- data.frame(c$p.value)
  colnames(y) = c(1,2,3)
  ttype  <- c(2,3,4,2,3,4,2,3,4)

  z <- cbind.data.frame(melt(y),ttype)
  z$pval <- cut(z$value, breaks = c(-1, 0.0001, 0.001, 0.01, 1))
  z$ttype <- factor(z$ttype, levels = c(4,3,2))
  ggplot(z, aes(variable, ttype, fill = pval)) + geom_tile(colour = "black") + scale_fill_brewer(palette = "Blues",direction=-1, limits=c("(-1,0.0001]","(0.0001,0.01]","(0.001,0.01]")) + theme_classic() + theme(text=element_text(size=20))    +  
    theme(text=element_text(size=12), panel.margin = unit(0, "cm"), 
          plot.margin = unit(c(0, 0, 0.01, 0), "cm"),axis.ticks = element_line(size = 0.1), 
          axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + 
    xlab(NULL) + ylab(NULL) 
}
```









## CG

### Shoot


```r
file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.shoot.CG.filtered.txt" 

uoutput_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/heatmap/shoot_CG_upstream.pdf"
goutput_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/heatmap/shoot_CG_genebody.pdf"
doutput_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/heatmap/shoot_CG_downstream.pdf"

upstream_pval(file_name)
```

```
## Error in upstream_pval(file_name): could not find function "posthoc.kruskal.nemenyi.test"
```

```r
ggsave(uoutput_name, width=1.8, height=1.5, dpi=600, units="cm") 
```

```
## Error in eval(expr, envir, enclos): could not find function "ggsave"
```

```r
genebody_pval(file_name)
```

```
## Error in genebody_pval(file_name): could not find function "posthoc.kruskal.nemenyi.test"
```

```r
ggsave(uoutput_name, width=1.8, height=1.5, dpi=600, units="cm")
```

```
## Error in eval(expr, envir, enclos): could not find function "ggsave"
```

```r
downstream_pval(file_name)
```

```
## Error in downstream_pval(file_name): could not find function "posthoc.kruskal.nemenyi.test"
```

```r
gsave(doutput_name, width=1.8, height=1.5, dpi=600, units="cm")
```

```
## Error in eval(expr, envir, enclos): could not find function "gsave"
```


### Root


```r
file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.root.CG.filtered.txt" 

uoutput_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/heatmap/root_CG_upstream.pdf"
goutput_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/heatmap/root_CG_genebody.pdf"
doutput_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/heatmap/root_CG_downstream.pdf"

rupstream_pval(file_name)
```

```
## Error in eval(expr, envir, enclos): could not find function "rupstream_pval"
```

```r
ggsave(uoutput_name, width=1.8, height=1.5, dpi=600, units="cm") 
```

```
## Error in eval(expr, envir, enclos): could not find function "ggsave"
```

```r
rgenebody_pval(file_name)
```

```
## Error in eval(expr, envir, enclos): could not find function "rgenebody_pval"
```

```r
ggsave(uoutput_name, width=1.8, height=1.5, dpi=600, units="cm")
```

```
## Error in eval(expr, envir, enclos): could not find function "ggsave"
```

```r
rdownstream_pval(file_name)
```

```
## Error in eval(expr, envir, enclos): could not find function "rdownstream_pval"
```

```r
gsave(doutput_name, width=1.8, height=1.5, dpi=600, units="cm")
```

```
## Error in eval(expr, envir, enclos): could not find function "gsave"
```

## Vascular


```r
file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.vascular.CpG.filtered.txt" 

uoutput_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/heatmap/vas_CG_upstream.pdf"
goutput_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/heatmap/vas_CG_genebody.pdf"
doutput_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/heatmap/vas_CG_downstream.pdf"

rupstream_pval(file_name)
```

```
## Error in eval(expr, envir, enclos): could not find function "rupstream_pval"
```

```r
ggsave(uoutput_name, width=1.8, height=1.5, dpi=600, units="cm") 
```

```
## Error in eval(expr, envir, enclos): could not find function "ggsave"
```

```r
rgenebody_pval(file_name)
```

```
## Error in eval(expr, envir, enclos): could not find function "rgenebody_pval"
```

```r
ggsave(uoutput_name, width=1.8, height=1.5, dpi=600, units="cm")
```

```
## Error in eval(expr, envir, enclos): could not find function "ggsave"
```

```r
rdownstream_pval(file_name)
```

```
## Error in eval(expr, envir, enclos): could not find function "rdownstream_pval"
```

```r
gsave(doutput_name, width=1.8, height=1.5, dpi=600, units="cm")
```

```
## Error in eval(expr, envir, enclos): could not find function "gsave"
```

## Nonvascular


```r
file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.nonvascular.CG.filtered.txt" 

uoutput_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/heatmap/non_CG_upstream.pdf"
goutput_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/heatmap/non_CG_genebody.pdf"
doutput_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/heatmap/non_CG_downstream.pdf"

rupstream_pval(file_name)
```

```
## Error in eval(expr, envir, enclos): could not find function "rupstream_pval"
```

```r
ggsave(uoutput_name, width=1.8, height=1.5, dpi=600, units="cm") 
```

```
## Error in eval(expr, envir, enclos): could not find function "ggsave"
```

```r
rgenebody_pval(file_name)
```

```
## Error in eval(expr, envir, enclos): could not find function "rgenebody_pval"
```

```r
ggsave(uoutput_name, width=1.8, height=1.5, dpi=600, units="cm")
```

```
## Error in eval(expr, envir, enclos): could not find function "ggsave"
```

```r
rdownstream_pval(file_name)
```

```
## Error in eval(expr, envir, enclos): could not find function "rdownstream_pval"
```

```r
gsave(doutput_name, width=1.8, height=1.5, dpi=600, units="cm")
```

```
## Error in eval(expr, envir, enclos): could not find function "gsave"
```



## CHG


## CHH


