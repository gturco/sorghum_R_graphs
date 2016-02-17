

```r
library(reshape)
library(ggplot2)
library(PMCMR)

all_pvals = function(file_name,level)
{
  x = read.table(file_name, col.names=c("pos","freq","mtype"))
  x$mtype <- tolower(x$mtype)
  x$mtype <- factor(x$mtype, levels = c("root","shoot","vascular","nonvascular"))
  
  down = x[which(x$pos > 1000),]
  up =  x[which(x$pos < 0),]
  genebody = x[which(x$pos >= 0 & x$pos <= 1000),]
  
  d = kruskal.test(freq ~ mtype, data = up) 
  u = kruskal.test(freq ~ mtype, data = down) 
  g = kruskal.test(freq ~ mtype, data = genebody) 
  l <- cbind(d$p.value,g$p.value,u$p.value)
  colnames(l) <- c("upstream","genebody","downstream")
  rownames(l) <- c("t")
  y <- data.frame(l) 
  a <- level
  z <- cbind.data.frame(melt(y),a)
  return(z)

}
```

CHH
=========


```r
chh1 = all_pvals("/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/Tissue/CHH/CHH_1.txt",1)
```

```
## Using  as id variables
```

```r
chh2 = all_pvals("/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/Tissue/CHH/CHH_2.txt",2)
```

```
## Using  as id variables
```

```r
chh3 = all_pvals("/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/Tissue/CHH/CHH_3.txt",3)
```

```
## Using  as id variables
```

```r
chh4 = all_pvals("/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/Tissue/CHH/CHH_4.txt",4)
```

```
## Using  as id variables
```

```r
z <- rbind(chh1,chh2,chh3,chh4)
z$pval <- cut(z$value, breaks = c(0, 0.0001, 0.001, 0.01, 1))
ggplot(z, aes(variable, a, fill = z$pval)) + geom_tile(colour = "black") + scale_fill_brewer(palette = "Blues",direction=-1, limits=c("(0,0.0001]","(0.0001,0.01]","(0.001,0.01]")) + theme_classic() + theme(text=element_text(size=20))    +  
  theme(text=element_text(size=12), panel.margin = unit(0, "cm"), 
        plot.margin = unit(c(0, 0, 0.01, 0), "cm"),axis.ticks = element_line(size = 0.1), 
        axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + 
  xlab(NULL) + ylab(NULL) 
```

```
## Warning: `axis.ticks.margin` is deprecated. Please set `margin` property of
## `axis.text` instead
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

CHG
==================


```r
chh1 = all_pvals("/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/Tissue/CHG/CHG_1.txt",1)
```

```
## Using  as id variables
```

```r
chh2 = all_pvals("/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/Tissue/CHG/CHG_2.txt",2)
```

```
## Using  as id variables
```

```r
chh3 = all_pvals("/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/Tissue/CHG/CHG_3.txt",3)
```

```
## Using  as id variables
```

```r
chh4 = all_pvals("/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/Tissue/CHG/CHG_4.txt",4)
```

```
## Using  as id variables
```

```r
z <- rbind(chh1,chh2,chh3,chh4)
z$pval <- cut(z$value, breaks = c(0, 0.0001, 0.001, 0.01, 1))
ggplot(z, aes(variable, a, fill = z$pval)) + geom_tile(colour = "black") + scale_fill_brewer(palette = "Blues",direction=-1, limits=c("(0,0.0001]","(0.0001,0.01]","(0.001,0.01]")) + theme_classic() + theme(text=element_text(size=20))    +  
  theme(text=element_text(size=12), panel.margin = unit(0, "cm"), 
        plot.margin = unit(c(0, 0, 0.01, 0), "cm"),axis.ticks = element_line(size = 0.1), 
        axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + 
  xlab(NULL) + ylab(NULL) 
```

```
## Warning: `axis.ticks.margin` is deprecated. Please set `margin` property of
## `axis.text` instead
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)
CG
==========


```r
chh1 = all_pvals("/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/Tissue/CG/CG_1.txt",1)
```

```
## Using  as id variables
```

```r
chh2 = all_pvals("/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/Tissue/CG/CG_2.txt",2)
```

```
## Using  as id variables
```

```r
chh3 = all_pvals("/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/Tissue/CG/CG_3.txt",3)
```

```
## Using  as id variables
```

```r
chh4 = all_pvals("/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/Tissue/CG/CG_4.txt",4)
```

```
## Using  as id variables
```

```r
z <- rbind(chh1,chh2,chh3,chh4)
z$pval <- cut(z$value, breaks = c(0, 0.0001, 0.001, 0.01, 1))
ggplot(z, aes(variable, a, fill = z$pval)) + geom_tile(colour = "black") + scale_fill_brewer(palette = "Blues",direction=-1, limits=c("(0,0.0001]","(0.0001,0.01]","(0.001,0.01]")) + theme_classic() + theme(text=element_text(size=20))    +  
  theme(text=element_text(size=12), panel.margin = unit(0, "cm"), 
        plot.margin = unit(c(0, 0, 0.01, 0), "cm"),axis.ticks = element_line(size = 0.1), 
        axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + 
  xlab(NULL) + ylab(NULL) 
```

```
## Warning: `axis.ticks.margin` is deprecated. Please set `margin` property of
## `axis.text` instead
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

