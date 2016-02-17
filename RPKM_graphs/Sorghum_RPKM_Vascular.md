Vascular
========================================================

This is an R Markdown document. Markdown is a simple formatting syntax for authoring web pages (click the **MD** toolbar button for help on Markdown).

When you click the **Knit HTML** button a web page will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

Vascular CHG 
============

```r
library(grid)
library(ggplot2)

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.vascular.CHG.filtered.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#b9d2e2","#6aacd4","#4e91a5","#1d66a3"))
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png) 



```r
file_name = "/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/vascular_CHG_TSS_all.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#b9d2e2","#6aacd4","#4e91a5","#1d66a3"))
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png) 




```r
library(PMCMR)
kruskal.test(freq ~ mtype, data = x) 
```

```
## 
## 	Kruskal-Wallis rank sum test
## 
## data:  freq by mtype
## Kruskal-Wallis chi-squared = 5.2541, df = 3, p-value = 0.1541
```

```r
posthoc.kruskal.nemenyi.test(freq ~ mtype, data=x, dist="Tukey")
```

```
## 
## 	Pairwise comparisons using Tukey and Kramer (Nemenyi) test	
##                    with Tukey-Dist approximation for independent samples 
## 
## data:  freq by mtype 
## 
##   1    2    3   
## 2 0.21 -    -   
## 3 0.20 1.00 -   
## 4 0.78 0.75 0.74
## 
## P value adjustment method: none
```

Vascular CG 
============

```r
library(grid)
library(ggplot2)

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.vascular.CpG.filtered.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#F97976","#F82F2B","#B92320","#783A39"))
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png) 



```r
file_name = "/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/vascular_CG_TSS_all.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#F97976","#F82F2B","#B92320","#783A39"))
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png) 




```r
library(PMCMR)
kruskal.test(freq ~ mtype, data = x) 
```

```
## 
## 	Kruskal-Wallis rank sum test
## 
## data:  freq by mtype
## Kruskal-Wallis chi-squared = 12.0335, df = 3, p-value = 0.007269
```

```r
posthoc.kruskal.nemenyi.test(freq ~ mtype, data=x, dist="Tukey")
```

```
## 
## 	Pairwise comparisons using Tukey and Kramer (Nemenyi) test	
##                    with Tukey-Dist approximation for independent samples 
## 
## data:  freq by mtype 
## 
##   1     2     3    
## 2 0.984 -     -    
## 3 0.176 0.337 -    
## 4 0.014 0.040 0.764
## 
## P value adjustment method: none
```
Vascular CHH 
============

```r
library(grid)
library(ggplot2)

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.vascular.CHH.filtered.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#53CA77","#15C048","#0E8130","#1A4026"))
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png) 



```r
file_name = "/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/vascular_CHH_TSS_all.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#53CA77","#15C048","#0E8130","#1A4026"))
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png) 




```r
library(PMCMR)
kruskal.test(freq ~ mtype, data = x) 
```

```
## 
## 	Kruskal-Wallis rank sum test
## 
## data:  freq by mtype
## Kruskal-Wallis chi-squared = 19.0259, df = 3, p-value = 0.00027
```

```r
posthoc.kruskal.nemenyi.test(freq ~ mtype, data=x, dist="Tukey")
```

```
## 
## 	Pairwise comparisons using Tukey and Kramer (Nemenyi) test	
##                    with Tukey-Dist approximation for independent samples 
## 
## data:  freq by mtype 
## 
##   1       2       3      
## 2 0.90789 -       -      
## 3 0.45557 0.85499 -      
## 4 0.00028 0.00383 0.04651
## 
## P value adjustment method: none
```
