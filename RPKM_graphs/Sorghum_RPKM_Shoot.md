SHOOT
========================================================

TODO:

Shoot CHG 
============

```r
library(grid)
library(ggplot2)

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.shoot.filtered.CHG.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + 
scale_color_manual(values = c("#b9d2e2","#6aacd4","#4e91a5","#1d66a3"))
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png) 



```r
file_name = "/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/shoot_CHG_TSS_all.txt"
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
## Kruskal-Wallis chi-squared = 6.3975, df = 3, p-value = 0.09379
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
## 2 0.221 -     -    
## 3 0.092 0.977 -    
## 4 0.735 0.809 0.562
## 
## P value adjustment method: none
```

Shoot CG 
============

```r
library(grid)
library(ggplot2)

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.shoot.CG.filtered.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#F97976","#F82F2B","#B92320","#783A39"))
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png) 



```r
file_name = "/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/shoot_CG_TSS_all.txt"
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
## Kruskal-Wallis chi-squared = 12.5666, df = 3, p-value = 0.005674
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
##   1      2      3     
## 2 0.0391 -      -     
## 3 0.0709 0.9959 -     
## 4 0.0053 0.9177 0.8203
## 
## P value adjustment method: none
```
Shoot CHH 
============

```r
library(grid)
library(ggplot2)

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.shoot.CHH.filtered.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#53CA77","#15C048","#0E8130","#1A4026"))
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png) 



```r
file_name = "/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/shoot_CHH_TSS_all.txt"
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
## Kruskal-Wallis chi-squared = 30.9519, df = 3, p-value = 8.701e-07
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
## 2 0.01736 -       -      
## 3 0.00015 0.58021 -      
## 4 9.3e-07 0.09639 0.72942
## 
## P value adjustment method: none
```
