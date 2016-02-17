## one way anova- when do you uses one way anova vs other types...

fit <- aov(y ~ A, data=mydataframe)
## eq is response ~ factor
### methylation ~ tissue_type
file_name = "/Users/gturco/Documents/Data/Sorg/anova/ALL_CHG_meth.txt"
x = read.table(file_name, col.names=c("freq","ctype","mtype"), sep="\t", header= FALSE)
x$mtype <- factor(x$mtype, levels = c("Root","Shoot","Vascular","Nonvascular"))

mRNA = x[x$ctype == "miRNA",]
cds = x[x$ctype == "CDs",]
te = x[x$ctype == "TE",]
td = x[x$ctype == "tandem duplicates",]
five = x[x$ctype == "5 promoter",]
three = x[x$ctype == "3 promoter",]

###CHG
kruskal.test(Ozone ~ Month, data = airquality) 
kruskal.test(freq ~ mtype, data=mRNA)
##NOTHING!!!

kruskal.test(freq ~ mtype, data=cds)
Vascular-Root         0.0080059093  0.0021662594  0.013845559 0.0024182

geom_errorbar(limits, position=dodge, width=0.25)


kruskal.test(freq ~ mtype, data=te)
posthoc.kruskal.nemenyi.test(freq ~ mtype, data=te, dist="Tukey")

diff         lwr         upr p adj
Shoot-Root            0.01578245  0.01086903  0.02069587     0
Vascular-Root        -0.01562093 -0.02082660 -0.01041526     0
Nonvascular-Root     -0.05940051 -0.06480964 -0.05399138     0
Vascular-Shoot       -0.03140338 -0.03663158 -0.02617518     0
Nonvascular-Shoot    -0.07518296 -0.08061378 -0.06975215     0
Nonvascular-Vascular -0.04377958 -0.04947616 -0.03808300     0

fit <- aov(freq ~ mtype, data=td)
summary(fit)
posthoc <- TukeyHSD(x=fit, 'mtype', conf.level=0.95)
Vascular-Root         0.014081086  0.004093830 0.024068342 0.0016652


fit <- aov(freq ~ mtype, data=five)
summary(fit)
posthoc <- TukeyHSD(x=fit, 'mtype', conf.level=0.95)
Vascular-Root        -0.019723284 -0.0265222999 -0.012924269 0.0000000
Nonvascular-Root     -0.016495460 -0.0235953770 -0.009395543 0.0000000
Vascular-Shoot       -0.025401862 -0.0322132839 -0.018590441 0.0000000
Nonvascular-Shoot    -0.022174038 -0.0292858362 -0.015062240 0.0000000


fit <- aov(freq ~ mtype, data=three)
summary(fit)
posthoc <- TukeyHSD(x=fit, 'mtype', conf.level=0.95)
Shoot-Root            0.0087721647  0.003016040  0.014528289 0.0005238
Vascular-Shoot       -0.0090669942 -0.014877799 -0.003256189 0.0003556
Nonvascular-Shoot    -0.0113143835 -0.017243924 -0.005384843 0.0000057

#CHH
#CG

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}



file_name = "/Users/gturco/Documents/Data/Sorg/anova/ALL_CHG_meth.txt"
x = read.table(file_name, col.names=c("freq","ctype","mtype"), sep="\t", header= FALSE)
x$mtype <- factor(x$mtype, levels = c("Root","Shoot","Vascular","Nonvascular"))
x$ctype <- factor(x$ctype, levels = c("CDs", "3 promoter","5 promoter","miRNA", "TE","tandem duplicates")) 
limits <- aes(ymax = resp + se, ymin=resp - se)
ggplot(data=x, aes(x = ctype, y = freq* 100 ,group=mtype, fill = mtype)) +  geom_bar(stat = "summary", fun.y = "mean",position = "dodge") + geom_errorbar(stat = "summary", fun.ymax, fun.ymin, position="dodge", width=0.25) + scale_fill_manual(values=c("#e1131e","#3c4a96","#009987","#71b436")) + theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("") + ylim(0,80)

ggsave("/Users/gturco/Documents/Data/Sorg/meth_bar/CG.pdf",  width=5.60, height=3.2, dpi=600, units="cm") 


library(PMCMR)

### CG
### SAVE WITH REGIONS

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.shoot.CG.filtered.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))

kruskal.test(freq ~ mtype, data = x) 
Kruskal-Wallis chi-squared = 24.4949, df = 3, p-value = 1.969e-05
posthoc.kruskal.nemenyi.test(freq ~ mtype, data=x, dist="Tukey")
1       2       3      
2 0.00533 -       -      
3 0.00056 0.93257 -      
4 3.4e-05 0.60283 0.91932

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.root.CG.filtered.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
kruskal.test(freq ~ mtype, data = x) 
Kruskal-Wallis chi-squared = 3.2469, df = 3, p-value = 0.3551


file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.nonvascular.CG.filtered.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
kruskal.test(freq ~ mtype, data = x) 
Kruskal-Wallis chi-squared = 21.0687, df = 3, p-value = 0.0001019
posthoc.kruskal.nemenyi.test(freq ~ mtype, data=x, dist="Tukey")
data:  freq by mtype 

1       2       3      
2 0.16145 -       -      
3 0.00203 0.44036 -      
4 0.00014 0.13841 0.91496

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.vascular.CpG.filtered.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
kruskal.test(freq ~ mtype, data = x) 
Kruskal-Wallis chi-squared = 28.8327, df = 3, p-value = 2.428e-06
posthoc.kruskal.nemenyi.test(freq ~ mtype, data=x, dist="Tukey")
1       2       3      
2 0.25833 -       -      
3 0.00043 0.14349 -      
4 6.1e-06 0.01205 0.79144


######## CHG
file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.shoot.filtered.CHG.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))

kruskal.test(freq ~ mtype, data = x) 
Kruskal-Wallis chi-squared = 4.114, df = 3, p-value = 0.2494
posthoc.kruskal.nemenyi.test(freq ~ mtype, data=x, dist="Tukey")


file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.root.filtered.CHG.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
kruskal.test(freq ~ mtype, data = x) 
Kruskal-Wallis chi-squared = 2.3017, df = 3, p-value = 0.5122

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.nonvascular.filtered.CHG.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
kruskal.test(freq ~ mtype, data = x) 
Kruskal-Wallis chi-squared = 11.1785, df = 3, p-value = 0.0108
posthoc.kruskal.nemenyi.test(freq ~ mtype, data=x, dist="Tukey")
data:  freq by mtype 

1     2     3    
2 0.897 -     -    
  3 0.635 0.233 -    
  4 0.076 0.010 0.616

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.vascular.CHG.filtered.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
kruskal.test(freq ~ mtype, data = x) 
Kruskal-Wallis chi-squared = 9.3911, df = 3, p-value = 0.02452
posthoc.kruskal.nemenyi.test(freq ~ mtype, data=x, dist="Tukey")
1    2    3   
2 0.33 -    -   
3 0.77 0.89 -   
4 0.62 0.02 0.13

###### CHH
file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.shoot.CHH.filtered.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))

kruskal.test(freq ~ mtype, data = x) 
Kruskal-Wallis chi-squared = 20.2697, df = 3, p-value = 0.0001492
posthoc.kruskal.nemenyi.test(freq ~ mtype, data=x, dist="Tukey")

1       2       3      
2 0.00014 -       -      
  3 0.00321 0.86278 -      
  4 0.03155 0.44390 0.89412

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.root.CHH.filtered.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
kruskal.test(freq ~ mtype, data = x) 
Kruskal-Wallis chi-squared = 5.2708, df = 3, p-value = 0.153

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.nonvascular.filtered.CHH.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
kruskal.test(freq ~ mtype, data = x) 
Kruskal-Wallis chi-squared = 6.2707, df = 3, p-value = 0.09916


file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.vascular.CHH.filtered.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
kruskal.test(freq ~ mtype, data = x) 
Kruskal-Wallis chi-squared = 6.3619, df = 3, p-value = 0.09527

#########TISSUE AND RPKM

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CHH_all_1.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))

kruskal.test(freq ~ mtype, data = x) 
posthoc.kruskal.nemenyi.test(freq ~ mtype, data=x, dist="Tukey")
