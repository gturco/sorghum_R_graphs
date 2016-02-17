library(grid)
library(ggplot2)



file_name = "/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/all.vascular.CHH.all.freq.txt"

x = read.table(file_name, col.names=c("percentage","freq"))
ggplot(data=x, aes(x = percentage, y = log(freq))) + geom_line()
ggsave("/Users/gturco/Documents/Data/Sorg/PCA/root_shoot.pdf", width=5.60, height=3.2, dpi=600, units="cm")  


+ geom_point() + scale_size_continuous(range = c(0.01, 1))  +
  theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) 

ggsave("/Users/gturco/Documents/Data/Sorg/PCA/root_shoot.pdf", width=5.60, height=3.2, dpi=600, units="cm")  




file_name = "/Users/gturco/Downloads/scatter/root_shoot.txt"
x = read.table(file_name, col.names=c("Root","Shoot","freq"))
ggplot(data=x, aes(x = Root, y = Shoot ,size=log(freq))) + geom_point() + scale_size_continuous(range = c(0.01, 1))  +
 theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) 

ggsave("/Users/gturco/Documents/Data/Sorg/PCA/root_shoot.pdf", width=5.60, height=3.2, dpi=600, units="cm")  


file_name = "/Users/gturco/Downloads/scatter/vas_non.txt"
x = read.table(file_name, col.names=c("Vascular","Nonvascular","freq"))
ggplot(data=x, aes(x = Vascular, y = Nonvascular ,size=log(freq))) + geom_point() + scale_size_continuous(range = c(0.01, 1))  +
  theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) 

ggsave("/Users/gturco/Documents/Data/Sorg/PCA/vas_non.pdf", width=5.60, height=3.2, dpi=600, units="cm")  


file_name = "/Users/gturco/Downloads/scatter/vas_root.txt"
x = read.table(file_name, col.names=c("Vascular","Root","freq"))
ggplot(data=x, aes(x = Vascular, y = Root ,size=log(freq))) + geom_point() + scale_size_continuous(range = c(0.01, 1))  +
  theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) 

ggsave("/Users/gturco/Documents/Data/Sorg/PCA/vas_root.pdf", width=5.60, height=3.2, dpi=600, units="cm")  





vas = cbind(c(0,1,3,6,7694),"Vascular")
non = cbind(c(0,0,2,4,6379),"Nonvascular")
root =cbind(c(0,2,5,10,14847),"Root")
shoot =cbind(c(0,2,4,9,22275),"Shoot")
x <- rbind(vas,non,root,shoot)
colnames(x) <- c("Coverage","Tissue")
x <- data.frame(x) 
x$Coverage = as.numeric(as.character(x$Coverage) )

ggplot(x, aes(Tissue,Coverage, fill=Tissue)) +  geom_boxplot() + ylim(0,10) + scale_fill_manual(values=c("#c6141c","#fccb0c","#009987","#3261a8")) +
  theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("Coverage") + scale_colour_manual(values=c("#c6141c","#fccb0c","#009987","#3261a8"))

ggsave("/Users/gturco/Documents/Data/Sorg/coverage/cg_box.pdf", width=5.60, height=3.2, dpi=600, units="cm")  


vas = cbind(c(0,1,3,7,7450),"Vascular")
non = cbind(c(0,1,2,5,6141),"Nonvascular")
root =cbind(c(0,2,6,12,13922),"Root")
shoot =cbind(c(0,2,5,11,21194),"Shoot")
x <- rbind(vas,non,root,shoot)
colnames(x) <- c("Coverage","Tissue")
x <- data.frame(x) 
x$Coverage = as.numeric(as.character(x$Coverage) )

ggplot(x, aes(Tissue,Coverage, fill=Tissue)) +  geom_boxplot() + ylim(0,10) + scale_fill_manual(values=c("#c6141c","#fccb0c","#009987","#3261a8")) +
  theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab(NULL) + scale_colour_manual(values=c("#c6141c","#fccb0c","#009987","#3261a8"))

ggsave("/Users/gturco/Documents/Data/Sorg/coverage/chg_box.pdf", width=5.60, height=3.2, dpi=600, units="cm")  


vas = cbind(c(0,2,4,8,7720),"Vascular")
non = cbind(c(0,1,3,6,6357),"Nonvascular")
root =cbind(c(0,3,7,14,15270),"Root")
shoot =cbind(c(0,3,7,12,22724),"Shoot")
x <- rbind(vas,non,root,shoot)
colnames(x) <- c("Coverage","Tissue")
x <- data.frame(x) 
x$Coverage = as.numeric(as.character(x$Coverage) )

ggplot(x, aes(Tissue,Coverage, fill=Tissue)) +  geom_boxplot() + ylim(0,10) + scale_fill_manual(values=c("#c6141c","#fccb0c","#009987","#3261a8")) +
  theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab(NULL) + scale_colour_manual(values=c("#c6141c","#fccb0c","#009987","#3261a8"))

ggsave("/Users/gturco/Documents/Data/Sorg/coverage/chh_box.pdf", width=5.60, height=3.2, dpi=600, units="cm")  


####
file_name = "/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/shoot_CG_TSS_all.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#F97976","#F82F2B","#B92320","#783A39"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/shoot_CG_content.pdf",  width=5.60, height=3.2, dpi=600, units="cm") 


file_name = "/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/root_CG_TSS_all.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#F97976","#F82F2B","#B92320","#783A39"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/root_CG_content.pdf",  width=5.60, height=3.2, dpi=600, units="cm") 


file_name = "/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/vascular_CG_TSS_all.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#F97976","#F82F2B","#B92320","#783A39"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/vascular_CG_content.pdf",  width=5.60, height=3.2, dpi=600, units="cm") 

file_name = "/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/nonvascular_CG_TSS_all.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#F97976","#F82F2B","#B92320","#783A39"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/nonvascular_CG_content.pdf",  width=5.60, height=3.2, dpi=600, units="cm") 


### CHG

file_name = "/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/vascular_CHG_TSS_all.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#b9d2e2","#6aacd4","#4e91a5","#1d66a3"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/vascular_CHG_content.pdf",  width=5.60, height=3.2, dpi=600, units="cm") 

file_name = "/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/nonvascular_CHG_TSS_all.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#b9d2e2","#6aacd4","#4e91a5","#1d66a3"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/nonvascular_CHG_content.pdf",  width=5.60, height=3.2, dpi=600, units="cm") 


file_name = "/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/root_CHG_TSS_all.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#b9d2e2","#6aacd4","#4e91a5","#1d66a3"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/root_CHG_content.pdf",  width=5.60, height=3.2, dpi=600, units="cm") 

file_name = "/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/shoot_CHG_TSS_all.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#b9d2e2","#6aacd4","#4e91a5","#1d66a3"))
ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/shoot_CHG_content.pdf",  width=5.60, height=3.2, dpi=600, units="cm") 

#### CHH

file_name = "/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/shoot_CHH_TSS_all.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#53CA77","#15C048","#0E8130","#1A4026"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/shoot_CHH_content.pdf",  width=5.60, height=3.2, dpi=600, units="cm") 

file_name = "/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/root_CHH_TSS_all.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#53CA77","#15C048","#0E8130","#1A4026"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/root_CHH_content.pdf",  width=5.60, height=3.2, dpi=600, units="cm") 


file_name = "/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/nonvascular_CHH_TSS_all.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#53CA77","#15C048","#0E8130","#1A4026"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/nonvascular_CHH_content.pdf",  width=5.60, height=3.2, dpi=600, units="cm") 

file_name = "/Users/gturco/Documents/Data/Sorg/11_25_2015/dist_data/TSS_cs_plots/vascular_CHH_TSS_all.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#53CA77","#15C048","#0E8130","#1A4026"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/vascular_CHH_content.pdf",  width=5.60, height=3.2, dpi=600, units="cm") 



