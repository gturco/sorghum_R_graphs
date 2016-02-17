shoot = cbind("shoot",1:4361)
root = cbind("root",1:4361)
vas = cbind("vas",1:4361)
nonvas = cbind("nonvas",1:4361)
df = rbind(shoot, root,vas,nonvas)

shoot = rbind(data.frame(genes= rep(1,1972)), data.frame(genes= rep(0,2389)),deparse.level = 0)
root = rbind(data.frame(genes=rep(0,1972)), data.frame(genes=rep(2,1995)), data.frame(genes=rep(0,394)))
vascular = rbind(data.frame(genes=rep(0,1804)), data.frame(genes=rep(3,168)), data.frame(genes=rep(0,1351)),data.frame(genes=rep(3,644)),data.frame(genes=rep(3,225)),data.frame(genes=rep(0,169)))
non_vascular=rbind(data.frame(genes=rep(0,1625)), data.frame(genes=rep(4,179)), data.frame(genes=rep(0,168)), data.frame(genes=rep(0,1093)),data.frame(genes=rep(4,258)),data.frame(genes=rep(0,644)),data.frame(genes=rep(0,225)),data.frame(genes=rep(4,169))) 
               
rnaseq= rbind(shoot,root,vascular,non_vascular)  
rnaseq_data=cbind(df,rnaseq)
colnames(rnaseq_data) <- c("tissue_type","genes","score")

rnaseq_data$genes <- factor(rnaseq_data$genes,levels = 1:4361)
rnaseq_data$tissue_type <- factor(rnaseq_data$tissue_type,levels = c("nonvas","vas","root","shoot"))
rnaseq_data$score <- factor(rnaseq_data$score,levels = c(0,1,2,3,4))
library(ggplot2)

col <- c("white", "#D55E00", "#0072B2", "#56B4E9", "#F0E442")
qplot(genes,tissue_type, fill=col[score], data=rnaseq_data, geom="tile") + scale_fill_identity()

col <- c("white", "#D55E00", "#0072B2", "#56B4E9", "#F0E442")