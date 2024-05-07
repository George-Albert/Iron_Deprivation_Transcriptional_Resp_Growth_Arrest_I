tab=read.table("outputs/data/Rodriguez/Fold_change_table.txt")
exp=read.table("outputs/data/Rodriguez/stats_tables/Iron_deprivation_effect_at_Exp.txt")
stat=read.table("outputs/data/Rodriguez/stats_tables/Iron_deprivation_effect_at_Stat.txt")
rodri=read.table("outputs/data/Rodriguez/rodriguez_clean.txt",header=TRUE)
rownames(tab)=tab$Name

## Primero JAC hizo una version preliminar, aquí estamos cargando sus inputs y lo rehecho para confirmar el match
dupes=rodri$RV_Number[which(duplicated(rodri$RV_Number))]

rodri=rodri[-which(rodri$RV_Number %in% dupes),]
rownames(rodri)=rodri$RV_Number

rodri=rodri[which(rownames(rodri) %in% rownames(exp)),]
rodri$Fold_Change=log2(rodri$Fold_Change)

tab=tab[which(rownames(tab) %in% rownames(rodri)),]

tab=tab[order(rownames(tab)),]
rodri=rodri[order(rownames(rodri)),]

length(which(rownames(tab)!=rownames(rodri)))
## Rodri es un perfect match con tab


exp=exp[which(rownames(exp) %in% rownames(tab)),]
stat=stat[which(rownames(stat) %in% rownames(tab)),]

tab=tab[order(rownames(tab)),]
exp=exp[order(rownames(exp)),]
stat=stat[order(rownames(stat)),]

tab$raw_exp=exp$log2FoldChange_raw
tab$shrunk_exp=exp$log2FoldChange_shrunken
tab$chisq_exp=exp$stat
tab$raw_stat=stat$log2FoldChange_raw
tab$shrunk_stat=stat$log2FoldChange_shrunken
tab$chisq_stat=stat$stat
tab$Gene=exp$Gene
tab$Gene[which(rownames(tab)=="Rv2377c")]="mbtH"
tab$Gene[which(rownames(tab)=="Rv3402c")]="Rv3402c"


exp$narrow_BH=p.adjust(exp$pvalue,method="BH")

up_by_us=rownames(exp)[which(exp$log2FoldChange_shrunken>0 & exp$BH<0.01)]
down_by_us=rownames(exp)[which(exp$log2FoldChange_shrunken<0 & exp$BH<0.01)]

tab$label=""
tab[up_by_us,"label"]=tab[up_by_us,"Gene"]
tab[down_by_us,"label"]=tab[down_by_us,"Gene"]


tab$Condition="Background"
tab[up_by_us,"Condition"]="Upregulated"
tab[down_by_us,"Condition"]="Downregulated"



library(ggplot2)
library(ggrepel)

pl=ggplot(tab)+
geom_point(aes(x=Fold_change,y=Fold_change_us,color=Condition))+
geom_text_repel(aes(x=Fold_change,y=Fold_change_us,label=label))+
theme_classic()+
theme(
axis.text.y   = element_text(size=14),
axis.text.x   = element_text(size=14),
axis.title.y  = element_text(size=14),
axis.title.x  = element_text(size=14),
#panel.background = element_blank(),
#panel.grid.major = element_blank(),
#panel.grid.minor = element_blank(),
#axis.line = element_line(colour = "black"),
panel.border = element_rect(colour = "black", fill=NA, size=1,linetype="solid"),
#legend.title=element_blank(),
legend.position="none",
#legend.text=element_text(size=14),
legend.key.size = unit(1, 'lines'))+xlab("log2FC Rodriguez et al.")+ylab("Raw log2FC Iron deprivation @ Exp5")+ggtitle("Pearson r=0.49  p=8.4E-18")+labs(color="")


test=cor.test(tab$Fold_change,tab$raw_exp)
test$p.value
#8.383047e-18
exp=read.table("outputs/Data/txt/Iron_deprivation_effect_at_Exp.txt")

plot(tab$Fold_change,tab$chisq_exp)

pdf("outputs/Figures/Fig_S7_scatter.pdf",height=6,width=6)
print(pl)
dev.off()

tab_1=data.frame(Effect="Rodriguez",logFC=abs(tab$Fold_change))
tab_2=data.frame(Effect="Iron at Exp",logFC=abs(tab$raw_exp))
tab_3=data.frame(Effect="Iron at Stat",logFC=abs(tab$raw_stat))
tabm=rbind(tab_1,tab_2,tab_3)

wilcox.test(c(101,103,123,10009,1008,18976),c(1,2,3,0,4,5),alt="greater")

wilcox.test(tab$raw_stat,tab$raw_exp,alt="greater")
wilcox.test(tab$raw_stat,tab$Fold_change,alt="greater")
