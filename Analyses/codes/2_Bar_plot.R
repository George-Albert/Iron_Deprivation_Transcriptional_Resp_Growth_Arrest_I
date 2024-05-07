
{
    library(edgeR)
    library(qvalue)
    library(limma)
    library(ggplot2)
    library(cowplot)
    library(reshape2)
    library(DESeq2)
    library(fdrtool)
    library(ashr)
    library(dplyr)
    library(locfdr)
    library(readxl)
    library(writexl)
    library(eulerr)

}

################################
### 1. Load & pretty up data ###
################################

{
    feature_data=read.table("inputs/processed/unfiltered/txts/feature_data.txt")

    Iron_deprivation_effect_at_Exp=read.table("outputs/Data/txt/Iron_deprivation_effect_at_Exp.txt")
    Iron_deprivation_effect_at_Stat=read.table("outputs/Data/txt/Iron_deprivation_effect_at_Stat.txt")
    Stat_effect_with_Fe=read.table("outputs/Data/txt/Stat_effect_with_Fe.txt")
    Stat_effect_without_Fe=read.table("outputs/Data/txt/Stat_effect_without_Fe.txt")
    Stat_Fedeprivation_interaction=read.table("outputs/Data/txt/Stat_Fedeprivation_interaction.txt")

    ## Check congruence of stat files:

    length(which(rownames(Iron_deprivation_effect_at_Exp)!=rownames(Iron_deprivation_effect_at_Stat)))
    length(which(rownames(Iron_deprivation_effect_at_Exp)!=rownames(Stat_effect_with_Fe)))
    length(which(rownames(Iron_deprivation_effect_at_Exp)!=rownames(Stat_effect_without_Fe)))
    length(which(rownames(Iron_deprivation_effect_at_Exp)!=rownames(Stat_Fedeprivation_interaction)))

    ## cleanup featuredata
    feature_data$Type[which(feature_data$Type==".")]="misc_RNA"
    feature_data=feature_data[which(feature_data$Gene_ID %in% rownames(Iron_deprivation_effect_at_Exp)),]
    rownames(feature_data)=feature_data$Gene_ID
    feature_data=feature_data[order(tolower(rownames(feature_data))),]
    length(which(rownames(Iron_deprivation_effect_at_Exp)!=rownames(feature_data)))

}
threshold=0.01
UP_Iron_at_Exp  =length(which(Iron_deprivation_effect_at_Exp$BH<threshold & Iron_deprivation_effect_at_Exp$log2FoldChange_shrunken>0))
DOWN_Iron_at_Exp=length(which(Iron_deprivation_effect_at_Exp$BH<threshold & Iron_deprivation_effect_at_Exp$log2FoldChange_shrunken<0))

UP_Iron_at_Stat  =length(which(Iron_deprivation_effect_at_Stat$BH<threshold & Iron_deprivation_effect_at_Stat$log2FoldChange_shrunken>0))
DOWN_Iron_at_Stat=length(which(Iron_deprivation_effect_at_Stat$BH<threshold & Iron_deprivation_effect_at_Stat$log2FoldChange_shrunken<0))

UP_arrest_w_Fe  =length(which(Stat_effect_with_Fe$BH<threshold & Stat_effect_with_Fe$log2FoldChange_shrunken>0))
DOWN_arrest_w_Fe=length(which(Stat_effect_with_Fe$BH<threshold & Stat_effect_with_Fe$log2FoldChange_shrunken<0))

UP_arrest_wout_Fe  =length(which(Stat_effect_without_Fe$BH<threshold & Stat_effect_without_Fe$log2FoldChange_shrunken>0))
DOWN_arrest_wout_Fe=length(which(Stat_effect_without_Fe$BH<threshold & Stat_effect_without_Fe$log2FoldChange_shrunken<0))

up_gene <- c(UP_Iron_at_Exp,UP_Iron_at_Stat,UP_arrest_w_Fe,UP_arrest_wout_Fe)
down_gene <- -c(DOWN_Iron_at_Exp,DOWN_Iron_at_Stat,DOWN_arrest_w_Fe,DOWN_arrest_wout_Fe)
Samples <- c("(i) Fe deprivation at exp","(ii) Fe deprivation at stat",
             "(iii) Growth arrest at Fe/+ ","(iv) Growth arrest at Fe/-")
reg_table <- data.frame(Effects=rep(Samples,2),n_genes=c(up_gene,down_gene))
reg_table <- reg_table[order(reg_table$Effects),]
reg_table$Regulated <- c(rep(c("Activation","Repression"),4))
reg_table$Regulated=factor(reg_table$Regulated,levels=c("Repression","Activation"))

frame <- data.frame(reg_table)
sample <- factor(frame$Effects)
n_gene <- frame$n_genes
# Sample_Effect <- reorder(sample,abs(n_gene))
breaks_values <- pretty(n_gene)

bar_pl <- ggplot(frame,aes(x=sample, y=n_gene,fill=Regulated))+
geom_bar(stat = "identity",width = 0.8, color = "black")+
coord_flip()+
scale_y_continuous(breaks = breaks_values,labels = abs(breaks_values))+
#scale_x_discrete(expand = c(0,2.5))+
theme_classic()+
geom_hline(yintercept=0)+
#theme(aspect.ratio = .9)+
labs(x="", y = "Number of DEGs")+
scale_fill_manual(name="Regulated Genes",values = c("#50486D","#FFA373"))+
theme(
axis.text.x   = element_text(size=14),
axis.text.y   = element_text(size=14),
axis.title.x  = element_text(size=14),
axis.title.y  = element_text(size=14),
axis.ticks.y = element_blank(),
#panel.background = element_blank(),
#panel.grid.major = element_blank(),
#panel.grid.minor = element_blank(),
#axis.line = element_line(colour = "black"),
panel.border = element_rect(colour = "black", fill=NA, size=1,linetype="solid"),
legend.title=element_blank(),
legend.position="top",
#legend.text=element_text(size=14),
legend.key.size = unit(1, 'lines'))



pdf(file = "outputs/Figures/Fig_3B_Bar_plot_minimal.pdf",width=5,height=2.5)
print(bar_pl)
dev.off()

















