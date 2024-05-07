#############################
### 0. Load dependencies. ###
#############################

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
  library(readxl) 
  library(xlsx)
  library(ggridges)
  library(ggeasy)
}


############################
### 1. Declare functions ###
############################
dcols=function(x){data.frame(colnames(x))}
my_name <- function(v1) {
  deparse(substitute(v1))
}


#                         =================================
#                         === Bar plot: regulated genes ===
#                         =================================

##### Create Regulated table of down and upregulated genes ######
{
  d1 <- -dim(IDEE_reg[which(IDEE_reg$Regulated == "DOWN"),])[1]
  up1 <- dim(IDEE_reg[which(IDEE_reg$Regulated == "UP"),])[1]
  d2 <- -dim(IDES_reg[which(IDES_reg$Regulated == "DOWN"),])[1]
  up2 <- dim(IDES_reg[which(IDES_reg$Regulated == "UP"),])[1]
  d3 <- -dim(SEWFe_reg[which(SEWFe_reg$Regulated == "DOWN"),])[1]
  up3 <- dim(SEWFe_reg[which(SEWFe_reg$Regulated == "UP"),])[1]
  d4 <- -dim(SEWOFe_reg[which(SEWOFe_reg$Regulated == "DOWN"),])[1]
  up4 <- dim(SEWOFe_reg[which(SEWOFe_reg$Regulated == "UP"),])[1]
  d5 <- -dim(SFeDI_reg[which(SFeDI_reg$Regulated == "DOWN"),])[1]
  up5 <- dim(SFeDI_reg[which(SFeDI_reg$Regulated == "UP"),])[1]
}

up_gene <- c(up1,up2,up3,up4,up5)
down_gene <- c(d1,d2,d3,d4,d5)
Samples <- c("Iron Deprivation effects at EXP","Iron Deprivation effects at STAT",
             "Stat_effect_with_Fe","Stat_effect_without_Fe","Stat_Fe_deprivation_interaction")
reg_table <- data.frame(Effects=rep(Samples,2),n_genes=c(up_gene,down_gene))
reg_table <- reg_table[order(reg_table$Effects),]
reg_table$Regulated <- c(rep("UP",5),rep("DOWN",5))

write.table(reg_table,file="plot_tables/reg_table.txt")

bar <- function(data,title){
  frame <- data.frame(data)
  sample <- factor(frame$Effects)
  n_gene <- frame$n_genes
  # Sample_Effect <- reorder(sample,abs(n_gene))
  breaks_values <- pretty(n_gene) 
  
  bar_pl <- ggplot(frame,aes(x=sample, y=n_gene,fill=Regulated))+
    geom_bar(stat = "identity",width = 0.8, color = "black")+coord_flip()+
    scale_y_continuous(breaks = breaks_values,labels = abs(breaks_values))+
    scale_x_discrete(expand = c(0,2.5))+
    theme_classic()+ ggtitle(title)+geom_hline(yintercept=0)+theme(aspect.ratio = .9)+
    labs(x="Samples", y = "Number of Genes")+theme(plot.title = element_text(hjust = 0.5))+
    scale_fill_manual(name="Regulated Genes",values = c("#FFA373","#50486D"))
  return(bar_pl)
}

bar_plot <- bar(reg_table,"Regulated Genes")

pdf(file = "Outputs_def/Figures/Bar_plot.pdf",width=8,height=5)
print(bar_plot)
dev.off()

















