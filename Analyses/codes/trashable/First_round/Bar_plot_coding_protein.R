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
#                         =================================
#                         === Bar plot: regulated genes ===
#                         =================================

##### Create Regulated protein table of down and upregulated genes ######
{
  d1 <- -dim(IDEE_protein[which(IDEE_protein$Regulated == "DOWN"),])[1]
  up1 <- dim(IDEE_protein[which(IDEE_protein$Regulated == "UP"),])[1]
  d2 <- -dim(IDES_protein[which(IDES_protein$Regulated == "DOWN"),])[1]
  up2 <- dim(IDES_protein[which(IDES_protein$Regulated == "UP"),])[1]
  d3 <- -dim(SEWFe_protein[which(SEWFe_protein$Regulated == "DOWN"),])[1]
  up3 <- dim(SEWFe_protein[which(SEWFe_protein$Regulated == "UP"),])[1]
  d4 <- -dim(SEWOFe_protein[which(SEWOFe_protein$Regulated == "DOWN"),])[1]
  up4 <- dim(SEWOFe_protein[which(SEWOFe_protein$Regulated == "UP"),])[1]
  d5 <- -dim(SFeDI_protein[which(SFeDI_protein$Regulated == "DOWN"),])[1]
  up5 <- dim(SFeDI_protein[which(SFeDI_protein$Regulated == "UP"),])[1]
}

up_gene <- c(up1,up2,up3,up4,up5)
down_gene <- c(d1,d2,d3,d4,d5)
Samples <- c("Iron Deprivation effects at EXP","Iron Deprivation effects at STAT",
             "Stat_effect_with_Fe","Stat_effect_without_Fe","Stat_Fe_deprivation_interaction")
reg_protein_table <- data.frame(Effects=rep(Samples,2),n_genes=c(up_gene,down_gene))
reg_protein_table <- reg_protein_table[order(reg_protein_table$Effects),]
reg_protein_table$Regulated <- c(rep(c("UP","DOWN"),(nrow(reg_protein_table)/2)))

write.table(reg_protein_table,file="Outputs_def/plot_tables/reg_protein_table.txt")

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

bar_plot <- bar(reg_protein_table,"Regulated Genes \n Protein Coding Genes")

pdf(file = "Outputs_def/Figures/Bar_plot_protein_genes.pdf",width=8,height=5)
print(bar_plot)
dev.off()

##### Create Regulated non protein coding table of down and upregulated genes ######


{
  d1 <- -dim(IDEE_no_protein[which(IDEE_no_protein$Regulated == "DOWN"),])[1]
  up1 <- dim(IDEE_no_protein[which(IDEE_no_protein$Regulated == "UP"),])[1]
  
  d2_1 <- -dim(IDES_no_protein[which(IDES_no_protein$Regulated == "DOWN" &
                                     IDES_no_protein$Type=="tRNA"),])[1]
  d2_2 <- -dim(IDES_no_protein[which(IDES_no_protein$Regulated == "DOWN" &
                                     IDES_no_protein$Type=="rRNA"),])[1]
  d2_3 <- -dim(IDES_no_protein[which(IDES_no_protein$Regulated == "DOWN" &
                                     IDES_no_protein$Type=="misc_RNA"),])[1]
  d2_4 <- -dim(IDES_no_protein[which(IDES_no_protein$Regulated == "DOWN" &
                                       IDES_no_protein$Type=="nc_RNA"),])[1]
  
  u2_1 <- dim(IDES_no_protein[which(IDES_no_protein$Regulated == "UP" &
                                       IDES_no_protein$Type=="tRNA"),])[1]
  u2_2 <- dim(IDES_no_protein[which(IDES_no_protein$Regulated == "UP" &
                                       IDES_no_protein$Type=="rRNA"),])[1]
  u2_3 <- dim(IDES_no_protein[which(IDES_no_protein$Regulated == "UP" &
                                       IDES_no_protein$Type=="misc_RNA"),])[1]
  u2_4 <- dim(IDES_no_protein[which(IDES_no_protein$Regulated == "UP" &
                                      IDES_no_protein$Type=="nc_RNA"),])[1]
  
  d3_1 <- -dim(SEWFe_no_protein[which(SEWFe_no_protein$Regulated == "DOWN" &
                                        SEWFe_no_protein$Type=="tRNA"),])[1]
  d3_2 <- -dim(SEWFe_no_protein[which(SEWFe_no_protein$Regulated == "DOWN" &
                                        SEWFe_no_protein$Type=="rRNA"),])[1]
  d3_3 <- -dim(SEWFe_no_protein[which(SEWFe_no_protein$Regulated == "DOWN" &
                                        SEWFe_no_protein$Type=="misc_RNA"),])[1]
  d3_4 <- -dim(SEWFe_no_protein[which(SEWFe_no_protein$Regulated == "DOWN" &
                                        SEWFe_no_protein$Type=="nc_RNA"),])[1]
  
  u3_1 <- dim(SEWFe_no_protein[which(SEWFe_no_protein$Regulated == "UP" &
                                        SEWFe_no_protein$Type=="tRNA"),])[1]
  u3_2 <- dim(SEWFe_no_protein[which(SEWFe_no_protein$Regulated == "UP" &
                                        SEWFe_no_protein$Type=="rRNA"),])[1]
  u3_3 <- dim(SEWFe_no_protein[which(SEWFe_no_protein$Regulated == "UP" &
                                        SEWFe_no_protein$Type=="misc_RNA"),])[1]
  u3_4 <- dim(SEWFe_no_protein[which(SEWFe_no_protein$Regulated == "UP" &
                                       SEWFe_no_protein$Type=="nc_RNA"),])[1]
  
  d4_1 <- -dim(SEWOFe_no_protein[which(SEWOFe_no_protein$Regulated == "DOWN" &
                                        SEWOFe_no_protein$Type=="tRNA"),])[1]
  d4_2 <- -dim(SEWOFe_no_protein[which(SEWOFe_no_protein$Regulated == "DOWN" &
                                        SEWOFe_no_protein$Type=="rRNA"),])[1]
  d4_3 <- -dim(SEWOFe_no_protein[which(SEWOFe_no_protein$Regulated == "DOWN" &
                                        SEWOFe_no_protein$Type=="misc_RNA"),])[1]
  d4_4 <- -dim(SEWOFe_no_protein[which(SEWOFe_no_protein$Regulated == "DOWN" &
                                         SEWOFe_no_protein$Type=="nc_RNA"),])[1]
  
  u4_1 <- dim(SEWOFe_no_protein[which(SEWOFe_no_protein$Regulated == "UP" &
                                        SEWOFe_no_protein$Type=="tRNA"),])[1]
  u4_2 <- dim(SEWOFe_no_protein[which(SEWOFe_no_protein$Regulated == "UP" &
                                        SEWOFe_no_protein$Type=="rRNA"),])[1]
  u4_3 <- dim(SEWOFe_no_protein[which(SEWOFe_no_protein$Regulated == "UP" &
                                        SEWOFe_no_protein$Type=="misc_RNA"),])[1]
  u4_4 <- dim(SEWOFe_no_protein[which(SEWOFe_no_protein$Regulated == "UP" &
                                        SEWOFe_no_protein$Type=="nc_RNA"),])[1]
  
  d5_1 <- -dim(SFeDI_no_protein[which(SFeDI_no_protein$Regulated == "DOWN" &
                                        SFeDI_no_protein$Type=="tRNA"),])[1]
  d5_2 <- -dim(SFeDI_no_protein[which(SFeDI_no_protein$Regulated == "DOWN" &
                                        SFeDI_no_protein$Type=="rRNA"),])[1]
  d5_3 <- -dim(SFeDI_no_protein[which(SFeDI_no_protein$Regulated == "DOWN" &
                                        SFeDI_no_protein$Type=="misc_RNA"),])[1]
  d5_4 <- -dim(SFeDI_no_protein[which(SFeDI_no_protein$Regulated == "DOWN" &
                                        SFeDI_no_protein$Type=="nc_RNA"),])[1]
  
  u5_1 <- dim(SFeDI_no_protein[which(SFeDI_no_protein$Regulated == "UP" &
                                        SFeDI_no_protein$Type=="tRNA"),])[1]
  u5_2 <- dim(SFeDI_no_protein[which(SFeDI_no_protein$Regulated == "UP" &
                                        SFeDI_no_protein$Type=="rRNA"),])[1]
  u5_3 <- dim(SFeDI_no_protein[which(SFeDI_no_protein$Regulated == "UP" &
                                        SFeDI_no_protein$Type=="misc_RNA"),])[1]
  u5_4 <- dim(SFeDI_no_protein[which(SFeDI_no_protein$Regulated == "UP" &
                                       SFeDI_no_protein$Type=="nc_RNA"),])[1]
}

up_gene <- c(up1,u2_1,u2_2,u2_3,u2_4,u3_1,u3_2,u3_3,u3_4,u4_1,u4_2,u4_3,u4_4,
             u5_1,u5_2,u5_3,u5_4)
down_gene <- c(d1,d2_1,d2_2,d2_3,d2_4,d3_1,d3_2,d3_3,d3_4,d4_1,d4_2,d4_3,d4_4,
               d5_1,d5_2,d5_3,d5_4)
Samples <- c(rep("Iron Deprivation effects at EXP",1),
             rep("Iron Deprivation effects at STAT",4),
             rep("Stat_effect_with_Fe",4),
             rep("Stat_effect_without_Fe",4),
             rep("Stat_Fe_deprivation_interaction",4))
{
reg_no_protein_table <- data.frame(Effects=rep(Samples,2),n_genes=c(up_gene,down_gene))
reg_no_protein_table <- reg_no_protein_table[order(reg_no_protein_table$Effects),]
reg_no_protein_table$Regulated <- c("UP","DOWN",rep(c(c(rep("UP",4),rep("DOWN",4))),4))
reg_no_protein_table$Type <- c("","",rep(c(c(c("tRNA","rRNA","misc_RNA","nc_RNA"),
                                               c("tRNA","rRNA","misc_RNA","nc_RNA"))),4))

reg_no_protein_table$Type[which(reg_no_protein_table$Regulated=="UP" &
                                  reg_no_protein_table$Type=="tRNA")] <- "tRNA_UP"
reg_no_protein_table$Type[which(reg_no_protein_table$Regulated=="DOWN" &
                                  reg_no_protein_table$Type=="tRNA")] <- "tRNA_DOWN"
reg_no_protein_table$Type[which(reg_no_protein_table$Regulated=="UP" &
                                  reg_no_protein_table$Type=="rRNA")] <- "rRNA_UP"
reg_no_protein_table$Type[which(reg_no_protein_table$Regulated=="DOWN" &
                                  reg_no_protein_table$Type=="rRNA")] <- "rRNA_DOWN"
reg_no_protein_table$Type[which(reg_no_protein_table$Regulated=="UP" &
                                  reg_no_protein_table$Type=="misc_RNA")] <- "misc_RNA_UP"
reg_no_protein_table$Type[which(reg_no_protein_table$Regulated=="DOWN" &
                                  reg_no_protein_table$Type=="misc_RNA")] <- "misc_RNA_DOWN"
reg_no_protein_table$Type[which(reg_no_protein_table$Regulated=="UP" &
                                  reg_no_protein_table$Type=="nc_RNA")] <- "nc_RNA_UP"
reg_no_protein_table$Type[which(reg_no_protein_table$Regulated=="DOWN" &
                                  reg_no_protein_table$Type=="nc_RNA")] <- "nc_RNA_DOWN"

reg_no_protein_table$Type <- as.factor(reg_no_protein_table$Type)

reg_no_protein_table$Type <-factor(reg_no_protein_table$Type,
                                   levels=c("", "tRNA_DOWN","rRNA_DOWN",
                                          "misc_RNA_DOWN","nc_RNA_DOWN", 
                                          "tRNA_UP","rRNA_UP","misc_RNA_UP",
                                          "nc_RNA_UP"))

}
write.table(reg_no_protein_table,file="Outputs_def/plot_tables/reg_no_protein_table.txt")

bar <- function(data,title){
  frame <- data.frame(data)
  sample <- factor(frame$Effects)
  n_gene <- frame$n_genes
  # Sample_Effect <- reorder(sample,abs(n_gene))
  breaks_values <- pretty(n_gene) 
  mycolors <- c("white","#FFD5C0","#ffbd9a","#FFA373","#FF7127","#72679B","#8076a5","#50486D","#2E2A3F")
  
  bar_pl <- ggplot(frame,aes(x=sample, y=n_gene,fill = Type))+
    geom_bar(stat = "identity",width = 0.8)+coord_flip()+
    scale_y_continuous(breaks = breaks_values,labels = abs(breaks_values))+
    scale_x_discrete(expand = c(0,2.5))+
    theme_classic()+ ggtitle(title)+geom_hline(yintercept=0)+theme(aspect.ratio = .9)+
    labs(x="Samples", y = "Number of Genes")+theme(plot.title = element_text(hjust = 0.5))+
    scale_fill_manual(name="Genes Type",values = mycolors)+
    theme(legend.key = element_rect(colour = 'white'))
    
    return(bar_pl)
}

bar_plot <- bar(reg_no_protein_table,"Regulated Genes \n Non Protein Coding Genes")

pdf(file = "Outputs_def/Figures/Bar_plot_no_protein_genes.pdf",width=8,height=5)
print(bar_plot)
dev.off()
















