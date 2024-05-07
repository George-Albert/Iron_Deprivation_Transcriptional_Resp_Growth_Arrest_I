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



#                         ====================================
#                         === Density Plot: each CONDITION ===
#                         ====================================
{
  dens1 <- data.frame(LogFC=IDEE_reg$log2FoldChange,samples= rep(Samples[1],
                                                                 nrow(IDEE_reg)))
  dens2 <- data.frame(LogFC=IDES_reg$log2FoldChange,samples= rep(Samples[2],
                                                                 nrow(IDES_reg)))
  dens3 <- data.frame(LogFC=SEWFe_reg$log2FoldChange,samples= rep(Samples[3],
                                                                  nrow(SEWFe_reg)))
  dens4 <- data.frame(LogFC=SEWOFe_reg$log2FoldChange,samples= rep(Samples[4],
                                                                   nrow(SEWOFe_reg)))
  dens5 <- data.frame(LogFC=SFeDI_reg$log2FoldChange,samples= rep(Samples[5],
                                                                  nrow(SFeDI_reg)))
  dens <- rbind(dens1,dens2,dens3,dens4,dens5)
  
}

write.table(dens,file="processed_all_JAC/dens_table.txt")

density <- function(data,title){
  
  frame <- data.frame(data)
  lfc <- frame$LogFC
  y <- factor(frame$samples)
  sam <- reorder(y,lfc)
  color <- c("#00AFBB", "#E7B800","#50486D", "#FC4E07","#FFA373")
  
  ggplot(frame, aes(x=lfc,y=sam,fill=y)) +
    geom_density_ridges(scale = 4, alpha = 0.6) + xlim(range(-5,5))+
    scale_fill_manual(values = color)+geom_vline(xintercept=0,lty="dotted")+
    labs(x="LogFC", y = "Density",title=title)+
    theme_ridges(font_size = 8,center_axis_labels = TRUE,grid = FALSE)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
}

pdf(file = "Outputs_def/plot_table/dens_plot.pdf",width=9,height=5)

density(dens,"Density Plot")

dev.off()