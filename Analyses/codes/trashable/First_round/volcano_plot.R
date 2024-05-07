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


#                         =====================
#                         ===  Volcano Plot ===
#                         =====================

# =================Create volcano table=======================================
{
  IDEE_all$Association <- "Not Sig."
  IDEE_all[which(IDEE_all$RV %in% IDEE_up$RV),]$Association<- "UP"
  IDEE_all[which(IDEE_all$RV %in% IDEE_down$RV),]$Association<- "DOWN"
  
  IDES_all$Association <- "Not Sig."
  IDES_all[which(IDES_all$RV %in% IDES_up$RV),]$Association<- "UP"
  IDES_all[which(IDES_all$RV %in% IDES_down$RV),]$Association<- "DOWN"
  
  SEWFe_all$Association <- "Not Sig."
  SEWFe_all[which(SEWFe_all$RV %in% SEWFe_up$RV),]$Association<- "UP"
  SEWFe_all[which(SEWFe_all$RV %in% SEWFe_down$RV),]$Association<- "DOWN"
  
  SEWOFe_all$Association <- "Not Sig."
  SEWOFe_all[which(SEWOFe_all$RV %in% SEWOFe_up$RV),]$Association<- "UP"
  SEWOFe_all[which(SEWOFe_all$RV %in% SEWOFe_down$RV),]$Association<- "DOWN"
  
  SFeDI_all$Association <- "Not Sig."
  SFeDI_all[which(SFeDI_all$RV %in% SFeDI_up$RV),]$Association<- "UP"
  SFeDI_all[which(SFeDI_all$RV %in% SFeDI_down$RV),]$Association<- "DOWN"
}

list_vol <- list(IDEE_all,IDES_all,SEWFe_all,SEWOFe_all,SFeDI_all)
mapply(write.table, x = list_vol, file = c("vol_IDEE.txt", "vol_IDES..txt",
                                             "vol_SEWFe.txt","vol_SEWOFe.txt",
                                             "vol_SFeDI.txt"))

volcan_plot <- function(data,y,xintercept,th,ymax,title){
  
  datos <- data.frame(data)
  x <- datos$log2FoldChange
  y <- -(log10(y))
  y[is.infinite(y)] <- NA
  
  sizes <- c("UP" = 0.5, "DOWN" = 0.5, "Not Sig." = 0.5) 
  alphas <- c("UP" = 1, "DOWN" = 1, "Not Sig." = 0.5)
  #FFA373
  ##50486D
  mycolors <- c("UP" ="green","Not Sig." ="grey","DOWN" ="magenta")
  
  ggplot(datos,aes(x=x,y=y,
                   size=Association,
                   alpha=Association,
                   fill=Association))+
    geom_point(shape = 21,colour = "black")+
    ylab("-log10(FDR)")+xlab("Log2FC")+theme_minimal()+
    scale_fill_manual(values=mycolors)+
    scale_size_manual(values = sizes,guide="none")+
    scale_alpha_manual(values = alphas, guide="none")+
    geom_vline(xintercept=c(-xintercept, xintercept),col="black",linetype = "dashed")+
    geom_hline(yintercept=-log10(th), col="black",linetype = "dashed")+
    scale_x_continuous(limits = c(-(max(x)), max(x)))+
    scale_y_continuous(limits = c(-1, ymax),expand = expansion(0))+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))
}


pdf(file = "Outputs_def/Figures/volcano_plot.pdf",width=6,height=5)

volcan_plot(IDEE_all,y=IDEE_all$BH,xintercept = 0.0,th=0.01,17,title =Samples[1])

volcan_plot(IDES_all,y=IDES_all$BH,xintercept = 0.0,th=0.01,240,title =Samples[2])

volcan_plot(SEWFe_all,y=SEWFe_all$BH,xintercept = 0.0,th=0.01,300,title =Samples[3])

volcan_plot(SEWOFe_all,y=SEWOFe_all$BH,xintercept = 0.0,th=0.01,300,title =Samples[4])

volcan_plot(SFeDI_all,y=SFeDI_all$BH,xintercept = 0.0,th=0.01,105,title =Samples[5])

dev.off()
