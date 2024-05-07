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
#                         ===  Scatter Plot ===
#                         =====================

# =================Create scatter table=======================================
SEWFe_path = "Outputs_def/plot_tables/volcano_table/vol_SEWFe.txt"
SEWOFe_path = "Outputs_def/plot_tables/volcano_table/vol_SEWOFe.txt"
SFeDI_path <- "Outputs_def/plot_tables/volcano_table/vol_SFeDI.txt"

SEWFe_all <- read.table(SEWFe_path)
SEWOFe_all<- read.table(SEWOFe_path)
SFeDI_all <- read.table(SFeDI_path)
stat_effect <- data.frame(RV1=SEWFe_all$RV,LogFC1=SEWFe_all$log2FoldChange,
                          FDR1=SEWFe_all$BH,Association1=SEWFe_all$Association,
                          LogFC2=SEWOFe_all$log2FoldChange,
                          FDR2=SEWOFe_all$BH,Association2=SEWOFe_all$Association,
                          LogFC_int=SFeDI_all$log2FoldChange,FDR_int=SFeDI_all$BH,
                          Association=SFeDI_all$Association)
#SEWFe_all[which(stat_effect$RV1 != stat_effect$RV2 ),]

up_g <- subset(stat_effect,Association1=="UP" & Association2=="UP")
down_g <- stat_effect[which(stat_effect$Association1=="DOWN" & stat_effect$Association2=="DOWN"),]


#1. Genes upreg a transitar a stat de manera independiente del Fe

stat_effect$Condition <- "Not Sig."
safe=stat_effect
up_ind_Fe <- subset(stat_effect,LogFC1>0 & FDR1<0.01 & LogFC2>0 & FDR2<0.01 &
                      FDR_int>0.2)

stat_effect[which(stat_effect$RV1 %in% up_ind_Fe$RV1),]$Condition <- "up_ind_Fe"

#2. Genes DOWNREG a transitar a stat de manera independiente del Fe

down_ind_Fe <- subset(stat_effect,LogFC1<0 & FDR1<0.01 & LogFC2<0 & FDR2<0.01 &
                        FDR_int >0.2)

stat_effect[which(stat_effect$RV1 %in% down_ind_Fe$RV1),]$Condition <- "down_ind_Fe"

#3. Genes Upreg. al transitar a stat de manera dependiente del hierro
#más upregulation without Fe

up_dep_WOFe <- subset(stat_effect,LogFC1>0 & LogFC2>0 & FDR2<0.01 & LogFC_int>0 &
                        FDR_int<0.01)

stat_effect[which(stat_effect$RV1 %in% up_dep_WOFe$RV1),]$Condition <- "up_dep_without_Fe"

#4. Genes Upreg al transitar a stat de manera dependiente del hierro
###más upregulation with Fe
up_dep_Fe <- subset(stat_effect,LogFC1>0 & LogFC2>0 & FDR1<0.01 & LogFC_int<0 &
                      FDR_int<0.01)

stat_effect[which(stat_effect$RV1 %in% up_dep_Fe$RV1),]$Condition <- "up_dep_with_Fe"

#5. Genes Downreg. al transitar a stat de manera dependiente del hierro
##más downregulation without Fe
##
down_dep_WOFe <- subset(stat_effect,LogFC1<0 & LogFC2<0 & FDR2<0.01 & LogFC_int<0 &
                          FDR_int<0.01)

stat_effect[which(stat_effect$RV1 %in% down_dep_WOFe$RV1),]$Condition <- "down_dep_without_Fe"

#6. Genes Downreg. al transitar a stat de manera dependiente del hierro
##más upregulation with Fe

down_dep_up_WFe <- subset(stat_effect,LogFC1<0 & FDR1<0.01 & LogFC2<0 & LogFC_int>0 &
                            FDR_int<0.01)

stat_effect[which(stat_effect$RV1 %in% down_dep_up_WFe$RV1),]$Condition <- "down_dep_up_without_Fe"

#7. Genes that are regulated with different signs in the transition to stat,
# depending on iron availability Up without iron, down with iron

down_up_WFe <- subset(stat_effect,LogFC1<0 & FDR1<0.01 & LogFC2>0 & FDR2<0.01 & 
                        LogFC_int>0 & FDR_int<0.01)

stat_effect[which(stat_effect$RV1 %in% down_up_WFe$RV1),]$Condition <- "down_up_reg_with_Fe"

#8. Genes that are regulated with different signs in the transition to stat, depending on iron availability:
#Up with iron, down without iron

up_WFe_down <- subset(stat_effect,LogFC1>0 & FDR1<0.01 & LogFC2<0 & FDR2<0.01 & 
                        LogFC_int<0 & FDR_int<0.01)

stat_effect[which(stat_effect$RV1 %in% up_WFe_down$RV1),]$Condition <- "up_WFe_down_without_Fe"

#=========================PLOT=================================================
scatter_pl <- function(data,x,y,up_g,down_g,title){
  mycolors <- c("#FFA373","magenta","blue","darkgrey","red","green","yellow")
  
  ggplot(data, aes(x=x, y=y,group=Condition,col=Condition))+
    geom_point(alpha=.5)+
    scale_color_manual(values=mycolors)+
    ggtitle(title)+theme_minimal_grid()+geom_abline(col="black")+
    geom_hline(yintercept=0, color="black", size=1) + 
    geom_vline(xintercept=0, color="black", size=1)+
    xlab("LogFC1: STAT Effects with Fe")+xlim(-15,15)+
    ylab("LogFC2: STAT Effects without Fe")+ylim(-12,12)+
    theme(plot.title = element_text(hjust = 0.5))
}



data=stat_effect
x = stat_effect$LogFC1
y=stat_effect$LogFC2

# require(stats)
# fit_lm<-lm(y ~ x, data = data)
# fit_lm

pdf(file = "Outputs_def/Figures/stat_eff_plot_1.pdf",width=8,height=5)

scatter_pl(data=data,x=x,y=y,up_ind_Fe, down_ind_Fe,"STAT Effects with and w/o Fe")

dev.off()
