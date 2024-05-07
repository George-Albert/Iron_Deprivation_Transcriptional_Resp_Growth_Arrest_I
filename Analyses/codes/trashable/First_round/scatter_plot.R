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

############################
### 2. Read data ###########
############################
{
  ## filename directory.
  path = "th_0.01_th_size_0/"
  
  file_1 = "Iron_deprivation_effect_at_Exp"
  file_2 = "Iron_deprivation_effect_at_Stat"
  file_3 = "Stat_effect_with_Fe"
  file_4 = "Stat_effect_without_Fe"
  file_5 = "Stat_Fedeprivation_interaction"
  
  all = "all.xlsx"
  down = "downreg.xlsx"
  up = "upreg.xlsx"
  value_up="UP"
  value_dow="DOWN"
  ### Iron Deprivation effects at EXP
  
  IDEE_all <- read.xlsx(paste0(path,"/",file_1,"/",all),1)
  IDEE_all=IDEE_all[order(IDEE_all$RV),]
  IDEE_down <- read.xlsx(paste0(path,"/",file_1,"/",down),1)
  IDEE_down=IDEE_down[order(IDEE_down$RV),]
  IDEE_up <- read.xlsx(paste0(path,"/",file_1,"/",up),1)
  IDEE_up=IDEE_up[order(IDEE_up$RV),]
  
  IDEE_reg <- rbind(IDEE_up,IDEE_down)
  IDEE_reg$Regulated <-c(rep(value_up,nrow(IDEE_up)),rep(value_dow,nrow(IDEE_down)))
  
  ### Iron Deprivation effects at STAT
  
  IDES_all <- read.xlsx(paste0(path,"/",file_2,"/",all),1)
  IDES_all=IDES_all[order(IDES_all$RV),]
  IDES_down <- read.xlsx(paste0(path,"/",file_2,"/",down),1)
  IDES_down=IDES_down[order(IDES_down$RV),]
  IDES_up <- read.xlsx(paste0(path,"/",file_2,"/",up),1)
  IDES_up=IDES_up[order(IDES_up$RV),]
  
  IDES_reg <- rbind(IDES_up,IDES_down)
  IDES_reg$Regulated <-c(rep(value_up,nrow(IDES_up)),rep(value_dow,nrow(IDES_down)))
  ### Stat_effect_with_Fe
  
  SEWFe_all <- read.xlsx(paste0(path,"/",file_3,"/",all),1)
  SEWFe_all=SEWFe_all[order(SEWFe_all$RV),]
  SEWFe_down <- read.xlsx(paste0(path,"/",file_3,"/",down),1)
  SEWFe_down=SEWFe_down[order(SEWFe_down$RV),]
  SEWFe_up <- read.xlsx(paste0(path,"/",file_3,"/",up),1)
  SEWFe_up=SEWFe_up[order(SEWFe_up$RV),]
  
  SEWFe_reg <- rbind(SEWFe_up,SEWFe_down)
  SEWFe_reg$Regulated <-c(rep(value_up,nrow(SEWFe_up)),rep(value_dow,nrow(SEWFe_down)))
  ### Stat_effect_without_Fe
  
  SEWOFe_all <- read.xlsx(paste0(path,"/",file_4,"/",all),1)
  SEWOFe_all=SEWOFe_all[order(SEWOFe_all$RV),]
  SEWOFe_down <- read.xlsx(paste0(path,"/",file_4,"/",down),1)
  SEWOFe_down=SEWOFe_down[order(SEWOFe_down$RV),]
  SEWOFe_up <- read.xlsx(paste0(path,"/",file_4,"/",up),1)
  SEWOFe_up=SEWOFe_up[order(SEWOFe_up$RV),]
  
  SEWOFe_reg <- rbind(SEWOFe_up,SEWOFe_down)
  SEWOFe_reg$Regulated <-c(rep(value_up,nrow(SEWOFe_up)),rep(value_dow,nrow(SEWOFe_down)))
  
  ### Stat_Fedeprivation_interaction
  
  SFeDI_all <- read.xlsx(paste0(path,"/",file_5,"/",all),1)
  SFeDI_down <- read.xlsx(paste0(path,"/",file_5,"/",down),1)
  SFeDI_up <- read.xlsx(paste0(path,"/",file_5,"/",up),1)
  
  SFeDI_reg <- rbind(SFeDI_up,SFeDI_down)
  SFeDI_reg$Regulated <-c(rep(value_up,nrow(SFeDI_up)),rep(value_dow,nrow(SFeDI_down)))
  ### Read the data.
  
  metadata = read.table("processed_all_JAC/metadata.txt")
  feature_data = read.table("processed_all_JAC/feature_data.txt")
  reads = read.table("processed_all_JAC/reads.txt")
  fpkms = read.table("processed_all_JAC/fpkm.txt")
  rownames(feature_data)=feature_data$Gene_ID
  
  metadata=metadata[which(metadata$Culture %in% c("C5","C6")),]
  metadata$short_setup=paste0(metadata$Culture,"_Fe_",metadata$Iron)
  metadata <- metadata[order(metadata$short_setup),]
  reads <- reads[,rownames(metadata)]
  fpkms <- fpkms[,rownames(metadata)]
  
  ## Check coherence:
  length(which(colnames(reads)!=rownames(metadata)))
  
  # Select geneset:
  filtrado=filter_OK(reads,fpkms,feature_data,threshold=2,metadata=metadata,column="short_setup")
  
  reads=filtrado[[1]]
  fpkms=filtrado[[2]]
  feature_data=filtrado[[3]]
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


#1. Genes upreg a transitar a stat de manera independiente del Fe:

stat_effect$Condition <- "Not Sig."
safe=stat_effect
up_ind_Fe <- subset(stat_effect,LogFC1>0 & FDR1<0.01 & LogFC2>0 & FDR2<0.01 &
                      FDR_int>0.2)

stat_effect[which(stat_effect$RV1 %in% up_ind_Fe$RV1),]$Condition <- "up_ind_Fe"

#2. Genes DOWNREG a transitar a stat de manera independiente del Fe:

down_ind_Fe <- subset(stat_effect,LogFC1<0 & FDR1<0.01 & LogFC2<0 & FDR2<0.01 &
                        FDR_int >0.2)

stat_effect[which(stat_effect$RV1 %in% down_ind_Fe$RV1),]$Condition <- "down_ind_Fe"




scatter_pl <- function(data,x,y,up_g,down_g,title){
  mycolors <- c("#FFA373","grey","#50486D")
  
  ggplot(data, aes(x=x, y=y,group=Condition,col=Condition)) +geom_point(alpha=0.4)+
    geom_point(alpha=0.4)+
    scale_color_manual(values=mycolors)+
    ggtitle(title)+theme_minimal_grid()+geom_abline(col="red")+
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

pdf(file = "Outputs_def/Figures/stat_eff_plot.pdf",width=6,height=5)

scatter_pl(data=data,x=x,y=y,up_ind_Fe, down_ind_Fe,"STAT Effects with and w/o Fe")

dev.off()