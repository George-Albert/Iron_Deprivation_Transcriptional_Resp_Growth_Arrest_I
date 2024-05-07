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
filter_OK=function(reads,fpkms,features,threshold=2,metadata,column){
  
  number_conditions=length(unique(metadata[,column]))
  medians=fpkms[,1:number_conditions]
  colnames(medians)=unique(metadata[,column])
  ## This is unnecessary, nbut just for clarifying:
  for(i in 1:number_conditions)
    medians[,i]=0
  for(i in 1:number_conditions)
  {
    set=which(metadata[,column]==colnames(medians)[i])
    chunk=fpkms[,set]
    medians[,i]=apply(chunk,1,median)
  }
  medians$max_median=apply(medians,1,max)
  genes_to_keep=rownames(medians)[which(medians$max_median>threshold)]
  length(genes_to_keep)
  
  filtered_reads=reads[genes_to_keep,]
  filtered_features=features[genes_to_keep,]
  filtered_fpkms=fpkms[genes_to_keep,]
  
  output=list(filtered_reads=filtered_reads,filtered_fpkms=filtered_fpkms,filtered_features=filtered_features)
  
  return(output)
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

#                         ====================================
#                         ===  Expression level: all genes ===
#                         ====================================
 
#====================Expression levels=========================================

y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
v=voom(y,design,plot=TRUE)
exp=v$E

dcols(exp)

exp1 <- data.frame(exp[which(rownames(exp) %in% IDEE_up$RV),])
exp1 <- exp1 - rowMeans(exp1)
exp1$RV <- rownames(exp1)

#resort columns
exp1 <- exp1 %>%
  select(RV, everything())

Fe_NO_exp <- data.frame(Exp_level=cbind(c(exp1[,2],exp1[,3])))
Fe_NO_exp$Condition <- "Fe_NO_EXP"
Fe_NO_exp$RV <- exp1[,1]

Fe_YES_exp <- data.frame(Exp_level=cbind(c(exp1[,4],exp1[,5])))
Fe_YES_exp$Condition <- "Fe_YES_EXP"
Fe_YES_exp$RV <- exp1[,1]

Fe_NO_stat <- data.frame(Exp_level=cbind(c(exp1[,6],exp1[,7])))
Fe_NO_stat$Condition <- "Fe_NO_STAT"
Fe_NO_stat$RV <- exp1[,1]

Fe_YES_stat <- data.frame(Exp_level=cbind(c(exp1[,8],exp1[,9])))
Fe_YES_stat$Condition <- "Fe_YES_STAT"
Fe_YES_stat$RV <- exp1[,1]

#Create table Espression, Condition and RV

DE_table <- rbind(Fe_NO_exp,Fe_YES_exp,Fe_NO_stat,Fe_YES_stat)

write.table(DE_table,file="Outputs_def/plot_tables/DE_table_all.txt")

exp_plt <- function(data,title){
  datos <- data.frame(data)
  y <- datos$Exp_level
  x <- datos$Condition
  group <- datos$RV
  mycolors <- c("brown","blue","green","darkorchid3","grey",
                "cyan1","red","yellow","orange")
  fill <- c("brown","blue","green","darkorchid3","grey",
            "cyan1","red","yellow","orange")
  # mycolors <- c("red","blue","red","blue")
  # fill <- c("red","blue",NA,NA)
  shape <- c(21,21,1,1)
  
  ggplot(datos,aes(x,y,
                   group=group,
                   col=group,
                   shape=Condition,
                   fill=group))+
    geom_jitter()+
    scale_color_manual(values=mycolors)+
    scale_shape_manual(values=shape,guide="none")+
    scale_fill_manual(values=fill,guide="none")+
    labs(color="Genes:")+
    xlim("Fe_NO_EXP","Fe_YES_EXP","Fe_NO_STAT","Fe_YES_STAT")+
    xlab("Conditions")+scale_y_continuous(breaks=seq(-4,10,2))+
    ylab("log2-counts per million (logCPM)")+theme_minimal()+
    ggtitle(title)+theme(plot.title = element_text(hjust = 0.5))
  
}

pdf(file = "Outputs_def/Figures/level_exp_rem_mean.pdf",width=8,height=5)
exp_plt(DE_table,"Differential expression levels")
dev.off()


#                         =====================================
#                         ===  Expression level: each genes ===
#                         =====================================
{
  Rv2377c <- data.frame(Rv2377c=exp[which(rownames(exp) %in% IDEE_up$RV[1]),])
  Rv2377c$Condition <- rownames(Rv2377c)
  safe=Rv2377c
  Rv2377c[Rv2377c =="C5_Fe_NO_1" | Rv2377c =="C5_Fe_NO_2"] <- "Fe_NO_EXP"
  Rv2377c[Rv2377c =="C5_Fe_YES_1" | Rv2377c =="C5_Fe_YES_2"] <- "Fe_YES_EXP"
  Rv2377c[Rv2377c =="C6_Fe_NO_1" | Rv2377c =="C6_Fe_NO_2"] <- "Fe_NO_STAT"
  Rv2377c[Rv2377c =="C6_Fe_YES_1" | Rv2377c =="C6_Fe_YES_2"] <- "Fe_YES_STAT"
  
  Rv2378c <- data.frame(Rv2378c=exp[which(rownames(exp) %in% IDEE_up$RV[2]),])
  Rv2378c$Condition <- rownames(Rv2378c)
  Rv2378c[Rv2378c =="C5_Fe_NO_1" | Rv2378c =="C5_Fe_NO_2"] <- "Fe_NO_EXP"
  Rv2378c[Rv2378c =="C5_Fe_YES_1" | Rv2378c =="C5_Fe_YES_2"] <- "Fe_YES_EXP"
  Rv2378c[Rv2378c =="C6_Fe_NO_1" | Rv2378c =="C6_Fe_NO_2"] <- "Fe_NO_STAT"
  Rv2378c[Rv2378c =="C6_Fe_YES_1" | Rv2378c =="C6_Fe_YES_2"] <- "Fe_YES_STAT"
  
  Rv2379c <- data.frame(Rv2379c=exp[which(rownames(exp) %in% IDEE_up$RV[3]),])
  Rv2379c$Condition <- rownames(Rv2379c)
  Rv2379c[Rv2379c =="C5_Fe_NO_1" | Rv2379c =="C5_Fe_NO_2"] <- "Fe_NO_EXP"
  Rv2379c[Rv2379c =="C5_Fe_YES_1" | Rv2379c =="C5_Fe_YES_2"] <- "Fe_YES_EXP"
  Rv2379c[Rv2379c =="C6_Fe_NO_1" | Rv2379c =="C6_Fe_NO_2"] <- "Fe_NO_STAT"
  Rv2379c[Rv2379c =="C6_Fe_YES_1" | Rv2379c =="C6_Fe_YES_2"] <- "Fe_YES_STAT"
  
  Rv2380c <- data.frame(Rv2380c=exp[which(rownames(exp) %in% IDEE_up$RV[4]),])
  Rv2380c$Condition <- rownames(Rv2380c)
  Rv2380c[Rv2380c =="C5_Fe_NO_1" | Rv2380c =="C5_Fe_NO_2"] <- "Fe_NO_EXP"
  Rv2380c[Rv2380c =="C5_Fe_YES_1" | Rv2380c =="C5_Fe_YES_2"] <- "Fe_YES_EXP"
  Rv2380c[Rv2380c =="C6_Fe_NO_1" | Rv2380c =="C6_Fe_NO_2"] <- "Fe_NO_STAT"
  Rv2380c[Rv2380c =="C6_Fe_YES_1" | Rv2380c =="C6_Fe_YES_2"] <- "Fe_YES_STAT"
  
  Rv2381c <- data.frame(Rv2381c=exp[which(rownames(exp) %in% IDEE_up$RV[5]),])
  Rv2381c$Condition <- rownames(Rv2381c)
  Rv2381c[Rv2381c =="C5_Fe_NO_1" | Rv2381c =="C5_Fe_NO_2"] <- "Fe_NO_EXP"
  Rv2381c[Rv2381c =="C5_Fe_YES_1" | Rv2381c =="C5_Fe_YES_2"] <- "Fe_YES_EXP"
  Rv2381c[Rv2381c =="C6_Fe_NO_1" | Rv2381c =="C6_Fe_NO_2"] <- "Fe_NO_STAT"
  Rv2381c[Rv2381c =="C6_Fe_YES_1" | Rv2381c =="C6_Fe_YES_2"] <- "Fe_YES_STAT"
  
  Rv2382c <- data.frame(Rv2382c=exp[which(rownames(exp) %in% IDEE_up$RV[6]),])
  Rv2382c$Condition <- rownames(Rv2382c)
  Rv2382c[Rv2382c =="C5_Fe_NO_1" | Rv2382c =="C5_Fe_NO_2"] <- "Fe_NO_EXP"
  Rv2382c[Rv2382c =="C5_Fe_YES_1" | Rv2382c =="C5_Fe_YES_2"] <- "Fe_YES_EXP"
  Rv2382c[Rv2382c =="C6_Fe_NO_1" | Rv2382c =="C6_Fe_NO_2"] <- "Fe_NO_STAT"
  Rv2382c[Rv2382c =="C6_Fe_YES_1" | Rv2382c =="C6_Fe_YES_2"] <- "Fe_YES_STAT"
  
  Rv2383c <- data.frame(Rv2383c=exp[which(rownames(exp) %in% IDEE_up$RV[7]),])
  Rv2383c$Condition <- rownames(Rv2383c)
  Rv2383c[Rv2383c =="C5_Fe_NO_1" | Rv2383c =="C5_Fe_NO_2"] <- "Fe_NO_EXP"
  Rv2383c[Rv2383c =="C5_Fe_YES_1" | Rv2383c =="C5_Fe_YES_2"] <- "Fe_YES_EXP"
  Rv2383c[Rv2383c =="C6_Fe_NO_1" | Rv2383c =="C6_Fe_NO_2"] <- "Fe_NO_STAT"
  Rv2383c[Rv2383c =="C6_Fe_YES_1" | Rv2383c =="C6_Fe_YES_2"] <- "Fe_YES_STAT"
  
  Rv2386c <- data.frame(Rv2386c=exp[which(rownames(exp) %in% IDEE_up$RV[8]),])
  Rv2386c$Condition <- rownames(Rv2386c)
  Rv2386c[Rv2386c =="C5_Fe_NO_1" | Rv2386c =="C5_Fe_NO_2"] <- "Fe_NO_EXP"
  Rv2386c[Rv2386c =="C5_Fe_YES_1" | Rv2386c =="C5_Fe_YES_2"] <- "Fe_YES_EXP"
  Rv2386c[Rv2386c =="C6_Fe_NO_1" | Rv2386c =="C6_Fe_NO_2"] <- "Fe_NO_STAT"
  Rv2386c[Rv2386c =="C6_Fe_YES_1" | Rv2386c =="C6_Fe_YES_2"] <- "Fe_YES_STAT"
  
  Rv3402c <- data.frame(Rv3402c=exp[which(rownames(exp) %in% IDEE_up$RV[9]),])
  Rv3402c$Condition <- rownames(Rv3402c)
  Rv3402c[Rv3402c =="C5_Fe_NO_1" | Rv3402c =="C5_Fe_NO_2"] <- "Fe_NO_EXP"
  Rv3402c[Rv3402c =="C5_Fe_YES_1" | Rv3402c =="C5_Fe_YES_2"] <- "Fe_YES_EXP"
  Rv3402c[Rv3402c =="C6_Fe_NO_1" | Rv3402c =="C6_Fe_NO_2"] <- "Fe_NO_STAT"
  Rv3402c[Rv3402c =="C6_Fe_YES_1" | Rv3402c =="C6_Fe_YES_2"] <- "Fe_YES_STAT"
}
list_genes <- list(Rv2377c,Rv2378c,Rv2379c,Rv2380c,Rv2381c,Rv2382c,Rv2383c,
                Rv2386c,Rv3402c)
#write.table((Rv2377c,),file="Outputs_def/plot_tables/DE_table_all.txt")

mapply(write.table, x = list_genes, file = c("DE_table_1.txt", "DE_table_2.txt",
                                    "DE_table_3.txt","DE_table_4.txt",
                                    "DE_table_5.txt","DE_table_6.txt",
                                    "DE_table_7.txt","DE_table_8.txt",
                                    "DE_table_9.txt"))
exp_plt_1 <- function(data,y,title){
  datos <- data.frame(data)
  
  x <- datos$Condition
  group <- datos$Condition
  mycolors <- c("red","blue","red","blue")
  fill <- c("red","blue",NA,NA)
  shape <- c(21,21,1,1)
  
  ggplot(datos,aes(x,y,
                   shape=Condition,
                   fill=Condition,
                   col=Condition))+
    geom_jitter()+
    scale_color_manual(values=mycolors)+
    scale_shape_manual(values=shape)+
    scale_fill_manual(values=fill)+
    xlim("Fe_NO_EXP","Fe_YES_EXP","Fe_NO_STAT","Fe_YES_STAT")+
    xlab("Conditions")+scale_y_continuous(breaks=seq(-4,10,2))+
    ylab("log2-counts per million (logCPM)")+theme_minimal()+
    guides(
      colour = guide_legend("Conditions"),
      fill = guide_legend("Conditions"),
      shape = guide_legend("Conditions"),
      
    )+
    ggtitle(title)+theme(plot.title = element_text(hjust = 0.5))
  
}

pdf(file = "Outputs_def/Figures/level_exp_each_gen.pdf",width=8,height=5)

exp_plt_1(Rv2377c,Rv2377c$Rv2377c,"Rv2377c")
exp_plt_1(Rv2378c,Rv2378c$Rv2378c,"Rv2378c")
exp_plt_1(Rv2379c,Rv2379c$Rv2379c,"Rv2379c")
exp_plt_1(Rv2380c,Rv2380c$Rv2380c,"Rv2380c")
exp_plt_1(Rv2381c,Rv2381c$Rv2381c,"Rv2381c")
exp_plt_1(Rv2382c,Rv2382c$Rv2382c,"Rv2382c")
exp_plt_1(Rv2383c,Rv2383c$Rv2383c,"Rv2383c")
exp_plt_1(Rv2386c,Rv2386c$Rv2386c,"Rv2386c")
exp_plt_1(Rv3402c,Rv3402c$Rv3402c,"Rv3402c")

dev.off()