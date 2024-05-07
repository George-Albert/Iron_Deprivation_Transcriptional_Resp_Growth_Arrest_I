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