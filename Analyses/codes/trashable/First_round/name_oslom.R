# Program to get the names and Rv of genes network runned by OSLOM

# *****************
# ***Dependecies***
# *****************

library(dplyr)
library(tidyverse)

# *****************
# ***functions*****
# *****************
name_gene <- function(path,tp_name,pr_data,Condition){
  # path: set the file path
  # tp_name: tp name of the hierarchical level, i.e (tp1, tp1, etc.)
  # pr_data: data of the protein interaction file
  # Condition: Condition we are working on, i.e (IDEE)
  # dir_0: directory of rv form tp file
  # dir_1: directory of names of genes from the tp file
  # dir_3:directory of rv and gene name from partition_level_0 file
  path <- path
  
  tp_path <- paste0(path,"/",tp_name)
  partition_path <- paste0(path,"/partitions_level_0")
  
  
  
  n_row <- length(count.fields(tp_path, sep = ""))
  n_col <- max(count.fields(tp_path, sep = ""))
  tp <- read.delim(tp_path,header=F, sep=" ",comment.char = "#",strip.white=T,
                   fill=T,quote = "",col.names=paste("column", 1:(n_col+1), sep="_"),
                   row.names=paste("module", (1:n_row)-1, sep="_"))
  tp <- tp[1:n_row,1:n_col]
  tp <- data.frame(t(tp))
  
  n_row <- length(count.fields(partition_path, sep = ""))
  n_col <- max(count.fields(partition_path, sep = ""))
  
  partition <- read.delim (partition_path,header=F, sep=" ",strip.white=T,
                           fill=T,quote = "",
                           col.names=paste("column", 1:(n_col+2),sep="_"))
  name_col <- paste(partition$column_1 [c(T,F)],partition$column_2[c(T,F)],sep = "_")
  
  partition <- read.delim (partition_path,header=F, sep=" ",strip.white=T,
                           comment.char = "#",fill=T,quote = "",
                           col.names=paste("column", 1:(n_col+2),sep="_"),
                           row.names=paste("module", (1:n_row)-1, sep="_"))
  partition <- partition[1:n_row,1:n_col]
  partition <- data.frame(t(partition))
  colnames(partition) <- name_col
  sel_col <- which(colnames(partition)=="#module_0")[10]
  partition <- partition[,sel_col:ncol(partition)]
  # partition <- partition[, !duplicated(colnames(partition),fromLast = T)]
  partition_1 <- data.frame(partition)
  
  tp_community_0 <- pr_data[which(pr_data$number_gene %in% tp$X10),]
  tp_community_1 <- pr_data[which(pr_data$number_gene %in% tp$X4),]
  rv_tp_community_0 <- tp_community_0$RV_gene
  name_tp_community_0 <- tp_community_0$Gene_name
  
  dir_0 <- paste0("rv_",Condition,"/")
  dir_1 <- paste0("name_",Condition,"/")
  dir_3 <- paste0("partition_",Condition,"/")
  
  for(i in 1:ncol(tp)) {       # for-loop over columns
    number <- which(pr_data$number_gene %in% tp[ ,i])
    # print(number)
    tp_community  <- pr_data[number,]
    rv_tp_community <- tp_community[,1]
    n_tp_community <- tp_community[,3]
    
    dir.create(file.path(dir_0,"rv"), recursive = TRUE,showWarnings = FALSE)
    dir.create(file.path(dir_1,"name"), recursive = TRUE,showWarnings = FALSE)
    
    file=paste0(dir_0,"rv/rv_",tp_name,"_community_",i-1,".txt")
    file1=paste0(dir_1,"name/name_",tp_name,"_community_",i-1,".txt")
    
    write.table(rv_tp_community,file=file,sep="\t",quote = FALSE,row.names = FALSE,
                col.names = FALSE)
    write.table(n_tp_community,file=file1,sep="\t",quote = FALSE,row.names = FALSE,
                col.names = FALSE)
  }
  
  for(i in 1:ncol(partition_1)) {       # for-loop over columns
    
    number_p <- which(pr_data$number_gene %in% partition_1[ ,i])
    # print(number)
    partition_0  <- pr_data[number_p,]
    rv_partition_community <- partition_0[,1]
    n_partition_community <- partition_0[,3]
    
    dir.create(file.path(dir_3,"rv"), recursive = TRUE,showWarnings = FALSE)
    dir.create(file.path(dir_3,"name"), recursive = TRUE,showWarnings = FALSE)
    
    dir.create(dir_3,showWarnings = FALSE)
    file=paste0(dir_3,"rv/rv_partition_community_",i-1,".txt")
    
    dir.create(dir_3,showWarnings = FALSE)
    file1=paste0(dir_3,"name/name_partition_community_",i-1,".txt")
    
    
    write.table(rv_partition_community,file=file,sep="\t",quote = FALSE,row.names = FALSE,
                col.names = FALSE)
    write.table(n_partition_community,file=file1,sep="\t",quote = FALSE,row.names = FALSE,
                col.names = FALSE)
  }
  mylist <- list(tp_community,rv_tp_community,n_tp_community)
  
  
  return(mylist)
}

tp_gene <- function(path,tp_name,pr_data,Condition){
  # path: set the file path
  # tp_name: tp name of the hierarchical level, i.e (tp1, tp1, etc.)
  # pr_data: data of the protein interaction file
  # Condition: Condition we are working on, i.e (IDEE)
  # dir_0: directory of rv form tp file
  # dir_1: directory of names of genes from the tp file
  # dir_3:directory of rv and gene name from partition_level_0 file
  path <- path
  
  tp_path <- paste0(path,"/",tp_name)
  partition_path <- paste0(path,"/partitions_level_0")
  
  
  
  n_row <- length(count.fields(tp_path, sep = ""))
  n_col <- max(count.fields(tp_path, sep = ""))
  tp <- read.delim(tp_path,header=F, sep=" ",comment.char = "#",strip.white=T,
                   fill=T,quote = "",col.names=paste("column", 1:(n_col+1), sep="_"),
                   row.names=paste("module", (1:n_row)-1, sep="_"))
  tp <- tp[1:n_row,1:n_col]
  tp <- data.frame(t(tp))
  
  tp_community_0 <- pr_data[which(pr_data$number_gene %in% tp$X10),]
  tp_community_1 <- pr_data[which(pr_data$number_gene %in% tp$X4),]
  rv_tp_community_0 <- tp_community_0$RV_gene
  name_tp_community_0 <- tp_community_0$Gene_name
  
  dir_0 <- paste0("rv_",Condition,"/")
  dir_1 <- paste0("name_",Condition,"/")
  
  for(i in 1:ncol(tp)) {       # for-loop over columns
    number <- which(pr_data$number_gene %in% tp[ ,i])
    # print(number)
    tp_community  <- pr_data[number,]
    rv_tp_community <- tp_community[,1]
    n_tp_community <- tp_community[,3]
    
    dir.create(file.path(dir_0,"rv"), recursive = TRUE,showWarnings = FALSE)
    dir.create(file.path(dir_1,"name"), recursive = TRUE,showWarnings = FALSE)
    
    file=paste0(dir_0,"rv/rv_",tp_name,"_community_",i-1,".txt")
    file1=paste0(dir_1,"name/name_",tp_name,"_community_",i-1,".txt")
    
    write.table(rv_tp_community,file=file,sep="\t",quote = FALSE,row.names = FALSE,
                col.names = FALSE)
    write.table(n_tp_community,file=file1,sep="\t",quote = FALSE,row.names = FALSE,
                col.names = FALSE)
  }
  mylist <- list(tp_community,rv_tp_community,n_tp_community)
  
  
  return(mylist)
}


path1 <- "Iron_deprivation_exp_oslom_runned/Iron_deprivation_effect_at_Exp_down.dat_oslo_files/"
path2 <- "Iron_deprivation_exp_oslom_runned/Iron_deprivation_effect_at_Exp_up.dat_oslo_files/"
path3 <- "Iron_deprivation_stat_oslom_runned/Iron_deprivation_effect_at_Stat_up.dat_oslo_files/"
path4 <- "Iron_deprivation_stat_oslom_runned/Iron_deprivation_effect_at_Stat_down.dat_oslo_files/"
path5 <- "Stat_effect_with_Fe_oslom_runned/Stat_effect_with_Fe_down.dat_oslo_files/"
path6 <- "Stat_effect_with_Fe_oslom_runned/Stat_effect_with_Fe_up.dat_oslo_files/"
path7 <- "Stat_effect_without_Fe_oslom_runned/Stat_effect_without_Fe_down.dat_oslo_files/"
path8 <- "Stat_effect_without_Fe_oslom_runned/Stat_effect_without_Fe_up.dat_oslo_files/"
path9 <- "Stat_Fedeprivation_interaction_runned/Stat_Fedeprivation_interaction_down.dat_oslo_files/"
path10 <- "Stat_Fedeprivation_interaction_runned/Stat_Fedeprivation_interaction_up.dat_oslo_files/"

mylist1 <- name_gene(path=path1,tp_name="tp",pr_data=IDEE_osl,Condition="IDEE_down")
mylist2 <- name_gene(path=path2,tp_name="tp",pr_data=IDEE_osl,Condition="IDEE_up")
mylist3 <- name_gene(path=path3,tp_name="tp",pr_data=IDEE_osl,Condition="IDES_up")
mylist4 <- name_gene(path=path4,tp_name="tp",pr_data=IDEE_osl,Condition="IDES_down")
mylist5 <- name_gene(path=path5,tp_name="tp",pr_data=IDEE_osl,Condition="SEWFe_down")
mylist6 <- name_gene(path=path6,tp_name="tp",pr_data=IDEE_osl,Condition="SEWFe_up")
mylist7 <- name_gene(path=path7,tp_name="tp",pr_data=IDEE_osl,Condition="SEWOFe_down")
mylist8 <- name_gene(path=path8,tp_name="tp",pr_data=IDEE_osl,Condition="SEWOFe_up")
mylist9 <- name_gene(path=path9,tp_name="tp",pr_data=IDEE_osl,Condition="SFeDI_down")
mylist10 <- name_gene(path=path10,tp_name="tp",pr_data=IDEE_osl,Condition="SFeDI_up")




mylist1_tp <- tp_gene(path=path1,tp_name="tp1",pr_data=IDEE_osl,Condition="IDEE_down")
mylist2_tp <- tp_gene(path=path2,tp_name="tp1",pr_data=IDEE_osl,Condition="IDEE_up")
mylist3_tp <- tp_gene(path=path3,tp_name="tp1",pr_data=IDEE_osl,Condition="IDES_up")
mylist4_tp <- tp_gene(path=path4,tp_name="tp1",pr_data=IDEE_osl,Condition="IDES_down")
mylist5_tp <- tp_gene(path=path5,tp_name="tp",pr_data=IDEE_osl,Condition="SEWFe_down")
mylist6_tp <- tp_gene(path=path6,tp_name="tp1",pr_data=IDEE_osl,Condition="SEWFe_up")
mylist7_tp <- tp_gene(path=path7,tp_name="tp1",pr_data=IDEE_osl,Condition="SEWOFe_down")
mylist8_tp <- tp_gene(path=path8,tp_name="tp1",pr_data=IDEE_osl,Condition="SEWOFe_up")
mylist9_tp <- tp_gene(path=path9,tp_name="tp1",pr_data=IDEE_osl,Condition="SFeDI_down")
mylist10_tp <- tp_gene(path=path10,tp_name="tp1",pr_data=IDEE_osl,Condition="SFeDI_up")
