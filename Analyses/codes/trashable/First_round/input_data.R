#############################
### 0. Load dependencies. ###
#############################

{ library(ashr)
  library(dplyr)
  library(readxl)
  library(writexl)
  library(readxl) 
  library(xlsx)
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
  path = "th_001_th_size_0_xlsx/" # Declare the file directory
  
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
  { 
    IDEE_all <- read.xlsx(paste0(path,"/",file_1,"/",all),1)
    IDEE_all=IDEE_all[order(IDEE_all$RV),]
    IDEE_down <- read.xlsx(paste0(path,"/",file_1,"/",down),1)
    IDEE_down=IDEE_down[order(IDEE_down$RV),]
    IDEE_up <- read.xlsx(paste0(path,"/",file_1,"/",up),1)
    IDEE_up=IDEE_up[order(IDEE_up$RV),]
    
    IDEE_reg <- rbind(IDEE_up,IDEE_down)
    IDEE_reg$Regulated <-c(rep(value_up,nrow(IDEE_up)),rep(value_dow,nrow(IDEE_down)))
    safe=IDEE_reg
    
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
    SFeDI_all=SFeDI_all[order(SFeDI_all$RV),]
    
    SFeDI_down <- read.xlsx(paste0(path,"/",file_5,"/",down),1)
    SFeDI_down=SFeDI_down[order(SFeDI_down$RV),]
    
    SFeDI_up <- read.xlsx(paste0(path,"/",file_5,"/",up),1)
    SFeDI_up=SFeDI_up[order(SFeDI_up$RV),]
    
    SFeDI_reg <- rbind(SFeDI_up,SFeDI_down)
    SFeDI_reg$Regulated <-c(rep(value_up,nrow(SFeDI_up)),rep(value_dow,nrow(SFeDI_down)))
  }
  ### Read the data.
  path1 <- "processed_all_JAC/metadata.txt"
  
  metadata = read.table(paste0(path1,"/metadata.txt"))
  feature_data = read.table(paste0(path1,"/feature_data.txt"))
  reads = read.table(paste0(path1,"/reads.txt"))
  fpkms = read.table(paste0(path1,"/fpkm.txt"))
  rownames(feature_data)=feature_data$Gene_ID
  
  sel_cul <- c("C5","C6")
  metadata=metadata[which(metadata$Culture %in% sel_cul),]
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

##########################Coding protein table##################################
{
unique(feature_data$Type)
prot_cod <- feature_data[which(feature_data$Type == "protein_coding"),]
print(nrow(prot_cod))
non_prot_cod <- feature_data[which(feature_data$Type != "protein_coding"),]
print(nrow(non_prot_cod))
prot_cod <- prot_cod[,c("Gene_ID","Type")]
non_prot_cod <- non_prot_cod[,c("Gene_ID","Type")]

rRNA <- feature_data[which(feature_data$Type == "rRNA"),]
print(nrow(rRNA))
tRNA <- feature_data[which(feature_data$Type == "tRNA"),]
print(nrow(tRNA))
ncRNA <- feature_data[which(feature_data$Type == "ncRNA"),]
print(nrow(ncRNA))
misc_RNA <- feature_data[which(feature_data$Type == "misc_RNA"|
                                 feature_data$Type == "."),]
print(nrow(misc_RNA))

IDEE_protein<-IDEE_reg[which(IDEE_reg$RV %in% row.names(prot_cod)),]
IDEE_protein$Type <-"protein_coding" 
IDEE_no_protein<-IDEE_reg[which(IDEE_reg$RV %in% row.names(non_prot_cod)),]

IDES_protein<-IDES_reg[which(IDES_reg$RV %in% prot_cod$Gene_ID),]
IDES_protein$Type <-"protein_coding" 
IDES_no_protein<-IDES_reg[which(IDES_reg$RV %in% row.names(non_prot_cod)),]
IDES_no_protein$Type <-"non_protein_coding" 
IDES_no_protein[which(IDES_no_protein$RV %in%rRNA$Gene_ID),]$Type <- "rRNA"
IDES_no_protein[which(IDES_no_protein$RV %in% tRNA$Gene_ID),]$Type <- "tRNA"
IDES_no_protein[which(IDES_no_protein$RV %in% ncRNA$Gene_ID),]$Type <- "nc_RNA"
IDES_no_protein[which(IDES_no_protein$RV %in% misc_RNA$Gene_ID),]$Type <- "misc_RNA"

SEWFe_protein<-SEWFe_reg[which(SEWFe_reg$RV %in% prot_cod$Gene_ID),]
SEWFe_protein$Type <-"protein_coding" 
SEWFe_no_protein<-SEWFe_reg[which(SEWFe_reg$RV %in% row.names(non_prot_cod)),]
SEWFe_no_protein$Type <-"non_protein_coding" 
SEWFe_no_protein[which(SEWFe_no_protein$RV %in%rRNA$Gene_ID),]$Type <- "rRNA"
SEWFe_no_protein[which(SEWFe_no_protein$RV %in% tRNA$Gene_ID),]$Type <- "tRNA"
SEWFe_no_protein[which(SEWFe_no_protein$RV %in% ncRNA$Gene_ID),]$Type <- "nc_RNA"
SEWFe_no_protein[which(SEWFe_no_protein$RV %in% misc_RNA$Gene_ID),]$Type <- "misc_RNA"

SEWOFe_protein<-SEWOFe_reg[which(SEWOFe_reg$RV %in% prot_cod$Gene_ID),]
SEWOFe_protein$Type <-"protein_coding" 
SEWOFe_no_protein<-SEWOFe_reg[which(SEWOFe_reg$RV %in% row.names(non_prot_cod)),]
SEWOFe_no_protein$Type <-"non_protein_coding" 
SEWOFe_no_protein[which(SEWOFe_no_protein$RV %in%rRNA$Gene_ID),]$Type <- "rRNA"
SEWOFe_no_protein[which(SEWOFe_no_protein$RV %in% tRNA$Gene_ID),]$Type <- "tRNA"
SEWOFe_no_protein[which(SEWOFe_no_protein$RV %in% ncRNA$Gene_ID),]$Type <- "nc_RNA"
SEWOFe_no_protein[which(SEWOFe_no_protein$RV %in% misc_RNA$Gene_ID),]$Type <- "misc_RNA"

SFeDI_protein<-SFeDI_reg[which(SFeDI_reg$RV %in% prot_cod$Gene_ID),]
SFeDI_protein$Type <-"protein_coding" 
SFeDI_no_protein<-SFeDI_reg[which(SFeDI_reg$RV %in% row.names(non_prot_cod)),]
SFeDI_no_protein$Type <-"non_protein_coding"
SFeDI_no_protein[which(SFeDI_no_protein$RV %in%rRNA$Gene_ID),]$Type <- "rRNA"
SFeDI_no_protein[which(SFeDI_no_protein$RV %in% tRNA$Gene_ID),]$Type <- "tRNA"
SFeDI_no_protein[which(SFeDI_no_protein$RV %in% ncRNA$Gene_ID),]$Type <- "nc_RNA"
SFeDI_no_protein[which(SFeDI_no_protein$RV %in% misc_RNA$Gene_ID),]$Type <- "misc_RNA"
}









