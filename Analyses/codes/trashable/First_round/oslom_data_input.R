library(xlsx)
library(dplyr)

#==========Read Verbose===============
#==============================Iron_deprivation_effect_at_EXPt========================

path_IDEE <- "networks/verbose/Iron_deprivation_effect_at_Exp.txt"
IDEE <- read.table(path_IDEE)

IDEE_up <- IDEE[c(3,4,11)]
IDEE_up <- IDEE_up[which(IDEE_up$weight>0),]

IDEE_down <- IDEE[c(3,4,11)]
IDEE_down <- IDEE_down[which(IDEE_down$weight<0),]
IDEE_down$weight <- abs(IDEE_down$weight)
which(IDEE_down$weight==0)

dir.create("Iron_deprivation_exp_oslom")
write.table(IDEE_up,file="Iron_deprivation_exp_oslom/Iron_deprivation_effect_at_Exp_up.dat",row.names = FALSE,
            col.names = FALSE)

write.table(IDEE_down,file="Iron_deprivation_exp_oslom/Iron_deprivation_effect_at_Exp_down.dat",row.names = FALSE,
            col.names = FALSE)
#==============================Iron_deprivation_effect_at_Stat========================

path_IDES <- "networks/verbose/Iron_deprivation_effect_at_Stat.txt"
IDES <- read.table(path_IDES)

IDES_up <- IDES[c(3,4,11)]
IDES_up <- IDES_up[which(IDES_up$weight>0),]

IDES_down <- IDES[c(3,4,11)]
IDES_down <- IDES_down[which(IDES_down$weight<0),]
IDES_down$weight <- abs(IDES_down$weight)
which(IDES_up$weight==0)

dir.create("Iron_deprivation_stat_oslom")
write.table(IDES_up,file="Iron_deprivation_stat_oslom/Iron_deprivation_effect_at_Stat_up.dat",row.names = FALSE,
            col.names = FALSE)

write.table(IDES_down,file="Iron_deprivation_stat_oslom/Iron_deprivation_effect_at_Stat_down.dat",row.names = FALSE,
            col.names = FALSE)

oslom_table <- function(path,var,var_up,var_down,dir,file_up,file_down){

  path <- path
  var <- read.table(path)
  
  var_up <- var[c(3,4,11)]
  var_up <- var_up[which(var_up$weight>0),]
  
  var_down <- var[c(3,4,11)]
  var_down <- var_down[which(var_down$weight<0),]
  var_down$weight <- abs(var_down$weight)
  
  dir.create(dir,showWarnings = FALSE)
  
  file_up <- paste0(dir,file_up)
  write.table(var_up,file=file_up,row.names = FALSE,col.names = FALSE)
  
  file_down <- paste0(dir,file_down)
  write.table(var_down,file=file_down,row.names = FALSE,col.names = FALSE)
  
  mylist <- list(var,var_up,var_down)
  return(mylist)
}
#==============================Stat_effect_with_Fe========================

path <- "networks/verbose/Stat_effect_with_Fe.txt"
dir <-  "Stat_effect_with_Fe_oslom/"
file_up <-  "Stat_effect_with_Fe_up.dat"
file_down <-  "Stat_effect_with_Fe_down.dat"

mylist <- oslom_table(path=path,var=SEWFe,var_up=SEWFe_up,var_down=SEWFe_down,
            dir=dir,file_up=file_up,
            file_down = file_down)

SEWFe <- data.frame(mylist[1])
SEWFe_up <- data.frame(mylist[2])
SEWFe_down <- data.frame(mylist[3])
#==============================Stat_effect_without_Fe========================
path <- "networks/verbose/Stat_effect_without_Fe.txt"
dir <-  "Stat_effect_without_Fe_oslom/"
file_up <-  "Stat_effect_without_Fe_up.dat"
file_down <-  "Stat_effect_without_Fe_down.dat"

mylist2 <- oslom_table(path=path,var=SEWOFe,var_up=SEWOFe_up,var_down=SEWOFe_down,
                      dir=dir,file_up=file_up,
                      file_down = file_down)

SEWFe <- data.frame(mylist2[1])
SEWFe_up <- data.frame(mylist2[2])
SEWFe_down <- data.frame(mylist2[3])

#==============================Stat_Fedeprivation_interaction=================
path <- "networks/verbose/Stat_Fedeprivation_interaction.txt"
dir <-  "Stat_Fedeprivation_interaction/"
file_up <-  "Stat_Fedeprivation_interaction_up.dat"
file_down <-  "Stat_Fedeprivation_interaction_down.dat"

mylist3 <- oslom_table(path=path,var=SFeI,var_up=SFeI_up,var_down=SFeI_down,
                       dir=dir,file_up=file_up,
                       file_down = file_down)

SEWFe <- data.frame(mylist3[1])
SEWFe_up <- data.frame(mylist3[2])
SEWFe_down <- data.frame(mylist3[3])

#======================************************=================================
#===============Getting tables of genes A and genes B===========================
#======================************************=================================

order_tab_osl <- function(data,name_table,dir,file) {
  var <- data[c(1,2,3,4,7,8,9,10)]
  var_osl_A <- var[c(1,3,5,7)]
  var_osl_B <- var[c(2,4,6,8)]
  
  names(var_osl_A) <- names(var_osl_B) 
  name_table <- rbind (var_osl_A,var_osl_B)
  colnames(name_table) <- c("RV_gene","number_gene","Gene_name","Annotation_gene")
  name_table <- name_table[!duplicated(name_table[,"number_gene"]),]
  
  dir.create(dir,showWarnings = FALSE)
  
  file <- paste0(dir,file)
  write.table(name_table,file=file)
  return(name_table)
    
}

dir <-  "order_oslom/"
file <-  "IDEE_osl.txt"

IDEE_osl <- order_tab_osl(IDEE,IDEE_osl,dir=dir,file=file)
# 
# file <-  "IDES_osl.txt"
# IDES_osl <- order_tab_osl(IDES,IDES_osl,dir=dir,file=file)
# 
# file <-  "SEWFe_osl.txt"
# SEWFe_osl <- order_tab_osl(SEWFe,SEWFe_osl,dir=dir,file=file)
# oslom_read <- read.table("Iron_deprivation_effect_at_Exp.dat")
# colnames(oslom_read) <- oslom_read[1,]
# oslom_read <- oslom_read[2:nrow(oslom_read),]
# oslom_read$number_geneA <- as.double(oslom_read$number_geneA)
# oslom_read$number_geneB <- as.double(oslom_read$number_geneB)
# oslom_read$censored_weight_down <- as.double(oslom_read$censored_weight_down)
# safe=oslom_read
# 
# oslom_read$censored_weight_down[which(oslom_read$censored_weight_down==0)] <- 1e-10
# oslom_read$censored_weight_down <- abs(oslom_read$censored_weight_down)
# 
# write.table(oslom_read,file="Iron_deprivation_effect_at_Exp.dat",row.names = FALSE,
#             col.names = FALSE)
