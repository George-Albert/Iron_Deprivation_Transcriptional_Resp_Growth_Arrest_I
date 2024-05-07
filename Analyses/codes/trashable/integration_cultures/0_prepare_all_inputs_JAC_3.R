### ############
### libraries###
### ############
{
library(readxl)
library(writexl)
library(dplyr)
}

### ############
### Functions###
### ############

dcols=function(x){data.frame(colnames(x))}
rm_word  <- function(df,column,start,stop,col_x) {
  
  if (column == "locus_tag" || column =="Locus_Tag") {
    
    RV_codes <- df[,column]
    # df$RV_codes <- RV_codes
    df <- data.frame(RV_codes = RV_codes, df)
    
    } else{
    
      # Check if all the genes have "gene-"
    beginings=substr(df[,column], start,stop)
    length(which(beginings=="gene-"))
    length(beginings)
    #Remove the "gene-" from Gene_ID column
    RV_codes=substr(df[,column], (stop+1),nchar(df[,column]))
    df[,column]=RV_codes
    colnames(df)[col_x] <- "RV_codes"
    
    }
  return(df)
}

#*******************************************************************************
#Exponential phase

#*C1 <- Glycerol, Dextrose, NO LCFA
#*C5 <- Glycerol, No Dextrose, NO LCFA
#*C7 <- Glycerol, Dextrose, LCFA
#*C11 <- NO Glycerol, NO Dextrose, LCFA
#**C14 <- NULL

#Stationary phase

#*C2 <- Glycerol, Dextrose, no LCFA
#*C6 <- Glycerol, No Dextrose, NO LCFA
#*C8 <- Glycerol, Dextrose, LCFA
#*C10 <- Glycerol, NO Dextrose, LCFA
#*C12 <- NO Glycerol, NO Dextrose, LCFA
#*C15 <- NULL

#long Stationary phase

#*C10 <- Glycerol, NO Dextrose, LCFA
#*C13 <- NO Glycerol, NO Dextrose, LCFA

#*******************************************************************************

### load a list of files

main_folder <- "inputs/raw"
my_files <- list.files(path = main_folder,pattern = "xlsx")
filenames <-substr(my_files,1,nchar(my_files)-5)
######==========================================================================
############################Load all files######################################
###### =========================================================================


# for(i in filenames){
#   filepath <- file.path(main_folder,paste0(i,".xlsx"))
#   assign(paste0("df",i), data.frame(read_excel(filepath)))
# }

all_df <- lapply(my_files,function(x){
  data.frame(read_excel(file.path(main_folder,x)))})

# names(all_df) <- gsub(".xlsx","",list.files(main_folder,full.names = FALSE))

q <- c(1:length(all_df)) #Amount of data frames we have

for(i in q){
  names(all_df)[i] <- paste0("df",i)
}

list2env(all_df, envir = .GlobalEnv)


###===============================================================================================================================
### Cleanup RV gene ids column, and make sure that the result is made of unique IDs, and that they match across input files.======
### ==============================================================================================================================
{
#C1_C2_rep1
one <- rm_word(df=df1, column="Feature_ID",start=1,stop=5,col_x=1)
#C1_C2_rep2
two <- rm_word(df=df2, column="locus_tag",start=1,stop=5,col_x=1)
colnames(two)[12:ncol(two)] <- colnames(one)[12:ncol(one)]   
#C11_C12_rC2_rep1
three <- rm_word(df=df3, column="locus_tag",start=1,stop=5,col_x=1)
#C11_C12_rC2_rep2
four <- rm_word(df=df4, column="locus_tag",start=1,stop=5,col_x=1)
#C14_C15_rC8_rep1
five <- rm_word(df=df5, column="locus_tag",start=1,stop=5,col_x=1)
#C14_C15_rC8_rep2
six <- rm_word(df=df6, column="locus_tag",start=1,stop=5,col_x=1)
#C5_C6_rep1
seven <- rm_word(df=df7, column="Gene_ID",start=1,stop=5,col_x=1)
#C5_C6_rep2
eight <- rm_word(df=df8, column="Gene_ID",start=1,stop=5,col_x=1)
#C7_C8_rep1
nine <- rm_word(df=df9, column="locus_tag",start=1,stop=5,col_x=1)
# #C7_C8_rep2
# ten <- rm_word(df=df10, column="locus_tag",start=1,stop=5,col_x=1)
}

### Order df
{
one <- one[order(one$RV_codes),]
two <- two[order(two$RV_codes),]
three <- three[order(three$RV_codes),]
four <- four[order(four$RV_codes),]
five <- five[order(five$RV_codes),]
six <- six[order(six$RV_codes),]
seven <- seven[order(seven$RV_codes),]
eight <- eight[order(eight$RV_codes),]
nine <- nine[order(nine$RV_codes),]
# ten <- ten[order(ten$RV_codes),]
}
#=================================================================================================================================
### Check the RV codes to be the sames.######################################################
#=================================================================================================================================
genes_out <- one[which(!one$RV_codes %in% three$RV_codes),]
genes_out$RV_codes

genes_out <- one[which(!one$RV_codes %in% four$RV_codes),]
genes_out$RV_codes

genes_out <- (one[which(!one$RV_codes %in% five$RV_codes),])
dim(genes_out)[1]
genes_out$RV_codes

genes_out <- one[which(!one$RV_codes %in% six$RV_codes),]
dim(genes_out)[1]
genes_out$RV_codes

### "Rv2280"  "Rv3476c" are removed in three (C11_C12_rC2_rep1),four (C14_C15_rC8_rep1),
### five (C11_C12_rC2_rep2),six (C14_C15_rC8_rep2)
add_gene <- function(df) {
  
  df[nrow(df) + 1,] = genes_out[1,]
  df[nrow(df) + 1,] = genes_out[2,]
  df[which(df$RV_codes=="Rv2280"),c(12:35)] <- "-"
  df[which(df$RV_codes=="Rv3476c"),c(12:35)] <- "-"
  df <- df[order(df$RV_codes),]
  return(df)
}

three <- add_gene(three)
four <- add_gene(four)
five <- add_gene(five)
six <- add_gene(six)


cul_list <- list(one,two,three,four,five,six,seven,eight,nine)

for(i in (2:length(cul_list))){
  
  print(colnames(cul_list[[i]][1]))
  print(paste0("ones is different at:",length(which(cul_list[[1]][1]!=cul_list[[i]][1]))))
  }

###feature data
feature_data=one[,c("RV_codes","Chrom","Start","End","Strand","Type","Locus_Tag",
                    "Entrez_Gene_ID","Gene_Symbol","Protein_ID", "Product")]
feature_data$length=abs(feature_data$End-feature_data$Start)

#=================================================================================================================================
### Change the names of the colnames by adding rep number.######################################################
#=================================================================================================================================
for (i in (12:ncol(one))) {
  
  colnames(one)[i] <- paste0(colnames(one)[i],"_rep1")
  
  colnames(two)[i] <- paste0(colnames(two)[i],"_rep2")
}


  colnames(three)[12:15] <- paste0(colnames(three)[12:15],"_rep3")
  colnames(three)[16:23] <- paste0(colnames(three)[16:23],"_rep1")
  colnames(three)[24:27] <- paste0(colnames(three)[24:27],"_rep3")
  colnames(three)[28:ncol(three)] <- paste0(colnames(three)[28:ncol(three)],"_rep1")
  
  colnames(four)[12:15] <- paste0(colnames(four)[12:15],"_rep4")
  colnames(four)[16:23] <- paste0(colnames(four)[16:23],"_rep2")
  colnames(four)[24:27] <- paste0(colnames(four)[24:27],"_rep4")
  colnames(four)[28:ncol(four)] <- paste0(colnames(four)[28:ncol(four)],"_rep2")
  
  colnames(five)[12:15] <- paste0(colnames(five)[12:15],"_rep1")
  colnames(five)[16:19] <- paste0(colnames(five)[16:19],"_rep2")
  colnames(five)[20:27] <- paste0(colnames(five)[20:27],"_rep1")
  colnames(five)[28:31] <- paste0(colnames(five)[28:31],"_rep2")
  colnames(five)[32:ncol(five)] <- paste0(colnames(five)[32:ncol(five)],"_rep1")
  
  colnames(six)[12:15] <- paste0(colnames(six)[12:15],"_rep2")
  colnames(six)[16:19] <- paste0(colnames(six)[16:19],"_rep3")
  colnames(six)[20:27] <- paste0(colnames(six)[20:27],"_rep2")
  colnames(six)[28:31] <- paste0(colnames(six)[28:31],"_rep3")
  colnames(six)[32:ncol(six)] <- paste0(colnames(six)[32:ncol(six)],"_rep2")

for (i in (11:ncol(seven))) {
  
  colnames(seven)[i] <- paste0(colnames(seven)[i],"_rep1")
  
  colnames(eight)[i] <- paste0(colnames(eight)[i],"_rep2")
}

for (i in (12:ncol(nine))) {
  
  colnames(nine)[i] <- paste0(colnames(nine)[i],"_rep1")
  
  # colnames(ten)[i] <- paste0(colnames(ten)[i],"_rep2")
}

### fpkm and read count matrix
fpkm=cbind(one[,c("X10_3MMmADN_F1C1_FPKM_rep1","X10_3MMmADNmF1C1_FPKM_rep1",
                  "X16_3MMmADN_F2C1_FPKM_rep1","X16_3MMmADNmF2C1_FPKM_rep1",
                  "X23_4MMmADN_F1C2_FPKM_rep1","X23_4MMmADN_F2C2_FPKM_rep1",
                  "X23_4MMmADNmF1C2_FPKM_rep1","X23_4MMmADNmF2C2_FPKM_rep1")],
           two[,c("X10_3MMmADN_F1C1_FPKM_rep2","X10_3MMmADNmF1C1_FPKM_rep2",
                  "X16_3MMmADN_F2C1_FPKM_rep2","X16_3MMmADNmF2C1_FPKM_rep2",
                  "X23_4MMmADN_F1C2_FPKM_rep2","X23_4MMmADN_F2C2_FPKM_rep2",
                  "X23_4MMmADNmF1C2_FPKM_rep2","X23_4MMmADNmF2C2_FPKM_rep2")],
           three[,c("X19_01_ADN_F1C2r_FPKM_rep3","X19_01_ADN_F2C2r_FPKM_rep3",
                    "X19_01_ADNmF1C2r_FPKM_rep3","X19_01_ADNmF2C2r_FPKM_rep3",
                    "X19_01_AN_F1C11_FPKM_rep1","X19_01_AN_F1C12_FPKM_rep1",
                    "X19_01_AN_F2C11_FPKM_rep1","X19_01_AN_F2C12_FPKM_rep1",
                    "X19_01_ANmF1C11_FPKM_rep1","X19_01_ANmF1C12_FPKM_rep1",
                    "X19_01_ANmF2C11_FPKM_rep1","X19_01_ANmF2C12_FPKM_rep1")],
           four[,c("X19_01_ADN_F1C2r_FPKM_rep4","X19_01_ADN_F2C2r_FPKM_rep4",
                    "X19_01_ADNmF1C2r_FPKM_rep4","X19_01_ADNmF2C2r_FPKM_rep4",
                    "X19_01_AN_F1C11_FPKM_rep2","X19_01_AN_F1C12_FPKM_rep2",
                    "X19_01_AN_F2C11_FPKM_rep2","X19_01_AN_F2C12_FPKM_rep2",
                    "X19_01_ANmF1C11_FPKM_rep2","X19_01_ANmF1C12_FPKM_rep2",
                    "X19_01_ANmF2C11_FPKM_rep2","X19_01_ANmF2C12_FPKM_rep2")],
           five[,c("X03_12_AN_F1C14_FPKM_rep1","X03_12_AN_F2C14_FPKM_rep1",
                   "X03_12_ANmF1C14_FPKM_rep1","X03_12_ANmF2C14_FPKM_rep1",
                   "X18_03_AN_F1C8r_FPKM_rep2","X18_03_AN_F2C8r_FPKM_rep2",
                   "X18_03_ANmF1C8r_FPKM_rep2","X18_03_ANmF2C8r_FPKM_rep2",
                   "X29_11_AN_F1C15_FPKM_rep1","X29_11_AN_F2C15_FPKM_rep1",
                   "X29_11_ANmF1C15_FPKM_rep1","X29_11_ANmF2C15_FPKM_rep1")],
           six[,c("X03_12_AN_F1C14_FPKM_rep2","X03_12_AN_F2C14_FPKM_rep2",
                   "X03_12_ANmF1C14_FPKM_rep2","X03_12_ANmF2C14_FPKM_rep2",
                   "X18_03_AN_F1C8r_FPKM_rep3","X18_03_AN_F2C8r_FPKM_rep3",
                   "X18_03_ANmF1C8r_FPKM_rep3","X18_03_ANmF2C8r_FPKM_rep3",
                   "X29_11_AN_F1C15_FPKM_rep2","X29_11_AN_F2C15_FPKM_rep2",
                   "X29_11_ANmF1C15_FPKM_rep2","X29_11_ANmF2C15_FPKM_rep2")],
           seven[,c("X19_4_MMmAN_F1C5_FPKM_rep1","X19_4_MMmAN_F2C5_FPKM_rep1",
                   "X19_4_MMmANmF1C5_FPKM_rep1","X19_4_MMmANmF2C5_FPKM_rep1",
                   "X8_4_MMmAN_F1C6_FPKM_rep1","X8_4_MMmAN_F2C6_FPKM_rep1",
                   "X8_4_MMmANmF1C6_FPKM_rep1","X8_4_MMmANmF2C6_FPKM_rep1")],
           eight[,c("X19_4_MMmAN_F1C5_FPKM_rep2","X19_4_MMmAN_F2C5_FPKM_rep2",
                  "X19_4_MMmANmF1C5_FPKM_rep2","X19_4_MMmANmF2C5_FPKM_rep2",
                  "X8_4_MMmAN_F1C6_FPKM_rep2","X8_4_MMmAN_F2C6_FPKM_rep2",
                  "X8_4_MMmANmF1C6_FPKM_rep2","X8_4_MMmANmF2C6_FPKM_rep2")],
           nine[,c("X11_10_ANmF1C7_FPKM_rep1","X11_10_ANmF2C7_FPKM_rep1",
                    "X11_10_AN_F1C7_FPKM_rep1","X11_10_AN_F2C7_FPKM_rep1",
                    "X19_10_ANmF1C8_FPKM_rep1","X19_10_ANmF2C8_FPKM_rep1",
                    "X19_10_AN_F1C8_FPKM_rep1","X19_10_AN_F2C8_FPKM_rep1")])

reads=cbind(one[,c("X10_3MMmADN_F1C1_Read_Count_rep1","X10_3MMmADNmF1C1_Read_Count_rep1",
                   "X16_3MMmADN_F2C1_Read_Count_rep1","X16_3MMmADNmF2C1_Read_Count_rep1",
                   "X23_4MMmADN_F1C2_Read_Count_rep1","X23_4MMmADN_F2C2_Read_Count_rep1",
                   "X23_4MMmADNmF1C2_Read_Count_rep1","X23_4MMmADNmF2C2_Read_Count_rep1")],
            two[,c("X10_3MMmADN_F1C1_Read_Count_rep2","X10_3MMmADNmF1C1_Read_Count_rep2",
                   "X16_3MMmADN_F2C1_Read_Count_rep2","X16_3MMmADNmF2C1_Read_Count_rep2",
                   "X23_4MMmADN_F1C2_Read_Count_rep2","X23_4MMmADN_F2C2_Read_Count_rep2",
                   "X23_4MMmADNmF1C2_Read_Count_rep2","X23_4MMmADNmF2C2_Read_Count_rep2")],
            three[,c("X19_01_ADN_F1C2r_Read_Count_rep3","X19_01_ADN_F2C2r_Read_Count_rep3",
                     "X19_01_ADNmF1C2r_Read_Count_rep3","X19_01_ADNmF2C2r_Read_Count_rep3",
                     "X19_01_AN_F1C11_Read_Count_rep1","X19_01_AN_F1C12_Read_Count_rep1",
                     "X19_01_AN_F2C11_Read_Count_rep1","X19_01_AN_F2C12_Read_Count_rep1",
                     "X19_01_ANmF1C11_Read_Count_rep1","X19_01_ANmF1C12_Read_Count_rep1",
                     "X19_01_ANmF2C11_Read_Count_rep1","X19_01_ANmF2C12_Read_Count_rep1")],
            four[,c("X19_01_ADN_F1C2r_Read_Count_rep4","X19_01_ADN_F2C2r_Read_Count_rep4",
                     "X19_01_ADNmF1C2r_Read_Count_rep4","X19_01_ADNmF2C2r_Read_Count_rep4",
                     "X19_01_AN_F1C11_Read_Count_rep2","X19_01_AN_F1C12_Read_Count_rep2",
                     "X19_01_AN_F2C11_Read_Count_rep2","X19_01_AN_F2C12_Read_Count_rep2",
                     "X19_01_ANmF1C11_Read_Count_rep2","X19_01_ANmF1C12_Read_Count_rep2",
                     "X19_01_ANmF2C11_Read_Count_rep2","X19_01_ANmF2C12_Read_Count_rep2")],
            five[,c("X03_12_AN_F1C14_Read_Count_rep1","X03_12_AN_F2C14_Read_Count_rep1",
                    "X03_12_ANmF1C14_Read_Count_rep1","X03_12_ANmF2C14_Read_Count_rep1",
                    "X18_03_AN_F1C8r_Read_Count_rep2","X18_03_AN_F2C8r_Read_Count_rep2",
                    "X18_03_ANmF1C8r_Read_Count_rep2","X18_03_ANmF2C8r_Read_Count_rep2",
                    "X29_11_AN_F1C15_Read_Count_rep1","X29_11_AN_F2C15_Read_Count_rep1",
                    "X29_11_ANmF1C15_Read_Count_rep1","X29_11_ANmF2C15_Read_Count_rep1")],
            six[,c("X03_12_AN_F1C14_Read_Count_rep2","X03_12_AN_F2C14_Read_Count_rep2",
                    "X03_12_ANmF1C14_Read_Count_rep2","X03_12_ANmF2C14_Read_Count_rep2",
                    "X18_03_AN_F1C8r_Read_Count_rep3","X18_03_AN_F2C8r_Read_Count_rep3",
                    "X18_03_ANmF1C8r_Read_Count_rep3","X18_03_ANmF2C8r_Read_Count_rep3",
                    "X29_11_AN_F1C15_Read_Count_rep2","X29_11_AN_F2C15_Read_Count_rep2",
                    "X29_11_ANmF1C15_Read_Count_rep2","X29_11_ANmF2C15_Read_Count_rep2")],
            seven[,c("X19_4_MMmAN_F1C5_Read_Count_rep1","X19_4_MMmAN_F2C5_Read_Count_rep1",
                    "X19_4_MMmANmF1C5_Read_Count_rep1","X19_4_MMmANmF2C5_Read_Count_rep1",
                    "X8_4_MMmAN_F1C6_Read_Count_rep1","X8_4_MMmAN_F2C6_Read_Count_rep1",
                    "X8_4_MMmANmF1C6_Read_Count_rep1","X8_4_MMmANmF2C6_Read_Count_rep1")],
            eight[,c("X19_4_MMmAN_F1C5_Read_Count_rep2","X19_4_MMmAN_F2C5_Read_Count_rep2",
                   "X19_4_MMmANmF1C5_Read_Count_rep2","X19_4_MMmANmF2C5_Read_Count_rep2",
                   "X8_4_MMmAN_F1C6_Read_Count_rep2","X8_4_MMmAN_F2C6_Read_Count_rep2",
                   "X8_4_MMmANmF1C6_Read_Count_rep2","X8_4_MMmANmF2C6_Read_Count_rep2")],
            nine[,c("X11_10_ANmF1C7_Read_Count_rep1","X11_10_ANmF2C7_Read_Count_rep1",
                     "X11_10_AN_F1C7_Read_Count_rep1","X11_10_AN_F2C7_Read_Count_rep1",
                     "X19_10_ANmF1C8_Read_Count_rep1","X19_10_ANmF2C8_Read_Count_rep1",
                     "X19_10_AN_F1C8_Read_Count_rep1","X19_10_AN_F2C8_Read_Count_rep1")])

rownames(reads)=feature_data$RV_codes
rownames(fpkm)=feature_data$RV_codes


metadata=data.frame(Original_ID=substr(colnames(reads),1,(nchar(colnames(reads)))))
# metadata[9:16,] <- metadata[1:8,] # Giving the same names than C1_C2 rep1 to C1_C2 rep2  
{
#Create Iron column
metadata$Iron <- NA

for (i in (1:nrow(metadata))) {
  
  if (grepl("mF", metadata[[1]][i], fixed = TRUE)) {
    
    metadata$Iron[i]= "YES"
    
  } else{
    metadata$Iron[i]= "NO"
  }
}

#Create Dextrose column
metadata$Dextrose <- NA

for (i in (1:nrow(metadata))) {
  
  if (grepl("C11", metadata[[1]][i], fixed = TRUE)||grepl("C12", metadata[[1]][i], fixed = TRUE)||
      grepl("C14", metadata[[1]][i], fixed = TRUE)||grepl("C15", metadata[[1]][i], fixed = TRUE)||
      grepl("C7", metadata[[1]][i], fixed = TRUE)||grepl("C8", metadata[[1]][i], fixed = TRUE)) {
    
    metadata$Dextrose[i]= "NO"
    
  } else{
    if (grepl("C1", metadata[[1]][i], fixed = TRUE)||grepl("C2", metadata[[1]][i], fixed = TRUE)){
      
      metadata$Dextrose[i]= "YES"
      
    } else{
      metadata$Dextrose[i]= "NO"
    }
    
  }
}

#Create Growth column

metadata$Growth <- NA

for (i in (1:nrow(metadata))) {
  
  if (grepl("C2", metadata[i,1], fixed = TRUE)||grepl("C6", metadata[i,1], fixed = TRUE)||
      grepl("C8", metadata[i,1], fixed = TRUE)||grepl("C12", metadata[i,1], fixed = TRUE)||
      grepl("C15", metadata[i,1], fixed = TRUE)) {
    
    metadata$Growth[i]= "STAT"
    
  } else{
    
    if (grepl("C1", metadata[i,1], fixed = TRUE)||grepl("C5", metadata[i,1], fixed = TRUE)||
        grepl("C7", metadata[i,1], fixed = TRUE)||grepl("C11", metadata[i,1], fixed = TRUE)||
        grepl("C14", metadata[i,1], fixed = TRUE)) {
      
      metadata$Growth[i]= "EXP"
      
    } else{
      metadata$Growth[i]= "STAT"
    }
    
    }
}


#Create LCFA column

metadata$LCFA <- NA

for (i in (1:nrow(metadata))) {
  if (grepl("C1_", metadata[i,1], fixed = TRUE)||grepl("C5", metadata[i,1], fixed = TRUE)||
      grepl("C2", metadata[i,1], fixed = TRUE)||grepl("C6", metadata[i,1], fixed = TRUE)||
      grepl("C14", metadata[i,1], fixed = TRUE)||grepl("C15", metadata[i,1], fixed = TRUE)){
    
    metadata$LCFA[i]= "NO"
    
    } else{
      
    if (grepl("C7", metadata[i,1], fixed = TRUE)||grepl("C8", metadata[i,1], fixed = TRUE)||
        grepl("C11", metadata[i,1], fixed = TRUE)||grepl("C12", metadata[i,1], fixed = TRUE)){
      
      metadata$LCFA[i]= "YES"
      
    } else{

      metadata$LCFA[i]= "YES"
    }
    }
}

#Create Glycerol column

metadata$Glycerol <- NA

for (i in (1:nrow(metadata))) {
  if (grepl("C11", metadata[i,1], fixed = TRUE)||grepl("C12", metadata[i,1], fixed = TRUE)||
      grepl("C14", metadata[i,1], fixed = TRUE)||grepl("C15", metadata[i,1], fixed = TRUE)){
    
    metadata$Glycerol[i]= "NO"
    
  } else{
    
    metadata$Glycerol[i]= "YES"
    
    }
}

#Create Experiment order column

metadata$Order_exp <- NA


for (i in (1:nrow(metadata))) {
  
  if (grepl("rep3", metadata[[1]][i], fixed = TRUE)||grepl("rep3", metadata[[1]][i], fixed = TRUE)) {
    
    metadata$Order_exp[i]= 3
    
  } else if (grepl("rep2", metadata[[1]][i], fixed = TRUE)||grepl("rep2", metadata[[1]][i], fixed = TRUE)){
    
    metadata$Order_exp[i]= 2
    
  } else if (grepl("rep1", metadata[[1]][i], fixed = TRUE)||grepl("rep1", metadata[[1]][i], fixed = TRUE)){
    
    metadata$Order_exp[i]= 1
  } else if (grepl("rep4", metadata[[1]][i], fixed = TRUE)||grepl("rep4", metadata[[1]][i], fixed = TRUE)){
    
    metadata$Order_exp[i]= 4
  }
  
}
#Create technical replicate column

metadata$Technical_replicate <- NA

for (i in (1:nrow(metadata))) {
  
  if (grepl("rep2", metadata[[1]][i], fixed = TRUE)||grepl("rep2", metadata[[1]][i], fixed = TRUE)) {
    
    metadata$Technical_replicate[i]= 1
    
  } else if (grepl("rep1", metadata[[1]][i], fixed = TRUE)||grepl("rep1", metadata[[1]][i], fixed = TRUE)){
    
    metadata$Technical_replicate[i]= 0
    
  } else if (grepl("rep3", metadata[[1]][i], fixed = TRUE)||grepl("rep3", metadata[[1]][i], fixed = TRUE)){
    
    metadata$Technical_replicate[i]= 0
  } else if (grepl("rep4", metadata[[1]][i], fixed = TRUE)||grepl("rep4", metadata[[1]][i], fixed = TRUE)){
    
    metadata$Technical_replicate[i]= 0
  }
  
}

#Create biological replicate column

metadata$Biological_replicate <- NA

for (i in (1:nrow(metadata))) {
  
  if (grepl("rep3", metadata[[1]][i], fixed = TRUE)||grepl("rep3", metadata[[1]][i], fixed = TRUE)) {
    
    metadata$Biological_replicate[i]= 1
    
  } else if (grepl("rep4", metadata[[1]][i], fixed = TRUE)||grepl("rep4", metadata[[1]][i], fixed = TRUE)){
    
    metadata$Biological_replicate[i]= 1
    
  } else if (grepl("rep1", metadata[[1]][i], fixed = TRUE)||grepl("rep1", metadata[[1]][i], fixed = TRUE)){
    
    metadata$Biological_replicate[i]= 0
  } else if (grepl("rep2", metadata[[1]][i], fixed = TRUE)||grepl("rep2", metadata[[1]][i], fixed = TRUE)){
    
    metadata$Biological_replicate[i]= 0
  }
  
}
metadata$Myc_strain <- "H737Rv"

for (i in (1:nrow(metadata))) {
  
  if (grepl("C11_Read_Count_rep1", metadata[[1]][i], fixed = TRUE)||grepl("C12_Read_Count_rep1", metadata[[1]][i], fixed = TRUE)||
      grepl("C14_Read_Count_rep1", metadata[[1]][i], fixed = TRUE)||grepl("C15_Read_Count_rep1", metadata[[1]][i], fixed = TRUE)||
      grepl("C1_Read_Count_rep2", metadata[[1]][i], fixed = TRUE)||grepl("C2_Read_Count_rep2", metadata[[1]][i], fixed = TRUE)) {
    
    metadata$Myc_strain[i]= "7199-99"
    
  }
}

}
metadata$Original_ID=substr(metadata$Original_ID, 2,(nchar(metadata$Original_ID)-16))

metadata$setup=paste0("Fe.",metadata$Iron,"_Dextrose.",metadata$Dextrose,"_Growth.",
                      metadata$Growth,"_LCFA.",metadata$LCFA,"_Glycerol.",metadata$Glycerol)

length(which(paste0("X",metadata$Original_ID,"_Read_Count_rep")==colnames(reads)))
length(which(paste0("X",metadata$Original_ID,"_FPKM_rep")==colnames(fpkm)))

# metadata$Culture=c(rep("C5",4),rep("C11",4),rep("C14",4),rep("C6",4),rep("C12",4),
#                    rep("C15",4),rep("C1",4),rep("C7",2),rep("C2",8),rep("C8",6),
#                    rep("C5",4),rep("C11",4),rep("C14",4),rep("C6",4),rep("C12",4),
#                    rep("C15",4),rep("C1",4),rep("C7",2),rep("C2",8),rep("C8",6))
reads=reads[,order(metadata$setup)]
fpkm=fpkm[,order(metadata$setup)]
metadata=metadata[order(metadata$setup),]

name_cul <- unlist(strsplit(metadata$Original_ID,c("F1","F2") ))
name_cul <- unlist(strsplit(name_cul,c("r") ))


name_cul <- name_cul[c(F,T)]
metadata$Culture <- name_cul

rownames(metadata)=paste0(metadata$Culture,"_Fe_",metadata$Iron,"_rep",
                          metadata$Order_exp,"_",rep(c(1,2),44))

metadata$Sample_ID <- rownames(metadata)

colnames(reads)=rownames(metadata)
colnames(fpkm)=rownames(metadata)

metadata=metadata[,c("Iron","Dextrose","Growth","LCFA","Glycerol","Culture",
                     "Sample_ID","Order_exp","Technical_replicate",
                     "Biological_replicate","Myc_strain","Original_ID",
                     "setup")]

write.table(metadata,"inputs/all_processed_data/metadata_not_ordered.txt")

# order_cul <- str_sort((metadata$Culture), numeric = TRUE)
metadata <- metadata[order(as.numeric(gsub("C","",metadata$Culture))),]

colnames(feature_data)[1]="Gene_ID"

feature_data = feature_data[,c("Gene_ID","Entrez_Gene_ID","Gene_Symbol","Protein_ID",
                             "Product","Chrom","Start","End","length","Strand","Type",
                             "Locus_Tag")]

dir.create("inputs/all_processed_data", showWarnings = F)

write_xlsx(metadata,"inputs/all_processed_data/metadata.xlsx")

write.table(metadata,"inputs/all_processed_data/metadata.txt")
write.table(feature_data,"inputs/all_processed_data/feature_data.txt")
write.table(reads,"inputs/all_processed_data/reads.txt")
write.table(fpkm,"inputs/all_processed_data/fpkm.txt")

