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
  
  # Check if all the genes have "gene-"
  beginings=substr(df[,column], start,stop)
  length(which(beginings=="gene-"))
  length(beginings)
  #Remove the "gene-" from Gene_ID column
  RV_codes=substr(df[,column], (stop+1),nchar(df[,column]))
  df[,column]=RV_codes
  colnames(df)[col_x] <- "RV_codes"
  return(df)
  
}

### Create a Data frame from an xlsx
main_folder <- "inputs/raw"

## JQ: OJO: La única razón para cargar df2 aquí es que lo necesitamos para los featuredata PERO LOS DATOS NO SE USAN: EL PAPER VA SOLO DE C5 y C6

df1=data.frame(read_excel(file.path(main_folder,"C5_C6_rep1.xlsx")))
df2=data.frame(read_excel(file.path(main_folder,"C1_C2_rep1.xlsx")))
# df3=data.frame(read_excel(file.path(main_folder,"C7_C8_rep1.xlsx")))

#=================================================================================================================================
### Cleanup RV gene ids column, and make sure that the result is made of unique IDs, and that they match across input files.======
### ==============================================================================================================================

one <- rm_word(df=df1, column="Gene_ID",start=1,stop=5,col_x=1)
two <- rm_word(df=df2, column="Feature_ID",start=1,stop=5,col_x=1)

# three$Gene_ID=three$locus_tag
# three <- three[order(three$Gene_ID),]

length(unique(one$RV_codes))
length(which(!one$RV_codes==two$RV_codes))
length(which(one$RV_codes==two$RV_codes))
# length(which(two$Gene_ID==three$Gene_ID))
# length(which(three$Gene_ID==one$Gene_ID))

#=================================================================================================================================
### Make sure that Chrom, Start, End, Strand Type, Locus_tag, Entrez_ID,    ######################################################
### Gene_Symbol & Protein_ID & product behave coherently across input files.######################################################
#=================================================================================================================================

for(i in c(2:ncol(one))){
  print((colnames(one)[i]))
  print(paste0("ones is different than two at:",length(which(one[,i]!=two[,i]))))
  }

## JQ: Notas tras observar esto, para antes del merging:
{
## Los RV_codes están OK, que es lo más importante,
## El nombre de,l chromosoma está ok en two:
## unique(two$Chrom) es "NC_000962.3", el código del assembly
## pero no en one, asi que:
one$Chrom=two$Chrom
## Las columnas conteniendo gene_symbol, protein_IDs, descriptions and biotypes se llaman diferente y están en columnas diferentes, así que, primero chequeo:

length(which(one$Gene_Symbol==two$Gene_Symbol))
# 1998. EN muchas, difieren, veamos qué pasa en estas casi dos mil:

one[which(one$Gene_Symbol!=two$Gene_Symbol),1:11]
two[which(one$Gene_Symbol!=two$Gene_Symbol),1:11]

## Haciendo eso veo que en todo caso, son símbolos que están ausentes en 1, y presentes en 2, así que plos feature data buenos, en esto, son los de two

# Lo mismo pasa con biotipo, descripción y protein_ID: más facil de ver, pues son columnas vacías en one:

unique(one$Protein_ID)
#[1] "."
unique(one$Description)
#[1] "."
unique(one$gene_biotype)
#[1] "."
## Lo cual indica que seguir como sigues aquí está bien ;)
}
###feature data ## JQ: Quito locus_tab: es redundante, y pongo el orden limpio ya aqui
feature_data=two[,c("RV_codes","Entrez_Gene_ID","Gene_Symbol","Protein_ID","Product","Type","Chrom","Strand","Start","End")]
feature_data$length=abs(feature_data$End-feature_data$Start)

### fpkm and read count matrix
# reads=cbind(one[,11:18],two[,12:19],three[,11:18])
## JQ: Aquí estás incluyendo más datos de la cuenta!
{
fpkm=one[,19:26]
reads=one[,11:18]
rownames(reads)=feature_data$RV_codes
rownames(fpkm)=feature_data$RV_codes


metadata=data.frame(Original_ID=colnames(reads))
metadata$Original_ID=substr(metadata$Original_ID, 2,(nchar(metadata$Original_ID)-11))
metadata$Iron=c("NO","NO","YES","YES","NO","NO","YES","YES")
metadata$Dextrose=c(rep("NO",8))
metadata$Growth=c(rep("EXP",4),rep("STAT",4))
metadata$LCFA=rep("NO",8)
metadata$setup=paste0("Fe.",metadata$Iron,"_Dextrose.",metadata$Dextrose,"_Growth.",metadata$Growth,"_LCFA.",metadata$LCFA)
metadata$Culture=c(rep("C5",4),rep("C6",4))
rownames(metadata)=paste0(metadata$Culture,"_Fe_",metadata$Iron,"_",rep(c(1,2),4))


reads=reads[,order(metadata$setup)]
fpkm=fpkm[,order(metadata$setup)]
metadata=metadata[order(metadata$setup),]

length(which(paste0("X",metadata$Original_ID,"_Read_Count")==colnames(reads)))
length(which(paste0("X",metadata$Original_ID,"_FPKM")==colnames(fpkm)))

colnames(reads)=rownames(metadata)
colnames(fpkm)=rownames(metadata)

metadata=metadata[,c(2,4,7,6,1,3,5)]
colnames(feature_data)[1]="Gene_ID"

##JQ Los saco todos en .txt y .xlsx. Ojo al recursive=TRUE
dir.create("inputs/processed/unfiltered/txts",recursive=TRUE)
dir.create("inputs/processed/unfiltered/xlsx",recursive=TRUE)

write_excel_with_rownames=function(tab,path,name_row_col){
    tab=cbind(tab,rownames(tab))
    tab=tab[,c(ncol(tab):(1:(ncol(tab)-1)))]
    colnames(tab)[1]=name_row_col
    write_xlsx(tab,path)
}

write_excel_with_rownames(metadata,"inputs/processed/unfiltered/xlsx/metadata.xlsx","Sample")
write_excel_with_rownames(feature_data,"inputs/processed/unfiltered/xlsx/feature_data.xlsx","Gene")
write_excel_with_rownames(reads,"inputs/processed/unfiltered/xlsx/reads.xlsx","Gene")
write_excel_with_rownames(fpkm,"inputs/processed/unfiltered/xlsx/fpkm.xlsx","Gene")

write.table(metadata,"inputs/processed/unfiltered/txts/metadata.txt")
write.table(feature_data,"inputs/processed/unfiltered/txts/feature_data.txt")
write.table(reads,"inputs/processed/unfiltered/txts/reads.txt")
write.table(fpkm,"inputs/processed/unfiltered/txts/fpkm.txt")

## JQ Resumen de ediciones aquí: https://capture.dropbox.com/z7Rqg8nuy4XAmvrb
