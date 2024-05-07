### Libraries

library(readxl)
library(writexl)
library(dplyr)

### Create a Data frame from an xlsx

one=data.frame(read_excel("inputs/raw/C5_C6_rep1.xlsx"))
two=data.frame(read_excel("inputs/raw/C1_C2_rep1.xlsx"))
three=data.frame(read_excel("inputs/raw/C7_C8_rep1.xlsx"))


### Cleanup RV gene ids column, and make sure that the result is made of unique IDs, and that they match across input files.

beginings=substr(one$Gene_ID, 1,5)
length(which(beginings=="gene-"))
length(beginings)
RV_codes=substr(one$Gene_ID, 6,nchar(one$Gene_ID))
one$Gene_ID=RV_codes

beginings=substr(two$Feature_ID, 1,5)
length(which(beginings=="gene-"))
length(beginings)
RV_codes=substr(two$Feature_ID, 6,nchar(two$Feature_ID))
two$Feature_ID=RV_codes

three$Gene_ID=three$locus_tag
three <- three[order(three$Gene_ID),]


length(unique(one$Gene_ID))
length(which(one$Gene_ID==two$Gene_ID))
length(which(two$Gene_ID==three$Gene_ID))
length(which(three$Gene_ID==one$Gene_ID))
### Make sure that Chrom, Start, End, Strand Type, Locus_tag, Entrez_ID,
### Gene_Symbol & Protein_ID & product behave coherently across input files.

for(i in c(2:11)){
  print(i)
  print(length(which(one[,i]!=three[,i])))}

###feature data
feature_data=two[,c(1:11)]
feature_data$length=abs(feature_data$End-feature_data$Start)

### read count matrix
reads=cbind(one[,11:18],two[,12:19],three[,11:18])
fpkm=cbind(one[,19:26],two[,20:27],three[,19:26])


rownames(reads)=feature_data$Feature_ID
rownames(fpkm)=feature_data$Feature_ID


metadata=data.frame(Original_ID=colnames(reads))
metadata$Original_ID=substr(metadata$Original_ID, 2,(nchar(metadata$Original_ID)-11))
metadata$Iron=c("NO","NO","YES","YES","NO","NO","YES","YES","NO","YES","NO","YES","NO","NO","YES","YES",
                "YES","YES","NO","NO","YES","YES","NO","NO")
metadata$Dextrose=c(rep("NO",8),rep("YES",8),rep("NO",8))

metadata$Growth=rep(c(rep("EXP",4),rep("STAT",4)),3)

metadata$LCFA=rep(c(rep("NO",16),rep("YES",8)),1)

metadata$setup=paste0("Fe.",metadata$Iron,"_Dextrose.",metadata$Dextrose,"_Growth.",metadata$Growth,"_LCFA.",metadata$LCFA)


reads=reads[,order(metadata$setup)]
fpkm=fpkm[,order(metadata$setup)]
metadata=metadata[order(metadata$setup),]

length(which(paste0("X",metadata$Original_ID,"_Read_Count")==colnames(reads)))
length(which(paste0("X",metadata$Original_ID,"_FPKM")==colnames(fpkm)))

metadata$Culture=rep(c("C5","C5","C7","C7","C6","C6","C8","C8","C1","C1","C2","C2"),2)
rownames(metadata)=paste0(metadata$Culture,"_Fe_",metadata$Iron,"_",rep(c(1,2),8))

colnames(reads)=rownames(metadata)
colnames(fpkm)=rownames(metadata)

metadata=metadata[,c(2,3,4,5,7,1,6)]
colnames(feature_data)[1]="Gene_ID"
feature_data=feature_data[,c(1,8,9,10,11,2,3,4,12,5,6,7)]

dir.create("inputs/all_processed")
write_xlsx(metadata,"inputs/processed/all_metadata.xlsx")

write.table(metadata,"inputs/all_processed/metadata.txt")
write.table(feature_data,"inputs/all_processed/feature_data.txt")
write.table(reads,"inputs/all_processed/reads.txt")
write.table(fpkm,"inputs/all_processed/fpkm.txt")

