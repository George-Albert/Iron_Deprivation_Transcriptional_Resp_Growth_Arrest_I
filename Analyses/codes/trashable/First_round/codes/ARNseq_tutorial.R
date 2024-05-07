library(edgeR)
library(qvalue)
library(limma)
library(qvalue)
library(ggplot2)
library(cowplot)
library(reshape2)
library(DESeq2)
library(fdrtool)
library(ashr)
library(dplyr)
library(locfdr)
library(Glimma)
library(RNAseq123)

### Read the data
metadata = read.table("inputs/processed/metadata.txt")
feature_data = read.table("inputs/processed/feature_data.txt")
reads = read.table("inputs/processed/reads.txt")
fpkms = read.table("inputs/processed/fpkm.txt")


### Delete C2 culture from all the tables
cultures = list("C1","C2","C5","C6")
### number of samples and genes
n_samples <- ncol(reads)
n_genes <- nrow (reads)
### remove last part of the ID culture
samplenames[1:8] <- substring(colnames(reads[,1:8]),1,8)
samplenames[9:ncol(reads)]<- substring(colnames(reads[,9:16]),1,9)
colnames(reads) <- samplenames 

### See if is there is some gene repeted3

feature_data <- feature_data[!duplicated(feature_data$Entrez_Gene_ID),]

