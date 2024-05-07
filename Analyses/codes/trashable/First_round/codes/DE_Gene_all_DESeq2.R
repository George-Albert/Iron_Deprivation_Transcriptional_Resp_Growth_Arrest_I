### load packages
{
  library(tidyverse)
  library(reshape2)
  library(ggrepel)
  library(cowplot)
  library(edgeR)
  library(limma)
  library(Rtsne)
  library(umap)
  library(qvalue)
  library(DESeq2)
  library(tximport)
  library(readxl)
  library(writexl)
  library(apeglm)
  library(biomaRt)
  library(ashr)
  library(mashr)
  library(RColorBrewer)
  
  options(width=1000)
}
### Read Data
metadata=read.table("inputs/processed_all_JAC/metadata.txt")
feature_data=read.table("inputs/processed_all_JAC/feature_data.txt")
reads=read.table("inputs/processed_all_JAC/reads.txt")
fpkms=read.table("inputs/processed_all_JAC/fpkm.txt")

value <- "EXP"
metadata<-metadata[which(metadata$Growth==value),]
metadata$short_setup=paste0(metadata$Culture,"_Fe_",metadata$Iron)
metadata <- metadata[order(metadata$short_setup),]
metadata <- metadata[order(metadata$Iron),]


reads <- reads[rownames(metadata)]
fpkms <- fpkms[rownames(metadata)]
### number of samples and genes
n_samples <- ncol(reads)
n_genes <- nrow (reads)
n_Fe_NO <- nrow(metadata[which(metadata$Iron=="NO"),])
n_Fe_Yes <- n_samples-n_Fe_NO

# ### remove last part of the ID culture
# samplenames <- substring(colnames(reads[,1:n_Fe_NO]),1,8)
# samplenames[(n_Fe_NO+1):n_samples]<- substring(colnames(reads[,(n_Fe_NO+1):n_Fe_Yes]),1,9)
# colnames(reads) <- samplenames 

all(colnames(reads) %in% rownames(metadata))
all(colnames(reads) == rownames(metadata))

### Create the design formula
design <- ~ short_setup

dds <- DESeqDataSetFromMatrix(countData = round(reads),colData = metadata,design= design)

dds <- DESeq(dds,test = "Wald")
res_nam <- resultsNames(dds) # lists the coefficients
### Iron_YES_vs_NO
res1 <- results(dds, name= res_nam[1], contrast = c("Iron","YES","NO"))
### Dextrose_YES_vs_NO
res2 <- results(dds, name=res_nam[2], contrast = c("Dextrose","YES","NO"))
### LCFA_YES_vs_NO
res3 <- results(dds, name=res_nam[3])

## This to make contrasts out of differences of coefficients
lista=list(add=c("condition_Cond1_vs_Control"),remove=c("condition_Cond2_vs_Control"))
res1_B=results(dds, contrast=lista)


## LogFC shrinkage (it does not affect the p values, or beta-standard errors)
res2_B <- lfcShrink(dds, contrast=lista, type="ashr")
res3_B <- lfcShrink(dds, coef="condition_PAM3CSK4_vs_NS", type="ashr")
