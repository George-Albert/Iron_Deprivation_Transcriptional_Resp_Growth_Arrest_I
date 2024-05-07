### Libraries

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
# library(tibble)

### Read the data
metadata = read.table("inputs/processed/metadata.txt")
feature_data = read.table("inputs/processed/feature_data.txt")
reads = read.table("inputs/processed/reads.txt")
fpkms = read.table("inputs/processed/fpkm.txt")

### number of samples and genes
n_samples <- ncol(reads)
n_genes <- nrow (reads)
n_Fe_NO <- 8
n_Fe_Yes <- n_samples
### remove last part of the ID culture
samplenames <- substring(colnames(reads[,1:n_Fe_NO]),1,n_Fe_NO)
samplenames[9:ncol(reads)]<- substring(colnames(reads[,(n_Fe_NO+1):n_Fe_Yes]),1,9)
colnames(reads) <- samplenames 


group <- factor(colnames(reads)[1:ncol(reads)])

y <- DGEList(counts = reads[1:ncol(reads)],group = group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)

### perform quasi-likelihood F-tests:

fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef = 2)
topTags(qlf)

### perform likelihood ratio tests:

# fit<-glmFit(y,design)
# lrt<-glmLRT(fit,coef=2)
# topTags(lrt)



