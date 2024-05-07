### Libraries

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

### Read the data
metadata = read.table("inputs/all_processed/metadata.txt")
feature_data = read.table("inputs/all_processed/feature_data.txt")
reads = read.table("inputs/all_processed/reads.txt")
fpkms = read.table("inputs/all_processed/fpkm.txt")


metadata$short_setup=paste0(metadata$Culture,"_Fe_",metadata$Iron)
### number of samples and genes
n_samples <- ncol(reads)
n_genes <- nrow (reads)
n_Fe_NO <- 12
n_Fe_Yes <- n_samples
### remove last part of the ID culture
samplenames <- substring(colnames(reads[,1:n_Fe_NO]),1,8)
samplenames[(n_Fe_NO+1):n_samples]<- substring(colnames(reads[,(n_Fe_NO+1):n_Fe_Yes]),1,9)
colnames(reads) <- samplenames 


group <- factor(colnames(reads)[1:ncol(reads)])

y <- DGEList(counts = reads[1:ncol(reads)],group = group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~group)
# y <- estimateDisp(y,design)

v=voom(y,design,plot=TRUE)


exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance[,c(1:5)])

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos[,1:7])

datos$Culture=factor(datos$Culture,levels=c("C5","C6","C7","C8","C1","C2"))
datos$short_setup=factor(datos$short_setup,levels=c("C5_Fe_YES","C5_Fe_NO",
                                                    "C6_Fe_YES","C6_Fe_NO",
                                                    "C7_Fe_YES","C7_Fe_NO",
                                                    "C8_Fe_YES","C8_Fe_NO",
                                                    "C1_Fe_YES","C1_Fe_NO",
                                                    "C2_Fe_YES","C2_Fe_NO"))
datos=datos[order(datos$short_setup),]
colores_base=c("red","forestgreen","dodgerblue","black","aquamarine2","chocolate")
fill_base=c("red",NA,"forestgreen",NA,"dodgerblue",NA,"black",NA,"aquamarine2",NA,
            "chocolate",NA)

###PCA plot PC1 vs PC2
pl=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1.5)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab("PC1: 50.3% variance ~growth phase")+ylab("PC2: 14.8% variance ~Dextrose, Iron (at stat phase only)")


dir.create("outputs_JAC")
pdf("outputs_JAC/pca_1_2.pdf",width=6,height=5)
print(pl)
dev.off()

###PCA plot PC1 vs PC3
pl=ggplot(datos)+geom_point(aes(x=PC1,y=PC3,color=Culture,fill=short_setup),shape=21,size=3,stroke=1.5)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab("PC1: 50.3% variance ~growth phase")+ylab("PC3: 14.8% variance ~Dextrose, Iron (at stat phase only)")

pdf("outputs_JAC/pca_1_3.pdf",width=6,height=5)
print(pl)
dev.off()

###PCA plot PC2 vs PC3

pl=ggplot(datos)+geom_point(aes(x=PC2,y=PC3,color=Culture,fill=short_setup)
                            ,shape=21,size=3,stroke=1.5)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab("PC2: 50.3% variance ~growth phase")+ylab("PC3: 14.8% variance ~Dextrose, Iron (at stat phase only)")

pdf("outputs_JAC/pca_2_3.pdf",width=6,height=5)
print(pl)
dev.off()
