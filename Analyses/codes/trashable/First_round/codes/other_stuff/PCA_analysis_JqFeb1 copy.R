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
metadata = read.table("inputs/processed_all_JAC/metadata.txt")
feature_data = read.table("inputs/processed_all_JAC/feature_data.txt")
reads = read.table("inputs/processed_all_JAC/reads.txt")
fpkms = read.table("inputs/processed_all_JAC/fpkm.txt")


metadata$short_setup=paste0(metadata$Culture,"_Fe_",metadata$Iron)
### number of samples and genes
n_samples <- ncol(reads)
n_genes <- nrow (reads)
n_Fe_NO <- length(which(metadata$Iron=="NO"))
n_Fe_Yes <- length(which(metadata$Iron=="YES"))
### remove last part of the ID culture: JQ Esto no conviene aquí: los colnames tienen que ser identificadores únicos de cada muetra: cuando tienes varias muestras que son replicados (con todos los atributos iguales) necesitas enumerarlas:
## samplenames <- substring(colnames(reads[,1:n_Fe_NO]),1,8)
## samplenames[(n_Fe_NO+1):n_samples]<- substring(colnames(reads[,(n_Fe_NO+1):n_Fe_Yes]),1,9)
## colnames(reads) <- samplenames



## Why did you attempted this?
## group <- factor(colnames(reads)[1:ncol(reads)])
#y <- DGEList(counts = reads[1:ncol(reads)],group = group)

y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
# y <- estimateDisp(y,design)

v=voom(y,design,plot=TRUE)


exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance[,c(1:5)])

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos[,1:7])

colores_base=c("forestgreen","black","red","dodgerblue","aquamarine2","chocolate")
fill_base=c("forestgreen",NA,"black",NA,"red",NA,"dodgerblue",NA,"aquamarine2",NA,
            "chocolate",NA)

###PCA plot PC1 vs PC2
pl1=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1.5)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab("PC1: 57.9% variance ~growth phase")+ylab("PC2: 18.5% variance ~Dextrose, Iron (at stat phase only)")


pdf("all.pdf",width=6,height=5)
print(pl1)
dev.off()

### Subset only C7 and C8
meta=metadata
counts=reads

metadata=metadata[which(metadata$Culture %in% c("C7","C8")),]
for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

reads=reads[,which(colnames(reads) %in% rownames(metadata))]

y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
# y <- estimateDisp(y,design)

v=voom(y,design,plot=TRUE)


exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance[,c(1:5)])

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos[,1:7])

datos$Culture=factor(datos$Culture,levels=c("C7","C8"))
datos$short_setup=factor(datos$short_setup,levels=c("C7_Fe_YES","C7_Fe_NO","C8_Fe_YES","C8_Fe_NO"))

colores_base=c("aquamarine2","chocolate")
fill_base=c("aquamarine2",NA,"chocolate",NA)

###PCA plot PC1 vs PC2
pl2=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1.5)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab("PC1: 54.4% variance ~growth phase")+ylab("PC2: 27.2% variance ~Dextrose, Iron (at stat phase only)")


## Remove C8 wout iron


metadata=meta[which(meta$short_setup %in% c("C7_Fe_YES","C7_Fe_NO","C8_Fe_YES")),]
for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

reads=counts[,which(colnames(counts) %in% rownames(metadata))]

y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
# y <- estimateDisp(y,design)

v=voom(y,design,plot=TRUE)


exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance[,c(1:5)])

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos)

datos$Culture=factor(datos$Culture,levels=c("C7","C8"))
datos$short_setup=factor(datos$short_setup,levels=c("C7_Fe_YES","C7_Fe_NO","C8_Fe_YES"))

colores_base=c("aquamarine2","chocolate")
fill_base=c("aquamarine2",NA,"chocolate")

###PCA plot PC1 vs PC2
pl3=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1.5)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab("PC1: 72.2% variance ~growth phase")+ylab("PC2: 12.6% variance ~Dextrose, Iron (at stat phase only)")



pdf("C7_C8FeYES_only.pdf",width=6,height=5)
print(pl3)
dev.off()


### Remove only C8 outlier:


metadata=meta[c(3,4,8,15,16,19,20),]
for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

reads=counts[,c(3,4,8,15,16,19,20)]

set_cult=c("C7","C8")
set=c("C7_Fe_YES","C7_Fe_NO","C8_Fe_YES","C8_Fe_NO")
# ("forestgreen","black","red","dodgerblue","aquamarine2","chocolate")
colores_base=c("aquamarine2","chocolate")
fill_base=c("aquamarine2",NA,"chocolate",NA)

y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
# y <- estimateDisp(y,design)

v=voom(y,design,plot=TRUE)


exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance)

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos)

datos$Culture=factor(datos$Culture,levels=set_cult)
datos$short_setup=factor(datos$short_setup,levels=set)


###PCA plot PC1 vs PC2
pl2b=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab(paste0("PC1: ",sum_pca[2,1]*100,"% variance ~growth phase"))+ylab(paste0("PC2: ",sum_pca[2,2]*100,"% variance ~growth phase"))
pdf("C7_C8_outlier_out.pdf",width=6,height=5)
print(pl2b)
dev.off()



### Remove only the other C8 sample:


metadata=meta[c(3,4,7,15,16,19,20),]
for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

reads=counts[,c(3,4,7,15,16,19,20)]

set_cult=c("C7","C8")
set=c("C7_Fe_YES","C7_Fe_NO","C8_Fe_YES","C8_Fe_NO")
# ("forestgreen","black","red","dodgerblue","aquamarine2","chocolate")
colores_base=c("aquamarine2","chocolate")
fill_base=c("aquamarine2",NA,"chocolate",NA)

y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
# y <- estimateDisp(y,design)

v=voom(y,design,plot=TRUE)


exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance)

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos)

datos$Culture=factor(datos$Culture,levels=set_cult)
datos$short_setup=factor(datos$short_setup,levels=set)


###PCA plot PC1 vs PC2
pl2c=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab(paste0("PC1: ",sum_pca[2,1]*100,"% variance ~growth phase"))+ylab(paste0("PC2: ",sum_pca[2,2]*100,"% variance ~growth phase"))
pdf("C7_C8_NONoutlier_out.pdf",width=6,height=5)
print(pl2c)
dev.off()






y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
# y <- estimateDisp(y,design)

v=voom(y,design,plot=TRUE)


exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance[,c(1:5)])

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos)

datos$Culture=factor(datos$Culture,levels=c("C7","C8"))
datos$short_setup=factor(datos$short_setup,levels=c("C7_Fe_YES","C7_Fe_NO","C8_Fe_YES","C8_Fe_NO"))

colores_base=c("aquamarine2","chocolate")
fill_base=c("aquamarine2",NA,"chocolate",NA)

###PCA plot PC1 vs PC2
pl2b=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1.5)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab("PC1: 57.7% variance ~growth phase")+ylab("PC2: 19.6% variance ~Dextrose, Iron (at stat phase only)")



pdf("C7_C8FeYES_only.pdf",width=6,height=5)
print(pl3)
dev.off()


## Remove C8 w iron


metadata=meta[which(meta$short_setup %in% c("C7_Fe_YES","C7_Fe_NO","C8_Fe_NO")),]
for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

reads=counts[,which(colnames(counts) %in% rownames(metadata))]

y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
# y <- estimateDisp(y,design)

v=voom(y,design,plot=TRUE)


exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance[,c(1:5)])

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos)

datos$Culture=factor(datos$Culture,levels=c("C7","C8"))
datos$short_setup=factor(datos$short_setup,levels=c("C7_Fe_YES","C7_Fe_NO","C8_Fe_NO"))

colores_base=c("aquamarine2","chocolate")
fill_base=c("aquamarine2",NA,NA)

###PCA plot PC1 vs PC2
pl4=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=0.5)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab("PC1: 61.7% variance ~growth phase")+ylab("PC2: 26.0% variance ~Dextrose, Iron (at stat phase only)")


pdf("C7_C8FeNO_only.pdf",width=6,height=5)
print(pl4)
dev.off()



## Only C7


metadata=meta[which(meta$short_setup %in% c("C7_Fe_YES","C7_Fe_NO")),]
for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

reads=counts[,which(colnames(counts) %in% rownames(metadata))]

y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
# y <- estimateDisp(y,design)

v=voom(y,design,plot=TRUE)


exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance)

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos)

datos$Culture=factor(datos$Culture,levels=c("C7"))
datos$short_setup=factor(datos$short_setup,levels=c("C7_Fe_YES","C7_Fe_NO"))

colores_base=c("aquamarine2")
fill_base=c("aquamarine2",NA)

###PCA plot PC1 vs PC2
pl5=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab("PC1: 53.0% variance ~growth phase")+ylab("PC2: 28.9% variance ~Dextrose, Iron (at stat phase only)")

pdf("C7.pdf",width=6,height=5)
print(pl5)
dev.off()


## Only C8


metadata=meta[which(meta$short_setup %in% c("C8_Fe_YES","C8_Fe_NO")),]
for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

reads=counts[,which(colnames(counts) %in% rownames(metadata))]

y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
# y <- estimateDisp(y,design)

v=voom(y,design,plot=TRUE)


exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance)

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos)

datos$Culture=factor(datos$Culture,levels=c("C8"))
datos$short_setup=factor(datos$short_setup,levels=c("C8_Fe_YES","C8_Fe_NO"))

colores_base=c("chocolate")
fill_base=c("chocolate",NA)

###PCA plot PC1 vs PC2
pl6=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab("PC1: 63.6% variance ~growth phase")+ylab("PC2: 24.4% variance ~Dextrose, Iron (at stat phase only)")


pdf("C8.pdf",width=6,height=5)
print(pl6)
dev.off()

## Only C2


metadata=meta[which(meta$short_setup %in% c("C2_Fe_YES","C2_Fe_NO")),]
for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

reads=counts[,which(colnames(counts) %in% rownames(metadata))]

y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
# y <- estimateDisp(y,design)

v=voom(y,design,plot=TRUE)


exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance)

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos)

datos$Culture=factor(datos$Culture,levels=c("C2"))
datos$short_setup=factor(datos$short_setup,levels=c("C2_Fe_YES","C2_Fe_NO"))
# ("forestgreen","black","red","dodgerblue","aquamarine2","chocolate")
colores_base=c("black")
fill_base=c("black",NA)

###PCA plot PC1 vs PC2
pl7=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab("PC1: 55.0% variance ~growth phase")+ylab("PC2: 31.4% variance ~Dextrose, Iron (at stat phase only)")

pdf("C2.pdf",width=6,height=5)
print(pl7)
dev.off()

## Only C1


metadata=meta[which(meta$short_setup %in% c("C1_Fe_YES","C1_Fe_NO")),]
for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

reads=counts[,which(colnames(counts) %in% rownames(metadata))]

y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
# y <- estimateDisp(y,design)

v=voom(y,design,plot=TRUE)


exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance)

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos)

datos$Culture=factor(datos$Culture,levels=c("C1"))
datos$short_setup=factor(datos$short_setup,levels=c("C1_Fe_YES","C1_Fe_NO"))
# ("forestgreen","black","red","dodgerblue","aquamarine2","chocolate")
colores_base=c("forestgreen")
fill_base=c("forestgreen",NA)

###PCA plot PC1 vs PC2
pl8=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab(paste0("PC1: ",sum_pca[2,1]*100,"% variance ~growth phase"))+ylab(paste0("PC2: ",sum_pca[2,2]*100,"% variance ~growth phase"))

pdf("C1.pdf",width=6,height=5)
print(pl8)
dev.off()


## Only C5


metadata=meta[which(meta$short_setup %in% c("C5_Fe_YES","C5_Fe_NO")),]
for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

reads=counts[,which(colnames(counts) %in% rownames(metadata))]

y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
# y <- estimateDisp(y,design)

v=voom(y,design,plot=TRUE)


exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance)

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos)

datos$Culture=factor(datos$Culture,levels=c("C5"))
datos$short_setup=factor(datos$short_setup,levels=c("C5_Fe_YES","C5_Fe_NO"))
# ("forestgreen","black","red","dodgerblue","aquamarine2","chocolate")
colores_base=c("red")
fill_base=c("red",NA)

###PCA plot PC1 vs PC2
pl9=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab(paste0("PC1: ",sum_pca[2,1]*100,"% variance ~growth phase"))+ylab(paste0("PC2: ",sum_pca[2,2]*100,"% variance ~growth phase"))



pdf("C5.pdf",width=6,height=5)
print(pl9)
dev.off()


## Only C6


metadata=meta[which(meta$short_setup %in% c("C6_Fe_YES","C6_Fe_NO")),]
for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

reads=counts[,which(colnames(counts) %in% rownames(metadata))]

y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
# y <- estimateDisp(y,design)

v=voom(y,design,plot=TRUE)


exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance)

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos)

datos$Culture=factor(datos$Culture,levels=c("C6"))
datos$short_setup=factor(datos$short_setup,levels=c("C6_Fe_YES","C6_Fe_NO"))
# ("forestgreen","black","red","dodgerblue","aquamarine2","chocolate")
colores_base=c("dodgerblue")
fill_base=c("dodgerblue",NA)

###PCA plot PC1 vs PC2
pl10=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab(paste0("PC1: ",sum_pca[2,1]*100,"% variance ~growth phase"))+ylab(paste0("PC2: ",sum_pca[2,2]*100,"% variance ~growth phase"))



pdf("C6.pdf",width=6,height=5)
print(pl10)
dev.off()


## Meaning sets: C1-c2    c5-c6  c7-c8   c5-c6-c7-c8

## Only C1-c2
set_cult=c("C1","C2")
set=c("C1_Fe_YES","C1_Fe_NO","C2_Fe_YES","C2_Fe_NO")
# ("forestgreen","black","red","dodgerblue","aquamarine2","chocolate")
colores_base=c("forestgreen","black")
fill_base=c("forestgreen",NA,"black",NA)

{

metadata=meta[which(meta$short_setup %in% set),]
for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

reads=counts[,which(colnames(counts) %in% rownames(metadata))]

y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
# y <- estimateDisp(y,design)

v=voom(y,design,plot=TRUE)


exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance)

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos)

datos$Culture=factor(datos$Culture,levels=set_cult)
datos$short_setup=factor(datos$short_setup,levels=set)
}

###PCA plot PC1 vs PC2
pl11=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab(paste0("PC1: ",sum_pca[2,1]*100,"% variance ~growth phase"))+ylab(paste0("PC2: ",sum_pca[2,2]*100,"% variance ~growth phase"))

pdf("C1_C2.pdf",width=6,height=5)
print(pl11)
dev.off()


{
    ## Only C1-c7
    set_cult=c("C1","C7")
    set=c("C1_Fe_YES","C1_Fe_NO","C7_Fe_YES","C7_Fe_NO")
    # ("forestgreen","black","red","dodgerblue","aquamarine2","chocolate")
    colores_base=c("forestgreen","aquamarine2")
    fill_base=c("forestgreen",NA,"aquamarine2",NA)

    {

    metadata=meta[which(meta$short_setup %in% set),]
    for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

    reads=counts[,which(colnames(counts) %in% rownames(metadata))]

    y <- DGEList(counts = reads)
    y <- calcNormFactors(y)
    design <- model.matrix(~short_setup,data=metadata)
    # y <- estimateDisp(y,design)

    v=voom(y,design,plot=TRUE)


    exp=v$E

    pca <-prcomp(t(exp),scale=T)
    sum_pca=data.frame(summary(pca)$importance)

    datos <- data.frame(pca$x)
    colnames(datos)=paste0("PC",c(1:ncol(datos)))
    length(which(rownames(datos)!=rownames(metadata)))
    datos=cbind(metadata,datos)

    datos$Culture=factor(datos$Culture,levels=set_cult)
    datos$short_setup=factor(datos$short_setup,levels=set)
    }

    ###PCA plot PC1 vs PC2
    pl11b=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab(paste0("PC1: ",sum_pca[2,1]*100,"% variance ~growth phase"))+ylab(paste0("PC2: ",sum_pca[2,2]*100,"% variance ~growth phase"))

    pdf("C1_C7.pdf",width=6,height=5)
    print(pl11b)
    dev.off()
    
}

{
    ## Only C5-c7
    set_cult=c("C5","C7")
    set=c("C5_Fe_YES","C5_Fe_NO","C7_Fe_YES","C7_Fe_NO")
    # ("forestgreen","black","red","dodgerblue","aquamarine2","chocolate")
    colores_base=c("red","aquamarine2")
    fill_base=c("red",NA,"aquamarine2",NA)

    {

    metadata=meta[which(meta$short_setup %in% set),]
    for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

    reads=counts[,which(colnames(counts) %in% rownames(metadata))]

    y <- DGEList(counts = reads)
    y <- calcNormFactors(y)
    design <- model.matrix(~short_setup,data=metadata)
    # y <- estimateDisp(y,design)

    v=voom(y,design,plot=TRUE)


    exp=v$E

    pca <-prcomp(t(exp),scale=T)
    sum_pca=data.frame(summary(pca)$importance)

    datos <- data.frame(pca$x)
    colnames(datos)=paste0("PC",c(1:ncol(datos)))
    length(which(rownames(datos)!=rownames(metadata)))
    datos=cbind(metadata,datos)

    datos$Culture=factor(datos$Culture,levels=set_cult)
    datos$short_setup=factor(datos$short_setup,levels=set)
    }

    ###PCA plot PC1 vs PC2
    pl11c=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab(paste0("PC1: ",sum_pca[2,1]*100,"% variance ~growth phase"))+ylab(paste0("PC2: ",sum_pca[2,2]*100,"% variance ~growth phase"))

    pdf("C5_C7.pdf",width=6,height=5)
    print(pl11c)
    dev.off()
    
}

{
    ## Only C5-c7
    set_cult=c("C1","C5")
    set=c("C1_Fe_YES","C1_Fe_NO", "C5_Fe_YES","C5_Fe_NO")
    # ("forestgreen","black","red","dodgerblue","aquamarine2","chocolate")
    colores_base=c("forestgreen","red")
    fill_base=c("forestgreen",NA,"red",NA)

    {

    metadata=meta[which(meta$short_setup %in% set),]
    for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

    reads=counts[,which(colnames(counts) %in% rownames(metadata))]

    y <- DGEList(counts = reads)
    y <- calcNormFactors(y)
    design <- model.matrix(~short_setup,data=metadata)
    # y <- estimateDisp(y,design)

    v=voom(y,design,plot=TRUE)


    exp=v$E

    pca <-prcomp(t(exp),scale=T)
    sum_pca=data.frame(summary(pca)$importance)

    datos <- data.frame(pca$x)
    colnames(datos)=paste0("PC",c(1:ncol(datos)))
    length(which(rownames(datos)!=rownames(metadata)))
    datos=cbind(metadata,datos)

    datos$Culture=factor(datos$Culture,levels=set_cult)
    datos$short_setup=factor(datos$short_setup,levels=set)
    }

    ###PCA plot PC1 vs PC2
    pl11e=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab(paste0("PC1: ",sum_pca[2,1]*100,"% variance ~growth phase"))+ylab(paste0("PC2: ",sum_pca[2,2]*100,"% variance ~growth phase"))

    pdf("C1_C5.pdf",width=6,height=5)
    print(pl11e)
    dev.off()
    
}


{
    ## Only C5-c7
    set_cult=c("C1","C5","C7")
    set=c("C1_Fe_YES","C1_Fe_NO", "C5_Fe_YES","C5_Fe_NO","C7_Fe_YES","C7_Fe_NO")
    # ("forestgreen","black","red","dodgerblue","aquamarine2","chocolate")
    colores_base=c("forestgreen","red","aquamarine2")
    fill_base=c("forestgreen",NA,"red",NA,"aquamarine2",NA)

    {

    metadata=meta[which(meta$short_setup %in% set),]
    for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

    reads=counts[,which(colnames(counts) %in% rownames(metadata))]

    y <- DGEList(counts = reads)
    y <- calcNormFactors(y)
    design <- model.matrix(~short_setup,data=metadata)
    # y <- estimateDisp(y,design)

    v=voom(y,design,plot=TRUE)


    exp=v$E

    pca <-prcomp(t(exp),scale=T)
    sum_pca=data.frame(summary(pca)$importance)

    datos <- data.frame(pca$x)
    colnames(datos)=paste0("PC",c(1:ncol(datos)))
    length(which(rownames(datos)!=rownames(metadata)))
    datos=cbind(metadata,datos)

    datos$Culture=factor(datos$Culture,levels=set_cult)
    datos$short_setup=factor(datos$short_setup,levels=set)
    }

    ###PCA plot PC1 vs PC2
    pl11d=ggplot(datos)+geom_point(aes(x=PC1,y=PC3,color=Culture,fill=short_setup),shape=21,size=3,stroke=1)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab(paste0("PC1: ",sum_pca[2,1]*100,"% variance ~growth phase"))+ylab(paste0("PC2: ",sum_pca[2,2]*100,"% variance ~growth phase"))

    pdf("C1_C5_C7.pdf",width=6,height=5)
    print(pl11d)
    dev.off()
    
}


## Only c5-c6
set_cult=c("C5","C6")
set=c("C5_Fe_YES","C5_Fe_NO","C6_Fe_YES","C6_Fe_NO")
# ("forestgreen","black","red","dodgerblue","aquamarine2","chocolate")
colores_base=c("red","dodgerblue")
fill_base=c("red",NA,"dodgerblue",NA)

{

metadata=meta[which(meta$short_setup %in% set),]
for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

reads=counts[,which(colnames(counts) %in% rownames(metadata))]

y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
# y <- estimateDisp(y,design)

v=voom(y,design,plot=TRUE)


exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance)

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos)

datos$Culture=factor(datos$Culture,levels=set_cult)
datos$short_setup=factor(datos$short_setup,levels=set)
}

###PCA plot PC1 vs PC2
pl12=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab(paste0("PC1: ",sum_pca[2,1]*100,"% variance ~growth phase"))+ylab(paste0("PC2: ",sum_pca[2,2]*100,"% variance ~growth phase"))


pdf("C5_C6.pdf",width=6,height=5)
print(pl12)
dev.off()


## Only c2-c6
set_cult=c("C2","C6")
set=c("C2_Fe_YES","C2_Fe_NO","C6_Fe_YES","C6_Fe_NO")
# ("forestgreen","black","red","dodgerblue","aquamarine2","chocolate")
colores_base=c("black","dodgerblue")
fill_base=c("black",NA,"dodgerblue",NA)

{

metadata=meta[which(meta$short_setup %in% set),]
for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

reads=counts[,which(colnames(counts) %in% rownames(metadata))]

y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
# y <- estimateDisp(y,design)

v=voom(y,design,plot=TRUE)


exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance)

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos)

datos$Culture=factor(datos$Culture,levels=set_cult)
datos$short_setup=factor(datos$short_setup,levels=set)
}

###PCA plot PC1 vs PC2
pl12b=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab(paste0("PC1: ",sum_pca[2,1]*100,"% variance ~growth phase"))+ylab(paste0("PC2: ",sum_pca[2,2]*100,"% variance ~growth phase"))


## Only c8-c6
set_cult=c("C6","C8")
set=c("C6_Fe_YES","C6_Fe_NO","C8_Fe_YES","C8_Fe_NO")
# ("forestgreen","black","red","dodgerblue","aquamarine2","chocolate")
colores_base=c("dodgerblue","chocolate")
fill_base=c("dodgerblue",NA,"chocolate",NA)

{

metadata=meta[which(meta$short_setup %in% set),]
for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

reads=counts[,which(colnames(counts) %in% rownames(metadata))]

y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
# y <- estimateDisp(y,design)

v=voom(y,design,plot=TRUE)


exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance)

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos)

datos$Culture=factor(datos$Culture,levels=set_cult)
datos$short_setup=factor(datos$short_setup,levels=set)
}

###PCA plot PC1 vs PC2
pl12c=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab(paste0("PC1: ",sum_pca[2,1]*100,"% variance ~growth phase"))+ylab(paste0("PC2: ",sum_pca[2,2]*100,"% variance ~growth phase"))

pdf("C6_C8.pdf",width=6,height=5)
print(pl12c)
dev.off()

## Only c7-c8
set_cult=c("C7","C8")
set=c("C7_Fe_YES","C7_Fe_NO","C8_Fe_YES","C8_Fe_NO")
# ("forestgreen","black","red","dodgerblue","aquamarine2","chocolate")
colores_base=c("aquamarine2","chocolate")
fill_base=c("aquamarine2",NA,"chocolate",NA)

{

metadata=meta[which(meta$short_setup %in% set),]
for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

reads=counts[,which(colnames(counts) %in% rownames(metadata))]

y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
# y <- estimateDisp(y,design)

v=voom(y,design,plot=TRUE)


exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance)

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos)

datos$Culture=factor(datos$Culture,levels=set_cult)
datos$short_setup=factor(datos$short_setup,levels=set)
}

###PCA plot PC1 vs PC2
pl13=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab(paste0("PC1: ",sum_pca[2,1]*100,"% variance ~growth phase"))+ylab(paste0("PC2: ",sum_pca[2,2]*100,"% variance ~growth phase"))
pdf("C7_C8.pdf",width=6,height=5)
print(pl13)
dev.off()



## Only c2-c8
set_cult=c("C2","C8")
set=c("C2_Fe_YES","C2_Fe_NO","C8_Fe_YES","C8_Fe_NO")
# ("forestgreen","black","red","dodgerblue","aquamarine2","chocolate")
colores_base=c("black","chocolate")
fill_base=c("black",NA,"chocolate",NA)

{

metadata=meta[which(meta$short_setup %in% set),]
for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

reads=counts[,which(colnames(counts) %in% rownames(metadata))]

y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
# y <- estimateDisp(y,design)

v=voom(y,design,plot=TRUE)


exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance)

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos)

datos$Culture=factor(datos$Culture,levels=set_cult)
datos$short_setup=factor(datos$short_setup,levels=set)
}

###PCA plot PC1 vs PC2
pl13b=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab(paste0("PC1: ",sum_pca[2,1]*100,"% variance ~growth phase"))+ylab(paste0("PC2: ",sum_pca[2,2]*100,"% variance ~growth phase"))
pdf("C2_C8.pdf",width=6,height=5)
print(pl13b)
dev.off()



## Only c5-c6-c7-c8
set_cult=c("C5","C6","C7","C8")
set=c("C5_Fe_YES","C5_Fe_NO","C6_Fe_YES","C6_Fe_NO","C7_Fe_YES","C7_Fe_NO","C8_Fe_YES","C8_Fe_NO")
# ("forestgreen","black","red","dodgerblue","aquamarine2","chocolate")
colores_base=c("red","dodgerblue","aquamarine2","chocolate")
fill_base=c("red",NA,"dodgerblue",NA,"aquamarine2",NA,"chocolate",NA)

{

metadata=meta[which(meta$short_setup %in% set),]
for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

reads=counts[,which(colnames(counts) %in% rownames(metadata))]

y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
# y <- estimateDisp(y,design)

v=voom(y,design,plot=TRUE)


exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance)

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos)

datos$Culture=factor(datos$Culture,levels=set_cult)
datos$short_setup=factor(datos$short_setup,levels=set)
}

###PCA plot PC1 vs PC2
pl14=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab(paste0("PC1: ",sum_pca[2,1]*100,"% variance ~growth phase"))+ylab(paste0("PC2: ",sum_pca[2,2]*100,"% variance ~growth phase"))

pdf("C5_C6_C7_C8.pdf",width=6,height=5)
print(pl14)
dev.off()


## Only c1-c2-c7-c8
set_cult=c("C1","C2","C7","C8")
set=c("C1_Fe_YES","C1_Fe_NO","C2_Fe_YES","C2_Fe_NO","C7_Fe_YES","C7_Fe_NO","C8_Fe_YES","C8_Fe_NO")
# ("forestgreen","black","red","dodgerblue","aquamarine2","chocolate")
colores_base=c("forestgreen","black","aquamarine2","chocolate")
fill_base=c("forestgreen",NA,"black",NA,"aquamarine2",NA,"chocolate",NA)

{

metadata=meta[which(meta$short_setup %in% set),]
for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

reads=counts[,which(colnames(counts) %in% rownames(metadata))]

y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
# y <- estimateDisp(y,design)

v=voom(y,design,plot=TRUE)


exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance)

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos)

datos$Culture=factor(datos$Culture,levels=set_cult)
datos$short_setup=factor(datos$short_setup,levels=set)
}

###PCA plot PC1 vs PC2
pl15=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab(paste0("PC1: ",sum_pca[2,1]*100,"% variance ~growth phase"))+ylab(paste0("PC2: ",sum_pca[2,2]*100,"% variance ~growth phase"))

pdf("C1_C2_C7_C8.pdf",width=6,height=5)
print(pl15)
dev.off()




metadata=meta[which(meta$short_setup %in% c("C6_Fe_YES","C6_Fe_NO")),]
for(i in 1:8){metadata[,i]=factor(metadata[,i],levels=unique(metadata[,i]))}

reads=counts[,which(colnames(counts) %in% rownames(metadata))]

y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
# y <- estimateDisp(y,design)

v=voom(y,design,plot=TRUE)


exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance)

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos)

datos$Culture=factor(datos$Culture,levels=c("C6"))
datos$short_setup=factor(datos$short_setup,levels=c("C6_Fe_YES","C6_Fe_NO"))
# ("forestgreen","black","red","dodgerblue","aquamarine2","chocolate")
colores_base=c("dodgerblue")
fill_base=c("dodgerblue",NA)

###PCA plot PC1 vs PC2
pl10=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab(paste0("PC1: ",sum_pca[2,1]*100,"% variance ~growth phase"))+ylab(paste0("PC2: ",sum_pca[2,2]*100,"% variance ~growth phase"))












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
