### Libraries
dcols=function(x){data.frame(colnames(x))}
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

## Select sampleset:
metadata=metadata[which(metadata$Culture %in% c("C1","C5","C6","C7")),]
reads=reads[,which(colnames(reads) %in% rownames(metadata))]
fpkms=fpkms[,which(colnames(fpkms) %in% rownames(metadata))]
length(which(colnames(reads)!=rownames(metadata)))

# Select geneset:
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
filtrado=filter_OK(reads,fpkms,feature_data,threshold=10,metadata=metadata,column="short_setup")

reads=filtrado[[1]]
fpkms=filtrado[[2]]
feature_data=filtrado[[3]]

## PCA.
y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
v=voom(y,design,plot=TRUE)
exp=v$E
pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance[,c(1:5)])

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos[,1:7])

colores_base=c("forestgreen","red","dodgerblue","aquamarine2")
fill_base=c("forestgreen",NA,"red",NA,"dodgerblue",NA,"aquamarine2",NA)


#colores_base=c("forestgreen","black","red","dodgerblue","aquamarine2","chocolate")
#fill_base=c("forestgreen",NA,"black",NA,"red",NA,"dodgerblue",NA,"aquamarine2",NA,"chocolate",NA)

###PCA plot PC1 vs PC2
pl1=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1.5)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+
xlab(paste0("PC1:",100*round(sum_pca$PC1[2],digits=3),"% variance explained"))+
ylab(paste0("PC2:",100*round(sum_pca$PC2[2],digits=3),"% variance explained"))

pdf("pca_C1567.pdf")
print(pl1)
dev.off()

## DE: limma.
{
dge <- DGEList(counts = reads)
dge <- calcNormFactors(dge)
design=model.matrix(~short_setup,data=metadata)
v=voom(dge,design,plot=TRUE)
fit <-lmFit(v,design)
fit <- eBayes(fit)

betas=data.frame(fit$coefficients)
ts=data.frame(fit$t)
ps=data.frame(fit$p.value)
BH=data.frame(ps)
ST=data.frame(ps)
errors <- sqrt(fit$s2.post) * fit$stdev.unscaled
dofs=data.frame(base=(fit$df.residual+fit$df.prior))

th=0.05
for(i in 1:ncol(BH))
{
    BH[,i]=p.adjust(ps[,i],method="BH")
    if(i==0){ST[,i]=NA}else{
    ST[,i]=qvalue(ps[,i])$qvalue}
    print(paste("i=",i,length(which(BH[,i]<th)),length(which(ST[,i]<th))))
}

DE_base=list(betas=betas,ts=ts,ps=ps,BH=BH,ST=ST,errors=errors,dofs=dofs)

## Get contrasts:

DE_object=DE_base
contrast=a
th=0.05

process_contrast=function(DE_object,contrast,th=0.05){

name <-deparse(substitute(contrast))
vec=rep(0,length(colnames(design)))
vec[abs(contrast)]=contrast/abs(contrast)
fit2 <- contrasts.fit(fit, vec)
fit2 <- eBayes(fit2)

b=data.frame(fit2$coefficients)
t=data.frame(fit2$t)
p=data.frame(fit2$p.value)
bh=data.frame(p)
st=data.frame(p)

e <- sqrt(fit2$s2.post) * fit2$stdev.unscaled
d=fit2$df.residual+fit2$df.prior

print(d[1])

for(i in 1:ncol(bh))
{
    bh[,i]=p.adjust(p[,i],method="BH")
    if(i==0){st[,i]=NA}else{
    st[,i]=qvalue(p[,i])$qvalue}
    print(paste("i=",i,name,length(which(bh[,i]<th)),length(which(st[,i]<th))))
}

betas_ext=cbind(DE_object$betas,new=b)
ts_ext=cbind(DE_object$ts,new=t)
ps_ext=cbind(DE_object$ps,new=p)
BH_ext=cbind(DE_object$BH,new=bh)
ST_ext=cbind(DE_object$ST,new=st)
errors_ext=cbind(DE_object$errors,new=e)
dofs_ext=cbind(DE_object$dofs,new=d)

colnames(betas_ext)[ncol(betas_ext)]=name
colnames(ts_ext)[ncol(ts_ext)]=name
colnames(ps_ext)[ncol(ps_ext)]=name
colnames(BH_ext)[ncol(BH_ext)]=name
colnames(ST_ext)[ncol(ST_ext)]=name
colnames(errors_ext)[ncol(errors_ext)]=name
colnames(dofs_ext)[ncol(dofs_ext)]=name

colnames(p)="P_value"
nbins=50
pl=ggplot(p)+geom_histogram(aes(x=P_value),bins=nbins,center = 1/(2*nbins))+ggtitle(name)

DE_output=list(betas=betas_ext,ts=ts_ext,ps=ps_ext,BH=BH_ext,ST=ST_ext,errors=errors_ext,dofs=dofs_ext,histogram=pl)

return(DE_output)
}


a=c(7,-8)
b=c(3,-4)
c=c(-2)
j=c(5,-6)
#lipids
d=c(8,-4)
f=c(3,-7)
# dextrose
e=c(2,-4)
g=c(-3)
#stat
h=c(6,-4)
i=c(5,-3)

## Iron x lipids
a_minus_b=c(7,-8,-3,4)
## iron x stat
j_minus_b=c(5,-6,-3,4)
## iron x dextrose
c_minus_b=c(-2,-3,4)

safe=DE_base

## Iron effects
DE_base=process_contrast(DE_object=DE_base,contrast=a)
pl_a=DE_base$histogram

DE_base=process_contrast(DE_object=DE_base,contrast=b)
pl_b=DE_base$histogram

DE_base=process_contrast(DE_object=DE_base,contrast=c)
pl_c=DE_base$histogram

DE_base=process_contrast(DE_object=DE_base,contrast=j)
pl_j=DE_base$histogram

## Lipids effects

DE_base=process_contrast(DE_object=DE_base,contrast=d)
pl_d=DE_base$histogram

DE_base=process_contrast(DE_object=DE_base,contrast=f)
pl_f=DE_base$histogram

## Dextrose effects

DE_base=process_contrast(DE_object=DE_base,contrast=e)
pl_e=DE_base$histogram

DE_base=process_contrast(DE_object=DE_base,contrast=g)
pl_g=DE_base$histogram

## Stat effects

DE_base=process_contrast(DE_object=DE_base,contrast=h)
pl_h=DE_base$histogram

DE_base=process_contrast(DE_object=DE_base,contrast=i)
pl_i=DE_base$histogram

# Iron x lipids

DE_base=process_contrast(DE_object=DE_base,contrast=a_minus_b)
pl_a_minus_b=DE_base$histogram

# # Iron x stat

DE_base=process_contrast(DE_object=DE_base,contrast=j_minus_b)
pl_j_minus_b=DE_base$histogram

# # Iron x dextrose

DE_base=process_contrast(DE_object=DE_base,contrast=c_minus_b)
pl_c_minus_b=DE_base$histogram
}

## DE: DESEQ2.

dds <- DESeqDataSetFromMatrix(countData = reads,colData = metadata,design= ~ short_setup)
dds <- DESeq(dds)
contrast=a_minus_b

process_contrast_DEseq=function(dds,contrast,th=0.05){
    
    name <-deparse(substitute(contrast))
    
    lista=list(add=resultsNames(dds)[contrast[which(contrast>0)]],remove=resultsNames(dds)[-contrast[which(contrast<0)]])
    result <- lfcShrink(dds, contrast=lista, type="ashr")

    result$BH=p.adjust(result$pvalue,method="BH")
    result$ST=qvalue(result$pvalue)$qvalue
    print(paste(name,length(which(result$BH<th)),length(which(result$ST<th))))
    nbins=50
    pl=ggplot(result)+geom_histogram(aes(x=pvalue),bins=nbins,center = 1/(2*nbins))+ggtitle(name)
    output=list(data=result,figure=pl)
    return(output)
    }

res_a=process_contrast_DEseq(dds,contrast=a)
res_b=process_contrast_DEseq(dds,contrast=b)
res_c=process_contrast_DEseq(dds,contrast=c)
res_j=process_contrast_DEseq(dds,contrast=j)

res_d=process_contrast_DEseq(dds,contrast=d)
res_f=process_contrast_DEseq(dds,contrast=f)

res_e=process_contrast_DEseq(dds,contrast=e)
res_g=process_contrast_DEseq(dds,contrast=g)

res_h=process_contrast_DEseq(dds,contrast=h)
res_i=process_contrast_DEseq(dds,contrast=i)

res_a_minus_b=process_contrast_DEseq(dds,contrast=a_minus_b)
res_j_minus_b=process_contrast_DEseq(dds,contrast=j_minus_b)
res_c_minus_b=process_contrast_DEseq(dds,contrast=c_minus_b)

## FALTA EDGER
