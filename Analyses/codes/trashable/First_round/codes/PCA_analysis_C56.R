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


### functions
dcols=function(x){data.frame(colnames(x))}

my_name <- function(v1) {
  deparse(substitute(v1))
}

histogram <- function(data,pvalue,title){
  frame <- data.frame(data)
  hist <- ggplot(frame)+geom_histogram(aes(x=pvalue),color="black",
                                       fill="lightblue",bins = (n_genes/(n_genes/50)),boundary=0)+
    scale_x_continuous(limits = c(0, 1)) +
    theme_minimal()+ ggtitle(title)+labs(x="PValues", y = "Counts")+
    theme(plot.title = element_text(hjust = 0.5))
  return(hist)
}

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

process_contrast_limma=function(DE_object,contrast,th=0.05,title){
  name <-deparse(substitute(contrast))
  vec=rep(0,length(colnames(design)))
  vec[abs(contrast)]=contrast/abs(contrast)
  fit2 <- contrasts.fit(fit, vec)
  fit2 <- eBayes(fit2)
  top=data.frame(topTable(fit2))
  b=data.frame(fit2$coefficients)
  t=data.frame(fit2$t)
  p=data.frame(fit2$p.value)
  bh=data.frame(p)
  st=data.frame(p)
  # top.table <- topTable(fit2, sort.by = "P", n = Inf)
  
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
  colnames(top)[ncol(top)]=name
  colnames(p)="P_value"
  nbins=50
  pl=histogram(p,p$P_value,title=title)
  DE_output=list(betas=betas_ext,ts=ts_ext,ps=ps_ext,BH=BH_ext,ST=ST_ext,
                 errors=errors_ext,dofs=dofs_ext,top=top,histogram=pl)
  
  return(DE_output)
}

process_contrast_DEseq=function(dds,contrast,th=0.05,title){
  name <-deparse(substitute(contrast))
  
  lista=list(add=resultsNames(dds)[contrast[which(contrast>0)]],remove=resultsNames(dds)[-contrast[which(contrast<0)]])
  result <- lfcShrink(dds, contrast=lista, type="ashr")
  
  result$BH=p.adjust(result$pvalue,method="BH")
  result$ST=qvalue(result$pvalue)$qvalue
  print(paste(name,length(which(result$BH<th)),length(which(result$ST<th))))
  nbins=50
  pl=histogram(result,result$pvalue,title)
  # pl=ggplot(result)+geom_histogram(aes(x=pvalue),bins=nbins,center = 1/(2*nbins))+ggtitle(name)
  output=list(data=result,figure=pl)
  return(output)
}


dir.create("histogramsC5_C6_JAC")
### Read the data
metadata = read.table("inputs/processed_all_JAC/metadata.txt")
feature_data = read.table("inputs/processed_all_JAC/feature_data.txt")
reads = read.table("inputs/processed_all_JAC/reads.txt")
fpkms = read.table("inputs/processed_all_JAC/fpkm.txt")


metadata$short_setup=paste0(metadata$Culture,"_Fe_",metadata$Iron)

metadata <- metadata[order(metadata$short_setup),]
metadata <- metadata[order(metadata$Iron),]

reads <- reads[rownames(metadata)]
fpkms <- fpkms[rownames(metadata)]
## Select sampleset:
metadata=metadata[which(metadata$Culture %in% c("C5","C6")),]
reads=reads[,which(colnames(reads) %in% rownames(metadata))]
fpkms=fpkms[,which(colnames(fpkms) %in% rownames(metadata))]
length(which(colnames(reads)!=rownames(metadata)))

# Select geneset:

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

colores_base=c("red","dodgerblue")
fill_base=c("red",NA,"dodgerblue",NA)


#colores_base=c("forestgreen","black","red","dodgerblue","aquamarine2","chocolate")
#fill_base=c("forestgreen",NA,"black",NA,"red",NA,"dodgerblue",NA,"aquamarine2",NA,"chocolate",NA)

###PCA plot PC1 vs PC2

pdf(file = "histogramsC5_C6_JAC/PCA_C5_C6.pdf",width=6,height=5)
pl1=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1.5)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+
xlab(paste0("PC1:",100*round(sum_pca$PC1[2],digits=3),"% variance explained"))+
ylab(paste0("PC2:",100*round(sum_pca$PC2[2],digits=3),"% variance explained"))
print(pl1)
dev.off()

## DE: limma.
design <- ~ short_setup
dge <- DGEList(counts = reads)
dge <- calcNormFactors(dge)
n_genes <- nrow(dge$counts)
design=model.matrix(design,data=metadata)
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
  # if(i==0){ST[,i]=NA}else{
  #   ST[,i]=qvalue(ps[,i])$qvalue}
  print(paste("i=",i,length(which(BH[,i]<th)),length(which(ST[,i]<th))))
}
DE_base=list(betas=betas,ts=ts,ps=ps,BH=BH,ST=ST,errors=errors,dofs=dofs)
  
  ## Get contrasts:
  
DE_object=DE_base
  # contrast=a
th=0.05
  
# Iron effects
a=c(-2)
b=c(3,-4)
  
    # Stat effects
c=c(4,-2)
d=c(3)
    
    ## Stat x iron
    
a_minus_b=c(-2,-3,4)

safe=DE_base
  
pdf("histogramsC5_C6_JAC/histograms_C5_C6_limma.pdf",width=6,height=5)
{
    ## Iron effects
DE_base=process_contrast_limma(DE_object=DE_base,contrast=a,title = "C5 vs -C5")
pl_a=DE_base$histogram
topC5vsC5 = DE_base$top
print(pl_a)
    
DE_base=process_contrast_limma(DE_object=DE_base,contrast=b,title = "C6 vs -C6")
pl_b=DE_base$histogram
topC6vsC6 = DE_base$top
print(pl_b)
    
    ## Stat effects
    
DE_base=process_contrast_limma(DE_object=DE_base,contrast=c,title = "C6 vs C5")
pl_c=DE_base$histogram
topC6vsC5 = DE_base$top
print(pl_c)
    
DE_base=process_contrast_limma(DE_object=DE_base,contrast=d,title = "-C6 vs -C5")
pl_d=DE_base$histogram
topC6vsC5_no = DE_base$top
print(pl_d)
    
    # Iron x Stat
    
DE_base=process_contrast_limma(DE_object=DE_base,contrast=a_minus_b,title = "C5_C5_Fe_No vs C6_C6_Fe_No")
pl_a_minus_b=DE_base$histogram
topC5_C5vsC6_C6 = DE_base$top
print(pl_a_minus_b)
}
  
dev.off()

top_genes <- cbind(topC5vsC5,topC6vsC5,topC6vsC5_no,topC6vsC6,topC5_C5vsC6_C6)

# top.table$Gene <- rownames(top.table)
# top.table <- top.table[,c("Gene", names(top.table)[1:6])]
# write.table(top.table, file = "table_C5_C6.txt", row.names = F, sep = "\t", quote = F)

## DE: DESEQ2.
{
dds <- DESeqDataSetFromMatrix(countData = reads,colData = metadata,design= ~ short_setup)
dds <- DESeq(dds)
  
dev.new()
pdf("histogramsC5_C6_JAC/histograms_C5_C6_DESeq2.pdf",width=6,height=5)
{
res_a=process_contrast_DEseq(dds,contrast=a,title = "C5 vs -C5")
pl_a <- res_a$figure
print(pl_a)
res_b=process_contrast_DEseq(dds,contrast=b,title = "C6 vs -C6")
pl_b <- res_b$figure
print(pl_b)
res_c=process_contrast_DEseq(dds,contrast=c,title = "C5 vs C6")
pl_c <- res_c$figure
print(pl_c)
res_j=process_contrast_DEseq(dds,contrast=d,title = "-C5 vs -C6_Fe_NO")
pl_d <- res_d$figure
print(pl_d)
    
res_a_minus_b=process_contrast_DEseq(dds,contrast=a_minus_b,title = "C5_C5_Fe_No vs C6_C6_Fe_No")
pl_a_minus_b<- res_a_minus_b$figure
print(pl_a_minus_b)
  
}
dev.off()
}

## EDGER
{
group <- factor(metadata$short_setup)
  
  ### 
{
dge <- DGEList(counts = reads,group = group)
min_count <- 20
keep <- filterByExpr(dge, min.count=min_count)
table(keep)
dge <- dge[keep,,keep.lib.sizes = FALSE]
n_genes <- nrow(dge$counts)

###Here, we retain only those genes that are represented at least 1cpm reads in 
#at least two samples (cpm=counts per million).
countsPerMillion <- cpm(dge)
summary(countsPerMillion)
countCheck <- countsPerMillion > 1
head(countCheck,5)
keep <- which(rowSums(countCheck) >= 2)
dge <- dge[keep,]
summary(cpm(dge)) 

dge <- calcNormFactors(dge)
}
  
  ### design matrix 
design=model.matrix(~0+group,data=metadata)
colnames(design) <- levels(group)
  ###Cox-Reid profile-adjusted likelihood (CR) method in estimating dispersions
  
dge<- estimateDisp(dge, design,robust=TRUE)
dge$common.dispersion
# dge <- estimateCommonDisp(dge)
# dge <- estimateTagwiseDisp(dge)
n_genes <- nrow(dge$counts)
  ###DE gene
  
fit <- glmQLFit(dge,design,robust = TRUE)
{
  ###                         Contrasts
  #=========================== Effect EXP - STAT =====================================
C6vsC5_Fe_NO <- makeContrasts(C6_Fe_NO-C5_Fe_NO, levels=design)
C6vsC5 <- makeContrasts(C6_Fe_YES-C5_Fe_YES, levels=design)
  #=========================== Effect of Fe =====================================
C5vsC5 <- makeContrasts(C5_Fe_NO-C5_Fe_YES, levels=design)
C6vsC6 <- makeContrasts(C6_Fe_NO-C6_Fe_YES, levels=design)
  #========================== IRON X STAT========================================
C6_C6_NOvsC5_C5_NO<- makeContrasts(C6vsC6-C5vsC5, levels=design)
}
  ########################## DE for each comparison ############################
path="DE_genes"
dir.create(path = path)


{
qlf.C6vsC5_Fe_NO <- glmQLFTest(fit, contrast=C6vsC5_Fe_NO)
tag_C6vsC5_Fe_NO <- topTags(qlf.C6vsC5_Fe_NO,n=(n_genes+1),adjust.method = "BH",sort.by = "PValue")
write.table(tag_C6vsC5_Fe_NO,file = paste0(path,"/",my_name(C6vsC5_Fe_NO),".txt"))

qlf.C6vsC5 <- glmQLFTest(fit, contrast=C6vsC5)
tag_C6vsC5 <- topTags(qlf.C6vsC5,n=(n_genes+1),adjust.method = "BH",sort.by = "PValue")
write.table(tag_C6vsC5,file = paste0(path,"/",my_name(C6vsC5),".txt"))
  
qlf.C5vsC5  <- glmQLFTest(fit, contrast=C5vsC5 )
tag_C5vsC5<- topTags(qlf.C5vsC5,n=(n_genes+1),adjust.method = "BH",sort.by = "PValue")
write.table(tag_C5vsC5,file = paste0(path,"/",my_name(C5vsC5),".txt"))
  
qlf.C6vsC6 <- glmQLFTest(fit, contrast=C6vsC6)
tag_C6vsC6 <- topTags(qlf.C6vsC6,n=(n_genes+1),adjust.method = "BH",sort.by = "PValue")
write.table(tag_C6vsC6,file = paste0(path,"/",my_name(C6vsC6),".txt"))

qlf.C6_C6_NOvsC5_C5_NO <- glmQLFTest(fit, contrast=C6_C6_NOvsC5_C5_NO)
tag_C6_C6_NOvsC5_C5_NO <- topTags(qlf.C6_C6_NOvsC5_C5_NO,n=(n_genes+1),adjust.method = "BH",
                                    sort.by = "PValue")
write.table(tag_C6_C6_NOvsC5_C5_NO,file = paste0(path,"/",my_name(C6_C6_NOvsC5_C5_NO),".txt"))

}



sum(tag_C6vsC5_Fe_NO$table$FDR < th)
sum(tag_C6vsC5$table$FDR < th)
sum(tag_C5vsC5$table$FDR < th)
sum(tag_C6vsC6$table$FDR < th)
sum(tag_C6_C6_NOvsC5_C5_NO$table$FDR < th)

### declare Pvalues 
p_C6vsC5_Fe_no <- qlf.C6vsC5_Fe_NO$table$PValue
p_C6vsC5 <- qlf.C6vsC5$table$PValue
p_C5vsC5 <- qlf.C5vsC5$table$PValue
p_C6vsC6 <- qlf.C6vsC6$table$PValue
p_C6_C6_novsC5_C5_no <- qlf.C6_C6_NOvsC5_C5_NO$table$PValue
  
  # Print the plot to a pdf file

  
  
pdf("histogramsC5_C6_JAC/histograms_C5_C6.pdf",width=6,height=5)
  
  ### Histogram Pvalue vs Count
{
hist_1 <- histogram(qlf.C6vsC5 ,p_C6vsC5,"C6vsC5")
print(hist_1)
hist_2 <- histogram(qlf.C6vsC5_Fe_NO ,p_C6vsC5_Fe_no,"-C6vs-C5_Fe_NO")
print(hist_2)
hist_3 <- histogram(qlf.C5vsC5,p_C5vsC5,"C5vs-C5")
print(hist_3)
hist_4 <- histogram(qlf.C6vsC6,p_C6vsC6,"C6vs-C6")
print(hist_4)
hist_4<- histogram(qlf.C6_C6_NOvsC5_C5_NO ,p_C6_C6_novsC5_C5_no,"C6_C6_NOvsC5_C5_NO")
print(hist_4)
}

dev.off()
}
