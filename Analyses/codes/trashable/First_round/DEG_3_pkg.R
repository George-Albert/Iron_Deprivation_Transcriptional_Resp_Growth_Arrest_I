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

### Read the data
meta_path = "inputs/processed_all_JAC/metadata.txt"
metadata = read.table(meta_path)

feat_path = "inputs/processed_all_JAC/feature_data.txt"
feature_data = read.table(feat_path)

read_path = "inputs/processed_all_JAC/reads.txt"
reads = read.table(read_path)

fpkms_path = "inputs/processed_all_JAC/fpkm.txt" 
fpkms = read.table(fpkms_path)

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

### number of samples and genes
n_samples <- ncol(reads)
n_genes <- nrow (reads)
n_Fe_NO <- nrow(metadata[which(metadata$Iron=="NO"),])
n_Fe_Yes <- n_samples-n_Fe_NO
### Threshold
th <- 0.05

group <- factor(metadata$short_setup)

### design matrix 
design=model.matrix(~0+group,data=metadata)
colnames(design) <- levels(group)

{
  #================================ Contrasts ==================================
  
  #=========================== Effect EXP - STAT ===============================
  C6vsC5_Fe_NO <- makeContrasts(C6_Fe_NO-C5_Fe_NO, levels=design)
  C6vsC5 <- makeContrasts(C6_Fe_YES-C5_Fe_YES, levels=design)
  #=========================== Effect of Fe ====================================
  C5vsC5 <- makeContrasts(C5_Fe_NO-C5_Fe_YES, levels=design)
  C6vsC6 <- makeContrasts(C6_Fe_NO-C6_Fe_YES, levels=design)
  #========================== IRON X STAT=======================================
  C6_C6_NOvsC5_C5_NO<- makeContrasts(C6vsC6-C5vsC5, levels=design)
}
###============================ EDGER===========================================

{
    dge <- DGEList(counts = reads,group = group)
    min_count <- 10
    keep <- filterByExpr(dge, min.count=min_count)
    table(keep)
    dge <- dge[keep,,keep.lib.sizes = FALSE]
   
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
  
###Cox-Reid profile-adjusted likelihood (CR) method in estimating dispersions
  
dge<- estimateDisp(dge, design,robust=TRUE)
dge$common.dispersion
  # dge <- estimateCommonDisp(dge)
  # dge <- estimateTagwiseDisp(dge)
n_genes <- nrow(dge$counts)
  ###DE gene
  
fit <- glmQLFit(dge,design,robust = TRUE)

######################### DE for each comparison ############################

path="DE_genes"
dir.create(path = path)
  
  
{
    qlf.C6vsC5_Fe_NO <- glmQLFTest(fit, contrast=C6vsC5_Fe_NO)
    tag_C6vsC5_Fe_NO <- topTags(qlf.C6vsC5_Fe_NO,n=(n_genes+1),adjust.method = "BH")
    write.table(tag_C6vsC5_Fe_NO,file = paste0(path,"/",my_name(C6vsC5_Fe_NO),".txt"))
    
    qlf.C6vsC5 <- glmQLFTest(fit, contrast=C6vsC5)
    tag_C6vsC5 <- topTags(qlf.C6vsC5,n=(n_genes+1),adjust.method = "BH")
    write.table(tag_C6vsC5,file = paste0(path,"/",my_name(C6vsC5),".txt"))
    
    qlf.C5vsC5  <- glmQLFTest(fit, contrast=C5vsC5 )
    tag_C5vsC5<- topTags(qlf.C5vsC5,n=(n_genes+1),adjust.method = "BH")
    write.table(tag_C5vsC5,file = paste0(path,"/",my_name(C5vsC5),".txt"))
    
    qlf.C6vsC6 <- glmQLFTest(fit, contrast=C6vsC6)
    tag_C6vsC6 <- topTags(qlf.C6vsC6,n=(n_genes+1),adjust.method = "BH")
    write.table(tag_C6vsC6,file = paste0(path,"/",my_name(C6vsC6),".txt"))
    
    qlf.C6_C6_NOvsC5_C5_NO <- glmQLFTest(fit, contrast=C6_C6_NOvsC5_C5_NO)
    tag_C6_C6_NOvsC5_C5_NO <- topTags(qlf.C6_C6_NOvsC5_C5_NO,n=(n_genes+1),adjust.method = "BH")
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

### ==========================LIMMA============================================
  dge1 <- DGEList(counts = reads)
  
  min_count <- 10
  keep <- filterByExpr(dge1, min.count=min_count)
  table(keep)
  dge1 <- dge1[keep,,keep.lib.sizes = FALSE]
  
  countsPerMillion <- cpm(dge1)
  summary(countsPerMillion)
  countCheck <- countsPerMillion > 1
  head(countCheck,5)
  keep <- which(rowSums(countCheck) >= 2)
  dge1 <- dge1[keep,]
  summary(cpm(dge1)) 
  
  n_genes <- nrow(dge1$counts)
  
  
  
  dge1 <- calcNormFactors(dge1)
  design <- model.matrix(~ 0 + group,data=metadata)
  colnames(design) <- levels(group)
  contrast.mat <- makeContrasts(C6vsC5_Fe_NO, C6vsC5, C6vsC6, C5vsC5, C6_C6_NOvsC5_C5_NO, levels = design)
  
  v=voom(dge1,design,plot=TRUE)
  
  fit <- lmFit(v, design = design)
  fit2 <- contrasts.fit(fit, contrast.mat)
  fit2 <- eBayes(fit2)
  
  deg_C6vsC5_NO <- topTable(fit2,coef = "C6vsC5_Fe_NO",number = n_genes,adjust.method = "BH", p.value = 0.05)
  write.table(deg_C6vsC5_NO,file = paste0(path,"/",my_name(C6vsC5_NO),"_limma.txt"))
  
  deg_C6vsC5 <- topTable(fit2,coef = "C6vsC5",number = n_genes,adjust.method = "BH",p.value = 0.05)
  write.table(deg_C6vsC5,file = paste0(path,"/",my_name(C6vsC5),"_limma.txt"))
  
  deg_C6vsC6 <- topTable(fit2,coef = "C6vsC6",number = n_genes,adjust.method = "BH",p.value = 0.05)
  write.table(deg_C6vsC6,file = paste0(path,"/",my_name(C6vsC6),"_limma.txt"))
  
  deg_C5vsC5 <- topTable(fit2,coef = "C5vsC5",number = n_genes,adjust.method = "BH",p.value = 0.05)
  write.table(deg_C5vsC5,file = paste0(path,"/",my_name(C5vsC5),"_limma.txt"))
  
  deg_C6_C6_NOvsC5_C5_NO <- topTable(fit2,coef = "C6_C6_NOvsC5_C5_NO",number = n_genes,adjust.method = "BH",p.value = 0.05)
  write.table(deg_C6_C6_NOvsC5_C5_NO,file = paste0(path,"/",my_name(C6_C6_NOvsC5_C5_NO),"limma.txt"))
  
  
  dim(deg_C6vsC5_NO)[1]
  dim(deg_C6vsC5)[1]
  dim(deg_C6vsC6)[1]
  dim(deg_C5vsC5)[1]
  dim(deg_C6_C6_NOvsC5_C5_NO)[1]
  
  pdf("histogramsC5_C6_JAC/histograms_C5_C6_limma_1.pdf",width=6,height=5)
  
  ### Histogram Pvalue vs Count
  {
    hist_1 <- histogram(fit2 ,fit2$p.value[,2],"C6vsC5")
    print(hist_1)
    hist_2 <- histogram(fit2 ,fit2$p.value[,1],"-C6vs-C5_Fe_NO")
    print(hist_2)
    hist_3 <- histogram(fit2,fit2$p.value[,4],"C5vs-C5")
    print(hist_3)
    hist_4 <- histogram(fit2,fit2$p.value[,3],"C6vs-C6")
    print(hist_4)
    hist_4<- histogram(fit2 ,fit2$p.value[,5],"C6_C6_NOvsC5_C5_NO")
    print(hist_4)
  }
  
  dev.off()

  
#=============================DESeq2=====================================
  
  
  design = ~ 0 + short_setup +
  colnames(design) <- levels(group)
  dds <- DESeqDataSetFromMatrix(countData = reads,colData = metadata,design= design)
  dds <- DESeq(dds,test = "Wald")
  resultsNames(dds)  
  
  
  
  
  
  
  
  
  
  
  