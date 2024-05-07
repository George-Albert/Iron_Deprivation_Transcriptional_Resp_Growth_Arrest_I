
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
{
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
}
### number of samples and genes
n_samples <- ncol(reads)
n_genes <- nrow (reads)
n_Fe_NO <- nrow(metadata[which(metadata$Iron=="NO"),])
n_Fe_Yes <- n_samples-n_Fe_NO

# ### remove last part of the ID culture
# samplenames <- substring(colnames(reads[,1:n_Fe_NO]),1,8)
# samplenames[(n_Fe_NO+1):n_samples]<- substring(colnames(reads[,(n_Fe_NO+1):n_Fe_Yes]),1,9)
# colnames(reads) <- samplenames 

group <- factor(metadata$short_setup)

### 
dge <- DGEList(counts = reads,group = group)
min_count <- 20
keep <- filterByExpr(dge, min.count=min_count)
table(keep)
dge <- dge[keep,,keep.lib.sizes = FALSE]
nrow(dge$counts)
dge <- calcNormFactors(dge)

### Examine the samples for outliers MDS plot
### Data exploration
plot_MDS <- function(dge,n_Fe_NO,n_Fe_Yes,bottom=5, left=4, top=4, right=2,n_samples){
  
  par(mar=c(bottom, left, top, right), xpd=TRUE)
  col <- hcl.colors(n_samples/2)
  col <- rep(col,2)
  plotMDS(dge, pch= c(rep(1,n_Fe_NO),rep(19,n_Fe_Yes)),col=col[dge$samples$group],
          cex=2)
  legend <- c(dge$samples$group[c(TRUE,FALSE)])
  legend("topright",inset=c(-0.8,0), legend=legend,pch=c(rep(1,3),rep(19,3)),
       col=col,cex=1.1)
}

plot_MDS(dge,n_Fe_NO,n_Fe_Yes,bottom=5, left=4, top=4, right=10,n_samples=n_samples)
### design matrix 
design=model.matrix(~0+group,data=metadata)
colnames(design) <- levels(group)
###Cox-Reid profile-adjusted likelihood (CR) method in estimating dispersions

dge<- estimateDisp(dge, design,robust=TRUE)
dge$common.dispersion
### Plot the dispersion

plotBCV(dge)

###DE gene

fit <- glmQLFit(dge,design,robust = TRUE)
### Visualize QL dispersions
plotQLDisp(fit)

###                         Contrasts
#============================ With Fe =====================================
### Contrast Effect of Dextrose 
C5vsC1 <- makeContrasts(C5_Fe_YES-C1_Fe_YES, levels=design)
### Contrast Effect of LCFA
C5vsC7 <- makeContrasts(C5_Fe_YES-C7_Fe_YES, levels=design)
#============================ Without Fe =====================================
### Contrast Effect of Dextrose 
C5vsC1_Fe_NO <- makeContrasts(C5_Fe_NO-C1_Fe_NO, levels=design)
### Contrast Effect of LCFA
C5vsC7_Fe_NO <- makeContrasts(C5_Fe_NO-C7_Fe_NO, levels=design) 
#=========================== Effect of Fe =====================================
C1vsC1 <- makeContrasts(C1_Fe_NO-C1_Fe_YES, levels=design)
C5vsC5 <- makeContrasts(C5_Fe_NO-C5_Fe_YES, levels=design)
C7vsC7 <- makeContrasts(C7_Fe_NO-C7_Fe_YES, levels=design)

########################## DE for each comparison ############################
qlf.C5vsC1 <- glmQLFTest(fit, contrast=C5vsC1)
tag_C5vsC1 <- topTags(qlf.C5vsC1,adjust.method = "BH",sort.by = "PValue")

qlf.C5vsC7 <- glmQLFTest(fit, contrast=C5vsC7)
tag_C5vsC7 <-topTags(qlf.C5vsC7,adjust.method = "BH",sort.by = "PValue")

qlf.C5vsC1_Fe_NO <- glmQLFTest(fit, contrast=C5vsC1_Fe_NO)
tag_C5vsC1_Fe_NO <- topTags(qlf.C5vsC1_Fe_NO,adjust.method = "BH",sort.by = "PValue")

qlf.C5vsC7_Fe_NO <- glmQLFTest(fit, contrast=C5vsC7_Fe_NO)
tag_C5vsC7_Fe_NO <- topTags(qlf.C5vsC7_Fe_NO,adjust.method = "BH",sort.by = "PValue")

qlf.C1vsC1 <- glmQLFTest(fit, contrast=C1vsC1 )
tag_C1vsC1<- topTags(qlf.C1vsC1,adjust.method = "BH",sort.by = "PValue")

qlf.C5vsC5  <- glmQLFTest(fit, contrast=C5vsC5 )
tag_C5vsC5<- topTags(qlf.C5vsC5,adjust.method = "BH",sort.by = "PValue")

qlf.C7vsC7 <- glmQLFTest(fit, contrast=C7vsC7)
tag_C7vsC7<- topTags(qlf.C7vsC7,adjust.method = "BH",sort.by = "PValue")

# p <- qlf.C1vsC1$table$PValue
# FDR <- p.adjust(p, method="BH")
# sum(FDR < 0.05)

summary_C5vsC1 <- summary(decideTests(qlf.C5vsC1,adjust.method = "BH"))
summary_C5vsC7 <- summary(decideTests(qlf.C5vsC7,adjust.method = "BH"))
summary_C5vsC1_Fe_NO <- summary(decideTests(qlf.C5vsC1_Fe_NO,adjust.method = "BH"))
summary_C5vsC7_Fe_NO <- summary(decideTests(qlf.C5vsC7_Fe_NO,adjust.method = "BH"))
summary_C1vsC1 <- summary(decideTests(qlf.C1vsC1,adjust.method = "BH"))
summary_C5vsC5 <- summary(decideTests(qlf.C5vsC5,adjust.method = "BH"))
summary_C7vsC7 <- summary(decideTests(qlf.C7vsC7,adjust.method = "BH"))
        
summary_test <- cbind(summary_C5vsC1,summary_C5vsC7,summary_C5vsC1_Fe_NO,
                      summary_C5vsC7_Fe_NO,summary_C1vsC1,summary_C5vsC5,
                      summary_C7vsC7)

### Plot MD
plotMD(qlf.C5vsC1)
plotMD(qlf.C5vsC7)
plotMD(qlf.C5vsC1_Fe_NO)
plotMD(qlf.C5vsC7_Fe_NO)
plotMD(qlf.C1vsC1)
plotMD(qlf.C5vsC5)
plotMD(qlf.C7vsC7)

histogram <- function(data,pvalue,title){
  frame <- data.frame(data)
  hist <- ggplot(frame)+geom_histogram(aes(x=pvalue),color="black", 
                                       fill="lightblue",bins=50)+
    theme_minimal()+ ggtitle(title)+theme(plot.title = element_text(hjust = 0.5))
  return(hist)
}

### declare Pvalues 
p_C5vsC1 <- qlf.C5vsC1$table$PValue
p_C5vsC7 <- qlf.C5vsC7$table$PValue
p_C5vsC1_no <- qlf.C5vsC1_Fe_NO$table$PValue
p_C5vsC7_no <- qlf.C5vsC7_Fe_NO$table$PValue
p_C1vsC1 <- qlf.C1vsC1$table$PValue
p_C5vsC5 <- qlf.C5vsC5$table$PValue
p_C7vsC7 <- qlf.C7vsC7$table$PValue

# Print the plot to a pdf file

dir.create("histograms_JAC")
pdf("histograms_JAC/histograms_all.pdf",width=6,height=5)

### Histogram Pvalue vs Count
hist_1 <- histogram(qlf.C5vsC1,p_C5vsC1,"C5vsC1")
print(hist_1)
hist_2 <- histogram(qlf.C5vsC7,p_C5vsC7,"C5vsC7")
print(hist_2)
hist_3 <- histogram(qlf.C5vsC1_Fe_NO,p_C5vsC1_no,"-C5vs-C1_Fe_NO")
print(hist_3)
hist_4 <- histogram(qlf.C5vsC7_Fe_NO,p_C5vsC7_no,"-C5vs-C7_Fe_NO")
print(hist_4)
hist_5 <- histogram(qlf.C1vsC1,p_C1vsC1,"C1vs-C1")
print(hist_5)
hist_6 <- histogram(qlf.C5vsC5,p_C5vsC5,"C5vs-C5")
print(hist_6)
hist_7 <- histogram(qlf.C7vsC7,p_C7vsC7,"C7vs-C7")
print(hist_7)

dev.off()

### Volcano Plots
volcan_plot <- function(data,x,y,ncol_sum,th_x,th_p,xintercept,title){
  
  datos <- data.frame(data)
  datos$Association <- "Not Sig."
  datos$Association[datos$logFC > th_x & datos$PValue < th_p] <- "UP"
  datos$Association[datos$logFC < (-th_x) & datos$PValue < th_p] <- "DOWN"
  
  mycolors <- c("blue", "black", "red")
  
ggplot(datos,aes(x=x,y=-log10(y),color=Association))+geom_point()+
    ylab("-log10(PValue)")+ theme_minimal()+geom_vline(xintercept=c(-xintercept, xintercept), 
                                                          col="black")+
  geom_hline(yintercept=-log10(0.05), col="black")+
  scale_color_manual(values=mycolors)+ggtitle(title)+
  theme(plot.title = element_text(hjust = 0.5))
}


dir.create("volcano_JAC")

pdf("histograms_JAC/volcano_all.pdf",width=6,height=5)

volcan_plot(qlf.C5vsC1,qlf.C5vsC7$table$logFC,qlf.C5vsC7$table$PValue,
            th_x = 0.5,th_p=0.05,xintercept = 0.5,title ="C5vsC1")
volcan_plot(qlf.C5vsC7,qlf.C5vsC7$table$logFC,qlf.C5vsC7$table$PValue,
            th_x = 0.5,th_p=0.05,xintercept = 0.5,title = "C5vsC7")

volcan_plot(qlf.C5vsC1_Fe_NO,qlf.C5vsC7$table$logFC,qlf.C5vsC7$table$PValue,
            th_x = 0.5,th_p=0.05,xintercept = 0.5,title ="C5vsC1_Fe_No")
volcan_plot(qlf.C5vsC7_Fe_NO,qlf.C5vsC7$table$logFC,qlf.C5vsC7$table$PValue,
            th_x = 0.5,th_p=0.05,xintercept = 0.5,title ="C5vsC7_Fe_No")

volcan_plot(qlf.C1vsC1,qlf.C5vsC7$table$logFC,qlf.C5vsC7$table$PValue,
            th_x = 0.5,th_p=0.05,xintercept = 0.5,title ="C1vsC1")
volcan_plot(qlf.C5vsC5,qlf.C5vsC7$table$logFC,qlf.C5vsC7$table$PValue,
            th_x = 0.5,th_p=0.05,xintercept = 0.5,title ="C5vsC5")
volcan_plot(qlf.C7vsC7,qlf.C5vsC7$table$logFC,qlf.C5vsC7$table$PValue,
            th_x = 0.5,th_p=0.05,xintercept = 0.5,title ="C7vsC7")
dev.off()
# abline(h=c(-1, 1), col="blue")

# v=voom(dge,design,plot=TRUE)
# fit <-lmFit(v,design)
# fit <- eBayes(fit)

# p <- qlf.C5vsC7$table$PValue
# hist(qvalue(p))
