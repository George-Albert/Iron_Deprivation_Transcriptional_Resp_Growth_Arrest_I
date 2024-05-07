#############################
### 0. Load dependencies. ###
#############################

{
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
  library(readxl)
  library(writexl)
  library(readxl) 
  library(xlsx)
  library(ggridges)
  library(ggeasy)
}


############################
### 1. Declare functions ###
############################
dcols=function(x){data.frame(colnames(x))}
my_name <- function(v1) {
  deparse(substitute(v1))
}


#################
### 2 Do PCA. ###
#################

{
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
  datos <- rename(datos,samples=short_setup)
}
  
  dir.create("Outputs_def/plot_tables/")
  write.table(datos,file="Outputs_def/plot_tables/PCA_table.txt")
  
  colores_base=c("red","dodgerblue")
  fill_base=c("red",NA,"dodgerblue",NA)
  shapes=c(21,1,21,1)
  pl1=ggplot(datos)+
    geom_point(aes(x=PC1,y=PC2,color=Culture,fill=samples,shape=samples),size=1.2,stroke=1.2)+
    scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+
    xlab(paste0("PC1:",100*round(sum_pca$PC1[2],digits=3),"% variance explained"))+
    ylab(paste0("PC2:",100*round(sum_pca$PC2[2],digits=3),"% variance explained"))+theme_classic()+
    scale_shape_manual(values=shapes)+
    theme(
      axis.text.y   = element_text(size=14),
      axis.text.x   = element_text(size=14),
      axis.title.y  = element_text(size=14),
      axis.title.x  = element_text(size=14),
      #panel.background = element_blank(),
      #panel.grid.major = element_blank(),
      #panel.grid.minor = element_blank(),
      #axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill=NA, size=1,linetype="solid"),
      #legend.title=element_blank(),
      #legend.position=c(0.8,0.8),
      #legend.text=element_text(size=14),
      legend.key.size = unit(1, 'lines'))
  
  dir.create("Outputs_def/Figures/",recursive=TRUE)
  pdf(file = "Outputs_def/Figures/PCA_C5_C6.pdf",width=6,height=5)
  print(pl1)
  dev.off()

