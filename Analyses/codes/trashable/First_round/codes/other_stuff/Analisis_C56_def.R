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
library(readxl)

### functions
dcols=function(x){data.frame(colnames(x))}
histogram <- function(data,pvalue,title){
  frame <- data.frame(data)
  n_genes=nrow(frame)
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
process_contrast_DEseq=function(dds,contrast,th=0.05,title){
  name <-deparse(substitute(contrast))
  
  lista=list(add=resultsNames(dds)[contrast[which(contrast>0)]],remove=resultsNames(dds)[-contrast[which(contrast<0)]])
  result <- lfcShrink(dds, contrast=lista, type="ashr")
  
  result$BH=p.adjust(result$pvalue,method="BH")
  print(paste(name,length(which(result$BH<th))))
  nbins=50
  pl=histogram(result,result$pvalue,title)
  # pl=ggplot(result)+geom_histogram(aes(x=pvalue),bins=nbins,center = 1/(2*nbins))+ggtitle(name)
  output=list(data=result,figure=pl)
  return(output)
}

dir.create("Outputs_def")
### Read the data. Jq: aqui había codigo redundante...
metadata = read.table("inputs/processed_all_JAC/metadata.txt")
feature_data = read.table("inputs/processed_all_JAC/feature_data.txt")
reads = read.table("inputs/processed_all_JAC/reads.txt")
fpkms = read.table("inputs/processed_all_JAC/fpkm.txt")
rownames(feature_data)=feature_data$Gene_ID

metadata=metadata[which(metadata$Culture %in% c("C5","C6")),]
metadata$short_setup=paste0(metadata$Culture,"_Fe_",metadata$Iron)
metadata <- metadata[order(metadata$short_setup),]
reads <- reads[,rownames(metadata)]
fpkms <- fpkms[,rownames(metadata)]

## Check coherence:
length(which(colnames(reads)!=rownames(metadata)))

# Select geneset:

filtrado=filter_OK(reads,fpkms,feature_data,threshold=2,metadata=metadata,column="short_setup")

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

pl1=ggplot(datos)+
geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1.5)+
scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+
xlab(paste0("PC1:",100*round(sum_pca$PC1[2],digits=3),"% variance explained"))+
ylab(paste0("PC2:",100*round(sum_pca$PC2[2],digits=3),"% variance explained"))

pdf(file = "Outputs_def/PCA_C5_C6.pdf",width=6,height=5)
print(pl1)
dev.off()

## Jq: Add PCs within culture.
## get contrasts
resultsNames(dds)
# Iron (Jq:) deprivation effects @ exp, and @ stat
a=c(-2)
b=c(3,-4)
# Stat effects w Fe and wout Fe
c=c(4,-2)
d=c(3)
## Stat x iron deprivation interaction
a_minus_b=c(-2,-3,4)

## DE: DESEQ2.
{
dds <- DESeqDataSetFromMatrix(countData = reads,colData = metadata,design= ~ short_setup)
dds <- DESeq(dds)

dev.new()
pdf("Outputs_def/histograms_C5_C6_DESeq2.pdf",width=6,height=5)
{
res_a=process_contrast_DEseq(dds,th=0.01,contrast=a,title = "Iron deprivation effect @ C5")
pl_a <- res_a$figure
print(pl_a)
res_b=process_contrast_DEseq(dds,th=0.01,contrast=b,title = "Iron deprivation effect @ C6")
pl_b <- res_b$figure
print(pl_b)
res_c=process_contrast_DEseq(dds,th=0.01,contrast=c,title = "Transition to stationary effects; with iron")
pl_c <- res_c$figure
print(pl_c)
res_d=process_contrast_DEseq(dds,th=0.01,contrast=d,title = "Transition to stationary effects; without iron")
pl_d <- res_d$figure
print(pl_d)
    
res_a_minus_b=process_contrast_DEseq(dds,contrast=a_minus_b,title = "Phase x iron interaction")
pl_a_minus_b<- res_a_minus_b$figure
print(pl_a_minus_b)
  
}
dev.off()
}


dir.create("Outputs_def/Stats_tables")
write.table(data.frame(res_a$data),"Outputs_def/Stats_tables/Iron_deprivation_effect_at_Exp.txt")
write.table(data.frame(res_b$data),"Outputs_def/Stats_tables/Iron_deprivation_effect_at_Stat.txt")
write.table(data.frame(res_c$data),"Outputs_def/Stats_tables/Stat_effect_with_Fe.txt")
write.table(data.frame(res_d$data),"Outputs_def/Stats_tables/Stat_effect_without_Fe.txt")
write.table(data.frame(res_a_minus_b$data),"Outputs_def/Stats_tables/Stat_Fedeprivation_interaction.txt")


how_much=function(tab,name){
    for(th in c(0.05,0.01))
    {
        for(th_size in c(0,0.2,0.5,1))
        {
            dir.create(paste0("Outputs_def/target_sets/th_",th,"_th_size_",th_size),recursive=TRUE)
            set=(which(abs(tab$log2FoldChange)>th_size & tab$BH<th))
            print(paste("th=",th,"th_size=",th_size,":",length(set)))
            target=rownames(tab)[set]
            sink(paste0("Outputs_def/target_sets/th_",th,"_th_size_",th_size,"/",name,".txt"))
            for(gen in target)
            cat(gen,"\n")
            sink()
        }
    }
}

how_much(res_a$data,name="Iron_deprivation_effect_at_Exp")
how_much(res_b$data,name="Iron_deprivation_effect_at_Stat")
how_much(res_c$data,name="Stat_effect_with_Fe")
how_much(res_d$data,name="Stat_effect_without_Fe")
how_much(res_a_minus_b$data,name="Stat_Fedeprivation_interaction")

library(readxl)
net_backbone=data.frame(read_xls("inputs/raw/ppi.xls"))
colnames(net_backbone)=c("RV_geneA","RV_geneB","Gene_name_A","Gene_name_B","Annotation_gene_A","Annotation_gene_B")

contrast=res_a$data

get_PPI_tab=function(contrast,name,net=net_backbone){
    contrast=data.frame(contrast)
    contrast$RV_geneA=rownames(contrast)
    contrast$RV_geneB=rownames(contrast)
    cosa=merge(net,contrast,by="RV_geneA")
    cosa=cosa[,c(1,2,3,4,5,6,8)]
    colnames(cosa)[c(2,7)]=c("RV_geneB","log2FC_gene_A")
    cosa=merge(cosa,contrast,by="RV_geneB")
    cosa=cosa[,c(2,1,3,4,5,6,7,9)]
    colnames(cosa)[c(1,8)]=c("RV_geneA","log2FC_gene_B")
    cosa$weight=cosa$log2FC_gene_A+cosa$log2FC_gene_B
    net_ok=cosa[,c(1,2,9)]
    dir.create("Outputs_def/networks/for_oslom",recursive=TRUE)
    dir.create("Outputs_def/networks/verbose",recursive=TRUE)
    write.table(cosa,paste0("Outputs_def/networks/verbose/",name,".txt"))
    write.table(net_ok,paste0("Outputs_def/networks/for_oslom/",name,".dat"),row.names=FALSE)
}

get_PPI_tab(contrast=res_a$data,name="Iron_deprivation_effect_at_Exp")
get_PPI_tab(contrast=res_b$data,name="Iron_deprivation_effect_at_Stat")
get_PPI_tab(contrast=res_c$data,name="Stat_effect_with_Fe")
get_PPI_tab(contrast=res_d$data,name="Stat_effect_without_Fe")
get_PPI_tab(contrast=res_a_minus_b$data,name="Stat_Fedeprivation_interaction")
