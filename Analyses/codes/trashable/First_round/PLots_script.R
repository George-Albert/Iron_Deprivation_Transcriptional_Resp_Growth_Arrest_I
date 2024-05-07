
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
bar <- function(data,title){
  frame <- data.frame(data)
  sample <- factor(frame$Effects)
  n_gene <- frame$n_genes
  Sample_Effect <- reorder(sample,abs(n_gene))
  breaks_values <- pretty(n_gene) 
  color <- c("#FFA373","#50486D")
  bar_pl <- ggplot(frame,aes(x=Sample_Effect, y=n_gene,fill=Regulated))+
    geom_bar(stat = "identity",width = 0.8, color = "black")+coord_flip()+
    scale_y_continuous(breaks = breaks_values,labels = abs(breaks_values))+
    scale_x_discrete(expand = c(0,2.5))+
    theme_classic()+ ggtitle(title)+geom_hline(yintercept=0)+theme(aspect.ratio = .9)+
    labs(x="Samples", y = "Number of Genes")+theme(plot.title = element_text(hjust = 0.5))+
    scale_fill_manual(name="Regulated Genes",values = color)
  return(bar_pl)
}
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
{
## filename directory.
path = "th_0.01_th_size_0/"

file_1 = "Iron_deprivation_effect_at_Exp"
file_2 = "Iron_deprivation_effect_at_Stat"
file_3 = "Stat_effect_with_Fe"
file_4 = "Stat_effect_without_Fe"
file_5 = "Stat_Fedeprivation_interaction"

all = "all.xlsx"
down = "downreg.xlsx"
up = "upreg.xlsx"
value_up="UP"
value_dow="DOWN"
### Iron Deprivation effects at EXP

IDEE_all <- read.xlsx(paste0(path,"/",file_1,"/",all),1)
IDEE_all=IDEE_all[order(IDEE_all$RV),]
IDEE_down <- read.xlsx(paste0(path,"/",file_1,"/",down),1)
IDEE_down=IDEE_down[order(IDEE_down$RV),]
IDEE_up <- read.xlsx(paste0(path,"/",file_1,"/",up),1)
IDEE_up=IDEE_up[order(IDEE_up$RV),]

IDEE_reg <- rbind(IDEE_up,IDEE_down)
IDEE_reg$Regulated <-c(rep(value_up,nrow(IDEE_up)),rep(value_dow,nrow(IDEE_down)))

### Iron Deprivation effects at STAT

IDES_all <- read.xlsx(paste0(path,"/",file_2,"/",all),1)
IDES_all=IDES_all[order(IDES_all$RV),]
IDES_down <- read.xlsx(paste0(path,"/",file_2,"/",down),1)
IDES_down=IDES_down[order(IDES_down$RV),]
IDES_up <- read.xlsx(paste0(path,"/",file_2,"/",up),1)
IDES_up=IDES_up[order(IDES_up$RV),]

IDES_reg <- rbind(IDES_up,IDES_down)
IDES_reg$Regulated <-c(rep(value_up,nrow(IDES_up)),rep(value_dow,nrow(IDES_down)))
### Stat_effect_with_Fe

SEWFe_all <- read.xlsx(paste0(path,"/",file_3,"/",all),1)
SEWFe_all=SEWFe_all[order(SEWFe_all$RV),]
SEWFe_down <- read.xlsx(paste0(path,"/",file_3,"/",down),1)
SEWFe_down=SEWFe_down[order(SEWFe_down$RV),]
SEWFe_up <- read.xlsx(paste0(path,"/",file_3,"/",up),1)
SEWFe_up=SEWFe_up[order(SEWFe_up$RV),]

SEWFe_reg <- rbind(SEWFe_up,SEWFe_down)
SEWFe_reg$Regulated <-c(rep(value_up,nrow(SEWFe_up)),rep(value_dow,nrow(SEWFe_down)))
### Stat_effect_without_Fe

SEWOFe_all <- read.xlsx(paste0(path,"/",file_4,"/",all),1)
SEWOFe_all=SEWOFe_all[order(SEWOFe_all$RV),]
SEWOFe_down <- read.xlsx(paste0(path,"/",file_4,"/",down),1)
SEWOFe_down=SEWOFe_down[order(SEWOFe_down$RV),]
SEWOFe_up <- read.xlsx(paste0(path,"/",file_4,"/",up),1)
SEWOFe_up=SEWOFe_up[order(SEWOFe_up$RV),]

SEWOFe_reg <- rbind(SEWOFe_up,SEWOFe_down)
SEWOFe_reg$Regulated <-c(rep(value_up,nrow(SEWOFe_up)),rep(value_dow,nrow(SEWOFe_down)))

### Stat_Fedeprivation_interaction

SFeDI_all <- read.xlsx(paste0(path,"/",file_5,"/",all),1)
SFeDI_down <- read.xlsx(paste0(path,"/",file_5,"/",down),1)
SFeDI_up <- read.xlsx(paste0(path,"/",file_5,"/",up),1)

SFeDI_reg <- rbind(SFeDI_up,SFeDI_down)
SFeDI_reg$Regulated <-c(rep(value_up,nrow(SFeDI_up)),rep(value_dow,nrow(SFeDI_down)))
### Read the data.

metadata = read.table("processed_all_JAC/metadata.txt")
feature_data = read.table("processed_all_JAC/feature_data.txt")
reads = read.table("processed_all_JAC/reads.txt")
fpkms = read.table("processed_all_JAC/fpkm.txt")
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
  
  
  colores_base=c("red","dodgerblue")
  fill_base=c("red",NA,"dodgerblue",NA)
  
  pl1=ggplot(datos)+
    geom_point(aes(x=PC1,y=PC2,color=Culture,fill=samples),shape=21,size=1.2,stroke=1.2)+
    scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+
    xlab(paste0("PC1:",100*round(sum_pca$PC1[2],digits=3),"% variance explained"))+
    ylab(paste0("PC2:",100*round(sum_pca$PC2[2],digits=3),"% variance explained"))+theme_classic()+
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
  
  ## JAC: Add PCs within culture.
  
}

#=======================barplot up and down==============================
##### Create Regulated table of down and upregulated genes ######
{
d1 <- -dim(IDEE_reg[which(IDEE_reg$Regulated == "DOWN"),])[1]
up1 <- dim(IDEE_reg[which(IDEE_reg$Regulated == "UP"),])[1]
d2 <- -dim(IDES_reg[which(IDES_reg$Regulated == "DOWN"),])[1]
up2 <- dim(IDES_reg[which(IDES_reg$Regulated == "UP"),])[1]
d3 <- -dim(SEWFe_reg[which(SEWFe_reg$Regulated == "DOWN"),])[1]
up3 <- dim(SEWFe_reg[which(SEWFe_reg$Regulated == "UP"),])[1]
d4 <- -dim(SEWOFe_reg[which(SEWOFe_reg$Regulated == "DOWN"),])[1]
up4 <- dim(SEWOFe_reg[which(SEWOFe_reg$Regulated == "UP"),])[1]
d5 <- -dim(SFeDI_reg[which(SFeDI_reg$Regulated == "DOWN"),])[1]
up5 <- dim(SFeDI_reg[which(SFeDI_reg$Regulated == "UP"),])[1]
}

up_gene <- c(up1,up2,up3,up4,up5)
down_gene <- c(d1,d2,d3,d4,d5)
Samples <- c("Iron Deprivation effects at EXP","Iron Deprivation effects at STAT",
            "Stat_effect_with_Fe","Stat_effect_without_Fe","Stat_Fe_deprivation_interaction")
reg_table <- data.frame(Effects=rep(Samples,2),n_genes=c(up_gene,down_gene))
reg_table <- reg_table[order(reg_table$Effects),]
reg_table$Regulated <- c(rep("UP",5),rep("DOWN",5))

write.table(reg_table,file="processed_all_JAC/reg_table.txt")



bar <- function(data,title){
  frame <- data.frame(data)
  sample <- factor(frame$Effects)
  n_gene <- frame$n_genes
  # Sample_Effect <- reorder(sample,abs(n_gene))
  breaks_values <- pretty(n_gene) 
  
  bar_pl <- ggplot(frame,aes(x=sample, y=n_gene,fill=Regulated))+
    geom_bar(stat = "identity",width = 0.8, color = "black")+coord_flip()+
    scale_y_continuous(breaks = breaks_values,labels = abs(breaks_values))+
    scale_x_discrete(expand = c(0,2.5))+
    theme_classic()+ ggtitle(title)+geom_hline(yintercept=0)+theme(aspect.ratio = .9)+
    labs(x="Samples", y = "Number of Genes")+theme(plot.title = element_text(hjust = 0.5))+
    scale_fill_manual(name="Regulated Genes",values = c("#FFA373","#50486D"))
  return(bar_pl)
}

bar_plot <- bar(reg_table,"Regulated Genes")

pdf(file = "Outputs_def/Figures/Bar_plot.pdf",width=8,height=5)
print(bar_plot)
dev.off()

#                       ===================
#                       === DENSITY PLOT===  
#                       ===================
{
dens1 <- data.frame(LogFC=IDEE_reg$log2FoldChange,samples= rep(Samples[1],
                                                               nrow(IDEE_reg)))
dens2 <- data.frame(LogFC=IDES_reg$log2FoldChange,samples= rep(Samples[2],
                                                               nrow(IDES_reg)))
dens3 <- data.frame(LogFC=SEWFe_reg$log2FoldChange,samples= rep(Samples[3],
                                                                nrow(SEWFe_reg)))
dens4 <- data.frame(LogFC=SEWOFe_reg$log2FoldChange,samples= rep(Samples[4],
                                                               nrow(SEWOFe_reg)))
dens5 <- data.frame(LogFC=SFeDI_reg$log2FoldChange,samples= rep(Samples[5],
                                                               nrow(SFeDI_reg)))
dens <- rbind(dens1,dens2,dens3,dens4,dens5)

}

write.table(dens,file="processed_all_JAC/dens_table.txt")

density <- function(data,title){
  
  frame <- data.frame(data)
  lfc <- frame$LogFC
  y <- factor(frame$samples)
  sam <- reorder(y,lfc)
  color <- c("#00AFBB", "#E7B800","#50486D", "#FC4E07","#FFA373")
  
  ggplot(frame, aes(x=lfc,y=sam,fill=y)) +
    geom_density_ridges(scale = 4, alpha = 0.6) + xlim(range(-5,5))+
    scale_fill_manual(values = color)+geom_vline(xintercept=0,lty="dotted")+
    labs(x="LogFC", y = "Density",title=title)+
    theme_ridges(font_size = 8,center_axis_labels = TRUE,grid = FALSE)+
    theme(plot.title = element_text(hjust = 0.5))
    
  
}

pdf(file = "Outputs_def/Figures/dens_plot.pdf",width=9,height=5)

density(dens,"Density Plot")

dev.off()


#=====================Normalized log expression values vs samples======================


stripchart(exp1~short_setup)
exp1[which(rownames(exp1) %in% IDEE_up$RV),]

#=============================Volcano Plot=====================================
{
IDEE_all$Association <- "Not Sig."
IDEE_all[which(IDEE_all$RV %in% IDEE_up$RV),]$Association<- "UP"
IDEE_all[which(IDEE_all$RV %in% IDEE_down$RV),]$Association<- "DOWN"

IDES_all$Association <- "Not Sig."
IDES_all[which(IDES_all$RV %in% IDES_up$RV),]$Association<- "UP"
IDES_all[which(IDES_all$RV %in% IDES_down$RV),]$Association<- "DOWN"

SEWFe_all$Association <- "Not Sig."
SEWFe_all[which(SEWFe_all$RV %in% SEWFe_up$RV),]$Association<- "UP"
SEWFe_all[which(SEWFe_all$RV %in% SEWFe_down$RV),]$Association<- "DOWN"

SEWOFe_all$Association <- "Not Sig."
SEWOFe_all[which(SEWOFe_all$RV %in% SEWOFe_up$RV),]$Association<- "UP"
SEWOFe_all[which(SEWOFe_all$RV %in% SEWOFe_down$RV),]$Association<- "DOWN"

SFeDI_all$Association <- "Not Sig."
SFeDI_all[which(SFeDI_all$RV %in% SFeDI_up$RV),]$Association<- "UP"
SFeDI_all[which(SFeDI_all$RV %in% SFeDI_down$RV),]$Association<- "DOWN"
}

volcan_plot <- function(data,y,xintercept,th,ymax,title){
  
  datos <- data.frame(data)
  x <- datos$log2FoldChange
  y <- -(log10(y))
  y[is.infinite(y)] <- NA
  
  sizes <- c("UP" = 0.5, "DOWN" = 0.5, "Not Sig." = 0.5) 
  alphas <- c("UP" = 1, "DOWN" = 1, "Not Sig." = 0.5)
  #FFA373
  ##50486D
  mycolors <- c("UP" ="green","Not Sig." ="grey","DOWN" ="magenta")
  
  ggplot(datos,aes(x=x,y=y,
                   size=Association,
                   alpha=Association,
                   fill=Association))+
    geom_point(shape = 21,colour = "black")+
    ylab("-log10(FDR)")+xlab("Log2FC")+theme_minimal()+
    scale_fill_manual(values=mycolors)+
    scale_size_manual(values = sizes,guide="none")+
    scale_alpha_manual(values = alphas, guide="none")+
    geom_vline(xintercept=c(-xintercept, xintercept),col="black",linetype = "dashed")+
    geom_hline(yintercept=-log10(th), col="black",linetype = "dashed")+
    scale_x_continuous(limits = c(-(max(x)), max(x)))+
    scale_y_continuous(limits = c(-1, ymax),expand = expansion(0))+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))
}


pdf(file = "Outputs_def/Figures/volcano_plot.pdf",width=6,height=5)

volcan_plot(IDEE_all,y=IDEE_all$BH,xintercept = 0.0,th=0.01,17,title =Samples[1])

volcan_plot(IDES_all,y=IDES_all$BH,xintercept = 0.0,th=0.01,240,title =Samples[2])

volcan_plot(SEWFe_all,y=SEWFe_all$BH,xintercept = 0.0,th=0.01,300,title =Samples[3])

volcan_plot(SEWOFe_all,y=SEWOFe_all$BH,xintercept = 0.0,th=0.01,300,title =Samples[4])

volcan_plot(SFeDI_all,y=SFeDI_all$BH,xintercept = 0.0,th=0.01,105,title =Samples[5])

dev.off()


#=======================Scatter Plot================================

stat_effect <- data.frame(RV1=SEWFe_all$RV,LogFC1=SEWFe_all$log2FoldChange,
                          Association1=SEWFe_all$Association,
                          RV2=SEWOFe_all$RV,LogFC2=SEWOFe_all$log2FoldChange,
                          Association2=SEWOFe_all$Association)

SEWFe_all[which(stat_effect$RV1 != stat_effect$RV1 ),]

up_g <- subset(stat_effect,Association1=="UP" & Association2=="UP")
down_g <- stat_effect[which(stat_effect$Association1=="DOWN" & stat_effect$Association2=="DOWN"),]




scatter_pl <- function(data,x,y,up_g,down_g,title){
  
  ggplot(data, aes(x=x, y=y)) +geom_point(alpha=0.4)+
    geom_point(data=up_g,aes(x=LogFC1, y=LogFC2),alpha=0.4,colour="#FFA373")+
    geom_point(data=down_g,aes(x=LogFC1, y=LogFC2),alpha=0.4,colour="#50486D")+
    ggtitle(title)+theme_minimal_grid()+geom_abline(col="red")+
    xlab("LogFC1: STAT Effects with Fe")+
    ylab("LogFC2: STAT Effects without Fe")+
    theme(plot.title = element_text(hjust = 0.5))
  
    
}

data=stat_effect
x = stat_effect$LogFC1
y=stat_effect$LogFC2

require(stats)
fit_lm<-lm(y ~ x, data = data)
fit_lm

pdf(file = "Outputs_def/Figures/stat_eff_plot.pdf",width=6,height=5)

scatter_pl(data=data,x=x,y=y,up_g, down_g,"STAT Effects with and w/o Fe")

dev.off()
#====================Expression levels=========================================

dcols(exp)

exp1 <- data.frame(exp[which(rownames(exp) %in% IDEE_up$RV),])
exp1 <- exp1 - rowMeans(exp1)
exp1$RV <- rownames(exp1)
exp1 <- exp1 %>%
   select(RV, everything())


gen_1 <- exp[which(rownames(exp) %in% IDEE_up$RV[1]),]

Fe_NO_exp <- data.frame(Exp_level=cbind(c(exp1[,2],exp1[,3])))
Fe_NO_exp$Condition <- "Fe_NO_EXP"
Fe_NO_exp$RV <- exp1[,1]

Fe_YES_exp <- data.frame(Exp_level=cbind(c(exp1[,4],exp1[,5])))
Fe_YES_exp$Condition <- "Fe_YES_EXP"
Fe_YES_exp$RV <- exp1[,1]

Fe_NO_stat <- data.frame(Exp_level=cbind(c(exp1[,6],exp1[,7])))
Fe_NO_stat$Condition <- "Fe_NO_STAT"
Fe_NO_stat$RV <- exp1[,1]

Fe_YES_stat <- data.frame(Exp_level=cbind(c(exp1[,8],exp1[,9])))
Fe_YES_stat$Condition <- "Fe_YES_STAT"
Fe_YES_stat$RV <- exp1[,1]

DE_table <- rbind(Fe_NO_exp,Fe_YES_exp,Fe_NO_stat,Fe_YES_stat)



exp_plt <- function(data,title){
  datos <- data.frame(data)
  y <- datos$Exp_level
  x <- datos$Condition
  group <- datos$RV
  mycolors <- c("brown","blue","green","darkorchid3","grey",
                "cyan1","red","yellow","orange")
  fill <- c("brown","blue","green","darkorchid3","grey",
            "cyan1","red","yellow","orange")
  # mycolors <- c("red","blue","red","blue")
  # fill <- c("red","blue",NA,NA)
  shape <- c(21,21,1,1)
  
  ggplot(datos,aes(x,y,
                   group=group,
                   col=group,
                   shape=Condition,
                   fill=group))+
    geom_jitter()+
    scale_color_manual(values=mycolors)+
    scale_shape_manual(values=shape,guide="none")+
    scale_fill_manual(values=fill,guide="none")+
    labs(color="Genes:")+
    xlim("Fe_NO_EXP","Fe_YES_EXP","Fe_NO_STAT","Fe_YES_STAT")+
    xlab("Conditions")+scale_y_continuous(breaks=seq(-4,10,2))+
    ylab("log2-counts per million (logCPM)")+theme_minimal()+
    ggtitle(title)+theme(plot.title = element_text(hjust = 0.5))
    
}

pdf(file = "Outputs_def/Figures/level_exp_rem_mean.pdf",width=8,height=5)
exp_plt(DE_table,"Differential expression levels")
dev.off()


#==================expression level for each gene==========================
{
Rv2377c <- data.frame(Rv2377c=exp[which(rownames(exp) %in% IDEE_up$RV[1]),])
Rv2377c$Condition <- rownames(Rv2377c)
safe=Rv2377c
Rv2377c[Rv2377c =="C5_Fe_NO_1" | Rv2377c =="C5_Fe_NO_2"] <- "Fe_NO_EXP"
Rv2377c[Rv2377c =="C5_Fe_YES_1" | Rv2377c =="C5_Fe_YES_2"] <- "Fe_YES_EXP"
Rv2377c[Rv2377c =="C6_Fe_NO_1" | Rv2377c =="C6_Fe_NO_2"] <- "Fe_NO_STAT"
Rv2377c[Rv2377c =="C6_Fe_YES_1" | Rv2377c =="C6_Fe_YES_2"] <- "Fe_YES_STAT"

Rv2378c <- data.frame(Rv2378c=exp[which(rownames(exp) %in% IDEE_up$RV[2]),])
Rv2378c$Condition <- rownames(Rv2378c)
Rv2378c[Rv2378c =="C5_Fe_NO_1" | Rv2378c =="C5_Fe_NO_2"] <- "Fe_NO_EXP"
Rv2378c[Rv2378c =="C5_Fe_YES_1" | Rv2378c =="C5_Fe_YES_2"] <- "Fe_YES_EXP"
Rv2378c[Rv2378c =="C6_Fe_NO_1" | Rv2378c =="C6_Fe_NO_2"] <- "Fe_NO_STAT"
Rv2378c[Rv2378c =="C6_Fe_YES_1" | Rv2378c =="C6_Fe_YES_2"] <- "Fe_YES_STAT"

Rv2379c <- data.frame(Rv2379c=exp[which(rownames(exp) %in% IDEE_up$RV[3]),])
Rv2379c$Condition <- rownames(Rv2379c)
Rv2379c[Rv2379c =="C5_Fe_NO_1" | Rv2379c =="C5_Fe_NO_2"] <- "Fe_NO_EXP"
Rv2379c[Rv2379c =="C5_Fe_YES_1" | Rv2379c =="C5_Fe_YES_2"] <- "Fe_YES_EXP"
Rv2379c[Rv2379c =="C6_Fe_NO_1" | Rv2379c =="C6_Fe_NO_2"] <- "Fe_NO_STAT"
Rv2379c[Rv2379c =="C6_Fe_YES_1" | Rv2379c =="C6_Fe_YES_2"] <- "Fe_YES_STAT"

Rv2380c <- data.frame(Rv2380c=exp[which(rownames(exp) %in% IDEE_up$RV[4]),])
Rv2380c$Condition <- rownames(Rv2380c)
Rv2380c[Rv2380c =="C5_Fe_NO_1" | Rv2380c =="C5_Fe_NO_2"] <- "Fe_NO_EXP"
Rv2380c[Rv2380c =="C5_Fe_YES_1" | Rv2380c =="C5_Fe_YES_2"] <- "Fe_YES_EXP"
Rv2380c[Rv2380c =="C6_Fe_NO_1" | Rv2380c =="C6_Fe_NO_2"] <- "Fe_NO_STAT"
Rv2380c[Rv2380c =="C6_Fe_YES_1" | Rv2380c =="C6_Fe_YES_2"] <- "Fe_YES_STAT"

Rv2381c <- data.frame(Rv2381c=exp[which(rownames(exp) %in% IDEE_up$RV[5]),])
Rv2381c$Condition <- rownames(Rv2381c)
Rv2381c[Rv2381c =="C5_Fe_NO_1" | Rv2381c =="C5_Fe_NO_2"] <- "Fe_NO_EXP"
Rv2381c[Rv2381c =="C5_Fe_YES_1" | Rv2381c =="C5_Fe_YES_2"] <- "Fe_YES_EXP"
Rv2381c[Rv2381c =="C6_Fe_NO_1" | Rv2381c =="C6_Fe_NO_2"] <- "Fe_NO_STAT"
Rv2381c[Rv2381c =="C6_Fe_YES_1" | Rv2381c =="C6_Fe_YES_2"] <- "Fe_YES_STAT"

Rv2382c <- data.frame(Rv2382c=exp[which(rownames(exp) %in% IDEE_up$RV[6]),])
Rv2382c$Condition <- rownames(Rv2382c)
Rv2382c[Rv2382c =="C5_Fe_NO_1" | Rv2382c =="C5_Fe_NO_2"] <- "Fe_NO_EXP"
Rv2382c[Rv2382c =="C5_Fe_YES_1" | Rv2382c =="C5_Fe_YES_2"] <- "Fe_YES_EXP"
Rv2382c[Rv2382c =="C6_Fe_NO_1" | Rv2382c =="C6_Fe_NO_2"] <- "Fe_NO_STAT"
Rv2382c[Rv2382c =="C6_Fe_YES_1" | Rv2382c =="C6_Fe_YES_2"] <- "Fe_YES_STAT"

Rv2383c <- data.frame(Rv2383c=exp[which(rownames(exp) %in% IDEE_up$RV[7]),])
Rv2383c$Condition <- rownames(Rv2383c)
Rv2383c[Rv2383c =="C5_Fe_NO_1" | Rv2383c =="C5_Fe_NO_2"] <- "Fe_NO_EXP"
Rv2383c[Rv2383c =="C5_Fe_YES_1" | Rv2383c =="C5_Fe_YES_2"] <- "Fe_YES_EXP"
Rv2383c[Rv2383c =="C6_Fe_NO_1" | Rv2383c =="C6_Fe_NO_2"] <- "Fe_NO_STAT"
Rv2383c[Rv2383c =="C6_Fe_YES_1" | Rv2383c =="C6_Fe_YES_2"] <- "Fe_YES_STAT"

Rv2386c <- data.frame(Rv2386c=exp[which(rownames(exp) %in% IDEE_up$RV[8]),])
Rv2386c$Condition <- rownames(Rv2386c)
Rv2386c[Rv2386c =="C5_Fe_NO_1" | Rv2386c =="C5_Fe_NO_2"] <- "Fe_NO_EXP"
Rv2386c[Rv2386c =="C5_Fe_YES_1" | Rv2386c =="C5_Fe_YES_2"] <- "Fe_YES_EXP"
Rv2386c[Rv2386c =="C6_Fe_NO_1" | Rv2386c =="C6_Fe_NO_2"] <- "Fe_NO_STAT"
Rv2386c[Rv2386c =="C6_Fe_YES_1" | Rv2386c =="C6_Fe_YES_2"] <- "Fe_YES_STAT"

Rv3402c <- data.frame(Rv3402c=exp[which(rownames(exp) %in% IDEE_up$RV[9]),])
Rv3402c$Condition <- rownames(Rv3402c)
Rv3402c[Rv3402c =="C5_Fe_NO_1" | Rv3402c =="C5_Fe_NO_2"] <- "Fe_NO_EXP"
Rv3402c[Rv3402c =="C5_Fe_YES_1" | Rv3402c =="C5_Fe_YES_2"] <- "Fe_YES_EXP"
Rv3402c[Rv3402c =="C6_Fe_NO_1" | Rv3402c =="C6_Fe_NO_2"] <- "Fe_NO_STAT"
Rv3402c[Rv3402c =="C6_Fe_YES_1" | Rv3402c =="C6_Fe_YES_2"] <- "Fe_YES_STAT"
}

exp_plt_1 <- function(data,y,title){
  datos <- data.frame(data)
  
  x <- datos$Condition
  group <- datos$Condition
  mycolors <- c("red","blue","red","blue")
  fill <- c("red","blue",NA,NA)
  shape <- c(21,21,1,1)
  
  ggplot(datos,aes(x,y,
                   shape=Condition,
                   fill=Condition,
                   col=Condition))+
    geom_jitter()+
    scale_color_manual(values=mycolors)+
    scale_shape_manual(values=shape)+
    scale_fill_manual(values=fill)+
    xlim("Fe_NO_EXP","Fe_YES_EXP","Fe_NO_STAT","Fe_YES_STAT")+
    xlab("Conditions")+scale_y_continuous(breaks=seq(-4,10,2))+
    ylab("log2-counts per million (logCPM)")+theme_minimal()+
    guides(
      colour = guide_legend("Conditions"),
      fill = guide_legend("Conditions"),
      shape = guide_legend("Conditions"),
      
    )+
    ggtitle(title)+theme(plot.title = element_text(hjust = 0.5))
  
}

pdf(file = "Outputs_def/Figures/level_exp_each_gen.pdf",width=8,height=5)

exp_plt_1(Rv2377c,Rv2377c$Rv2377c,"Rv2377c")
exp_plt_1(Rv2378c,Rv2378c$Rv2378c,"Rv2378c")
exp_plt_1(Rv2379c,Rv2379c$Rv2379c,"Rv2379c")
exp_plt_1(Rv2380c,Rv2380c$Rv2380c,"Rv2380c")
exp_plt_1(Rv2381c,Rv2381c$Rv2381c,"Rv2381c")
exp_plt_1(Rv2382c,Rv2382c$Rv2382c,"Rv2382c")
exp_plt_1(Rv2383c,Rv2383c$Rv2383c,"Rv2383c")
exp_plt_1(Rv2386c,Rv2386c$Rv2386c,"Rv2386c")
exp_plt_1(Rv3402c,Rv3402c$Rv3402c,"Rv3402c")

dev.off()

