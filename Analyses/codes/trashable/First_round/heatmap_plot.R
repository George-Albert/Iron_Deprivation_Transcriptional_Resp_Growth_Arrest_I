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
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(gplots)
  library(circlize)
  library(magick)
}


############################
### 1. Declare functions ###
############################
dcols=function(x){data.frame(colnames(x))}
my_name <- function(v1) {
  deparse(substitute(v1))
}

calc_ht_size = function(ht, unit = "inch") {
  pdf(NULL)
  ht = draw(ht)
  w = ComplexHeatmap:::width(ht)
  w = convertX(w, unit, valueOnly = TRUE)
  h = ComplexHeatmap:::height(ht)
  h = convertY(h, unit, valueOnly = TRUE)
  dev.off()
  
  c(w, h)
  }

# reads <- reads[which(rownames(reads) %in% non_prot_cod$Gene_ID),]

#===============================
#=====Non Prot Coding Genes=====
#===============================
y <- DGEList(counts = reads)
y <- calcNormFactors(y)
design <- model.matrix(~short_setup,data=metadata)
v=voom(y,design,plot=TRUE)
exp=v$E
exp <- exp - rowMeans(exp)
# exp1 <- data.frame(exp[which(rownames(exp) %in% IDEE_protein$RV),])
# exp1$RV <- rownames(exp1)
exp1 <- exp[which(rownames(exp) %in% non_prot_cod$Gene_ID),]
exp2 <- exp[which(rownames(exp) %in% rRNA$Gene_ID),]
exp3 <- exp[which(rownames(exp) %in% tRNA$Gene_ID),]
exp4 <- exp[which(rownames(exp) %in% ncRNA$Gene_ID),]
exp5 <- exp[which(rownames(exp) %in% rownames(misc_RNA)),]



col<- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
myColours <-  c('red3', "white", 'royalblue')

my.matrix <- as.matrix(exp1)
my.matrix2 <- as.matrix(exp2)
my.matrix3<- as.matrix(exp3)
my.matrix4 <- as.matrix(exp4)
my.matrix5 <- as.matrix(exp5)

# gene_info <- data.frame(RV=rownames(exp1))
# heatmap(my.matrix)
# heatmap.2(my.matrix)

col_fun = colorRamp2(c(-3, 0, 3), c("#FFA373", "white","#50486D"))
col_fun(seq(-5, 5))

pdf(file = "Outputs_def/Figures/heatmap_non_prot_gen.pdf",width=5,height=14)

#==========================NON protein coding genes=============================
{
  col_fun = colorRamp2(c(round(min(exp1)), 0, round(max(exp1))), c("#FFA373",
                                                                   "white",
                                                                   "#50486D"))
  h1=Heatmap(my.matrix,name="Expression Level",
             width = ncol(my.matrix)*unit(5, "mm"), 
             height = nrow(my.matrix)*unit(5, "mm"),
             col = col_fun,column_order = sort(colnames(my.matrix)),
             row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 8),
             column_title = "Conditions",
             column_title_side = "bottom",
             show_column_dend = T,
             column_dend_side = "top",
             row_title = "Non Coding Protein Genes",
             show_row_names = T ,
             show_row_dend = T,
             row_dend_side = "right")
  draw(h1)
  size = calc_ht_size(h1)
  size
  
  # col_fun = colorRamp2(c(round(min(exp2)), 0, round(max(exp2))), c("#FFA373", "white","#50486D"))
  h2=Heatmap(my.matrix2,name="Expression Level",
             width = ncol(my.matrix2)*unit(5, "mm"), 
             height = nrow(my.matrix2)*unit(5, "mm"),
             col = col_fun,column_order = sort(colnames(my.matrix)),
             row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 8),
             column_title = "Conditions",
             column_title_side = "bottom",
             show_column_dend = T,
             column_dend_side = "top",
             row_title = "rRNA Genes",
             show_row_names = T ,
             show_row_dend = T,
             row_dend_side = "right")
  draw(h2)
  
  # col_fun = colorRamp2(c(round(min(exp3)), 0, round(max(exp3))), c("#FFA373", "white","#50486D"))
  h3=Heatmap(my.matrix3,name="Expression Level",
             width = ncol(my.matrix3)*unit(5, "mm"), 
             height = nrow(my.matrix3)*unit(5, "mm"),
             col = col_fun,column_order = sort(colnames(my.matrix)),
             row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 8),
             column_title = "Conditions",
             column_title_side = "bottom",
             show_column_dend = T,
             column_dend_side = "top",
             row_title = "tRNA Genes",
             show_row_names = T ,
             show_row_dend = T,
             row_dend_side = "right")
  
  draw(h3)
  
  # col_fun = colorRamp2(c(round(min(exp4)), 0, round(max(exp4))), c("#FFA373", "white","#50486D"))
  h4=Heatmap(my.matrix4,name="Expression Level",
             width = ncol(my.matrix4)*unit(5, "mm"), 
             height = nrow(my.matrix4)*unit(5, "mm"),
             col = col_fun,column_order = sort(colnames(my.matrix)),
             row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 8),
             column_title = "Conditions",
             column_title_side = "bottom",
             show_column_dend = T,
             column_dend_side = "top",
             row_title = "nc_RNA Genes",
             show_row_names = T,
             show_row_dend = T,
             row_dend_side = "right")
  draw(h4)
  
  # col_fun = colorRamp2(c(round(min(exp5)), 0, round(max(exp5))), c("#FFA373", "white","#50486D"))
  h5=Heatmap(my.matrix5,name="Expression Level",
             width = ncol(my.matrix5)*unit(5, "mm"), 
             height = nrow(my.matrix5)*unit(5, "mm"),
             col = col_fun,column_order = sort(colnames(my.matrix)),
             row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 8),
             column_title = "Conditions",
             column_title_side = "bottom",
             show_column_dend = T,
             column_dend_side = "top",
             row_title = "misc_RNA Genes",
             show_row_names = T,
             show_row_dend = T,
             row_dend_side = "right")
 
  draw(h5)
}
dev.off()


#============================MYOBACTIN HEATMAP===============================

exp6 <- exp[which(rownames(exp) %in% IDEE_reg$RV),]
my.matrix6 <- as.matrix(exp6)

pdf(file = "Outputs_def/Figures/heatmap_mycobactin_gen.pdf",width=4,height=4)

col_fun = colorRamp2(c(round(min(exp6)), 0, round(max(exp6))), c("#FFA373",
                                                                 "white",
                                                                 "#50486D"))
h6=Heatmap(my.matrix6,name="Expression Level",
           width = ncol(my.matrix6)*unit(5, "mm"), 
           height = nrow(my.matrix6)*unit(5, "mm"),
           col = col_fun,column_order = sort(colnames(my.matrix6)),
           row_names_gp = gpar(fontsize = 8),
           column_names_gp = gpar(fontsize = 8),
           column_title = "Conditions",
           column_title_side = "bottom",
           show_column_dend = T,
           column_dend_side = "top",
           row_title = "Mycobactin Genes",
           show_row_names = T ,
           show_row_dend = T,
           row_dend_side = "right")
draw(h6)
size = calc_ht_size(h6)
size
dev.off()
