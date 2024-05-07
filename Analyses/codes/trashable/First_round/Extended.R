


raw_data <- results(dds)

beta <- raw_data$log2FoldChange
SE <- raw_data$lfcSE
z <- (beta/SE)

raw_data$Z <- z


process_contrast_DEseq_extended=function(dds,contrast,th=0.05,title){
  name <-deparse(substitute(contrast))
  
  lista=list(add=resultsNames(dds)[contrast[which(contrast>0)]],remove=resultsNames(dds)[-contrast[which(contrast<0)]])
  result <- lfcShrink(dds, contrast=lista, type="ashr")
  
  beta <- result$log2FoldChange
  SE <- result$lfcSE
  z <- (beta/SE)
  
  result$Z <- z
  # names(result)[names(result) == 'Z'] <- paste0("Z_",siglas)
  
  
  result$BH=p.adjust(result$pvalue,method="BH")
  print(paste(name,length(which(result$BH<th))))
  nbins=50
  pl=histogram(result,result$pvalue,title)
  # pl=ggplot(result)+geom_histogram(aes(x=pvalue),bins=nbins,center = 1/(2*nbins))+ggtitle(name)
  output=list(data=result,figure=pl)
  return(output)
}

get_PPI_tab_extended=function(raw_data,contrast,name,net=net_backbone){
  raw=data.frame(raw_data)
  raw$RV_geneA=rownames(raw)
  raw$RV_geneB=rownames(raw)
  raw$number_geneA=1:nrow(raw)
  raw$number_geneB=1:nrow(raw)
  cosa=merge(net,raw,by="RV_geneA")
  cosa=cosa[,c("RV_geneA","RV_geneB.x",
               "number_geneA",
               "Gene_name_A","Gene_name_B",
               "Annotation_gene_A","Annotation_gene_B",
               "log2FoldChange","lfcSE")]
  colnames(cosa)[c(2,8,9)] = c("RV_geneB","log2FC_gene_A_raw","lfcSE_A")
  cosa=merge(cosa,raw,by="RV_geneB")
  cosa=cosa[c("RV_geneA.x","RV_geneB","number_geneA.x","number_geneB",
              "log2FC_gene_A_raw","lfcSE_A",
              "log2FoldChange","lfcSE",
              "Gene_name_A","Gene_name_B",
              "Annotation_gene_A","Annotation_gene_B")]
  colnames(cosa)[c(1,3,7,8)]=c("RV_geneA","number_geneA",
                               "log2FC_gene_B_raw","lfcSE_B")
  
  safe <- cosa
  contrast <- data.frame(contrast)
  log2FoldChange <- contrast$log2FoldChange
  contrast <- cbind(RV=rownames(contrast),log2FoldChange=log2FoldChange)
  colnames(contrast) <- c("RV_geneA","log2FoldChange")
  cosa <- merge(cosa,contrast,by="RV_geneA")
  colnames(cosa)[13] <- "log2FC_geneA_sh"
  colnames(contrast) <- c("RV_geneB","log2FoldChange")
  cosa=merge(cosa,contrast,by="RV_geneB")
  colnames(cosa)[14] <- "log2FC_geneB_sh"
  cosa$log2FC_geneA_sh <- as.numeric(cosa$log2FC_geneA_sh)
  cosa$log2FC_geneB_sh <- as.numeric(cosa$log2FC_geneB_sh)
  
  cosa$weight=(cosa$log2FC_geneA_sh+cosa$log2FC_geneB_sh)
  
  beta1 <- cosa$log2FC_gene_A_raw
  beta2 <- cosa$log2FC_gene_B_raw
  SE1 <- cosa$lfcSE_A
  SE2 <- cosa$lfcSE_B
  SD <- sqrt(SE1^2 +SE2^2 )
  
  cosa$Z <- (beta1+beta2)/SD
  
  
  # cosa$weight=(cosa$log2FC_gene_A+cosa$log2FC_gene_B)
  cosa$abs_weight=abs(cosa$log2FC_gene_A+cosa$log2FC_gene_B)
  cosa$censored_weight=pmax(0,(cosa$log2FC_gene_A+cosa$log2FC_gene_B))
  cosa$censored_weight_down=pmin(0,(cosa$log2FC_gene_A+cosa$log2FC_gene_B))
  
  # plink
  cosa$plink <- pnorm(cosa$Z)
  cosa$FDR <- p.adjust(cosa$plink,method = "BH")
  
  net_raw=cosa[,c(3,4,15)]
  net_abs=cosa[,c(3,4,17)]
  net_censored=cosa[,c(3,4,18)]
  net_censored_down=cosa[,c(3,4,19)]
  
  dir.create("Outputs_def/networks/for_oslom/raw_extended",recursive=TRUE)
  dir.create("Outputs_def/networks/for_oslom/abs_extended",recursive=TRUE)
  dir.create("Outputs_def/networks/for_oslom/censored_extended",recursive=TRUE)
  dir.create("Outputs_def/networks/for_oslom/censored_down_extended",recursive=TRUE)
  dir.create("Outputs_def/networks/verbose_extended",recursive=TRUE)
  
  write.table(cosa,paste0("Outputs_def/networks/verbose/",name,".txt"))
  write.table(net_raw,paste0("Outputs_def/networks/for_oslom/raw_extended/",name,".dat"),row.names=FALSE)
  write.table(net_abs,paste0("Outputs_def/networks/for_oslom/abs_extended/",name,".dat"),row.names=FALSE)
  write.table(net_censored,paste0("Outputs_def/networks/for_oslom/censored_extended/",name,".dat"),row.names=FALSE)
  write.table(net_censored_down,paste0("Outputs_def/networks/for_oslom/censored_down_extended/",name,".dat"),row.names=FALSE)
}

