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
}

############################
### 1. Declare functions ###
############################

{
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
    write_results=function(tab,name,fdata=feature_data,th=0.01,th_size=0){
        if(length(which(rownames(fdata)!=rownames(tab)))>0){
            print("tablas no congruentes")
            return(0)}
        
        dir.create(paste0("Outputs_def/Stats_tables/th_",th,"_th_size_",th_size,"/",name),showWarnings = FALSE,recursive=TRUE)
        tab$t=tab$log2FoldChange/tab$lfcSE
        tab$RV=rownames(tab)
        tab$Gene=fdata$Gene_Symbol

        tab=tab[order(-tab$t),c(8,9,1,2,3,4,6)]
        
        upreg=tab[which(tab$BH<th & tab$log2FoldChange>th_size),]
        downreg=tab[which(tab$BH<th & tab$log2FoldChange<(-th_size)),]
        downreg=downreg[order(downreg$BH),]
        
        write_xlsx(tab,paste0("Outputs_def/Stats_tables/th_",th,"_th_size_",th_size,"/",name,"/all.xlsx"))
        write_xlsx(upreg,paste0("Outputs_def/Stats_tables/th_",th,"_th_size_",th_size,"/",name,"/upreg.xlsx"))
        write_xlsx(downreg,paste0("Outputs_def/Stats_tables/th_",th,"_th_size_",th_size,"/",name,"/downreg.xlsx"))

            }
    target_genesets=function(tab,name){
        for(th in c(0.05,0.01))
        {
            for(th_size in c(0,0.2,0.5,1))
            {
                dir.create(paste0("Outputs_def/target_sets/th_",th,"_th_size_",th_size,"/",name),showWarnings = FALSE,recursive=TRUE)
                
                set=(which(abs(tab$log2FoldChange)>th_size & tab$BH<th))
                print(paste("All: th=",th,"th_size=",th_size,":",length(set)))
                target=rownames(tab)[set]
                sink(paste0("Outputs_def/target_sets/th_",th,"_th_size_",th_size,"/",name,"/all.txt"))
                for(gen in target)
                cat(gen,"\n")
                sink()
                
                set=(which((tab$log2FoldChange)>th_size & tab$BH<th))
                print(paste("Up: th=",th,"th_size=",th_size,":",length(set)))
                target=rownames(tab)[set]
                sink(paste0("Outputs_def/target_sets/th_",th,"_th_size_",th_size,"/",name,"/upreg.txt"))
                for(gen in target)
                cat(gen,"\n")
                sink()
                
                set=(which((tab$log2FoldChange)<(-th_size) & tab$BH<th))
                print(paste("Up: th=",th,"th_size=",th_size,":",length(set)))
                target=rownames(tab)[set]
                sink(paste0("Outputs_def/target_sets/th_",th,"_th_size_",th_size,"/",name,"/downreg.txt"))
                for(gen in target)
                cat(gen,"\n")
                sink()
            }
        }
    }
    get_PPI_tab=function(contrast,name,net=net_backbone){
        contrast=data.frame(contrast)
        contrast$RV_geneA=rownames(contrast)
        contrast$RV_geneB=rownames(contrast)
        contrast$number_geneA=1:nrow(contrast)
        contrast$number_geneB=1:nrow(contrast)
        cosa=merge(net,contrast,by="RV_geneA")
        cosa=cosa[,c(1,2,14,3,4,5,6,8)]
        colnames(cosa)[c(2,8)]=c("RV_geneB","log2FC_gene_A")
        cosa=merge(cosa,contrast,by="RV_geneB")
        cosa=cosa[,c(2,1,3,16,8,10,4,5,6,7)]
        colnames(cosa)[c(1,3,4,6)]=c("RV_geneA","number_geneA","number_geneB","log2FC_gene_B")
        cosa$weight=(cosa$log2FC_gene_A+cosa$log2FC_gene_B)
        cosa$abs_weight=abs(cosa$log2FC_gene_A+cosa$log2FC_gene_B)
        cosa$censored_weight=pmax(0,(cosa$log2FC_gene_A+cosa$log2FC_gene_B))
        cosa$censored_weight_down=pmin(0,(cosa$log2FC_gene_A+cosa$log2FC_gene_B))

        net_raw=cosa[,c(3,4,11)]
        net_abs=cosa[,c(3,4,12)]
        net_censored=cosa[,c(3,4,13)]
        net_censored_down=cosa[,c(3,4,14)]

        dir.create("Outputs_def/networks/for_oslom/raw",recursive=TRUE)
        dir.create("Outputs_def/networks/for_oslom/abs",recursive=TRUE)
        dir.create("Outputs_def/networks/for_oslom/censored",recursive=TRUE)
        dir.create("Outputs_def/networks/for_oslom/censored_down",recursive=TRUE)
        dir.create("Outputs_def/networks/verbose",recursive=TRUE)
        
        write.table(cosa,paste0("Outputs_def/networks/verbose/",name,".txt"))
        write.table(net_raw,paste0("Outputs_def/networks/for_oslom/raw/",name,".dat"),row.names=FALSE)
        write.table(net_abs,paste0("Outputs_def/networks/for_oslom/abs/",name,".dat"),row.names=FALSE)
        write.table(net_censored,paste0("Outputs_def/networks/for_oslom/censored/",name,".dat"),row.names=FALSE)
        write.table(net_censored_down,paste0("Outputs_def/networks/for_oslom/censored_down/",name,".dat"),row.names=FALSE)
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
      
      # safe <- cosa
      contrast <- data.frame(contrast)
      log2FoldChange <- contrast$log2FoldChange
      contrast1 <- cbind(RV=rownames(contrast),log2FoldChange=log2FoldChange)
      contrast <- data.frame(contrast1)
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
      cosa$abs_weight=abs(cosa$log2FC_geneA_sh+cosa$log2FC_geneB_sh)
      cosa$censored_weight=pmax(0,(cosa$log2FC_geneA_sh+cosa$log2FC_geneB_sh))
      cosa$censored_weight_down=pmin(0,(cosa$log2FC_geneA_sh+cosa$log2FC_geneB_sh))
      
      # plink
      cosa$plink <- 2*pnorm(abs(cosa$Z),lower.tail = F)
      cosa$FDR <- p.adjust(cosa$plink,method = "BH")
      cosa <- cosa[cosa$FDR<0.05,]
      print(dim(cosa)[1])
      
      net_raw=cosa[,c(3,4,15)]
      net_abs=cosa[,c(3,4,17)]
      net_censored=cosa[,c(3,4,18)]
      net_censored_down=cosa[,c(3,4,19)]
      
      dir.create("Outputs_def/networks/for_oslom/raw_extended",recursive=TRUE,showWarnings = F)
      dir.create("Outputs_def/networks/for_oslom/abs_extended",recursive=TRUE,showWarnings = F)
      dir.create("Outputs_def/networks/for_oslom/censored_extended",recursive=TRUE,showWarnings = F)
      dir.create("Outputs_def/networks/for_oslom/censored_down_extended",recursive=TRUE,showWarnings = F)
      dir.create("Outputs_def/networks/verbose_extended",recursive=TRUE,showWarnings = F)
      
      write.table(cosa,paste0("Outputs_def/networks/verbose_extended/",name,".txt"))
      write.table(net_raw,paste0("Outputs_def/networks/for_oslom/raw_extended/",name,".dat"),row.names=FALSE)
      write.table(net_abs,paste0("Outputs_def/networks/for_oslom/abs_extended/",name,".dat"),row.names=FALSE)
      write.table(net_censored,paste0("Outputs_def/networks/for_oslom/censored_extended/",name,".dat"),row.names=FALSE)
      write.table(net_censored_down,paste0("Outputs_def/networks/for_oslom/censored_down_extended/",name,".dat"),row.names=FALSE)
    }
    
}

##################################################################################
### 2. Load data and filter desired cultures and sufficiently expressed genes. ###
##################################################################################

{
    ## Create output directory.
    dir.create("Outputs_def", showWarnings = F)
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

    colores_base=c("red","dodgerblue")
    fill_base=c("red",NA,"dodgerblue",NA)
    
    pl1=ggplot(datos)+
    geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=1.2,stroke=1.2)+
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

################
### 3 Do DE. ###
################

{
    ## get contrasts
    # Iron deprivation effects @ exp, and @ stat
    a=c(-2)
    b=c(3,-4)
    # Stat effects w Fe and wout Fe
    c=c(4,-2)
    d=c(3)
    ## Stat x iron deprivation interaction
    a_minus_b=c(-2,-3,4)

    dds <- DESeqDataSetFromMatrix(countData = reads,colData = metadata,design= ~ short_setup)
    dds <- DESeq(dds)
    
    
    raw_data <- results(dds)
    # raw_data <- data.frame(raw_data)
    # raw_data <- cbind(RV=rownames(raw_data),raw_data)
   
    
    
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
    
    
    dir.create("Outputs_def/Figures/signal/histograms/",recursive=TRUE)
    pdf("Outputs_def/Figures/signal/histograms/Iron_effect_at_exp.pdf",width=6,height=5)
    print(pl_a)
    dev.off()

    pdf("Outputs_def/Figures/signal/histograms/Iron_effect_at_stat.pdf",width=6,height=5)
    print(pl_b)
    dev.off()

    pdf("Outputs_def/Figures/signal/histograms/Stat_effect_w_Fe.pdf",width=6,height=5)
    print(pl_c)
    dev.off()

    pdf("Outputs_def/Figures/signal/histograms/Stat_effect_wout_Fe.pdf",width=6,height=5)
    print(pl_d)
    dev.off()

    pdf("Outputs_def/Figures/signal/histograms/interacts.pdf",width=6,height=5)
    print(pl_a_minus_b)
    dev.off()
    
    ## JAC: make 5 volcano plots:
    ## x axis: log2FoldChange; y axis (-log10(BH))
    ## Color either gray (for all genes), or firebrick4 (for example), for DE genes: fdr<0.01 & |log2FC|>0.5
    ## Add, with ggrepel, labels for the genes farther away from the origin in each plot. Use gene_names, not RV codes.
    dir.create("Outputs_def/Figures/signal/histograms/",recursive=TRUE)
}

#######################
### 3 Write results ###
#######################

{
    write_results(tab=data.frame(res_a$data),name="Iron_deprivation_effect_at_Exp")
    write_results(tab=data.frame(res_b$data),name="Iron_deprivation_effect_at_Stat")
    write_results(tab=data.frame(res_c$data),name="Stat_effect_with_Fe")
    write_results(tab=data.frame(res_d$data),name="Stat_effect_without_Fe")
    write_results(tab=data.frame(res_a_minus_b$data),name="Stat_Fedeprivation_interaction")

    target_genesets(res_a$data,name="Iron_deprivation_effect_at_Exp")
    target_genesets(res_b$data,name="Iron_deprivation_effect_at_Stat")
    target_genesets(res_c$data,name="Stat_effect_with_Fe")
    target_genesets(res_d$data,name="Stat_effect_without_Fe")
    target_genesets(res_a_minus_b$data,name="Stat_Fedeprivation_interaction")
}

#####################################
### 4 Build PPI weighted networks ###
#####################################

{
    net_backbone=data.frame(read_xls("inputs/raw/ppi.xls"))
    colnames(net_backbone)=c("RV_geneA","RV_geneB","Gene_name_A","Gene_name_B","Annotation_gene_A","Annotation_gene_B")

    get_PPI_tab(contrast=res_a$data,name="Iron_deprivation_effect_at_Exp")
    get_PPI_tab(contrast=res_b$data,name="Iron_deprivation_effect_at_Stat")
    get_PPI_tab(contrast=res_c$data,name="Stat_effect_with_Fe")
    get_PPI_tab(contrast=res_d$data,name="Stat_effect_without_Fe")
    get_PPI_tab(contrast=res_a_minus_b$data,name="Stat_Fe_deprivation_interaction")

    get_PPI_tab_extended(raw_data=raw_data1,contrast=res_a,name="Iron_deprivation_effect_at_Exp")
    get_PPI_tab_extended(raw_data=raw_data2,contrast=res_b,name="Iron_deprivation_effect_at_Stat")
    get_PPI_tab_extended(raw_data=raw_data3,contrast=res_c,name="Stat_effect_with_Fe")
    get_PPI_tab_extended(raw_data=raw_data4,contrast=res_d,name="Stat_effect_without_Fe")
    get_PPI_tab_extended(raw_data=raw_data5,contrast=res_a_minus_b,name="Stat_Fe_deprivation_interaction")
    
    
    
    
}
