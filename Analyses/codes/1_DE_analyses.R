#############################
### 0. Load dependencies. ###
#############################

{
    library(limma)
    library(edgeR)
    library(qvalue)
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
    library(eulerr)
    library(igraph)

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
      means=fpkms[,1:number_conditions]
      colnames(means)=unique(metadata[,column])
      ## This is unnecessary, nbut just for clarifying:
      for(i in 1:number_conditions)
        means[,i]=0
      for(i in 1:number_conditions)
      {
        set=which(metadata[,column]==colnames(means)[i])
        chunk=fpkms[,set]
        means[,i]=apply(chunk,1,mean)
      }
      means$max_mean=apply(means,1,max)
      genes_to_keep=rownames(means)[which(means$max_mean>threshold)]
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
      result_a <- results(dds, contrast=lista)
      colnames(result_a)[c(2,3)]=paste0(colnames(result_a)[c(2,3)],"_raw")
      result_b <- lfcShrink(dds, contrast=lista, type="ashr")
      colnames(result_b)[c(2,3)]=paste0(colnames(result_b)[c(2,3)],"_shrunken")

      result=cbind(result_a,result_b[,c(2,3)])
      result$BH=p.adjust(result$pvalue,method="BH")
      result=result[,c(1,2,3,7,8,4,5,9)]
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
        
        dir.create(paste0("outputs/Data/xlsx/th_",th,"_th_size_",th_size,"/",name),showWarnings = FALSE,recursive=TRUE)
        dir.create(paste0("outputs/Data/txt/th_",th,"_th_size_",th_size,"/"),showWarnings = FALSE,recursive=TRUE)

        #tab$t=tab$log2FoldChange_raw/tab$lfcSE_raw
        tab$RV=rownames(tab)
        tab$Gene=fdata$Gene_Symbol
        tab=tab[order(tolower(tab$RV)),c(9,10,1:8)]
        write.table(tab,paste0("outputs/Data/txt/",name,".txt"))

        tab=tab[order(-tab$stat),]
        
        upreg=tab[which(tab$BH<th & tab$log2FoldChange_shrunken>th_size),]
        downreg=tab[which(tab$BH<th & tab$log2FoldChange_shrunken<(-th_size)),]
        downreg=downreg[order(downreg$BH),]
        
        write_xlsx(tab,paste0("outputs/Data/xlsx/th_",th,"_th_size_",th_size,"/",name,"/all.xlsx"))
        write_xlsx(upreg,paste0("outputs/Data/xlsx/th_",th,"_th_size_",th_size,"/",name,"/upreg.xlsx"))
        write_xlsx(downreg,paste0("outputs/Data/xlsx/th_",th,"_th_size_",th_size,"/",name,"/downreg.xlsx"))
            }
    target_genesets=function(tab,name){
        for(th in c(0.05,0.01))
        {
            for(th_size in c(0,0.2,0.5,1))
            {
                dir.create(paste0("outputs/Data/target_sets/th_",th,"_th_size_",th_size,"/",name),showWarnings = FALSE,recursive=TRUE)
                
                set=(which(abs(tab$log2FoldChange_shrunken)>th_size & tab$BH<th))
                print(paste("All: th=",th,"th_size=",th_size,":",length(set)))
                target=rownames(tab)[set]
                sink(paste0("outputs/Data/target_sets/th_",th,"_th_size_",th_size,"/",name,"/all.txt"))
                for(gen in target)
                cat(gen,"\n")
                sink()
                
                set=(which((tab$log2FoldChange_shrunken)>th_size & tab$BH<th))
                print(paste("Up: th=",th,"th_size=",th_size,":",length(set)))
                target=rownames(tab)[set]
                sink(paste0("outputs/Data/target_sets/th_",th,"_th_size_",th_size,"/",name,"/upreg.txt"))
                for(gen in target)
                cat(gen,"\n")
                sink()
                
                set=(which((tab$log2FoldChange_shrunken)<(-th_size) & tab$BH<th))
                print(paste("Up: th=",th,"th_size=",th_size,":",length(set)))
                target=rownames(tab)[set]
                sink(paste0("outputs/Data/target_sets/th_",th,"_th_size_",th_size,"/",name,"/downreg.txt"))
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
}

##################################################################################
### 2. Load data and filter desired cultures and sufficiently expressed genes. ###
##################################################################################

{
    ## Create output directory.
    dir.create("outputs")
    ### Read the data. Jq: aqui había codigo redundante...
    metadata = read.table("inputs/processed/unfiltered/txts/metadata.txt")
    feature_data = read.table("inputs/processed/unfiltered/txts/feature_data.txt")
    reads = read.table("inputs/processed/unfiltered/txts/reads.txt")
    fpkms = read.table("inputs/processed/unfiltered/txts/fpkm.txt")
    rownames(feature_data)=feature_data$Gene_ID

    ## JQ: Esto ya no sería necesario, pero no importa
    metadata=metadata[which(metadata$Culture %in% c("C5","C6")),]
    metadata$short_setup=paste0(metadata$Culture,"_Fe_",metadata$Iron)
    metadata <- metadata[order(metadata$short_setup),]
    reads <- reads[,rownames(metadata)]
    fpkms <- fpkms[,rownames(metadata)]
    
    ## Check coherence:
    length(which(colnames(reads)!=rownames(metadata)))

    # Select geneset:
    filtrado=filter_OK(reads,fpkms,feature_data,threshold=2,metadata=metadata,column="short_setup")

    reads=filtrado[["filtered_reads"]]
    fpkms=filtrado[["filtered_fpkm"]]
    feature_data=filtrado[["filtered_features"]]
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
    sum_pca=data.frame(summary(pca)$importance)

    datos <- data.frame(pca$x)
    colnames(datos)=paste0("PC",c(1:ncol(datos)))
    length(which(rownames(datos)!=rownames(metadata)))
    datos=cbind(metadata,datos)

    colores_base=c("red","dodgerblue")
    fill_base=c("red",NA,"dodgerblue",NA)
    
    pl1=ggplot(datos)+
      geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=1.5,stroke=1.2)+
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
        legend.position="none",
        #legend.text=element_text(size=14),
        legend.key.size = unit(1, 'lines'))

    dir.create("outputs/Figures/",recursive=TRUE)
    pdf(file = "outputs/Figures/Fig_3A_PCA_C5_C6.pdf",width=4,height=4)
    print(pl1)
    dev.off()
    
    datos$culture_number=rep(0,8)
    datos$culture_number[which(metadata$Culture=="C6")]=1
    
    datos$iron_number=rep(0,8)
    datos$iron_number[which(metadata$Iron=="NO")]=1
    
    ## JQ: Vamos a ir imprimiendo los stats que se van sacando y que se mencionan luego en el texto.
    dir.create("outputs/Stats/",recursive=TRUE)

    sink("outputs/Stats/PC1_vs_culture_correlation_all_samples.txt")
    print("cor.test(datos$PC1,datos$culture_number)")
    cor.test(datos$PC1,datos$culture_number)
    sink()
    
    #t = 13.119, df = 6, p-value = 1.21e-05
    #alternative hypothesis: true correlation is not equal to 0
    #95 percent confidence interval:
    # 0.9057638 0.9970361
    #sample estimates:
          cor
    #0.9830122
    
    
    sink("outputs/Stats/PC2_vs_iron_correlation_stat_samples.txt")
    print("cor.test(datos$PC2[which(datos$Culture==\"C6\")],datos$iron_number[which(datos$Culture==\"C6\")])")
    cor.test(datos$PC2[which(datos$Culture=="C6")],datos$iron_number[which(datos$Culture=="C6")])
    sink()
    #t = 48.202, df = 2, p-value = 0.0004301
    #alternative hypothesis: true correlation is not equal to 0
    #95 percent confidence interval:
    #0.9785515 0.9999915
    #sample estimates:
    #      cor
    #0.9995699
    
    cors=rep(0,8)
    ps=rep(0,8)
    for(i in 9:16)
    {
        test=cor.test(datos[which(datos$Culture=="C5"),i],datos$iron_number[which(datos$Culture=="C5")])
        cors[i-8]=test$estimate
        ps[i-8]=test$p.value
    }
    
    cors
    #[1]  0.7040703  0.5709442  0.2677408  0.1600443 -0.2623416  0.8868217 -0.3868404 -0.9016488
    ps
    #[1] 0.29592974 0.42905582 0.73225917 0.83995571 0.73765839 0.11317825 0.61315963 0.09835119
    # Nothing significant to save or report.
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
    ## Stat x iron deprivation interaction. JQ: Ojo el nombre es chungo: no es a_minus_b sino al reves. No lo cambio para no liar
    a_minus_b=c(2,3,-4)
    

    dds <- DESeqDataSetFromMatrix(countData = reads,colData = metadata,design= ~ short_setup)
    dds <- DESeq(dds)
    resultsNames(dds)
   
   ## JQ: THis is to obtain the normalized counts for the heatmap.
   
   get_voom_like_normalized_data=function(desq){
     cuentas=counts(desq)
     normalization_coeficients=dds@colData$sizeFactor
     depths=apply(cuentas,2,sum)+1
     cuentas=cuentas+0.5
     
     normalized_counts=cuentas
     for(i in 1:ncol(normalized_counts))
     normalized_counts[,i]=normalized_counts[,i]/(depths[i]*normalization_coeficients[i]/1E6)
     
     normalized_counts=log2(normalized_counts)
     return(normalized_counts)
    }
   normalized_dds=get_voom_like_normalized_data(dds)
   
   medias=normalized_dds[,1:2]
   
   set_no=which(metadata$Iron=="NO")
   set_yes=which(metadata$Iron=="YES")
   
   medias[,1]=apply(normalized_dds[,set_no],1,mean)
   medias[,2]=apply(normalized_dds[,set_yes],1,mean)
   
   for(i in set_no)
   normalized_dds[,i]=normalized_dds[,i]-medias[,1]

   for(i in set_yes)
   normalized_dds[,i]=normalized_dds[,i]-medias[,2]
   
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

    ##JQ These get saved,but got unpublished
    dir.create("outputs/Figures/unpublished/signal/histograms/",recursive=TRUE)
    pdf("outputs/Figures/unpublished/signal/histograms/Iron_effect_at_exp.pdf",width=6,height=5)
    print(pl_a)
    dev.off()

    pdf("outputs/Figures/unpublished/signal/histograms/Iron_effect_at_stat.pdf",width=6,height=5)
    print(pl_b)
    dev.off()

    pdf("outputs/Figures/unpublished/signal/histograms/Stat_effect_w_Fe.pdf",width=6,height=5)
    print(pl_c)
    dev.off()

    pdf("outputs/Figures/unpublished/signal/histograms/Stat_effect_wout_Fe.pdf",width=6,height=5)
    print(pl_d)
    dev.off()

    pdf("outputs/Figures/unpublished/signal/histograms/interacts.pdf",width=6,height=5)
    print(pl_a_minus_b)
    dev.off()
    
    ## Possible other vizs:
    ## JAC: make 5 volcano plots:
    ## x axis: log2FoldChange; y axis (-log10(BH))
    ## Color either gray (for all genes), or firebrick4 (for example), for DE genes: fdr<0.01 & |log2FC|>0.5
    ## Add, with ggrepel, labels for the genes farther away from the origin in each plot. Use gene_names, not RV codes.
}

#######################
### 4 Write results ###
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
    
    ## Write global, summary table with 5 sheets.
    
    Iron_deprivation_effect_at_Exp= cbind(rownames(res_a$data),        data.frame(res_a$data)[,c(2,4,8)])
    Iron_deprivation_effect_at_Stat=cbind(rownames(res_b$data),        data.frame(res_b$data)[,c(2,4,8)])
    Stat_effect_with_Fe=            cbind(rownames(res_c$data),        data.frame(res_c$data)[,c(2,4,8)])
    Stat_effect_without_Fe=         cbind(rownames(res_d$data),        data.frame(res_d$data)[,c(2,4,8)])
    Stat_Fedeprivation_interaction= cbind(rownames(res_a_minus_b$data),data.frame(res_a_minus_b$data)[,c(2,4,8)])
    
    colnames(Iron_deprivation_effect_at_Exp)=c("Gene_ID","log2FC_raw","log2FC_shrunk","Fdr_BH")
    colnames(Iron_deprivation_effect_at_Stat)=c("Gene_ID","log2FC_raw","log2FC_shrunk","Fdr_BH")
    colnames(Stat_effect_with_Fe)=c("Gene_ID","log2FC_raw","log2FC_shrunk","Fdr_BH")
    colnames(Stat_effect_without_Fe)=c("Gene_ID","log2FC_raw","log2FC_shrunk","Fdr_BH")
    colnames(Stat_Fedeprivation_interaction)=c("Gene_ID","log2FC_raw","log2FC_shrunk","Fdr_BH")

    
    DE_results_summary=list(
    "Iron_deprivation_effect_at_Exp" = Iron_deprivation_effect_at_Exp,
    "Iron_deprivation_effect_at_Stat" = Iron_deprivation_effect_at_Stat,
    "Stat_effect_with_Fe" = Stat_effect_with_Fe,
    "Stat_effect_without_Fe" = Stat_effect_without_Fe,
    "Stat_Fedeprivation_interaction" = Stat_Fedeprivation_interaction
    )
    
    write_xlsx(DE_results_summary,"outputs/Data/xlsx/DE_results_summary.xlsx")
    
    ## Write Supplementary Tables
    
    Iron_deprivation_effect_at_Exp= cbind(rownames(res_a$data),        data.frame(res_a$data)[,c(1,4,6,7,8)])
    Iron_deprivation_effect_at_Stat=cbind(rownames(res_b$data),        data.frame(res_b$data)[,c(1,4,6,7,8)])
    Stat_effect_with_Fe=            cbind(rownames(res_c$data),        data.frame(res_c$data)[,c(1,4,6,7,8)])
    Stat_effect_without_Fe=         cbind(rownames(res_d$data),        data.frame(res_d$data)[,c(1,4,6,7,8)])
    Stat_Fedeprivation_interaction= cbind(rownames(res_a_minus_b$data),data.frame(res_a_minus_b$data)[,c(1,4,6,7,8)])
    
    length(which(Stat_effect_with_Fe$log2FoldChange_shrunken>0 & Stat_effect_with_Fe$BH<0.01))
    #[1] 1144
    length(which(Stat_effect_with_Fe$log2FoldChange_shrunken<0 & Stat_effect_with_Fe$BH<0.01))
    #[1] 1083
    length(which(Stat_effect_without_Fe$log2FoldChange_shrunken>0 & Stat_effect_without_Fe$BH<0.01))
    #[1] 1252
    length(which(Stat_effect_without_Fe$log2FoldChange_shrunken<0 & Stat_effect_without_Fe$BH<0.01))
    #[1] 1149
    
    
    length(which(Stat_effect_with_Fe$BH<0.01 & (feature_data$Type %in% "protein_coding")))
    # 2200
    length(which(Stat_effect_without_Fe$BH<0.01 & (feature_data$Type %in% "protein_coding")))
    # 2369

    
    length(which(Stat_effect_with_Fe$BH<0.01 & (!feature_data$Type %in% "protein_coding")))
    # 27
    length(which(Stat_effect_without_Fe$BH<0.01 & (!feature_data$Type %in% "protein_coding")))
    # 32
    
    
    ## JQ. Estos tres números me cambian con respevcto a lo que había.
    length(which(Stat_effect_with_Fe$log2FoldChange_shrunken<0 & Stat_effect_with_Fe$BH<0.01 & (feature_data$Type %in% "protein_coding")))
    #[1] 1081 (antes 1083)
    length(which(Stat_effect_without_Fe$log2FoldChange_shrunken>0 & Stat_effect_without_Fe$BH<0.01 & (feature_data$Type %in% "protein_coding")))
    #[1] 1221 (antes 1252)
    length(which(Stat_effect_without_Fe$log2FoldChange_shrunken<0 & Stat_effect_without_Fe$BH<0.01 & (feature_data$Type %in% "protein_coding")))
    #[1] 1148 (antes 1149)

    
    colnames(Iron_deprivation_effect_at_Exp)=c("RV_number","BaseMean","log2FC","Statistics","P_value","FDR")
    colnames(Iron_deprivation_effect_at_Stat)=c("RV_number","BaseMean","log2FC","Statistics","P_value","FDR")
    colnames(Stat_effect_with_Fe)=c("RV_number","BaseMean","log2FC","Statistics","P_value","FDR")
    colnames(Stat_effect_without_Fe)=c("RV_number","BaseMean","log2FC","Statistics","P_value","FDR")
    colnames(Stat_Fedeprivation_interaction)=c("RV_number","BaseMean","log2FC","Statistics","P_value","FDR")
   
    
    add_names=function(caso){
        caso$Gene_name=feature_data$Gene_Symbol
        return(caso[,c(1,7,2:6)])
    }
    Iron_deprivation_effect_at_Exp=add_names(Iron_deprivation_effect_at_Exp)
    Iron_deprivation_effect_at_Stat=add_names(Iron_deprivation_effect_at_Stat)
    Stat_effect_with_Fe=add_names(Stat_effect_with_Fe)
    Stat_effect_without_Fe=add_names(Stat_effect_without_Fe)
    Stat_Fedeprivation_interaction=add_names(Stat_Fedeprivation_interaction)


    
    Table_S2=list(
    "Iron_deprivation_effect_at_Exp" = Iron_deprivation_effect_at_Exp,
    "Iron_deprivation_effect_at_Stat" = Iron_deprivation_effect_at_Stat,
    "Stat_effect_with_Fe" = Stat_effect_with_Fe,
    "Stat_effect_without_Fe" = Stat_effect_without_Fe,
    "Stat_Fedeprivation_interaction" = Stat_Fedeprivation_interaction
    )
    
    write_xlsx(Table_S2,"outputs/Data/xlsx/Table_S4.xlsx")

    Table_S1=cbind(feature_data[,c(1,3,5,6)],reads)
    colnames(Table_S1)=c("RV_number","Gene_name","Description","Biotype","Exp5-Fe_1","Exp5-Fe_2","Exp5+Fe_1","Exp5+Fe_2","Stat6-Fe_1","Stat6-Fe_2","Stat6+Fe_1","Stat6+Fe_2")
    write_xlsx(Table_S1,"outputs/Data/xlsx/Table_S3.xlsx")
    
}

#########################################################
### 5 Retrieve assorted visualizations and statistics ###
#########################################################

{
    iron_atexp=res_a$data
    iron_atstat=res_b$data
    stat_withiron=res_c$data
    stat_withoutiron=res_d$data
    int=res_a_minus_b$data

    ## Check coherence & count genes post QC

    length(which(rownames(stat_withiron) == rownames(stat_withoutiron)))
    # 3899
    length(which(rownames(stat_withiron) == rownames(iron_atstat)))
    # 3899
    length(which(rownames(stat_withiron) == rownames(iron_atexp)))
    # 3899
    length(which(rownames(stat_withiron) == rownames(int)))


    threshold=0.01
    threshold_int=0.05

    # Fraction of DEGs
    100*length(which(iron_atexp$BH<threshold))/nrow(iron_atexp)
    # 0.256476
    100*length(which(iron_atstat$BH<threshold))/nrow(iron_atstat)
    # 20.44114
    100*length(which(stat_withiron$BH<threshold))/nrow(stat_withiron)
    # 57.11721
    100*length(which(stat_withoutiron$BH<threshold))/nrow(stat_withoutiron)
    # 61.57989

    ## Fractions & #s of interactions.
    100*length(which(int$BH<threshold))/nrow(int)
    # 12.38779
    100*length(which(int$BH<threshold_int))/nrow(int)
    # 18.31239

    length(which(int$BH<threshold))
    # 483
    length(which(int$BH<threshold_int))
    # 714

    ## Scatter and Venn diagrams

    resumen=data.frame(
    logFC_with=stat_withiron$log2FoldChange_shrunken,
    logFC_without=stat_withoutiron$log2FoldChange_shrunken,
    logFC_int=int$log2FoldChange_shrunken,
    BH_with=stat_withiron$BH,
    BH_without=stat_withoutiron$BH,
    BH_int=int$BH)
    rownames(resumen)=rownames(stat_withiron)

    resumen$label="Background"
    resumen$label[which((resumen$BH_with<threshold & resumen$BH_without<threshold) & resumen$logFC_with*resumen$logFC_without>0)]="Coherent response"
    resumen$label[which((resumen$BH_with<threshold & resumen$BH_without<threshold) & resumen$logFC_with*resumen$logFC_without<0)]="Flipped response"
    resumen=resumen[order(resumen$label),]

    DEGs_1=rownames(resumen)[which(resumen$BH_with<threshold)]
    DEGs_2=rownames(resumen)[which(resumen$BH_without<threshold)]

    Genes_list=list(Iron_rich=DEGs_1,Iron_deprived=DEGs_2)

    #colores=c("#B8141F","#A3C9F5","#FF7A14")

    names(Genes_list)=c("Fe/+","Fe/-")
    pdf("outputs/Figures/Fig_S1_Euler_diagram.pdf",width=4,height=3)
    plot(euler(Genes_list), quantities = list(fontsize=10),n=2000L, fills = list(fill = c("#A3C9F5","#FF7A14"), alpha = 0.5))
    dev.off()

    a=length(which(resumen$BH_with<threshold & resumen$BH_without<threshold))
    b=length(which(resumen$BH_with>threshold & resumen$BH_without<threshold))
    c=length(which(resumen$BH_with<threshold & resumen$BH_without>threshold))
    d=length(which(resumen$BH_with>threshold & resumen$BH_without>threshold))


    tab <- matrix(c(a,b,c,d), nrow = 2,
                               dimnames =
                        list(c("DEG_with", "No_DEG_with"),
                             c("DEG_without", "No_DEG_without")))
    test=fisher.test(tab, alternative = "greater")
    test$p.value
    # 1.442693e-286

    pl_scatter=ggplot(resumen)+ geom_point(aes(x=logFC_with,y=logFC_without,color=label),alpha=0.8)+scale_colour_manual(values=c("grey20","forestgreen","firebrick4"))+
    geom_abline(slope=1,intercept=0,color="forestgreen")+xlab("Response to growth arrest in iron-rich cultures (log2FC)")+ylab("Response to growth arrest in iron-deprived cultures (log2FC)")+
    theme(legend.position=c(0.8,0.2),legend.title=element_blank())

    pdf("outputs/Figures/Fig_3C_scatter_def.pdf",width=5,height=5)
    print(pl_scatter)
    dev.off()

    pl_scatter=ggplot(resumen)+ geom_point(aes(x=logFC_with,y=logFC_without,color=label),alpha=0.8)+scale_colour_manual(values=c("grey20","forestgreen","firebrick4"))+
    geom_abline(slope=1,intercept=0,color="forestgreen")+xlab("Response to growth arrest in iron-rich cultures (log2FC)")+ylab("Response to growth arrest in iron-deprived cultures (log2FC)")+
    theme(legend.position="none",legend.title=element_blank())

    pdf("outputs/Figures/Fig_3C_scatter_def_bg.pdf",width=5,height=5)
    print(pl_scatter)
    dev.off()

    test=cor.test(resumen$logFC_with,resumen$logFC_without)
    test$estimate
    #0.8784119
    test$p.value
    #0

    length(which(resumen$label=="Coherent response"))
    # 1874
    length(which(resumen$label=="Coherent response"))/(nrow(resumen)-length(which(resumen$label=="Background")))
    # 0.9847609
    length(which(resumen$BH_int<0.05))
    # 714
    length(which(resumen$BH_without<threshold))
    # 2401
    length(which(resumen$BH_with<threshold & resumen$BH_without<threshold))
    # 1903

    threshold_int=0.05
    threshold=0.01
    resumen$label_int="Background"
    resumen$label_int[which( resumen$logFC_with>0 & resumen$logFC_without>0)]="UP"
    resumen$label_int[which( resumen$logFC_with>0 & resumen$logFC_without>0 & resumen$BH_int<threshold_int & resumen$logFC_int>0)]="UP_More_without_Fe"
    resumen$label_int[which( resumen$logFC_with>0 & resumen$logFC_without>0 & resumen$BH_int<threshold_int & resumen$logFC_int<0)]="UP_More_with_Fe"

    resumen$label_int[which( resumen$logFC_with<0 & resumen$logFC_without<0)]="DOWN"
    resumen$label_int[which( resumen$logFC_with<0 & resumen$logFC_without<0 & resumen$BH_int<threshold_int & resumen$logFC_int>0)]="DOWN_More_with_Fe"
    resumen$label_int[which( resumen$logFC_with<0 & resumen$logFC_without<0 & resumen$BH_int<threshold_int & resumen$logFC_int<0)]="DOWN_More_without_Fe"

    resumen$label_int[which(resumen$logFC_with<0 & resumen$logFC_without>0 & resumen$BH_int<threshold_int)]="FLIPPED_int_down_up_weak"
    resumen$label_int[which(resumen$logFC_with>0 & resumen$logFC_without<0 & resumen$BH_int<threshold_int)]="FLIPPED_int_up_down_weak"
    resumen$label_int[which(resumen$BH_with<threshold & resumen$BH_without<threshold & resumen$logFC_with<0 & resumen$logFC_without>0 & resumen$BH_int<threshold_int)]="FLIPPED_int_down_up_strong"
    resumen$label_int[which(resumen$BH_with<threshold & resumen$BH_without<threshold & resumen$logFC_with>0 & resumen$logFC_without<0 & resumen$BH_int<threshold_int)]="FLIPPED_int_up_down_strong"


    data.frame(summary(factor(resumen$label_int)))
    # summary.factor.resumen.label_int..
    # Background                                                378
    # DOWN                                                     1480
    # DOWN_More_with_Fe                                          27
    # DOWN_More_without_Fe                                      175
    # FLIPPED_int_down_up_strong                                 17
    # FLIPPED_int_down_up_weak                                   85
    # FLIPPED_int_up_down_strong                                 12
    # FLIPPED_int_up_down_weak                                   72
    # UP                                                       1327
    # UP_More_with_Fe                                            72
    # UP_More_without_Fe                                        254

    ## Assym. distributions

    up_chunk=resumen[which(resumen$label_int %in% c("UP_More_with_Fe","UP_More_without_Fe")),]
    up_chunk$int_rebuilt=abs(up_chunk$logFC_without)-abs(up_chunk$logFC_with)
    up_chunk$label_int="Upregulated_genes"

    down_chunk=resumen[which(resumen$label_int %in% c("DOWN_More_with_Fe","DOWN_More_without_Fe")),]
    down_chunk$int_rebuilt=abs(down_chunk$logFC_without)-abs(down_chunk$logFC_with)
    down_chunk$label_int="Downregulated_genes"

    tab=rbind(up_chunk,down_chunk)

    pl_1=ggplot(up_chunk)+geom_density(aes(x=int_rebuilt),fill="#50486D")+xlim(-5,5)+xlab("")+theme(legend.position="none")+geom_vline(xintercept=0)
    pl_2=ggplot(down_chunk)+geom_density(aes(x=int_rebuilt),fill="#FFA373")+xlim(-5,5)+xlab("Difference of absolute effect sizes for growth arrest")+theme(legend.position="none")+geom_vline(xintercept=0)

    pl_both=plot_grid(pl_1,pl_2,ncol=1,align="h")

    pdf("outputs/Figures/Fig_4A_coherent_interactions.pdf",width=4,height=3)
    print(pl_both)
    dev.off()


    ## Para definir UP_More_with_Fe UP_More_without_Fe DOWN_More_without_Fe and DOWN_More_with_Fe threshold aquí pòdria ser igualmente 0,05, o 0.01: lo dejo a 1 (no threshold) porque así es todo más fácil de explicar.

    #summary(factor(resumen$label_int))
    dir.create("outputs/Data/target_sets/specific_sets",recursive=TRUE)

    set=rownames(resumen)[which(resumen$label_int=="UP_More_without_Fe")]
    sink("outputs/Data/target_sets/specific_sets/up_more_wout_Fe.txt")
    for(gen in set)
    cat(gen,"\n")
    sink()
    
    normalized_counts_EU=normalized_dds[which(rownames(normalized_dds) %in% set),]

    set=rownames(resumen)[which(resumen$label_int=="DOWN_More_without_Fe")]
    sink("outputs/Data/target_sets/specific_sets/down_more_wout_Fe.txt")
    for(gen in set)
    cat(gen,"\n")
    sink()
    
    normalized_counts_ED=normalized_dds[which(rownames(normalized_dds) %in% set),]


    set=rownames(resumen)[which(resumen$label_int=="UP_More_with_Fe")]
    sink("outputs/Data/target_sets/specific_sets/up_more_with_Fe.txt")
    for(gen in set)
    cat(gen,"\n")
    sink()

    set=rownames(resumen)[which(resumen$label_int=="DOWN_More_with_Fe")]
    sink("outputs/Data/target_sets/specific_sets/down_more_with_Fe.txt")
    for(gen in set)
    cat(gen,"\n")
    sink()
    
    
    load_and_pretty_up_go=function(path){
    go=data.frame(read_xls(path))
    go=go[which(!duplicated(go$ID)),c(1,2,4,5,10,11,12)]
    colnames(go)=c("GO_Term_ID","GO_Term_Name","P_value","Fdr","Percent_of_genes_in_target","Number_of_Genes_in_target","Genes_in_target")
    go=go[order(go$Fdr),]
    return(go)
    }
    
    go1=load_and_pretty_up_go("outputs/Data/ClueGO_results/sets_to_visualize_as_networks/Iron_at_stat_up_0.05/Iron_at_stat_up_0.05.xls")
    go2=load_and_pretty_up_go("outputs/Data/ClueGO_results/sets_to_visualize_as_networks/Iron_at_stat_down_0.05/Iron_at_stat_down_0.05.xls")
    go3=load_and_pretty_up_go("outputs/Data/ClueGO_results/sets_to_visualize_as_networks/UP_more_withoutFE0.05/UP_more_withoutFe_0.05.xls")
    go4=load_and_pretty_up_go("outputs/Data/ClueGO_results/sets_to_visualize_as_networks/Down_more_withoutFe0.05/Down_More_withoutF0.05.xls")
    
    guide=Iron_deprivation_effect_at_Stat
    
    guide$Iron_up=0
    guide$Iron_up[which(guide$FDR<0.01 & guide$log2FC>0)]=1
    
    guide$Iron_down=0
    guide$Iron_down[which(guide$FDR<0.01 & guide$log2FC<0)]=1
    
    set_enhanced_up=rownames(resumen)[which(resumen$label_int=="UP_More_without_Fe")]
    set_enhanced_down=rownames(resumen)[which(resumen$label_int=="DOWN_More_without_Fe")]

    guide$Enhanced_growth_arrest_up=0
    guide$Enhanced_growth_arrest_up[which(rownames(guide) %in% set_enhanced_up)]=1

    guide$Enhanced_growth_arrest_down=0
    guide$Enhanced_growth_arrest_down[which(rownames(guide) %in% set_enhanced_down)]=1
    
    summary(factor(guide$Iron_up))
    #   0    1
    #3447  452
    summary(factor(guide$Iron_down))
    #   0    1
    #3554  345
    summary(factor(guide$Enhanced_growth_arrest_up))
    #   0    1
    #3645  254
    summary(factor(guide$Enhanced_growth_arrest_down))
    #   0    1
    #3724  175
    
    guide=guide[,c(1,2,8:11)]

    Table_S4=list(
    "Readme"=guide,
    "Iron_up_@Stat6" = go1,
    "Iron_down_@Stat6" = go2,
    "Enhanced_growtharrest_up_@-Fe" = go3,
    "Enhanced_growtharrest_down_@-Fe" = go4
    )
    
    write_xlsx(Table_S4,"outputs/Data/xlsx/Table_S4.xlsx")

    

    threshold=0.01
    hits_iron_atexp=iron_atexp[which(iron_atexp$BH<threshold),]
    hits_iron_atexp$test="Fe deprivation at Exp5"

    hits_iron_atstat=iron_atstat[which(iron_atstat$BH<threshold),]
    hits_iron_atstat$test="Fe deprivation at Stat6"

    hits_stat_withiron=stat_withiron[which(stat_withiron$BH<threshold),]
    hits_stat_withiron$test="Growth arrest at +Fe"

    hits_stat_withoutiron=stat_withoutiron[which(stat_withoutiron$BH<threshold),]
    hits_stat_withoutiron$test="Growth arrest at -Fe"


    table=data.frame(rbind(hits_iron_atstat,hits_stat_withiron,hits_stat_withoutiron))

    colors_with=c("#A3C9F5","#1461B8")

    colors_without=c("mistyrose","red","firebrick4")
    colors_without=c("#FEE2CD","#FFA866","#FF7A14")
    
    colores=c("#B8141F","#A3C9F5","#FF7A14")
    #geom_segment(aes(x=x,xend=x,y=0,yend=logFC),color="#77B12F",size=1,data=mbt_exp)+
    #geom_segment(aes(x=x,xend=x,y=0,yend=logFC),color="#B8141F",size=1,data=mbt_stat)+

    
    
    pl=ggplot(table)+
    geom_density(aes(x=abs(log2FoldChange_shrunken),fill=test),alpha=0.7)+
    scale_fill_manual(values=colores)+
    xlim(0,6)+
    theme(legend.position=c(0.7,0.7))+theme_classic()+geom_vline(xintercept=0)+theme(
    axis.text.x   = element_text(size=14),
    axis.text.y   = element_text(size=14),
    axis.title.x  = element_text(size=14),
    axis.title.y  = element_text(size=14),
    #axis.ticks.y = element_blank(),
    legend.background = element_blank(),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1,linetype="solid"),
    legend.title=element_blank(),
    legend.position=c(0.7,0.8),
    #legend.text=element_text(size=14),
    legend.key.size = unit(1, 'lines'))+geom_hline(yintercept=0)+xlab("|log2(FC)|")+ylab("Density")

    pdf("outputs/Figures/Unpublished/abs_densities.pdf",width=4.5,height=2.5)
    print(pl)
    dev.off()

    test=wilcox.test(abs(hits_stat_withoutiron$log2FoldChange_shrunken),abs(hits_iron_atstat$log2FoldChange_shrunken), alternative = "greater")
    test$p.value
    # 4.329263e-80
    test=wilcox.test(abs(hits_stat_withoutiron$log2FoldChange_shrunken),abs(hits_stat_withiron$log2FoldChange_shrunken), alternative = "greater")
    test$p.value
    # 3.644891e-27
    test=wilcox.test(abs(hits_stat_withiron$log2FoldChange_shrunken),abs(hits_iron_atstat$log2FoldChange_shrunken), alternative = "greater")
    test$p.value
    # 4.306936e-37
}

## JQ: writing the normalized counts tables for plotting the heatmaps.

write.table(normalized_counts_EU,"outputs/Data/txt/normalized_counts_EU.txt")
write.table(normalized_counts_ED,"outputs/Data/txt/normalized_counts_ED.txt")



### Check specific genes within resumen

resumen=resumen[order(rownames(resumen)),]
feature_data=feature_data[order(rownames(feature_data)),]
length(which(rownames(resumen)!=rownames(feature_data)))
resumen$Gene_name=feature_data$Gene_Symbol
resumen$RV_number=rownames(resumen)

check=resumen[,c(10,9,7,8)]
# Busco UP_MORE_WOUT_FE y DOWN_MORE_WOUT_FE

genes_black=c("Rv0476","Rv0211","Rv2221","Rv2443","Rv3436c","Rv0132c","Rv3842c")
genes_blue=c("Rv1127","Rv2241","Rv1092","Rv0896","Rv1475","Rv3339c","Rv3432c","Rv1731","Rv1731","Rv1098c","Rv1837","Rv1130","Rv1131")

check[genes_black,]
check[genes_blue,]
## Rv0211 is a typo

