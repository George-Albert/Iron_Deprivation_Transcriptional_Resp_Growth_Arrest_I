# michael's version
# https://support.bioconductor.org/p/91218/

tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

options(width=10000)
{
    ul=function(tab,n=5)
    {
        rs=min(nrow(tab),n)
        cs=min(ncol(tab),n)
        return(tab[1:rs,1:cs])
    }
    ur=function(tab,n=5)
    {
        rs=min(nrow(tab),n)
        cs=min(ncol(tab),n)
        return(tab[1:rs,(ncol(tab)-cs):ncol(tab)])
    }
    dcols=function(tab){data.frame(colnames(tab))}
    drows=function(tab){data.frame(rownames(tab))}
}

    library(edgeR)
    library(qvalue)
    library(limma)
    library(qvalue)
    library(ggplot2)
    library(cowplot)
    library(reshape2)
    library(ashr)
    library(DESeq2)
    library(fdrtool)

    
    metadata=read.table("inputs/processed/metadata.txt")
    feature_data=read.table("inputs/processed/feature_data.txt")
    reads=read.table("inputs/processed/reads.txt")
    fpkms=read.table("inputs/processed/fpkm.txt")
    
    metadata$short_setup=paste0(metadata$Culture,"_Fe_",metadata$Iron)

###
### Jq: Subset muestras de interés-
###



    medias_within_group=fpkms[,c(1:8)]
    colnames(medias_within_group)=unique(metadata$short_setup)
    min_reads_detected=medias_within_group
    for(i in 1:ncol(medias_within_group))
    {
        chunk=fpkms[,which(metadata$short_setup %in% colnames(medias_within_group)[i])]
        chunk_reads=reads[,which(metadata$short_setup %in% colnames(medias_within_group)[i])]
        medias_within_group[,i]=apply(chunk,1,mean)
        min_reads_detected[,i]=apply(chunk_reads,1,min)
    }
    
    min_reads_detected$grand_min=apply(min_reads_detected,1,min)
    medias_within_group$grand_min=apply(medias_within_group,1,min)

    min_reads_detected$grand_max=apply(min_reads_detected,1,max)
    medias_within_group$grand_max=apply(medias_within_group,1,max)
    
    keep=which(medias_within_group$grand_max>10)
    
    reads=reads[keep,]
    min_reads_detected=min_reads_detected[keep,]
    
    ## Here, with a criterium of being above 10 fpkms in at least one setup, I only trash 7% of genes, and the genes kept are above a mean # of reads=9 in at least one of the setups.
    
    metadata$Iron=factor(metadata$Iron,levels=c("YES","NO"))
    metadata$short_setup=factor(metadata$short_setup)
    metadata$short_setup=relevel(metadata$short_setup,ref="C5_Fe_YES")
 
##
##  Modelización usando voom+limma.
##
    
    dge <- DGEList(counts = reads)
    ## TMM (trimmed M-means) descrito en el papel de edger
    dge <- calcNormFactors(dge)
    design=model.matrix(~short_setup,data=metadata)
    v=voom(dge,design,plot=TRUE)

    fit <-lmFit(v,design)
    fit <- eBayes(fit)
    
    ## Dextrose effects (Yes minus NO)
    
    ## @ exp phase, with iron: C1_Fe_YES - C5_Fe_YES (reference)
    c1=c(3)
    ## @ exp phase, without iron: C1_Fe_NO - C5_Fe_NO
    c2=c(2,-6)
    ## @ stat phase, with iron: C2_Fe_YES - C6_Fe_Yes
    c3=c(5,-8)
    ## @ stat phase, without iron: C2_Fe_NO - C6_Fe_NO
    c4=c(4,-7)
    
    
    ## Growth Phase effects (Stat minus exp)
    
    ## Without desxtrose, with iron: C6_Fe_YES - C5_Fe_YES (reference)
    c5=c(8)
    ## Without desxtrose, without iron: C6_Fe_NO - C5_Fe_NO
    c6=c(7,-6)
    ## With desxtrose, with iron: C2_Fe_YES - C1_Fe_Yes
    c7=c(5,-3)
    ## With desxtrose, without iron: C2_Fe_NO - C1_Fe_NO
    c8=c(4,-2)

    ## Iron deprovation effects (No minus Yes)
    
    ## @ exp phase, without dextrose: C5_Fe_NO - C5_Fe_YES (reference)
    c9=c(6)
    ## @ exp phase, with dextrose: C1_Fe_NO - C1_Fe_YES
    c10=c(2,-3)
    ## @ stat phase, without dextrose: C6_Fe_NO - C6_Fe_YES
    c11=c(7,-8)
    ## @ stat phase, with dextrose: C2_Fe_NO - C2_Fe_YES
    c12=c(4,-5)
    
    ### First, I build effect sizes, normalized t-statistics, pvalues and fdr matrices from contrasts1, 5 and 9 (accessible without contrast.fit.
    
    
    

    betas=data.frame(fit$coefficients)
    ts=data.frame(fit$t)
    ps=data.frame(fit$p.value)
    fdrs=data.frame(ps)
    ash_fdrs=data.frame(ps)
    ash_betas=data.frame(betas)

    errors <- sqrt(fit$s2.post) * fit$stdev.unscaled
    dofs=data.frame(base=(fit$df.residual+fit$df.prior))

    for(i in 1:ncol(fdrs))
    {
        fdrs[,i]=p.adjust(ps[,i],method="BH")
        ash_object=ash(betas[,i],errors[,i],df=dofs[1,1])
        ash_fdrs[,i]=ash_object$result$lfsr
        ash_betas[,i]=ash_object$result$PosteriorMean
        print(paste(colnames(fdrs)[i],length(which(fdrs[,i]<0.1 & abs(betas[,i])>0)),length(which(fdrs[,i]<0.05 & abs(betas[,i])>0))))
    }
    
    betas=betas[,c(3,8,6)]
    ts=ts[,c(3,8,6)]
    ps=ps[,c(3,8,6)]
    fdrs=fdrs[,c(3,8,6)]
    ash_fdrs=ash_fdrs[,c(3,8,6)]
    ash_betas=ash_betas[,c(3,8,6)]

    errors=errors[,c(3,8,6)]

    colnames(betas)=c("c1","c5","c9")
    colnames(ts)=c("c1","c5","c9")
    colnames(ps)=c("c1","c5","c9")
    colnames(fdrs)=c("c1","c5","c9")
    colnames(ash_fdrs)=c("c1","c5","c9")
    colnames(ash_betas)=c("c1","c5","c9")

    colnames(errors)=c("c1","c5","c9")

    
    DE_results=list(betas=betas,ash_betas=ash_betas,ts=ts,ps=ps,fdrs=fdrs,ash_fdrs=ash_fdrs,errors=errors,dofs=dofs)


    
    process_contrast=function(DE_object,contrast){

    name <-deparse(substitute(contrast))
    vec=rep(0,length(colnames(design)))
    vec[abs(contrast)]=contrast/abs(contrast)
    fit2 <- contrasts.fit(fit, vec)
    fit2 <- eBayes(fit2)
    
    b=data.frame(fit2$coefficients)
    t=data.frame(fit2$t)
    p=data.frame(fit2$p.value)
    f=data.frame(p)
    af=data.frame(p)
    ab=data.frame(b)

    e <- sqrt(fit2$s2.post) * fit2$stdev.unscaled
    d=fit2$df.residual+fit2$df.prior
    
    print(d[1])

    for(i in 1:ncol(f))
    {
        f[,i]=p.adjust(p[,i],method="BH")
        ash_object_local=ash(b[,i],e[,i],df=d[1])
        af[,i]=ash_object_local$result$lfsr
        ab[,i]=ash_object_local$result$PosteriorMean
        print(paste(colnames(f)[i],length(which(f[,i]<0.1 & abs(b[,i])>0)),length(which(f[,i]<0.05 & abs(b[,i])>0))))
    }
    
    betas_ext=cbind(DE_object$betas,new=b)
    ts_ext=cbind(DE_object$ts,new=t)
    ps_ext=cbind(DE_object$ps,new=p)
    fdrs_ext=cbind(DE_object$fdrs,new=f)
    ash_fdrs_ext=cbind(DE_object$ash_fdrs,new=af)
    ash_betas_ext=cbind(DE_object$ash_betas,new=ab)

    errors_ext=cbind(DE_object$errors,new=e)
    dofs_ext=cbind(DE_object$dofs,new=d)

    colnames(betas_ext)[ncol(betas_ext)]=name
    colnames(ts_ext)[ncol(ts_ext)]=name
    colnames(ps_ext)[ncol(ps_ext)]=name
    colnames(fdrs_ext)[ncol(fdrs_ext)]=name
    colnames(ash_fdrs_ext)[ncol(ash_fdrs_ext)]=name
    colnames(ash_betas_ext)[ncol(ash_betas_ext)]=name

    colnames(errors_ext)[ncol(errors_ext)]=name
    colnames(dofs_ext)[ncol(dofs_ext)]=name

    DE_output=list(betas=betas_ext,ash_betas=ash_betas_ext,ts=ts_ext,ps=ps_ext,fdrs=fdrs_ext,ash_fdrs=ash_fdrs_ext,errors=errors_ext,dofs=dofs_ext)
    
    return(DE_output)
    }

DE_results=process_contrast(DE_object=DE_results,contrast=c2)
DE_results=process_contrast(DE_object=DE_results,contrast=c3)
DE_results=process_contrast(DE_object=DE_results,contrast=c4)

DE_results=process_contrast(DE_object=DE_results,contrast=c6)
DE_results=process_contrast(DE_object=DE_results,contrast=c7)
DE_results=process_contrast(DE_object=DE_results,contrast=c8)

DE_results=process_contrast(DE_object=DE_results,contrast=c10)
DE_results=process_contrast(DE_object=DE_results,contrast=c11)
DE_results=process_contrast(DE_object=DE_results,contrast=c12)

DE_results$betas=DE_results$betas[,order( as.numeric(substr(colnames(DE_results$betas),2,nchar(colnames(DE_results$betas)))))]
DE_results$ts=DE_results$ts[,order( as.numeric(substr(colnames(DE_results$ts),2,nchar(colnames(DE_results$ts)))))]
DE_results$ps=DE_results$ps[,order( as.numeric(substr(colnames(DE_results$ps),2,nchar(colnames(DE_results$ps)))))]
DE_results$fdrs=DE_results$fdrs[,order( as.numeric(substr(colnames(DE_results$fdrs),2,nchar(colnames(DE_results$fdrs)))))]
DE_results$ash_fdrs=DE_results$ash_fdrs[,order( as.numeric(substr(colnames(DE_results$ash_fdrs),2,nchar(colnames(DE_results$ash_fdrs)))))]
DE_results$ash_betas=DE_results$ash_betas[,order( as.numeric(substr(colnames(DE_results$ash_betas),2,nchar(colnames(DE_results$ash_betas)))))]

DE_results$errors=DE_results$errors[,order( as.numeric(substr(colnames(DE_results$errors),2,nchar(colnames(DE_results$errors)))))]
## The degrees of freedom are just a number for all contrasts and genes...





get_signal_figures=function(x,y,xalt,yalt,z,index,th_x,th_y,path,nombre_file){
    
    datos=data.frame(logFC=x[,index],logFDR=y[,index],ash_logFC=xalt[,index],log_ash_lFDR=yalt[,index],p=z[,index])
    datos$Association="N.S."
    datos$Association[which(abs(datos$logFC)>th_x & datos$logFDR>th_y)]="DEG"
    datos$Association=factor(datos$Association,levels=c("N.S.","DEG"))

    
    datos$Association_ash="N.S."
    datos$Association_ash[which(abs(datos$ash_logFC)>th_x & datos$log_ash_lFDR>th_y)]="DEG"
    datos$Association_ash=factor(datos$Association_ash,levels=c("N.S.","DEG"))

    
    pl_volcano=ggplot(datos)+geom_point(aes(x=logFC,y=logFDR,color=Association))+ylab("-log10(FDR)")
    pl_volcano_ash=ggplot(datos)+geom_point(aes(x=ash_logFC,y=log_ash_lFDR,color=Association_ash))+ylab("-log10(ash_loc_FDR)")

    pl_histogram=ggplot(datos)+geom_histogram(aes(x=p),bins=50)
    pl_signal=plot_grid(pl_histogram,pl_volcano,pl_volcano_ash,ncol=3,align="h")
    
    system(paste0("mkdir -p ",path))
    
    pdf(paste0(path,"/",nombre_file,".pdf"),width=12,height=4)
    print(pl_signal)
    dev.off()
    
}




writer_set=function(set,path,nombre_file){
    
    sink(paste0(path,nombre_file,".txt"))
    for(gen in set)
    cat(gen,"\n")
    sink()
}

    
get_signal_barplot=function(hits,path,nombre_file){
    hitss=data.frame(t(hits))
    hitss$down=-hitss$down
    hitss=t(hitss)
    hits_melt=melt(as.matrix(hitss[c(1:2),]),id=0)
    pl=ggplot(hits_melt, aes(fill=Var1, y=value, x=Var2)) +     geom_bar(position="stack", stat="identity")+coord_flip()
    
    system(paste0("mkdir -p ",path))
    
    pdf(paste0(path,"/",nombre_file,".pdf"),width=5,height=5)
    print(pl)
    dev.off()
}


signal_portrait=function(DE_object,th_fdr,th_beta,name){

hits=DE_results$betas[1:3,]
rownames(hits)=c("up","down","all")
ash_hits=hits

betas=DE_object$betas
ts=DE_object$ts
ps=DE_object$ps
fdrs=DE_object$fdrs
ash_fdrs=DE_object$ash_fdrs
ash_betas=DE_object$ash_betas

errors=DE_object$errors

for(i in 1:12)
{
    
    upreg_genes=rownames(fdrs)[which(fdrs[,i]<th_fdr & betas[,i]>th_beta)]
    downreg_genes=rownames(fdrs)[which(fdrs[,i]<th_fdr & betas[,i]<(-th_beta))]
    
    ash_upreg_genes=rownames(ash_fdrs)[which(ash_fdrs[,i]<th_fdr & ash_betas[,i]>th_beta)]
    ash_downreg_genes=rownames(ash_fdrs)[which(ash_fdrs[,i]<th_fdr & ash_betas[,i]<(-th_beta))]
    
    hits[,i]=c(length(upreg_genes),length(downreg_genes),length(upreg_genes)+length(downreg_genes))
    ash_hits[,i]=c(length(ash_upreg_genes),length(ash_downreg_genes),length(ash_upreg_genes)+length(ash_downreg_genes))

    lista_i=list(up=upreg_genes,down=downreg_genes,all=c(upreg_genes,downreg_genes))
    ash_lista_i=list(up=ash_upreg_genes,down=ash_downreg_genes,all=c(ash_upreg_genes,ash_downreg_genes))

    if(i==1)
    {
        targets=list(contrast=lista_i)
        names(targets)[1]="c1"
        
        ash_targets=list(contrast=ash_lista_i)
        names(ash_targets)[1]="c1"

    }else{
        targets[[i]]=lista_i
        names(targets)[i]=paste0("c",i)
        
        ash_targets[[i]]=ash_lista_i
        names(ash_targets)[i]=paste0("c",i)

    }
}

for(i in 1:length(targets))
{
    system(paste0("mkdir -p outputs/DE_portrait/",name,"/target_sets/BH/",names(targets)[i]))
    system(paste0("mkdir -p outputs/DE_portrait/",name,"/target_sets/ash/",names(targets)[i]))

    writer_set(targets[[i]]$up,path=paste0("outputs/DE_portrait/",name,"/target_sets/BH/",names(targets)[i]),nombre_file="/up")
    writer_set(targets[[i]]$down,path=paste0("outputs/DE_portrait/",name,"/target_sets/BH/",names(targets)[i]),nombre_file="/down")
    writer_set(targets[[i]]$all,path=paste0("outputs/DE_portrait/",name,"/target_sets/BH/",names(targets)[i]),nombre_file="/all")
    
    writer_set(ash_targets[[i]]$up,path=paste0("outputs/DE_portrait/",name,"/target_sets/ash/",names(targets)[i]),nombre_file="/up")
    writer_set(ash_targets[[i]]$down,path=paste0("outputs/DE_portrait/",name,"/target_sets/ash/",names(targets)[i]),nombre_file="/down")
    writer_set(ash_targets[[i]]$all,path=paste0("outputs/DE_portrait/",name,"/target_sets/ash/",names(targets)[i]),nombre_file="/all")

    get_signal_figures(x=betas,y=-log10(fdrs),xalt=ash_betas,yalt=-log10(ash_fdrs),z=ps,i,th_x=th_beta,th_y=-log10(th_fdr),path=paste0("outputs/DE_portrait/",name,"/signal_plots/"),nombre_file=names(targets)[i])
    
}

get_signal_barplot(hits=hits,path=paste0("outputs/DE_portrait/",name,"/barplots/"),nombre_file="BH")
get_signal_barplot(hits=ash_hits,path=paste0("outputs/DE_portrait/",name,"/barplots/"),nombre_file="ash")

output=list(targets_BH=targets,ash_targets=ash_targets,hits=hits,ash_hits=ash_hits)

return(output)
}

conservative=signal_portrait(DE_object=DE_results,th_fdr=0.01,th_beta=1,name="conservative")
standard=signal_portrait(DE_object=DE_results,th_fdr=0.05,th_beta=0,name="standard")
liberal=signal_portrait(DE_object=DE_results,th_fdr=0.1,th_beta=0,name="liberal")

###
### PCA to check what is going on with the contrast 9 and 10.
###

exp=v$E

pca <-prcomp(t(exp),scale=T)
sum_pca=data.frame(summary(pca)$importance[,c(1:15)])

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))
datos=cbind(metadata,datos[,1:6])

datos$Culture=factor(datos$Culture,levels=c("C5","C1","C6","C2"))
datos$short_setup=factor(datos$short_setup,levels=c("C5_Fe_YES","C5_Fe_NO","C1_Fe_YES","C1_Fe_NO","C6_Fe_YES","C6_Fe_NO","C2_Fe_YES","C2_Fe_NO"))
datos=datos[order(datos$short_setup),]
colores_base=c("red","forestgreen","dodgerblue","black")
fill_base=c("red",NA,"forestgreen",NA,"dodgerblue",NA,"black",NA)

pl=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Culture,fill=short_setup),shape=21,size=3,stroke=1.5)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab("PC1: 50.3% variance ~growth phase")+ylab("PC1: 14.8% variance ~Dextrose, Iron (at stat phase only)")

system("mkdir -p outputs/summary_figures")
pdf("outputs/summary_figures/pca.pdf",width=6,height=4.5)
print(pl)
dev.off()

pl=ggplot(datos)+geom_point(aes(x=PC1,y=PC3,color=Culture,fill=short_setup),shape=21,size=3,stroke=1.5)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab("PC1: 50.3% variance ~growth phase")+ylab("PC1: 14.8% variance ~Dextrose, Iron (at stat phase only)")




system("mkdir -p outputs/summary_figures")
pdf("outputs/summary_figures/pca_1_3.pdf",width=6,height=4.5)
print(pl)
dev.off()


pl=ggplot(datos)+geom_point(aes(x=PC2,y=PC3,color=Culture,fill=short_setup),shape=21,size=3,stroke=1.5)+scale_colour_manual(values=colores_base)+scale_fill_manual(values=fill_base)+xlab("PC1: 50.3% variance ~growth phase")+ylab("PC1: 14.8% variance ~Dextrose, Iron (at stat phase only)")




system("mkdir -p outputs/summary_figures")
pdf("outputs/summary_figures/pca_1_3.pdf",width=6,height=4.5)
print(pl)
dev.off()

##
## betas distributions
##

betas_melt=melt(as.matrix(DE_results$betas),id=0)

pl=ggplot(betas_melt)+geom_density(aes(x=value,fill=Var2,color=Var2),alpha=0)+ylim(0,4)
pdf("outputs/summary_figures/beta_densities.pdf",width=10,height=6)
print(pl)
dev.off()


betas_melt=melt(as.matrix(DE_results$ash_betas),id=0)

pl=ggplot(betas_melt)+geom_density(aes(x=value,fill=Var2,color=Var2),alpha=0)+ylim(0,4)
pdf("outputs/summary_figures/ash_beta_densities.pdf",width=10,height=6)
print(pl)
dev.off()

###
### Beta distributions.
###

chunk=reads[,c(9,10,1,2)]
meta_chunk=metadata[which(rownames(metadata)%in%colnames(chunk)),]

chunk=chunk[,order(colnames(chunk))]
meta_chunk=meta_chunk[order(rownames(meta_chunk)),]



dds <- DESeqDataSetFromMatrix(countData = chunk,
                              colData = meta_chunk,
                              design= ~ Iron)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- data.frame(results(dds, name="Iron_NO_vs_YES"))

testo <- fdrtool(res$stat, statistic= "normal", plot = T)
pvals_corrected=data.frame(p=testo$pval,q=testo$qval,lfdr=testo$lfdr)

pl=ggplot(pvals_corrected)+geom_histogram(aes(x=p),bins=50)
length(which(testo$lfdr<0.01 & abs(res$log2FoldChange)>1))

##
chunk=reads[,c(1:4)]
meta_chunk=metadata[which(rownames(metadata)%in%colnames(chunk)),]

chunk=chunk[,order(colnames(chunk))]
meta_chunk=meta_chunk[order(rownames(meta_chunk)),]



dds <- DESeqDataSetFromMatrix(countData = chunk,
                              colData = meta_chunk,
                              design= ~ Growth)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- data.frame(results(dds, name="Growth_STAT_vs_EXP"))

testo <- fdrtool(res$stat, statistic= "normal", plot = T)
pvals_corrected=data.frame(p=testo$pval,q=testo$qval,lfdr=testo$lfdr)

pl=ggplot(pvals_corrected)+geom_histogram(aes(x=p),bins=50)
length(which(testo$q<0.01 & abs(res$log2FoldChange)>1))




names(DE_results)




    chunk=reads[,c(9,10,1,2)]
    meta_chunk=metadata[which(rownames(metadata)%in%colnames(chunk)),]

    chunk=chunk[,order(colnames(chunk))]
    meta_chunk=meta_chunk[order(rownames(meta_chunk)),]



    dds <- DESeqDataSetFromMatrix(countData = chunk,
                                  colData = meta_chunk,
                                  design= ~ Iron)
    dds <- DESeq(dds)
    resultsNames(dds) # lists the coefficients
    res <- data.frame(results(dds, name="Iron_NO_vs_YES"))



    testo <- fdrtool(res$stat, statistic= "normal", plot = T)
    good_bh=p.adjust(testo$pval, method = "BH")
    length(which(good_bh<0.01))

    for(i in 1:3){
    testo <- fdrtool((ts[,i]-median(ts[,i])), statistic= "normal", plot = T)
    good_bh=p.adjust(testo$pval, method = "BH")
    print(paste(i,length(which(good_bh<0.05 & abs(betas[,i])>0.5))))
    }

    testo <- fdrtool(qnorm(ps[,3]), statistic= "normal", plot = T)
    good_bh=p.adjust(testo$pval, method = "BH")
    length(which(good_bh<0.01))
    
    matt=ash(betahat=betas[,3],sebetahat=errors[,3],df=dofs[1,1])

    for(i in 1:3)
    {
    testo <- locfdr(ts[,1])
    good_bh=p.adjust(testo$pval, method = "BH")
    print(paste(i,length(which(good_bh<0.05 & abs(betas[,i])>0.5))))
    }



    things <- fdrtool(Zs, statistic= "normal", plot = T)


    things$param[1, "sd"]

    # or to shrink log fold changes association with condition:
    res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")



    dge_chunk <- DGEList(counts = chunk)
    dge_chunk <- calcNormFactors(dge_chunk)
    design_chunk=model.matrix(~Iron,data=meta_chunk)
    v_chunk=voom(dge_chunk,design_chunk,plot=TRUE)

    fit_chunk <-lmFit(v_chunk,design_chunk)
    fit_chunk <- eBayes(fit_chunk)

    betas_chunk=data.frame(fit_chunk$coefficients)
    ts_chunk=data.frame(fit_chunk$t)
    ps_chunk=data.frame(fit_chunk$p.value)
    fdrs_chunk=data.frame(ps_chunk)

    for(i in 1:ncol(fdrs_chunk))
    {
        fdrs_chunk[,i]=p.adjust(ps_chunk[,i],method="BH")
        print(paste(colnames(fdrs_chunk)[i],length(which(fdrs_chunk[,i]<0.1 & abs(betas_chunk[,i])>0)),length(which(fdrs_chunk[,i]<0.05 & abs(betas_chunk[,i])>0))))
    }

    pl_chunk=ggplot(ps_chunk)+geom_histogram(aes(x=IronNO),bins=50)

