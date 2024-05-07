########################################
########################################
###### D1: Human immune-pop Quach ######
########################################
########################################

############################
### 0. Load dependencies ###
############################

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

    options(width=1000)
}

##########################################################################################
## 1. Fetch raw data from Quach et al. and build clean and messed up deliverable tables ##
##########################################################################################


###################################################
## 9. Differential expression. 1: only treatment ##
###################################################

{
### Pongamos que tienes solo una covariable (condition), una tabla de metadata 
  
  # y una matriz de reads, entonces:

### voom+limma DE

{

    dge <- DGEList(counts = reads)
    dge <- calcNormFactors(dge)
    design=model.matrix(~condition,data=metadata)
    v=voom(dge,design,plot=TRUE)

    fit <-lmFit(v,design)
    fit <- eBayes(fit)

    betas=data.frame(fit$coefficients)
    ts=data.frame(fit$t)
    ps=data.frame(fit$p.value)
    BH=data.frame(ps)
    errors <- sqrt(fit$s2.post) * fit$stdev.unscaled

    th=0.05
    for(i in 1:ncol(BF))
    {
        BH[,i]=p.adjust(ps[,i],method="BH")
        print(paste("i=",i,length(which(BH[,i]<th))))
    }
    
    ### What about comparisons between other pairs.

    ## You can cahnge the reference level of the Condition factor,
    ## (DEVELOP)
    ## Or you can use the contrast.fit function
    ## (DEVELOP)
    
    ## We can add batch effects.
}

### DESeq2

{

    dds <- DESeqDataSetFromMatrix(countData = round(reads),colData = metadata,design= ~ condition)
    dds <- DESeq(dds)
    resultsNames(dds) # lists the coefficients
    res1 <- results(dds, name="condition_Cond1_vs_Control")
    res2 <- results(dds, name="condition_Cond2_vs_Control")
    
    ## This to make contrasts out of differences of coeficients
    lista=list(add=c("condition_Cond1_vs_Control"),remove=c("condition_Cond2_vs_Control"))
    res1_B=results(dds, contrast=lista)
    
    
    ## LogFC shrinkage (it does not affect the p values, or beta-standard errors)
    res2_B <- lfcShrink(dds, contrast=lista, type="ashr")
    res3_B <- lfcShrink(dds, coef="condition_PAM3CSK4_vs_NS", type="ashr")

    
    ## Some random exam
}

### EDGER

{
    # Para el trto del offset y demás (pero lo de abajo debería estar simplemente bien)(
    # https://rdrr.io/github/mikelove/tximeta/src/R/dgelist.R
        
    y <- DGEList(reads)
    y <- calcNormFactors(y)
    design <- model.matrix(~condition,data=metadata)
    y <- estimateDisp(y,design)
    
    #To perform quasi-likelihood F-tests:
    fit <- glmQLFit(y,design)
    
    ## Ojo: toptags da log2FCs; fit$coeffs; natural logs!! https://support.bioconductor.org/p/50273/
    
    betas_edger=log2(exp(fit$coefficients))
    
    
    # filtering using the design information
    
}

}
