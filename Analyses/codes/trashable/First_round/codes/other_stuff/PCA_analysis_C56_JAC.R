fcShrink
function (dds, coef, contrast, res, type = c("apeglm", "ashr",
    "normal"), lfcThreshold = 0, svalue = FALSE, returnList = FALSE,
    format = c("DataFrame", "GRanges", "GRangesList"), apeAdapt = TRUE,
    apeMethod = "nbinomCR", parallel = FALSE, BPPARAM = bpparam(),
    quiet = FALSE, ...)
{
    stopifnot(is(dds, "DESeqDataSet"))
    if (!missing(res)) {
        if (!is(res, "DESeqResults"))
            stop("res should be a DESeqResults object, for GRanges output use 'format'")
    }
    type <- match.arg(type, choices = c("apeglm", "ashr", "normal"))
    format <- match.arg(format, choices = c("DataFrame", "GRanges","GRangesList"))
    if (length(resultsNames(dds)) == 0) {
        if (type != "apeglm" | (type == "apeglm" & apeAdapt)) {
            stop("first run DESeq() before running lfcShrink()")
        }
    }
    betaPrior <- attr(dds, "betaPrior")
    if (!is.null(betaPrior) && betaPrior) {
        stop("lfcShrink() should be used downstream of DESeq() with betaPrior=FALSE (the default)")
    }
    stopifnot(length(lfcThreshold) == 1 && lfcThreshold >= 0)
    resultsNamesDDS <- resultsNames(dds)
    if (type == "apeglm" & !apeAdapt & length(resultsNames(dds)) ==
        0) {
        resultsNamesDDS <- colnames(model.matrix(design(dds),
            data = colData(dds)))
    }
    if (!missing(coef)) {
        if (is.numeric(coef)) {
            stopifnot(coef <= length(resultsNamesDDS))
            coefAlpha <- resultsNamesDDS[coef]
            coefNum <- coef
        }
        else if (is.character(coef)) {
            stopifnot(coef %in% resultsNamesDDS)
            coefNum <- which(resultsNamesDDS == coef)
            coefAlpha <- coef
        }
    }
    if (missing(res)) {
        if (!missing(coef)) {
            res <- results(dds, name = coefAlpha)
        }
        else if (!missing(contrast)) {
            if (type == "normal" & is.numeric(contrast)) {
                stop("for type='normal' and numeric contrast, user must provide 'res' object")
            }
            res <- results(dds, contrast = contrast)
        }
        else {
            stop("one of coef or contrast required if 'res' is missing")
        }
    }
    if (type %in% c("normal", "apeglm")) {
        if (is.null(dispersions(dds))) {
            stop("type='normal' and 'apeglm' require dispersion estimates, first call estimateDispersions()")
        }
        stopifnot(all(rownames(dds) == rownames(res)))
        if (parallel) {
            nworkers <- getNworkers(BPPARAM)
            parallelIdx <- factor(sort(rep(seq_len(nworkers),
                length.out = nrow(dds))))
        }
    }
    if (type == "normal") {
        if (missing(coef) & missing(contrast)) {
            stop("type='normal' requires either 'coef' or 'contrast' be specified.")
        }
        if (!quiet)
            message("using 'normal' for LFC shrinkage, the Normal prior from Love et al (2014).\n\nNote that type='apeglm' and type='ashr' have shown to have less bias than type='normal'.\nSee ?lfcShrink for more details on shrinkage type, and the DESeq2 vignette.\nReference: https://doi.org/10.1093/bioinformatics/bty895")
        if (is(design(dds), "formula")) {
            if (attr(dds, "modelMatrixType") == "user-supplied") {
                if (!missing(contrast)) {
                  stop("user-supplied design matrix supports shrinkage only with 'coef'")
                }
                modelMatrix <- attr(dds, "modelMatrix")
            }
            else {
                termsOrder <- attr(terms.formula(design(dds)),
                  "order")
                interactionPresent <- any(termsOrder > 1)
                if (interactionPresent) {
                  stop("LFC shrinkage type='normal' not implemented for designs with interactions")
                }
                modelMatrix <- NULL
            }
        }
        else if (is(design(dds), "matrix")) {
            if (!missing(contrast)) {
                stop("user-supplied design matrix supports shrinkage only with 'coef'")
            }
            modelMatrix <- design(dds)
        }
        stopifnot(missing(coef) | missing(contrast))
        betaCols <- grep("log2 fold change \\(MLE\\)", mcols(mcols(dds))$description)
        stopifnot(length(betaCols) > 0)
        if (!any(grepl("MLE_", names(mcols(dds))[betaCols]))) {
            names(mcols(dds))[betaCols] <- paste0("MLE_", names(mcols(dds))[betaCols])
        }
        if (missing(contrast)) {
            modelMatrixType <- "standard"
        }
        else {
            modelMatrixType <- "expanded"
            expMM <- makeExpandedModelMatrix(dds)
            resNames <- colnames(expMM)
            if (is(contrast, "character")) {
                stopifnot(length(contrast) == 3)
                stopifnot(contrast[1] %in% names(colData(dds)))
                stopifnot(is(colData(dds)[[contrast[1]]], "factor"))
                stopifnot(all(contrast[2:3] %in% levels(colData(dds)[[contrast[1]]])))
            }
            else {
                message("resultsNames from the expanded model matrix to be referenced by 'contrast':")
                message(paste0("'", paste(resNames, collapse = "', '"),
                  "'"))
            }
            contrast <- checkContrast(contrast, resNames)
        }
        attr(dds, "modelMatrixType") <- modelMatrixType
        betaPriorVar <- estimateBetaPriorVar(dds, modelMatrix = modelMatrix)
        stopifnot(length(betaPriorVar) > 0)
        if (!parallel) {
            dds.shr <- nbinomWaldTest(dds, betaPrior = TRUE,
                betaPriorVar = betaPriorVar, modelMatrix = modelMatrix,
                modelMatrixType = modelMatrixType, quiet = TRUE)
        }
        else {
            dds.shr <- do.call(rbind, bplapply(levels(parallelIdx),
                function(l) {
                  nbinomWaldTest(dds[parallelIdx == l, , drop = FALSE],
                    betaPrior = TRUE, betaPriorVar = betaPriorVar,
                    modelMatrix = modelMatrix, modelMatrixType = modelMatrixType,
                    quiet = TRUE)
                }, BPPARAM = BPPARAM))
        }
        if (missing(contrast)) {
            res.shr <- results(dds.shr, name = coefAlpha, lfcThreshold = lfcThreshold)
        }
        else {
            res.shr <- results(dds.shr, contrast = contrast,
                lfcThreshold = lfcThreshold, parallel = parallel,
                BPPARAM = BPPARAM)
        }
        if (lfcThreshold > 0) {
            change.cols <- c("log2FoldChange", "lfcSE", "stat",
                "pvalue", "padj")
        }
        else {
            change.cols <- c("log2FoldChange", "lfcSE")
        }
        for (column in change.cols) {
            res[[column]] <- res.shr[[column]]
        }
        mcols(res, use.names = TRUE)[change.cols, "description"] <- mcols(res.shr,
            use.names = TRUE)[change.cols, "description"]
        deseq2.version <- packageVersion("DESeq2")
        metadata(res)[["lfcThreshold"]] <- lfcThreshold
        priorInfo(res) <- list(type = "normal", package = "DESeq2",
            version = deseq2.version, betaPriorVar = betaPriorVar)
        res <- resultsFormatSwitch(object = dds, res = res, format = format)
        return(res)
    }
    else if (type == "apeglm") {
        if (!requireNamespace("apeglm", quietly = TRUE)) {
            stop("type='apeglm' requires installing the Bioconductor package 'apeglm'")
        }
        if (!missing(contrast)) {
            stop("type='apeglm' shrinkage only for use with 'coef'")
        }
        stopifnot(!missing(coef))
        if (apeAdapt) {
            incomingCoef <- gsub(" ", "_", sub("log2 fold change \\(MLE\\): ",
                "", mcols(res)[2, 2]))
            if (coefAlpha != incomingCoef) {
                stop("'coef' should specify same coefficient as in results 'res'")
            }
        }
        if (!quiet)
            message("using 'apeglm' for LFC shrinkage. If used in published research, please cite:\n    Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for\n    sequence count data: removing the noise and preserving large differences.\n    Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895")
        Y <- counts(dds)
        modelMatrixType <- attr(dds, "modelMatrixType")
        if (!is.null(modelMatrixType) && modelMatrixType == "user-supplied") {
            design <- attr(dds, "modelMatrix")
        }
        else {
            design <- model.matrix(design(dds), data = colData(dds))
        }
        disps <- dispersions(dds)
        if (is.null(normalizationFactors(dds))) {
            offset <- matrix(log(sizeFactors(dds)), nrow = nrow(dds),
                ncol = ncol(dds), byrow = TRUE)
        }
        else {
            offset <- log(normalizationFactors(dds))
        }
        if ("weights" %in% assayNames(dds)) {
            weights <- assays(dds)[["weights"]]
        }
        else {
            weights <- matrix(1, nrow = nrow(dds), ncol = ncol(dds))
        }
        if (apeAdapt) {
            mle <- log(2) * cbind(res$log2FoldChange, res$lfcSE)
        }
        else {
            mle <- NULL
        }
        if (apeMethod == "general") {
            log.lik <- apeglm::logLikNB
        }
        else {
            log.lik <- NULL
        }
        if (lfcThreshold > 0) {
            message(paste0("computing FSOS 'false sign or small' s-values (T=",
                round(lfcThreshold, 3), ")"))
            svalue <- TRUE
            apeT <- log(2) * lfcThreshold
        }
        else {
            apeT <- NULL
        }
        if (!parallel) {
            fit <- apeglm::apeglm(Y = Y, x = design, log.lik = log.lik,
                param = disps, coef = coefNum, mle = mle, threshold = apeT,
                weights = weights, offset = offset, method = apeMethod,
                ...)
        }
        else {
            fitList <- bplapply(levels(parallelIdx), function(l) {
                idx <- parallelIdx == l
                apeglm::apeglm(Y = Y[idx, , drop = FALSE], x = design,
                  log.lik = log.lik, param = disps[idx], coef = coefNum,
                  mle = mle, threshold = apeT, weights = weights[idx,
                    , drop = FALSE], offset = offset[idx, , drop = FALSE],
                  method = apeMethod, ...)
            })
            fit <- list()
            ape.cols <- c("map", "sd", "fsr", "svalue", "interval",
                "diag")
            if (lfcThreshold > 0) {
                ape.cols <- c(ape.cols, "thresh")
            }
            for (param in ape.cols) {
                fit[[param]] <- do.call(rbind, lapply(fitList,
                  `[[`, param))
            }
            fit$prior.control <- fitList[[1]]$prior.control
            fit$svalue <- apeglm::svalue(fit$fsr[, 1])
        }
        stopifnot(nrow(fit$map) == nrow(dds))
        conv <- fit$diag[, "conv"]
        if (!all(conv[!is.na(conv)] == 0)) {
            message("Some rows did not converge in finding the MAP")
        }
        res$log2FoldChange <- log2(exp(1)) * fit$map[, coefNum]
        res$lfcSE <- log2(exp(1)) * fit$sd[, coefNum]
        mcols(res)$description[2] <- sub("MLE", "MAP", mcols(res)$description[2])
        mcols(res)$description[3] <- sub("standard error", "posterior SD",
            mcols(res)$description[3])
        if (svalue) {
            coefAlphaSpaces <- gsub("_", " ", coefAlpha)
            res <- res[, 1:3]
            if (lfcThreshold > 0) {
                res$svalue <- as.numeric(apeglm::svalue(fit$thresh))
                description <- paste0("FSOS s-value (T=", round(lfcThreshold,
                  3), "): ", coefAlphaSpaces)
            }
            else {
                res$svalue <- as.numeric(fit$svalue)
                description <- paste0("s-value: ", coefAlphaSpaces)
            }
            mcols(res)[4, ] <- DataFrame(type = "results", description = description)
        }
        else {
            res <- res[, c(1:3, 5:6)]
        }
        metadata(res)[["lfcThreshold"]] <- lfcThreshold
        priorInfo(res) <- list(type = "apeglm", package = "apeglm",
            version = packageVersion("apeglm"), prior.control = fit$prior.control)
        res <- resultsFormatSwitch(object = dds, res = res, format = format)
        if (returnList) {
            return(list(res = res, fit = fit))
        }
        else {
            return(res)
        }
    }
    else if (type == "ashr") {
        if (!requireNamespace("ashr", quietly = TRUE)) {
            stop("type='ashr' requires installing the CRAN package 'ashr'")
        }
        if (!quiet)
            message("using 'ashr' for LFC shrinkage. If used in published research, please cite:\n    Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.\n    https://doi.org/10.1093/biostatistics/kxw041")
        if (lfcThreshold > 0)
            message("lfcThreshold is not used by type='ashr'")
        betahat <- res$log2FoldChange
        sebetahat <- res$lfcSE
        fit <- ashr::ash(betahat, sebetahat, mixcompdist = "normal",
            method = "shrink", ...)
        res$log2FoldChange <- fit$result$PosteriorMean
        res$lfcSE <- fit$result$PosteriorSD
        mcols(res)$description[2] <- sub("MLE", "MMSE", mcols(res)$description[2])
        mcols(res)$description[3] <- sub("standard error", "posterior SD",
            mcols(res)$description[3])
        if (svalue) {
            coefAlphaSpaces <- sub(".*p-value: ", "", mcols(res)$description[5])
            res <- res[, 1:3]
            res$svalue <- fit$result$svalue
            mcols(res)[4, ] <- DataFrame(type = "results", description = paste("s-value:",
                coefAlphaSpaces))
        }
        else {
            res <- res[, c(1:3, 5:6)]
        }
        metadata(res)[["lfcThreshold"]] <- lfcThreshold
        priorInfo(res) <- list(type = "ashr", package = "ashr",
            version = packageVersion("ashr"), fitted_g = fit$fitted_g)
        res <- resultsFormatSwitch(object = dds, res = res, format = format)
        if (returnList) {
            return(list(res = res, fit = fit))
        }
        else {
            return(res)
        }
    }
}


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
process_contrast_limma=function(DE_object,contrast,th=0.05,title){
  name <-deparse(substitute(contrast))
  vec=rep(0,length(colnames(design)))
  vec[abs(contrast)]=contrast/abs(contrast)
  fit2 <- contrasts.fit(fit, vec)
  fit2 <- eBayes(fit2)
  b=data.frame(fit2$coefficients)
  t=data.frame(fit2$t)
  p=data.frame(fit2$p.value)
  bh=data.frame(p)
  # top.table <- topTable(fit2, sort.by = "P", n = Inf)
  
  e <- sqrt(fit2$s2.post) * fit2$stdev.unscaled
  d=fit2$df.residual+fit2$df.prior
  
  print(d[1])
  
  for(i in 1:ncol(bh))
  {
    bh[,i]=p.adjust(p[,i],method="BH")
    print(paste("i=",i,name,length(which(bh[,i]<th))))
    
  }
  ## Jq: ojo esto saca los DEGs sólo del ultimo contraste si se corren varios.
  top=data.frame(topTable(fit2),number=length(which(bh[,i]<th)))

  betas_ext=cbind(DE_object$betas,new=b)
  ts_ext=cbind(DE_object$ts,new=t)
  ps_ext=cbind(DE_object$ps,new=p)
  BH_ext=cbind(DE_object$BH,new=bh)
  errors_ext=cbind(DE_object$errors,new=e)
  dofs_ext=cbind(DE_object$dofs,new=d)
  
  colnames(betas_ext)[ncol(betas_ext)]=name
  colnames(ts_ext)[ncol(ts_ext)]=name
  colnames(ps_ext)[ncol(ps_ext)]=name
  colnames(BH_ext)[ncol(BH_ext)]=name
  colnames(errors_ext)[ncol(errors_ext)]=name
  colnames(dofs_ext)[ncol(dofs_ext)]=name
  colnames(top)[ncol(top)]=name
  colnames(p)="P_value"
  nbins=50
  pl=histogram(p,p$P_value,title=title)
  DE_output=list(betas=betas_ext,ts=ts_ext,ps=ps_ext,BH=BH_ext,
                 errors=errors_ext,dofs=dofs_ext,top=top,histogram=pl)
  
  return(DE_output)
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

dir.create("histogramsC5_C6_JAC")
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

pdf(file = "histogramsC5_C6_JAC/PCA_C5_C6.pdf",width=6,height=5)
print(pl1)
dev.off()

## Jq: Add PCs within culture.

## DE: limma. Jq: again, economy of code:
dge <- DGEList(counts = reads)
dge <- calcNormFactors(dge)
design=model.matrix(~ short_setup,data=metadata)
v=voom(dge,design,plot=TRUE)
fit <-lmFit(v,design)
fit <- eBayes(fit)
  
betas=data.frame(fit$coefficients)
ts=data.frame(fit$t)
ps=data.frame(fit$p.value)
BH=data.frame(ps)
errors <- sqrt(fit$s2.post) * fit$stdev.unscaled
dofs=data.frame(base=(fit$df.residual+fit$df.prior))
  
th=0.01
for(i in 1:ncol(BH))
{
  BH[,i]=p.adjust(ps[,i],method="BH")
  print(paste("i=",i,length(which(BH[,i]<th))))
}
DE_base=list(betas=betas,ts=ts,ps=ps,BH=BH,errors=errors,dofs=dofs)
  
  ## Get contrasts:
DE_object=DE_base
  # contrast=a
th=0.05
  
# Iron (Jq:) deprivation effects
a=c(-2)
b=c(3,-4)
  
    # Stat effects
c=c(4,-2)
d=c(3)
    
    ## Stat x iron
    
a_minus_b=c(-2,-3,4)

safe=DE_base
  
pdf("histogramsC5_C6_JAC/histograms_C5_C6_limma.pdf",width=6,height=5)
{
    ## Iron effects
DE_base=process_contrast_limma(DE_object=DE_base,th=0.01,contrast=a,title = "Iron deprivation effect @ C5")
pl_a=DE_base$histogram
topC5vsC5 = DE_base$top
print(pl_a)
    
DE_base=process_contrast_limma(DE_object=DE_base,th=0.01,contrast=b,title = "Iron deprivation effect @ C6")
pl_b=DE_base$histogram
topC6vsC6 = DE_base$top
print(pl_b)
    
    ## Stat effects
    
DE_base=process_contrast_limma(DE_object=DE_base,th=0.01,contrast=c,title = "Transition to stationary effects; with iron")
pl_c=DE_base$histogram
topC6vsC5 = DE_base$top
print(pl_c)
    
DE_base=process_contrast_limma(DE_object=DE_base,th=0.01,contrast=d,title = "Transition to stationary effects; without iron")
pl_d=DE_base$histogram
topC6vsC5_no = DE_base$top
print(pl_d)
    
    # Iron x Stat
    
DE_base=process_contrast_limma(DE_object=DE_base,th=0.01,contrast=a_minus_b,title = "Phase x iron interaction")
pl_a_minus_b=DE_base$histogram
topC5_C5vsC6_C6 = DE_base$top
print(pl_a_minus_b)
}
dev.off()

## Jq: reiterativo. top_genes <- cbind(topC5vsC5,topC6vsC5,topC6vsC5_no,topC6vsC6,topC5_C5vsC6_C6)
## Tal como lo saco tienen distintas rows, los 10 primeros no nos da ninguna info, todo, ya está en la lista DE_base
# top.table$Gene <- rownames(top.table)
# top.table <- top.table[,c("Gene", names(top.table)[1:6])]
# write.table(top.table, file = "table_C5_C6.txt", row.names = F, sep = "\t", quote = F)

## DE: DESEQ2.
{
dds <- DESeqDataSetFromMatrix(countData = reads,colData = metadata,design= ~ short_setup)
dds <- DESeq(dds)
result <- results(dds, name="short_setup_C6_Fe_NO_vs_C5_Fe_NO")

  res.ash <- lfcShrink(dds=dds, coef=2, type="ashr")
  result <- lfcShrink(dds, contrast=lista, type="ashr")

dev.new()
pdf("histogramsC5_C6_JAC/histograms_C5_C6_DESeq2.pdf",width=6,height=5)
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

## EDGER
{
    metadata$short_setup <- factor(metadata$short_setup)
    ## Aquí todo esto que haces es reiterativo, pues ya habíamos filtrado low exp.
    
    y <- DGEList(reads)
    y <- calcNormFactors(y)
    
    design=model.matrix(~0+short_setup,data=metadata)
    colnames(design) <- levels(metadata$short_setup)
    
    #To perform quasi-likelihood F-tests:
    y <- estimateDisp(y,design)
    fit <- glmQLFit(y,design,robust=TRUE)
    
    ## Ojo: toptags da log2FCs; fit$coeffs; natural logs!! https://support.bioconductor.org/p/50273/
    ## betas_edger=log2(exp(fit$coefficients))
    ## plot(betas$conditionLPS,betas_edger[,"conditionLPS"])

{
  ###                         Contrasts
  #=========================== Effect EXP - STAT =====================================
C6vsC5_Fe_NO <- makeContrasts(C6_Fe_NO-C5_Fe_NO, levels=design)
C6vsC5 <- makeContrasts(C6_Fe_YES-C5_Fe_YES, levels=design)
  #=========================== Effect of Fe =====================================
C5vsC5 <- makeContrasts(C5_Fe_NO-C5_Fe_YES, levels=design)
C6vsC6 <- makeContrasts(C6_Fe_NO-C6_Fe_YES, levels=design)
  #========================== IRON X STAT========================================
C6_C6_NOvsC5_C5_NO<- makeContrasts(C6vsC6-C5vsC5, levels=design)
}
  ########################## DE for each comparison ############################
{
qlf.C6vsC5_Fe_NO <- glmQLFTest(fit, contrast=C6vsC5_Fe_NO)
tag_C6vsC5_Fe_NO <- topTags(qlf.C6vsC5_Fe_NO,n=nrow(reads),adjust.method = "BH",sort.by = "none")
  
qlf.C6vsC5 <- glmQLFTest(fit, contrast=C6vsC5)
tag_C6vsC5 <- topTags(qlf.C6vsC5,n=nrow(reads),adjust.method = "BH",sort.by = "none")
  
qlf.C5vsC5  <- glmQLFTest(fit, contrast=C5vsC5 )
tag_C5vsC5<- topTags(qlf.C5vsC5,n=nrow(reads),adjust.method = "BH",sort.by = "none")
  
qlf.C6vsC6 <- glmQLFTest(fit, contrast=C6vsC6)
tag_C6vsC6 <- topTags(qlf.C6vsC6,n=nrow(reads),adjust.method = "BH",sort.by = "none")
  
qlf.C6_C6_NOvsC5_C5_NO <- glmQLFTest(fit, contrast=C6_C6_NOvsC5_C5_NO)
tag_C6_C6_NOvsC5_C5_NO <- topTags(qlf.C6_C6_NOvsC5_C5_NO,n=nrow(reads),adjust.method = "BH",
                                    sort.by = "none")
}

th=0.01
sum(tag_C6vsC5_Fe_NO$table$FDR < th)
sum(tag_C6vsC5$table$FDR < th)
sum(tag_C5vsC5$table$FDR < th)
sum(tag_C6vsC6$table$FDR < th)
sum(tag_C6_C6_NOvsC5_C5_NO$table$FDR < th)

### declare Pvalues 
p_C6vsC5_Fe_no <- qlf.C6vsC5_Fe_NO$table$PValue
p_C6vsC5 <- qlf.C6vsC5$table$PValue
p_C5vsC5 <- qlf.C5vsC5$table$PValue
p_C6vsC6 <- qlf.C6vsC6$table$PValue
p_C6_C6_novsC5_C5_no <- qlf.C6_C6_NOvsC5_C5_NO$table$PValue
  
  # Print the plot to a pdf file

  
  
pdf("histogramsC5_C6_JAC/histograms_C5_C6_edgeR.pdf",width=6,height=5)
  
  ### Histogram Pvalue vs Count
{
hist_1 <- histogram(qlf.C6vsC5 ,p_C6vsC5,"Transition to stationary effects; with iron")
print(hist_1)
hist_2 <- histogram(qlf.C6vsC5_Fe_NO ,p_C6vsC5_Fe_no,"Transition to stationary effects; without iron")
print(hist_2)
hist_3 <- histogram(qlf.C5vsC5,p_C5vsC5,"Iron deprivation effect @ C5")
print(hist_3)
hist_4 <- histogram(qlf.C6vsC6,p_C6vsC6,"Iron deprivation effect @ C6")
print(hist_4)
hist_4<- histogram(qlf.C6_C6_NOvsC5_C5_NO ,p_C6_C6_novsC5_C5_no,"Phase x iron interaction")
print(hist_4)
}

dev.off()
}


tag_C6vsC5_Fe_NO
tag_C6vsC5
tag_C5vsC5
tag_C6vsC6
tag_C6_C6_NOvsC5_C5_NO


res_a
res_b
res_c
res_d
res_a_minus_b

DE_base$betas

deseq=res_b$data$log2FoldChange
edger=tag_C6vsC6$table$logFC
limma=DE_base$betas$b

x=DE_base$betas$b
y=tag_C6vsC6$table$logFC
z=res_b$data$log2FoldChange
labx="limma"
laby="EDGER"
labz="DESEQ"
library(cowplot)

compare_signal_betas=function(x,y,z,labx="limma",laby="EDGER",labz="DESEQ",title,limite=3)
{
tab=data.frame(x=x,y=y,z=z)
tabx=data.frame(stats=x,case=labx)
taby=data.frame(stats=y,case=laby)
tabz=data.frame(stats=z,case=labz)

tabm=rbind(tabx,taby,tabz)
plxy=ggplot(tab)+geom_abline(slope=1,intercept=0,color="red")+geom_point(aes(x=x,y=y))+xlab(labx)+ylab(laby)
plxz=ggplot(tab)+geom_abline(slope=1,intercept=0,color="red")+geom_point(aes(x=x,y=z))+xlab(labx)+ylab(labz)
plyz=ggplot(tab)+geom_abline(slope=1,intercept=0,color="red")+geom_point(aes(x=y,y=z))+xlab(laby)+ylab(labz)

pl_densities=ggplot(tabm)+geom_density(aes(x=stats,fill=case),alpha=0.2)+xlim(-limite,limite)+ggtitle(title)

pl_all=plot_grid(pl_densities,plxy,plxz,plyz,ncol=2,align="hv")
return(pl_all)
}


plot_1=compare_signal_betas(x=DE_base$betas$b,y=tag_C5vsC5$table$logFC,z=res_a$data$log2FoldChange,title="Fe deprivation effect @exp")
pdf("histogramsC5_C6_JAC/methods_comparison_Fe_effect_at_exp.pdf",width=10,height=10)
print(plot_1)
dev.off()

plot_2=compare_signal_betas(x=DE_base$betas$b,y=tag_C6vsC6$table$logFC,z=res_b$data$log2FoldChange,title="Fe deprivation effect @stat")
pdf("histogramsC5_C6_JAC/methods_comparison_Fe_effect_at_stat.pdf",width=10,height=10)
print(plot_2)
dev.off()

plot_3=compare_signal_betas(x=DE_base$betas$c,y=tag_C6vsC5$table$logFC,z=res_c$data$log2FoldChange,title="Stat effect w. Fe",limite=5)
pdf("histogramsC5_C6_JAC/methods_comparison_stat_ffect_w_Fe.pdf",width=10,height=10)
print(plot_3)
dev.off()

plot_4=compare_signal_betas(x=DE_base$betas$d,y=tag_C6vsC5_Fe_NO$table$logFC,z=res_d$data$log2FoldChange,title="Stat effect wout. Fe",limite=5)
pdf("histogramsC5_C6_JAC/methods_comparison_stat_effect_wout_Fe.pdf",width=10,height=10)
print(plot_4)
dev.off()

plot_5=compare_signal_betas(x=DE_base$betas$a_minus_b,y=-tag_C6_C6_NOvsC5_C5_NO$table$logFC,z=res_a_minus_b$data$log2FoldChange,title="Stat X Fe Interaction")
pdf("histogramsC5_C6_JAC/methods_comparison_stat_Fe_interaction.pdf",width=10,height=10)
print(plot_5)
dev.off()


## Jq: Conclusion: DESeq nos está diciendo la verdad: cuando el modelo está over-dispersed, limma y edger dan effect sizes muy grandes y en gran medida aleatorios (lo cual no es del todo mentir, puesto que no dan p valores significativos) De cualquier modo, DEseq2 lo hace mejor, porque los hunde (shrinkage) a casi cero. Nos quedamos con DEseq.

dir.create("Stats_tables")
write.table(data.frame(res_a$data),"Stats_tables/Iron_deprivation_effect_at_Exp.txt")
write.table(data.frame(res_b$data),"Stats_tables/Iron_deprivation_effect_at_Stat.txt")
write.table(data.frame(res_c$data),"Stats_tables/Stat_effect_with_Fe.txt")
write.table(data.frame(res_d$data),"Stats_tables/Stat_effect_without_Fe.txt")
write.table(data.frame(res_a_minus_b$data),"Stats_tables/Stat_Fedeprivation_interaction.txt")


how_much=function(tab,name){
    for(th in c(0.05,0.01))
    {
        for(th_size in c(0,0.2,0.5,1))
        {
            dir.create(paste0("target_sets/th_",th,"_th_size_",th_size),recursive=TRUE)
            set=(which(abs(tab$log2FoldChange)>th_size & tab$BH<th))
            print(paste("th=",th,"th_size=",th_size,":",length(set)))
            target=rownames(tab)[set]
            sink(paste0("target_sets/th_",th,"_th_size_",th_size,"/",name,".txt"))
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


