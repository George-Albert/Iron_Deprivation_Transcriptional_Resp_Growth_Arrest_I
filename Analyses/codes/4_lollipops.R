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
    library(eulerr)

}

################################
### 1. Load & pretty up data ###
################################

{
    feature_data=read.table("inputs/processed/unfiltered/txts/feature_data.txt")

    Iron_deprivation_effect_at_Exp=read.table("outputs/Data/txt/Iron_deprivation_effect_at_Exp.txt")
    Iron_deprivation_effect_at_Stat=read.table("outputs/Data/txt/Iron_deprivation_effect_at_Stat.txt")
    Stat_effect_with_Fe=read.table("outputs/Data/txt/Stat_effect_with_Fe.txt")
    Stat_effect_without_Fe=read.table("outputs/Data/txt/Stat_effect_without_Fe.txt")
    Stat_Fedeprivation_interaction=read.table("outputs/Data/txt/Stat_Fedeprivation_interaction.txt")

    ## Check congruence of stat files:

    length(which(rownames(Iron_deprivation_effect_at_Exp)!=rownames(Iron_deprivation_effect_at_Stat)))
    length(which(rownames(Iron_deprivation_effect_at_Exp)!=rownames(Stat_effect_with_Fe)))
    length(which(rownames(Iron_deprivation_effect_at_Exp)!=rownames(Stat_effect_without_Fe)))
    length(which(rownames(Iron_deprivation_effect_at_Exp)!=rownames(Stat_Fedeprivation_interaction)))

    ## cleanup featuredata
    feature_data$Type[which(feature_data$Type==".")]="misc_RNA"
    feature_data=feature_data[which(feature_data$Gene_ID %in% rownames(Iron_deprivation_effect_at_Exp)),]
    rownames(feature_data)=feature_data$Gene_ID
    feature_data=feature_data[order(tolower(rownames(feature_data))),]
    length(which(rownames(Iron_deprivation_effect_at_Exp)!=rownames(feature_data)))

}

#######################################
### 2. Flipping genes lollipop plot ###
#######################################

{
    ### 1. Select genes
    threshold=0.01
    delta=0.01
    flipping=rownames(Stat_effect_without_Fe)[which(Stat_effect_without_Fe$BH<threshold & Stat_effect_with_Fe$BH<threshold & Stat_effect_without_Fe$log2FoldChange_shrunken*Stat_effect_with_Fe$log2FoldChange_shrunken<0)]

    ### Prepare tables to plot
    flip_w=Stat_effect_with_Fe[flipping,c(2,6,10)]
    flip_wout=Stat_effect_without_Fe[flipping,c(2,6,10)]

    flip_w=flip_w[order(flip_wout$log2FoldChange_shrunken),]
    flip_wout=flip_wout[order(flip_wout$log2FoldChange_shrunken),]

    length(which(rownames(flip_w)!=rownames(flip_wout)))

    flip_w$x=c(1:nrow(flip_w))-delta
    flip_wout$x=c(1:nrow(flip_w))+delta

    colnames(flip_w)=c("Gene","logFC","BH","x")
    colnames(flip_wout)=c("Gene","logFC","BH","x")

    colors_with=c("#A3C9F5","#1461B8")

    colors_without=c("mistyrose","red","firebrick4")
    colors_without=c("#FEE2CD","#FFA866","#FF7A14")

    values_without=(c(2,15,115)-2)/(115-2)
    values_with=(c(2,15)-2)/13

    flip_w$max=pmax(flip_w$logFC,flip_wout$logFC)
    flip_w$Gene[which(flip_w$Gene==".")]=rownames(flip_w)[which(flip_w$Gene==".")]

    row_1=c(NA,1,0.01,30)
    row_2=c(NA,0.5,1E-60,30)
    row_3=c(NA,-1,1E-120,30)

    flip_wout=rbind(flip_wout,row_1,row_2,row_3)

    row_1=c(NA,1,0.01,30)
    row_2=c(NA,0.5,1E-15,30)

    flip_w=rbind(flip_w,row_1,row_2)

    
    flip_w$shift=nchar(flip_w$Gene)/10
    ## This is tinkered manually
    ajustes=c(0.05,-0.05,0,0,0,0,0,0.1,-0.05,0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.05,0.15,0.2,0.15,0.1,0.2,0.2)
    flip_w$shift=flip_w$shift+ajustes

    
    ## This is problematic later for the mycobactins, where there are non-significant cases that only with size cannot be distinghished.
    ## Here I could preserve it, but later I would need to introduce a slightly vis for the mbts, so I will create that vis directly here.
    flip_plot=ggplot()+
    geom_segment(aes(x=x,xend=x,y=0,yend=logFC),color="dodgerblue",data=flip_w)+
    geom_segment(aes(x=x,xend=x,y=0,yend=logFC),color="firebrick4",data=flip_wout)+
    geom_point(aes(x=x,y=logFC,size=-log10(BH)),color="dodgerblue",data=flip_w)+
    geom_point(aes(x=x,y=logFC,size=-log10(BH)),color="firebrick4",data=flip_wout)+
    scale_size(range = c(4,10))+coord_flip()+geom_text(aes(x=x,y=max+shift,label=Gene),data=flip_w)+ylim(-3,5)+ylab("logFC response to growth arrest")+
    theme(
    axis.text.x   = element_text(size=14),
    axis.text.y   = element_blank(),
    axis.title.x  = element_text(size=14),
    axis.title.y  = element_blank(),
    axis.ticks.y = element_blank(),
    #panel.background = element_blank(),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1,linetype="solid"),
    #legend.title=element_blank(),
    legend.position="top",
    #legend.text=element_text(size=14),
    legend.key.size = unit(1, 'lines'))+geom_hline(yintercept=0)

    ## which is this one:
    flip_plot_ok=ggplot()+
    geom_segment(aes(x=x,xend=x,y=0,yend=logFC),color="#1461B8",data=flip_w)+
    geom_segment(aes(x=x,xend=x,y=0,yend=logFC),color="#FF7A14",data=flip_wout)+
    geom_point(aes(x=x,y=logFC,color=-log10(BH),size=-log10(BH)),data=flip_w)+scale_colour_gradientn(colours = colors_with,values =values_with, breaks = c(2,8,15))+
    geom_point(aes(x=x,y=logFC,fill=-log10(BH),size=-log10(BH)),color="grey90",shape=21,data=flip_wout)+scale_fill_gradientn(colours = colors_without,values =values_without,breaks = c(2,60,120))+
    geom_point(aes(x=x,y=logFC,size=-log10(BH)),fill=NA,color="grey90",shape=21,data=flip_w)+
    scale_size(range = c(4,10))+coord_flip()+geom_text(aes(x=x,y=max+shift,label=Gene),data=flip_w)+ylim(-3,5)+ylab("logFC response to growth arrest")+
    theme(
    axis.text.x   = element_text(size=14),
    axis.text.y   = element_blank(),
    axis.title.x  = element_text(size=14),
    axis.title.y  = element_blank(),
    axis.ticks.y = element_blank(),
    #panel.background = element_blank(),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1,linetype="solid"),
    #legend.title=element_blank(),
    legend.position="right",
    #legend.text=element_text(size=14),
    legend.key.size = unit(1, 'lines'))+geom_hline(yintercept=0)
    
    pdf("outputs/Figures/Unpublished/flipping_genes_stat_FCs.pdf",width=7,height=8)
    print(flip_plot_ok)
    dev.off()
    ## It will require manual tinkering.

}

####################################
### 3. Mycobactins lollipop plot ###
####################################

{
    
    feature_data$Gene_Symbol[which(rownames(feature_data)=="Rv2377c")]="mbtH"
    mycobactins=paste0("mbt",c("A","B","C","D","E","F","G","H","I","J","K","L","M","N"))
    
    mycobactins_data=feature_data[which(feature_data$Gene_Symbol %in% mycobactins),]
    mycobactins_RV=rownames(mycobactins_data)


    threshold=0.01
    delta=0.05

    ### Prepare tables to plot
    mbt_exp=Iron_deprivation_effect_at_Exp[mycobactins_RV,c(2,6,10)]
    mbt_stat=Iron_deprivation_effect_at_Stat[mycobactins_RV,c(2,6,10)]

    mbt_stat=mbt_stat[order(mbt_exp$BH),]
    mbt_exp=mbt_exp[order(mbt_exp$BH),]
    
    mbt_stat$Gene[which(rownames(mbt_stat)=="Rv2377c")]="mbtH"
    mbt_exp$Gene[which(rownames(mbt_exp)=="Rv2377c")]="mbtH"

    length(which(rownames(mbt_exp)!=rownames(mbt_stat)))

    mbt_exp$x=c(1:nrow(mbt_exp))-delta
    mbt_stat$x=c(1:nrow(mbt_exp))+delta

    colnames(mbt_exp)=c("Gene","logFC","BH","x")
    colnames(mbt_stat)=c("Gene","logFC","BH","x")

    colors_with=c("grey20","grey95","#E7F7D4","#77B12F")
    colors_without=c("grey20","grey95","#F5A3A9","#ED3B47","#B8141F")

    values_without=c(0,1.99,2.01,120,240)/240
    values_with=c(0,1.99,2.01,16)/16

    mbt_stat$max=pmax(mbt_exp$logFC,mbt_stat$logFC)

    row_1=c(NA,1,1,15)
    row_2=c(NA,1,0.01,15)
    row_3=c(NA,1,1E-120,15)
    row_4=c(NA,1.5,1E-240,15)

    mbt_stat=rbind(mbt_stat,row_1,row_2,row_3,row_4)
  
    row_1=c(NA,1,1,15)
    row_2=c(NA,1,0.01,15)
    row_3=c(NA,1,1E-16,15)

    mbt_exp=rbind(mbt_exp,row_1,row_2,row_3)
    

    mbt_stat$shift=nchar(mbt_stat$Gene)/4
    ## This is tinkered manually
    #ajustes=c(0.05,-0.05,0,0,0,0,0,0.1,-0.05,0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.05,0.15,0.2,0.15,0.1,0.2,0.2)
    #mbt_exp$shift=mbt_exp$shift+ajustes
    
    mbt_plot=ggplot()+
    geom_segment(aes(x=x,xend=x,y=0,yend=logFC),color="#77B12F",size=1,data=mbt_exp)+
    geom_segment(aes(x=x,xend=x,y=0,yend=logFC),color="#B8141F",size=1,data=mbt_stat)+
    geom_point(aes(x=x,y=logFC,color=-log10(BH),size=-log10(BH)),data=mbt_exp)+
    scale_colour_gradientn(colours = colors_with,values =values_with, breaks = c(0,2,16))+
    geom_point(aes(x=x,y=logFC,fill=-log10(BH),size=-log10(BH)),color="grey80",stroke=1,shape=21,data=mbt_stat)+
    scale_fill_gradientn(colours = colors_without,values =values_without,breaks = c(0,2,120,240))+
    geom_point(aes(x=x,y=logFC,size=-log10(BH)),fill=NA,color="grey80",stroke=1,shape=21,data=mbt_exp)+
    scale_size(range = c(4,10))+
    coord_flip()+
    geom_text(aes(x=x,y=max+shift,label=Gene),data=mbt_stat)+ylim(-0.1,8.5)+ylab("logFC response to Fe deprivation")+
    theme(
    axis.text.x   = element_text(size=14),
    axis.text.y   = element_blank(),
    axis.title.x  = element_text(size=14),
    axis.title.y  = element_blank(),
    axis.ticks.y = element_blank(),
    #panel.background = element_blank(),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1,linetype="solid"),
    #legend.title=element_blank(),
    legend.position="right",
    #legend.text=element_text(size=14),
    legend.key.size = unit(1, 'lines'))+geom_hline(yintercept=0)

    pdf("outputs/Figures/Fig_S3_mycobactins_Fe_deprivation_FCs.pdf",width=5,height=6)
    print(mbt_plot)
    dev.off()
    ## It will require manual tinkering.

}

############################################################
### 4. Stat responses vs tuberculist FGSEA lollipop plot ###
############################################################

{
    with=read.delim("outputs/Data/Tuberculist/SEWFe/fgsea_output_tables/SEWFe.txt")[,1:7]
    wout=read.delim("outputs/Data/Tuberculist/SEWOFe/fgsea_output_tables/SEWOFe.txt")[,1:7]

    rownames(with)=with$pathway
    rownames(wout)=wout$pathway
    
    with=with[order(rownames(with)),]
    wout=wout[order(rownames(wout)),]

    length(which(rownames(with)!=rownames(wout)))
    
    with=with[,c(6,7,3)]
    wout=wout[,c(6,7,3)]
    
    colnames(with)=paste0(colnames(with),"_w")
    colnames(wout)=paste0(colnames(wout),"_wout")
    
    all=cbind(with,wout)

    all$pathway=rownames(all)
    
    delta=0.1
    all$x_w=c(1:nrow(all))-delta
    all$x_wout=c(1:nrow(all))+delta
    all$x=c(1:nrow(all))
    all$max=pmax(all$NES_w,all$NES_wout)
    all$shift=nchar(all$pathway)/10

    
    colors_with=c("grey20","grey90","#A3C9F5","#1461B8")

    colors_without=c("grey20","grey90","#FEE2CD","#FF7A14")

    values=c(0,-log10(0.051),-log10(0.05),6)/6
    
    all$y=4.9

    
    all=all[order(-all$NES_w),]
    cosa=c(0,100,1,0,100,1,"test",9.9,10.1,10,1,0.4,4.9)
    cosa2=c(0,100,1E-6,0,100,1E-6,"test",10.9,11.1,11,1,0.4,4.9)
    all=rbind(all,cosa)
    all=rbind(all,cosa2)
    
    dim(all)
    for(i in c(1:6,8:ncol(all)))
    all[,i]=as.numeric(all[,i])

    tuberculist=ggplot(all)+
    geom_segment(aes(x=x_w,xend=x_w,y=0,yend=NES_w),color="#1461B8")+
    geom_segment(aes(x=x_wout,xend=x_wout,y=0,yend=NES_wout),color="#FF7A14")+
    geom_point(aes(x=x_w,y=NES_w,color=-log10(padj_w),size=size_w))+scale_colour_gradientn(colours = colors_with,values =values, breaks = c(0,-log10(0.05),6))+
    geom_point(aes(x=x_wout,y=NES_wout,fill=-log10(padj_wout),size=size_wout),color="grey90",shape=21)+scale_fill_gradientn(colours = colors_without,values =values,breaks = c(0,-log10(0.05),6))+
    geom_point(aes(x=x_w,y=NES_w,size=size_w),fill=NA,color="grey90",shape=21)+
    scale_size(range = c(2,6))+
    coord_flip()+
    #geom_text(aes(x=x,y=y,label=pathway),hjust = 0)+ylim(-3,3)+ylab("NES among responses to growth arrest")+
    theme(
    axis.text.x   = element_text(size=14),
    axis.text.y   = element_blank(),
    axis.title.x  = element_text(size=14),
    axis.title.y  = element_blank(),
    axis.ticks.y = element_blank(),
    #panel.background = element_blank(),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1,linetype="solid"),
    #legend.title=element_blank(),
    legend.position="right",
    #legend.text=element_text(size=14),
    legend.key.size = unit(1, 'lines'))+geom_hline(yintercept=0)+ylab("NES")
    
    
    names=ggplot(all)+
    #geom_segment(aes(x=x_w,xend=x_w,y=0,yend=NES_w),color="#1461B8")+
    #geom_segment(aes(x=x_wout,xend=x_wout,y=0,yend=NES_wout),color="#FF7A14")+
    #geom_point(aes(x=x_w,y=NES_w,color=-log10(padj_w),size=size_w))+scale_colour_gradientn(colours = colors_with,values =values_with, breaks = c(2,8,15))+
    #geom_point(aes(x=x_wout,y=NES_wout,fill=-log10(padj_wout),size=size_wout),color="grey90",shape=21)+scale_fill_gradientn(colours = colors_without,values =values_without,breaks = c(2,60,120))+
    #geom_point(aes(x=x_w,y=NES_w,size=size_w),fill=NA,color="grey90",shape=21)+
    #scale_size(range = c(2,6))+
    geom_text(aes(x=x,y=y,label=pathway),hjust = 1)+
    ylim(-0.1,5)+
    ylab("NES among responses to growth arrest")+
    theme_classic()+    coord_flip()+
    theme(
    axis.text.x   = element_blank(),
    axis.text.y   = element_blank(),
    axis.title.x  = element_blank(),
    axis.title.y  = element_blank(),
    axis.ticks.y = element_blank(),
    #panel.background = element_blank(),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    )
    
    pl_tub=plot_grid(names,tuberculist,ncol=2,align="h",axis="tb")
    
    
    pdf("outputs/Figures/Fig_3D_tuberculist_fgsea.pdf",width=7,height=4)
    print(pl_tub)
    dev.off()

}

###################################################################
### 5. Growth arrest responses for assorted genes lollipop plot ###
###################################################################

{

    data=data.frame(RV=Stat_effect_with_Fe$RV,
    Gene=Stat_effect_with_Fe$Gene,
    LogFC_with_Fe=Stat_effect_with_Fe$log2FoldChange_shrunken,
    Fdr_with_Fe=Stat_effect_with_Fe$BH,
    LogFC_without_Fe=Stat_effect_without_Fe$log2FoldChange_shrunken,
    Fdr_without_Fe=Stat_effect_without_Fe$BH,
    LogFC_Interaction=Stat_Fedeprivation_interaction$log2FoldChange_shrunken,
    Fdr_Interaction=Stat_Fedeprivation_interaction$BH)
    rownames(data)=data$RV
    
    atp_genes=c(paste0("Rv130",c(4:9)),"Rv1310","Rv1311")

    Stat_effect_with_Fe[atp_genes,]
    Stat_effect_without_Fe[atp_genes,]
    Iron_deprivation_effect_at_Exp[atp_genes,]
    Iron_deprivation_effect_at_Stat[atp_genes,]
    Stat_Fedeprivation_interaction[atp_genes,]

    get_lolliplot=function(geneset,min_fdr_wout=1E-24,min_fdr_w=1E-9,delta=0.04,shift=0.2,slope_vs=0.25,intercept=2,deltax=1.2){
        
        
        
    name=deparse(substitute(geneset))

    
    local=data[unique(c(which(rownames(data) %in% geneset),which(data$Gene %in% geneset))),]
    
    ausencias=geneset[which((!geneset %in% rownames(data))& (!geneset %in% data$Gene))]
    if(length(ausencias)>0){
        print(paste("In",name,"faltan genes:\n"))
        print(ausencias)
    }
    
    repeticiones=geneset[which(duplicated(geneset))]
    if(length(repeticiones)>0){
        print(paste("In",name,"hay dupes:\n"))
        print(repeticiones)
    }

    
    actual_min_with=min(local$Fdr_with_Fe[which(local$Fdr_with_Fe>0)])
    actual_min_without=min(local$Fdr_without_Fe[which(local$Fdr_without_Fe>0)])

    actual_min_both=min(actual_min_with,actual_min_without)
    local$Fdr_with_Fe[which(local$Fdr_with_Fe==0)]=actual_min_both/10
    local$Fdr_without_Fe[which(local$Fdr_without_Fe==0)]=actual_min_both/10

    
    min_logfdr_with=-log10(1)
    max_logfdr_with=max(-log10(min_fdr_w),-log10(actual_min_with))
    min_logfdr_without=-log10(1)
    max_logfdr_without=max(-log10(min_fdr_wout),-log10(actual_min_without))
    
    minima_x=min(c(local$LogFC_with_Fe,local$LogFC_without_Fe))
    maxima_x=max(c(local$LogFC_with_Fe,local$LogFC_without_Fe))
    
    limite_fv=max(abs(c(minima_x,maxima_x,1.6)))+deltax

    
    local$Int="N.S."
    local$Int[which(local$Fdr_Interaction<0.05 & (local$LogFC_without_Fe-local$LogFC_with_Fe)>0)]="Up"
    local$Int[which(local$Fdr_Interaction<0.05 & (local$LogFC_without_Fe-local$LogFC_with_Fe)<0)]="Down"
    local$Int
    ### Prepare tables to plot
    

    
    flip_w=local[,c(2,3,4,9)]
    flip_wout=local[,c(2,5,6,9)]

    flip_w=flip_w[order(flip_wout$LogFC_without_Fe),]
    flip_wout=flip_wout[order(flip_wout$LogFC_without_Fe),]

    length(which(rownames(flip_w)!=rownames(flip_wout)))

    flip_w$x=c(1:nrow(flip_w))-delta
    flip_wout$x=c(1:nrow(flip_w))+delta

    colnames(flip_w)=c("Gene","logFC","BH","Label","x")
    colnames(flip_wout)=c("Gene","logFC","BH","Label","x")
    
    colors_with=c("grey20","grey95","#A3C9F5","#1461B8")
    colors_without=c("grey20","grey95","#FEE2CD","#FF7A14")

    
    values_without=c(0,1.99,2.01,max_logfdr_without)/max_logfdr_without
    values_with=c(0,1.99,2.01,max_logfdr_with)/max_logfdr_with
    
    #values_without=(c(min_logfdr_without,(min_logfdr_without+max_logfdr_without)/2,max_logfdr_without)-min_logfdr_without)/(max_logfdr_without-min_logfdr_without)
    #values_with=(c(min_logfdr_with,(min_logfdr_with+max_logfdr_with)/2,max_logfdr_with)-min_logfdr_with)/(max_logfdr_with-min_logfdr_with)

    flip_w$Gene[which(flip_w$Gene==".")]=rownames(flip_w)[which(flip_w$Gene==".")]

    row_1=c(NA,-0.5,1,"Up",nrow(flip_wout)+1)
    row_2=c(NA,-1, 10^-max_logfdr_without,"Down",nrow(flip_wout)+1)
    row_3=c(NA,-1.5,10^-max_logfdr_without,"N.S.",nrow(flip_wout)+1)

    flip_wout=rbind(flip_wout,row_1,row_2,row_3)

    row_1=c(NA,0.5,1,"Up",nrow(flip_w)+1)
    row_2=c(NA,1,10^-max_logfdr_with,"Down",nrow(flip_w)+1)
    row_3=c(NA,1.5,10^-max_logfdr_with,"N.S.",nrow(flip_w)+1)

    flip_w=rbind(flip_w,row_1,row_2,row_3)
    for(i in c(2,3,5)){
    flip_w[,i]=as.numeric(flip_w[,i])
    flip_wout[,i]=as.numeric(flip_wout[,i])
    }
    
    flip_w$max=pmax(rep(0,nrow(flip_w)),pmax(flip_w$logFC,flip_wout$logFC))
    flip_w$shift=nchar(flip_w$Gene)/10+shift
    flip_w$Label=factor(flip_w$Label,levels=c("N.S.","Up","Down"))
    color_labels=c("black","forestgreen","firebrick")
    
    ## This is tinkered manually
    #ajustes=c(0.05,-0.05,0,0,0,0,0,0.1,-0.05,0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.05,0.15,0.2,0.15,0.1,0.2,0.2)
    #flip_w$shift=flip_w$shift+ajustes
    flip_w$RV_Gene=paste0(rownames(flip_w),"_",flip_w$Gene)
    
    flip_plot_ok=ggplot()+
    geom_segment(aes(x=x,xend=x,y=0,yend=logFC),color="#1461B8",data=flip_w)+
    geom_segment(aes(x=x,xend=x,y=0,yend=logFC),color="#FF7A14",data=flip_wout)+
    geom_point(aes(x=x,y=logFC,color=-log10(BH),size=-log10(BH)),data=flip_w)+
    scale_colour_gradientn(colours = colors_with,values =values_with, breaks = c(0,2,max_logfdr_with))+
    geom_point(aes(x=x,y=logFC,fill=-log10(BH),size=-log10(BH)),color="grey90",shape=21,data=flip_wout)+
    scale_fill_gradientn(colours = colors_without,values =values_without,breaks = c(0,2,max_logfdr_without))+
    geom_point(aes(x=x,y=logFC,size=-log10(BH)),fill=NA,color="grey90",shape=21,data=flip_w)+
    scale_size(range = c(4,10))+
    coord_flip()+   new_scale_color() +
    geom_text(aes(x=x+delta,y=max+shift,label=RV_Gene,color=Label),data=flip_w)+scale_colour_manual(values = color_labels)+
    ylim(-limite_fv,limite_fv)+
    ylab("logFC response to growth arrest")+
    theme(
    axis.text.x   = element_text(size=14),
    axis.text.y   = element_blank(),
    axis.title.x  = element_text(size=14),
    axis.title.y  = element_blank(),
    axis.ticks.y = element_blank(),
    #panel.background = element_blank(),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1,linetype="solid"),
    #legend.title=element_blank(),
    legend.position="right",
    #legend.text=element_text(size=14),
    legend.key.size = unit(1, 'lines'))+geom_hline(yintercept=0)+ggtitle(paste0(name,". Growth arrests: Orange_-Fe, blue:+Fe\n","Interactions: Black:N.S.\n","Red: growth arrest upreg. enhanced in -Fe\n","Green: growth arrest downreg. enhanced in -Fe"))
    
    dir.create(paste0("outputs/Figures/Unpublished/lollipop_plot_genesets/",name),recursive=TRUE)
    write.table(local,paste0("outputs/Figures/Unpublished/lollipop_plot_genesets/",name,"/stats.txt"))
    
    pdf(paste0("outputs/Figures/Unpublished/lollipop_plot_genesets/",name,"/plot.pdf"),width=8,height=nrow(local)*slope_vs+intercept)
    print(flip_plot_ok)
    dev.off()
    return(list(table=local,pl=flip_plot_ok))
}

    
    oxphos=c("Rv1360","Rv1695","Rv2858c","Rv3145","Rv3146","Rv3147","Rv3148","Rv3149","Rv3150","Rv3151",
    "Rv3152","Rv3153","Rv3154","Rv3155","Rv3156","Rv3157","Rv3158","Rv3559c","Rv1260","Rv1304","Rv1305","Rv1306","Rv1307","Rv1308",
    "Rv1309","Rv1310","Rv1311","Rv1846c","Rv0527","Rv0529","Rv0764c","Rv0766c","Rv0136","Rv3518c","Rv0886")
    
    ## Rv0766 is Rv0766c; y faltaban tres mencionados en el word: "Rv3157","Rv3158","Rv3155"
    
    
    oxphos_vis=get_lolliplot(geneset=oxphos,delta=0.04,shift=1.4,slope_vs=0.25,intercept=2,deltax=3)

    carbon_metabolism=c("Rv1872c","Rv0694","Rv0951","Rv0952","Rv0427c","Rv0428c","Rv0429c","Rv1240","Rv2589","Rv0211")
    ## RV0211 is Rv0211

    carbon_metabolism_vis=get_lolliplot(geneset=carbon_metabolism,delta=0.04,shift=0.6,slope_vs=0.25,intercept=4,deltax=1.5)


TCA_cycle=c("Rv1127c","Rv2241","Rv1092c","Rv0896","Rv1475c","Rv3339c","Rv1098c","Rv1837c","Rv1130","Rv1131","Rv1731","Rv3432c","Rv2221c","Rv0467","Rv0211","Rv2443","Rv3436c","Rv0132c","Rv3842c","Rv1438","Rv1122")

# missannotations
# "Rv1127" "Rv1092" "Rv1475" "Rv3339" "Rv1837" "Rv3432" "Rv2221"
#"Rv1127c" "Rv1092c" "Rv1475c" "Rv3339c" "Rv1837c" "Rv3432c" "Rv2221c"

#TCA_cycle_vis=get_lolliplot(geneset=TCA_cycle,delta=0.04,shift=0.6,slope_vs=0.25,intercept=4,deltax=1.5)

cupper_detoxification=c("Rv0967","Rv0968","Rv0969","Rv0970","Rv0190","Rv2963","Rv0846c","Rv0847","Rv0186A")
## Rv0846 was Rv0846c
cadmium_detoxification=c("Rv1991A","Rv1991c","Rv1992c","Rv1993c","Rv1994c","Rv1995","Rv2640c","Rv2641","Rv2642","Rv2643")
## "Rv1992","Rv1993" were "Rv1992c","Rv1993c"
zinc_detoxification=c("Rv2358","Rv2359",paste0("Rv028",c(2:9)),"Rv0290","Rv0291","Rv0292","ctpC","Rv3269","cmtR")
niquel_cobalt_detoxification=c("Rv0827c","Rv2025c","Rv3744","Rv3743c","Rv3742c")
other_metal_detoxification_channels=c("Rv1239c","Rv2025c")




cupper_detoxification_vis=get_lolliplot(geneset=cupper_detoxification,delta=0.04,shift=2.5,slope_vs=0.25,intercept=4,deltax=5)
cadmium_detoxification_vis=get_lolliplot(geneset=cadmium_detoxification,delta=0.04,shift=4,slope_vs=0.25,intercept=4,deltax=7)
zinc_detoxification_vis=get_lolliplot(geneset=zinc_detoxification,delta=0.04,shift=4,slope_vs=0.25,intercept=4,deltax=7)
niquel_cobalt_detoxification_vis=get_lolliplot(geneset=niquel_cobalt_detoxification,delta=0.04,shift=3,slope_vs=0.25,intercept=4,deltax=6)
other_metal_detoxification_channels_vis=get_lolliplot(geneset=other_metal_detoxification_channels,delta=0.04,shift=3,slope_vs=0.25,intercept=6,deltax=6)

pl_metales=plot_grid(cupper_detoxification_vis[[2]],
cadmium_detoxification_vis[[2]],
zinc_detoxification_vis[[2]],
niquel_cobalt_detoxification_vis[[2]],
other_metal_detoxification_channels_vis[[2]],ncol=2)

dir.create("outputs/Figures/Unpublished/sup_lollipops/",recursive=TRUE)
pdf("outputs/Figures/Unpublished/sup_lollipops/metals.pdf",height=20,width=14)
print(pl_metales)
dev.off()

alternative_glycerol=c("Rv1438","Rv1122","Rv3842c")

alternative_glycerol_vis=get_lolliplot(geneset=alternative_glycerol,delta=0.04,shift=3,slope_vs=0.25,intercept=6,deltax=6)

toxin_anti_toxin_systems=c("Rv0598c","Rv0300","Rv1991A","Rv1991c")
stable_RNAs=c("rrl","rnpB","ssr","Rv2280","Rv3476c")

toxin_anti_toxin_systems_vis=get_lolliplot(geneset=toxin_anti_toxin_systems,delta=0.04,shift=3,slope_vs=0.25,intercept=6,deltax=6)
stable_RNAs_vis=get_lolliplot(geneset=stable_RNAs,delta=0.04,shift=3,slope_vs=0.25,intercept=6,deltax=6)

pl_RNAs_TAs=plot_grid(toxin_anti_toxin_systems_vis[[2]],
stable_RNAs_vis[[2]],
ncol=1)

dir.create("outputs/Figures/Unpublished/sup_lollipops/",recursive=TRUE)
pdf("outputs/Figures/Unpublished/sup_lollipops/RNAs_TAs.pdf",height=15,width=7)
print(pl_RNAs_TAs)
dev.off()


### Parto el TCA en trozos


TCA_cycle_strict=c("Rv0896","Rv1475c","Rv3339c","Rv1098c")
acetyl_coA=c("Rv1127c","Rv2241","Rv1092c")
methylcitrate=c("Rv1130","Rv1131")
glyoxylate_shunt=c("Rv0467","Rv1837c")
GABA_synthesis=c("Rv1731","Rv3432c")
glutamine_pathways=c("Rv2221c","Rv2443","Rv3436c")
alternative_glycerol_routes=c("Rv0132c","Rv3842c","Rv1438","Rv1122")

TCA_cycle_strict_vis=get_lolliplot(geneset=TCA_cycle_strict,delta=0.04,shift=3,slope_vs=0.25,intercept=6,deltax=6)
acetyl_coA_vis=get_lolliplot(geneset=acetyl_coA,delta=0.04,shift=3,slope_vs=0.25,intercept=6,deltax=6)
methylcitrate_vis=get_lolliplot(geneset=methylcitrate,delta=0.04,shift=3,slope_vs=0.25,intercept=6,deltax=6)
glyoxylate_shunt_vis=get_lolliplot(geneset=glyoxylate_shunt,delta=0.04,shift=3,slope_vs=0.25,intercept=6,deltax=6)
GABA_synthesis_vis=get_lolliplot(geneset=GABA_synthesis,delta=0.04,shift=3,slope_vs=0.25,intercept=6,deltax=6)
glutamine_pathways_vis=get_lolliplot(geneset=glutamine_pathways,delta=0.04,shift=3,slope_vs=0.25,intercept=6,deltax=6)
alternative_glycerol_routes_vis=get_lolliplot(geneset=alternative_glycerol_routes,delta=0.04,shift=3,slope_vs=0.25,intercept=6,deltax=6)

pl_TCA=plot_grid(TCA_cycle_strict_vis[[2]],
acetyl_coA_vis[[2]],
methylcitrate_vis[[2]],
glyoxylate_shunt_vis[[2]],
GABA_synthesis_vis[[2]],
glutamine_pathways_vis[[2]],
alternative_glycerol_routes_vis[[2]],ncol=3)

dir.create("outputs/Figures/Unpublished/sup_lollipops/",recursive=TRUE)
pdf("outputs/Figures/Unpublished/sup_lollipops/TCA.pdf",height=20,width=21)
print(pl_TCA)
dev.off()



### Oxphos by parts

atp_synthases=c("Rv1304","Rv1305","Rv1306","Rv1307","Rv1308","Rv1309","Rv1310","Rv1311")
NADH_dehydrogenase=c("Rv3145","Rv3146","Rv3147","Rv3148","Rv3149","Rv3150","Rv3151","Rv3152","Rv3153","Rv3154","Rv3155","Rv3156","Rv3157","Rv3158")
cytochrome_P450=c("Rv0764c","Rv0766c","Rv0136","Rv3518c")
cytochrome_C=c("Rv0527","Rv0529")
oxphos_others=c("Rv1846c","Rv0886","Rv3559c","Rv1360","Rv1695","Rv1260","Rv2858c")

atp_synthases_vis=get_lolliplot(geneset=atp_synthases,delta=0.04,shift=3,slope_vs=0.25,intercept=6,deltax=6)
NADH_dehydrogenase_vis=get_lolliplot(geneset=NADH_dehydrogenase,delta=0.04,shift=3,slope_vs=0.25,intercept=6,deltax=6)
cytochrome_P450_vis=get_lolliplot(geneset=cytochrome_P450,delta=0.04,shift=3,slope_vs=0.25,intercept=6,deltax=6)
cytochrome_C_vis=get_lolliplot(geneset=cytochrome_C,delta=0.04,shift=3,slope_vs=0.25,intercept=6,deltax=6)
oxphos_others_vis=get_lolliplot(geneset=oxphos_others,delta=0.04,shift=3,slope_vs=0.25,intercept=6,deltax=6)


pl_oxphos=plot_grid(atp_synthases_vis[[2]],
NADH_dehydrogenase_vis[[2]],
cytochrome_P450_vis[[2]],
cytochrome_C_vis[[2]],
oxphos_others_vis[[2]],
ncol=2)

dir.create("outputs/Figures/Unpublished/sup_lollipops/",recursive=TRUE)
pdf("outputs/Figures/Unpublished/sup_lollipops/oxphos.pdf",height=20,width=14)
print(pl_oxphos)
dev.off()




PDIM_synthases=c("Rv2937","Rv2942","Rv2953")
PDIM_synthases_vis=get_lolliplot(geneset=PDIM_synthases,delta=0.04,shift=3,slope_vs=0.25,intercept=6,deltax=6)

}


## JQ: No conocia la brujería esta del newscalecolor, ;)
