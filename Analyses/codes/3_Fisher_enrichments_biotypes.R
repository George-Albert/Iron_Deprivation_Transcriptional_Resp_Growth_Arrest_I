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

#########################################################################################################################################
### 2. (DO NOT RUN) Test mashr before going further: this does not produce shrunken logFC that preserve linearity either, so, not run ###
#########################################################################################################################################

if(FALSE){
library(mashr)
betas=data.frame(
Iron_deprivation_effect_at_Exp=Iron_deprivation_effect_at_Exp$log2FoldChange_raw,
Iron_deprivation_effect_at_Stat=Iron_deprivation_effect_at_Stat$log2FoldChange_raw,
Stat_effect_with_Fe=Stat_effect_with_Fe$log2FoldChange_raw,
Stat_effect_without_Fe=Stat_effect_without_Fe$log2FoldChange_raw,
Stat_Fedeprivation_interaction=Stat_Fedeprivation_interaction$log2FoldChange_raw)

rownames(betas)=rownames(Iron_deprivation_effect_at_Exp)


SEs=data.frame(
Iron_deprivation_effect_at_Exp=Iron_deprivation_effect_at_Exp$lfcSE_raw,
Iron_deprivation_effect_at_Stat=Iron_deprivation_effect_at_Stat$lfcSE_raw,
Stat_effect_with_Fe=Stat_effect_with_Fe$lfcSE_raw,
Stat_effect_without_Fe=Stat_effect_without_Fe$lfcSE_raw,
Stat_Fedeprivation_interaction=Stat_Fedeprivation_interaction$lfcSE_raw)

rownames(SEs)=rownames(Iron_deprivation_effect_at_Exp)
mash_data = mash_set_data(as.matrix(betas), as.matrix(SEs))

U.c = cov_canonical(mash_data)
print(names(U.c))
m.c = mash(mash_data, U.c)
lfsr=get_lfsr(m.c)
lfdr=get_lfdr(m.c)
shrunken_betas=get_pm(m.c)
shrunken_SEs=get_psd(m.c)
}

#####################################
### 3. Get FET enrichments (Ugly) ###
#####################################

{
    get_Fisher_column=function(biotype,stat,direction,threshold=0.01,feat=feature_data){
    
    genes_in_type=feat$Gene_ID[which(feat$Type==biotype)]
    genes_in_tgt=rownames(stat)[which(direction*stat$log2FoldChange_shrunken>0 & stat$BH<threshold)]
    
    yes_yes=length(which(genes_in_type %in% genes_in_tgt))
    yes_no=length(genes_in_type)-yes_yes
    no_yes=length(genes_in_tgt)-yes_yes
    no_no=nrow(stat)-yes_no-no_yes-yes_yes
    
    contingency_table <- data.frame(
        "In_biotype" = c(yes_yes, yes_no),
        "Out_biotype" = c(no_yes,no_no),
        row.names = c("In_target", "Out_target"),
        stringsAsFactors = FALSE
    )

    test=fisher.test(contingency_table,alternative="two.sided")
    test$p.value
    test$estimate
    output=c(test$estimate,test$conf.int[1:2],test$p.value,yes_yes,yes_no,no_yes,no_no)
    return(output)
}
    get_Fisher_table=function(stat,direction,threshold=0.01,feat=feature_data){
        types=unique(feat$Type)
        
        result=(data.frame(Type=types,Estimate=NA,Low_CI=NA,Top_CI=NA,p=NA,TYPEyes_TGTyes=NA,TYPEyes_TGTno=NA,TYPEno_TGTyes=NA,TYPEno_TGTno=NA))
        rownames(result)=result$Type
        result=data.frame(t(result[,2:ncol(result)]))
        for(type in types){
            print(type)
            result[,type]=get_Fisher_column(biotype=type,stat=stat,direction=direction,threshold=0.01,feat=feature_data)
        }
        
        result=data.frame(t(result))
        result$Bonferroni=pmin(1,ncol(result)*result$p)
        
        return(result)
    }

    CT_Iron_dep_at_Exp_UP  =get_Fisher_table(stat=Iron_deprivation_effect_at_Exp,direction=1)
    CT_Iron_dep_at_Exp_DOWN=get_Fisher_table(stat=Iron_deprivation_effect_at_Exp,direction=-1)

    CT_Iron_dep_at_Stat_UP  =get_Fisher_table(stat=Iron_deprivation_effect_at_Stat,direction=1)
    CT_Iron_dep_at_Stat_DOWN=get_Fisher_table(stat=Iron_deprivation_effect_at_Stat,direction=-1)

    CT_Stat_effect_with_Fe_UP  =get_Fisher_table(stat=Stat_effect_with_Fe,direction=1)
    CT_Stat_effect_with_Fe_DOWN=get_Fisher_table(stat=Stat_effect_with_Fe,direction=-1)

    CT_Stat_effect_without_Fe_UP  =get_Fisher_table(stat=Stat_effect_without_Fe,direction=1)
    CT_Stat_effect_without_Fe_DOWN=get_Fisher_table(stat=Stat_effect_without_Fe,direction=-1)

    CT_interaction_UP  =get_Fisher_table(stat=Stat_Fedeprivation_interaction,direction=1)
    CT_interaction_DOWN=get_Fisher_table(stat=Stat_Fedeprivation_interaction,direction=-1)

    ## Mirando lo que sale: conclusion: único resultado significativo: tRNAs y ncRNAs are depleted among gens downregulated upon transition to stat.
    ## The rest of categories seem equally depleted but are too small to give significant p-values/Bonferroni FDRS.
    ## Conclusión: dar el OR con solo dos biotypos (prot-coding y no-prot-coding) y solo para stat effects, y finalmente visualizar con un Volcano o con un density plot.
}

##########################################################
### 4. Get summarize FET enrichments (Relevant:REPORT) ###
##########################################################

{
    # First: summarize biotypes:
    feature_data$Type[which(feature_data$Type!="protein_coding")]="Non_protein_coding"
    
    # Now, check enrichments:
    CT_Iron_dep_at_Exp_UP  =get_Fisher_table(stat=Iron_deprivation_effect_at_Exp,direction=1)
    CT_Iron_dep_at_Exp_DOWN=get_Fisher_table(stat=Iron_deprivation_effect_at_Exp,direction=-1)

    CT_Iron_dep_at_Stat_UP  =get_Fisher_table(stat=Iron_deprivation_effect_at_Stat,direction=1)
    CT_Iron_dep_at_Stat_DOWN=get_Fisher_table(stat=Iron_deprivation_effect_at_Stat,direction=-1)

    CT_Stat_effect_with_Fe_UP  =get_Fisher_table(stat=Stat_effect_with_Fe,direction=1)
    CT_Stat_effect_with_Fe_DOWN=get_Fisher_table(stat=Stat_effect_with_Fe,direction=-1)

    CT_Stat_effect_without_Fe_UP  =get_Fisher_table(stat=Stat_effect_without_Fe,direction=1)
    CT_Stat_effect_without_Fe_DOWN=get_Fisher_table(stat=Stat_effect_without_Fe,direction=-1)

    CT_interaction_UP  =get_Fisher_table(stat=Stat_Fedeprivation_interaction,direction=1)
    CT_interaction_DOWN=get_Fisher_table(stat=Stat_Fedeprivation_interaction,direction=-1)

    ### non-protein coding genes are marginally enriched among genes upreg. in Stat transition either with (OR=1.64 CI: 0.942-2.81, p=0.0668) or wout Fe (OR=2.14 CI: 1.25-3.66, p=3.64E-3)
    ### non-protein coding genes are strongly depleted among genes downreg. in Stat transition either with (OR=8.50E-2 CI: 0.0100-0.322, p=9.07E-7) or wout Fe (OR=0.0384 CI: 9.60E-4-0.223, p=1.79E-8)
}

#####################################
### 5. Visualize: volcanos (Ugly) ###
#####################################

{
    volcano_plotter=function(stat,feat=feature_data,threshold=0.01){
        stat$Type=feat$Type
        
        stat$Label="Background"
        stat$Label[which(stat$BH<threshold)]="DE"
        stat$Label[which(stat$Type=="Non_protein_coding")]="Non_PC"
        stat$Label[which(stat$BH<threshold & stat$Type=="Non_protein_coding")]="DE_Non_PC"

        stat=stat[order(stat$Label),]
        
        pl=ggplot(stat)+geom_point(aes(x=log2FoldChange_shrunken,y=-log10(BH),color=Label))+
        theme_classic()+
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
        legend.position="top",
        #legend.text=element_text(size=14),
        legend.key.size = unit(1, 'lines'))+geom_hline(yintercept=threshold)
        
        return(pl)
    }

    volcano_Stat_effect_without_Fe=volcano_plotter(stat=Stat_effect_without_Fe)
    volcano_Stat_effect_with_Fe=volcano_plotter(stat=Stat_effect_with_Fe)

    ## Un poco feotes: demasiado streched in Y: go with density plots:
}


#################################################
### 5. Visualize: densities (Relevant:REPORT) ###
#################################################


density_plotter=function(stat,feat=feature_data,threshold=0.01,culture){
    stat$Type=feat$Type
    stat=stat[which(stat$BH<threshold),]
    
    pc=stat$log2FoldChange_shrunken[which(stat$Type=="protein_coding")]
    npc=stat$log2FoldChange_shrunken[which(stat$Type=="Non_protein_coding")]
    test=wilcox.test(npc,pc,alt="greater")
    print(paste0("Wilcox rank test p=",test$p.value," (alt: npc>pc)"))

    
    pl=ggplot(stat)+geom_density(aes(x=log2FoldChange_shrunken,fill=Type),alpha=0.8)+
    xlim(-10,10)+
    xlab(paste0("logFC Stat effect in ",culture))+
    theme(legend.position="none")+
    geom_vline(xintercept=0)+
    theme_classic()+
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
    legend.position="top",
    #legend.text=element_text(size=14),
    legend.key.size = unit(1, 'lines'))
    
    return(pl)
}

densities_Stat_effect_without_Fe=density_plotter(stat=Stat_effect_without_Fe,culture="Fe/-")
# "Wilcox rank test p=2.82513408605732e-11 (alt: npc>pc)"
densities_Stat_effect_with_Fe=density_plotter(stat=Stat_effect_with_Fe,culture="Fe/+")
# "Wilcox rank test p=2.15908812585328e-10 (alt: npc>pc)"
pl_densities=plot_grid(densities_Stat_effect_without_Fe,densities_Stat_effect_with_Fe,nrow=2,align="hv")

pdf("outputs/Figures/Fig_S2_densities_stats_effects_per_biotype.pdf",width=6,height=6)
print(pl_densities)
dev.off()
