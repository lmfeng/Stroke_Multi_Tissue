# Load required libraries
library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(venn)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
###################differential expressed genes
# neut_tmp=subset(neut,tissue%in%c('CBM','TBM')&celltype%in%c("BDNeu_a", "neutLinkedmacro", "immNeu", "PreNeu", "BDNeu_c", "matNeu", "BDNeu_b"))
# Idents(neut_tmp)=neut_tmp$state
####1. CBM: CTRL vs EXP和TBM:CTRL vs EXP###################
setwd('/home/liyixiang/Fig2')
gexpr <- readRDS('gexpr_cluster1.rds')
DimPlot(gexpr,group.by = 'Cluster1')
#1) CBM: MCAO vs Sham
gexpr_cbm=subset(gexpr,Cluster1%in%c("Neutrophils")&tissue%in%c('CBM'))
Idents(gexpr_cbm)=gexpr_cbm$state
table(gexpr$Cluster1)
gexpr_cbm_state_degs=FindMarkers(gexpr_cbm,ident.1='MCAO',ident.2='Sham')
gexpr_cbm_state_degs_sig=gexpr_cbm_state_degs#[neut_cbm_state_degs$avg_log2FC>0.1&neut_cbm_state_degs$p_val_adj<1e-2,] #&neut_cbm_state_degs$p_val_adj<1e-2

####2)FBM: CTRL vs EXP
gexpr_fbm=subset(gexpr,Cluster1%in%c("Neutrophils")&tissue%in%c('FBM'))
Idents(gexpr_fbm)=gexpr_fbm$state

gexpr_fbm_state_degs=FindMarkers(gexpr_fbm,ident.1='MCAO',ident.2='Sham')
gexpr_fbm_state_degs_sig=gexpr_fbm_state_degs#[neut_tbm_state_degs$avg_log2FC>0.1&neut_tbm_state_degs$p_val_adj<1e-2,] #&neut_tbm_state_degs$p_val_adj<1e-2
# table(neut_tbm_state_degs_sig$cluster)
# neut_tbm_state_degs_sig_list=split(neut_tbm_state_degs_sig$gene,neut_tbm_state_degs_sig$cluster)

#SFIG3) Neut_Cd14 vs Neut_Gngt2
gexpr_cd14=subset(gexpr,Cluster%in%c("Neut_Cd14",'Neut_Gngt2'))
Idents(gexpr_cd14)=gexpr_cd14$Cluster
table(gexpr$Cluster)
gexpr_cd14_state_degs=FindMarkers(gexpr_cd14,ident.1='Neut_Cd14',ident.2 = 'Neut_Gngt2')
gexpr_cd14_state_degs_sig=gexpr_cd14_state_degs#[neut_cbm_state_degs$avg_log2FC>0.1&neut_cbm_state_degs$p_val_adj<1e-2,] #&neut_cbm_state_degs$p_val_adj<1e-2
gexpr_gngt2_state_degs=FindMarkers(gexpr,ident.1='Neut_Gngt2')
gexpr_gngt2_state_degs_sig=gexpr_gngt2_state_degs#[neut_cbm_state_degs$avg_log2FC>0.1&neut_cbm_state_degs$p_val_adj<1e-2,] #&neut_cbm_state_degs$p_val_adj<1e-2


########venn
library(venn)
degs_list=list('Higher in Neut_Cd14'=rownames(gexpr_cd14_state_degs_sig)[gexpr_cd14_state_degs_sig$avg_log2FC>0 & gexpr_cd14_state_degs_sig$p_val_adj<0.05],
               'Higher in Neut_Gngt2'=rownames(gexpr_cd14_state_degs_sig)[gexpr_cd14_state_degs_sig$avg_log2FC<0 & gexpr_cd14_state_degs_sig$p_val_adj<0.05]
)

degs_list=list('CBM:Higher in MCAO'=rownames(gexpr_cbm_state_degs_sig)[gexpr_cbm_state_degs_sig$avg_log2FC>0],
               'CBM:Higher in Sham'=rownames(gexpr_cbm_state_degs_sig)[gexpr_cbm_state_degs_sig$avg_log2FC<0],
               'FBM:Higher in Sham'=rownames(gexpr_fbm_state_degs_sig)[gexpr_fbm_state_degs_sig$avg_log2FC<0],
               'FBM:Higher in MCAO'=rownames(gexpr_fbm_state_degs_sig)[gexpr_fbm_state_degs_sig$avg_log2FC>0]
               
)
CBM_Higher_MCAO <- gexpr_cbm_state_degs_sig[gexpr_cbm_state_degs_sig$avg_log2FC>0.25&gexpr_cbm_state_degs_sig$p_val_adj<0.05,]
CBM_Higher_Sham<- gexpr_cbm_state_degs_sig[gexpr_cbm_state_degs_sig$avg_log2FC<0.25 &gexpr_cbm_state_degs_sig$p_val_adj<0.05,]
CBM_Higher_MCAO_top100 <- CBM_Higher_MCAO%>%top_n(n =100,wt=avg_log2FC)
CBM_Higher_Sham_top100 <- CBM_Higher_Sham%>%top_n(n =100,wt=abs(avg_log2FC))

FBM_Higher_MCAO <- gexpr_fbm_state_degs_sig[gexpr_fbm_state_degs_sig$avg_log2FC>0.25&gexpr_fbm_state_degs_sig$p_val_adj<0.05,]
FBM_Higher_Sham<- gexpr_fbm_state_degs_sig[gexpr_fbm_state_degs_sig$avg_log2FC<0.25 &gexpr_fbm_state_degs_sig$p_val_adj<0.05,]
FBM_Higher_MCAO_top100 <- FBM_Higher_MCAO%>%top_n(n =100,wt=avg_log2FC)
FBM_Higher_Sham_top100 <- FBM_Higher_Sham%>%top_n(n =100,wt=abs(avg_log2FC))

degs_list=list('CBM:Higher in MCAO'=rownames(CBM_Higher_MCAO_top100),
               'CBM:Higher in Sham'=rownames(CBM_Higher_Sham_top100),
               'FBM:Higher in Sham'=rownames(FBM_Higher_MCAO_top100),
               'FBM:Higher in MCAO'=rownames(FBM_Higher_Sham_top100)
               
)
degs_list=list('CBM:Higher in MCAO'=rownames(CBM_Higher_MCAO),
               'CBM:Higher in Sham'=rownames(CBM_Higher_Sham),
               'FBM:Higher in Sham'=rownames(FBM_Higher_Sham),
               'FBM:Higher in MCAO'=rownames(FBM_Higher_MCAO)
               
)

pdf("Neut_venn.pdf",width = 8,height = 8)
venn::venn(degs_list,zcolor = c("#e3323c",'#5487c7','#5fb85e','#f58628'),ilcs = 0.9,opacity = 1)
dev.off()
degs_list2=list(both_up=intersect(degs_list[["CBM_up"]],degs_list[["TBM_up"]]),
                both_down=intersect(degs_list[["CBM_down"]],degs_list[["TBM_down"]]),
                CBM_up=setdiff(degs_list[["CBM_up"]],c(degs_list[["TBM_up"]],degs_list[["TBM_down"]])),
                CBM_down=setdiff(degs_list[["CBM_down"]],c(degs_list[["TBM_up"]],degs_list[["TBM_down"]])),
                TBM_down=setdiff(degs_list[["TBM_down"]],c(degs_list[["CBM_up"]],degs_list[["CBM_down"]])),
                TBM_up=setdiff(degs_list[["TBM_up"]],c(degs_list[["CBM_up"]],degs_list[["CBM_down"]]))
)

#degs_df=melt(degs_list2)
#write.csv(degs_df,'venn_degs_list2.csv',row.names = F)

#######################pathway enrichment

library(clusterProfiler)
library(enrichplot)
library(msigdbr)

##HALLMARK
m_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% dplyr::select(gs_name, gene_symbol)
m_t2g$gs_name<-gsub('HALLMARK_','',m_t2g$gs_name)

x <-compareCluster(degs_list, 'enricher', TERM2GENE = m_t2g)
dotplot(x, showCategory = 50, font.size = 10)+RotatedAxis()
ggsave("Both_CBM_TBM'_HALLMARK.pdf",width=20,height=9)

##GO BP
gobp = msigdbr(species = "Mus musculus",
               category = "C5",
               subcategory = 'BP') %>% dplyr::select(gs_name, gene_symbol)
gobp$gs_name=gsub('GOBP_','',gobp$gs_name)
x <- compareCluster(degs_list, 'enricher', TERM2GENE = gobp)
x2=x@compareClusterResult[!duplicated(x@compareClusterResult$ID),]
#saveRDS(x,'Both_CBM_TBM_GO.rds')
x3=x
x3@compareClusterResult=x2
dotplot(x3, showCategory = 10, font.size = 10, label_format = 100)
ggsave("cnet_Both_cd14_gngt2_GO_dotplot.pdf",width=12,height=16)
write.csv(x2,"全部亚型Both_CBM_TBM_GOunique_res.csv")      
p=cnetplot(x3)
p
ggsave(p,filename="cnet_Both_cd14_gngt2_GOunique.pdf",width=10,height=12)

####KEGG
kegg = msigdbr(species = "Mus musculus",
               category = "C2",
               subcategory = 'KEGG') %>% dplyr::select(gs_name, gene_symbol)
x <-compareCluster(degs_list, 'enricher', TERM2GENE = kegg)
dotplot(x, showCategory = 10, font.size = 10)
ggsave("Both_CBM_TBM_KEGG.pdf",width=20,height=9)
#####CC
gobp = msigdbr(species = "Mus musculus",
               category = "C5",
               subcategory = 'CC') %>% dplyr::select(gs_name, gene_symbol)
gobp$gs_name=gsub('GOCC_','',gobp$gs_name)
x <- compareCluster(degs_list, 'enricher', TERM2GENE = gobp)
x2=x@compareClusterResult[!duplicated(x@compareClusterResult$ID),]
#saveRDS(x,'Both_CBM_TBM_GO.rds')
x3=x
x3@compareClusterResult=x2
dotplot(x3, showCategory = 10, font.size = 10, label_format = 100)
ggsave("cnet_Both_cd14_gngt2_GO_dotplot.pdf",width=12,height=16)
