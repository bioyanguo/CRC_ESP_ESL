
# Loading env ---------------------------------------------------------------------
if(T){
  setwd("/home/lyg/lyg/CRC")
  set.seed(123456)
  library(igraph)
  library(reshape2)
  library(WGCNA)
  library(fdrtool)
  library(ggplot2)
  library(patchwork)
  library(Seurat)#4.0
  library(tidyverse)
  library(ggpubr)
  library(SeuratData)
  library(cowplot)
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(clustree)
  library(scDblFinder)
  library(BiocParallel)
  library(SingleCellExperiment)
  library(harmony)
  library(scater)
  library(CellChat)
  library(ggplot2)
  library(ggalluvial)
  library(NMF)
  library(monocle3)
  library(survival)
  library(survminer)
  library(dplyr)
  library(glmnet)
  library(GGally)
  library(rms)
  library(survivalROC)
  library(plotROC)
  library(grplasso)
  library(stringr)
  library(SCENIC)
  library(AUCell)
  library(GSEABase)
  library(DT)
  library(plotly)
  library(GEOquery)
  library(doMC);library(doRNG)
  library(ComplexHeatmap)
  library(pheatmap)
  library(FactoMineR)
  library(factoextra)
  library(infercnv)
  library(nichenetr)
  library(ggrepel)
  library(corrplot)
  library(GSVA)
  library(qusage)
  library(Hmisc)
  library(copykat)
  library(immunedeconv)
  library(fdrtool)
  library(WGCNA)
  library(reticulate)
}


#scRNA-seq GSE132465 process---------------------------------------------------------------

counts=fread("1rawdata/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt",data.table = T)
counts = column_to_rownames(counts,var = "Index")
counts[1:3,1:3];dim(counts);gc()
write.csv(rownames(counts),"counts_name_orig.csv")
library(clusterProfiler)
gene_id_scRNA<-bitr(rownames(counts),fromType = "SYMBOL", toType = c("ENSEMBL","ENTREZID"),OrgDb = org.Hs.eg.db)

nofound = sESPs$gene[which(! sESPs$gene %in% rownames(counts))]
nofound2 = pncRNAs$gene[which(!pncRNAs$gene %in% rownames(counts))]

gene_id<-bitr(nofound,fromType = "SYMBOL", toType = c("ENSEMBL","ENTREZID"),OrgDb = org.Hs.eg.db)
#repleace gene that change name
# "UTP11"="UTP11L"   "SPOUT1"="C9orf114" 
# "ATP5F1B"="ATP5B"  "RTF2"="RTFDC1"   "INTS11"="CPSF3L" 
# "YJU2"="CCDC94"    "RACK1"="GNB2L1"   "MTREX"="SKIV2L2" 
# "ATP5MF"="ATP5J2"   "BUD23"= "WBSCR22" 
gene_id<-bitr(nofound2,fromType = "SYMBOL", toType = c("ENSEMBL","ENTREZID"),OrgDb = org.Hs.eg.db)
gene_id = merge(gene_id,gene_id_scRNA,by = "ENSEMBL")


rownames(counts)[grep("GPRASP2",rownames(counts))]


meta=fread("1rawdata/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt")
meta = column_to_rownames(meta,var = "Index")
meta[1:3,1:3]
scRNA = CreateSeuratObject(counts = counts,meta.data = meta)
saveRDS(scRNA,"1rawdata/scRNA.RDS")

rm(counts,meta);gc()
##remove doublet by scDblFinder
Idents(scRNA)="orig.ident"
scRNA = as.SingleCellExperiment(scRNA)
scRNA <- scDblFinder(scRNA, samples="orig.ident"
);gc()
scRNA=as.Seurat(scRNA)
table(scRNA$scDblFinder.class)
scRNA@meta.data %>%
  ggplot(aes(x=orig.ident,fill=scDblFinder.class))+
  geom_bar()+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Ncell")

Idents(scRNA)="scDblFinder.class"
scRNA = subset(scRNA,idents="singlet")
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m = which(rownames(scRNA@assays$RNA) %in% HB.genes)
HB.genes = rownames(scRNA@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes)
scRNA[["percent.ribo"]] <- PercentageFeatureSet(scRNA, pattern = "^RP[SL]")
scRNA[["percent.HSP"]] <- PercentageFeatureSet(scRNA, pattern = "^HSP")
# Add number of genes per UMI for each cell to metadata
scRNA$log10GenesPerUMI <- log10(scRNA$nFeature_RNA) / log10(scRNA$nCount_RNA);gc()

##QC
# cell with >1,000 unique molecular identifier (UMI) counts;
# >200 genes and <6,000 genes;
# and <20% of mitochondrial gene expression in UMI counts.

length(which(scRNA$nFeature_RNA > 500 & scRNA$percent.mt < 20 & 
               scRNA$nCount_RNA > 500 & scRNA$percent.HB<0.1 & 
               scRNA$log10GenesPerUMI >0.78))
cell = rownames(scRNA@meta.data)[which(scRNA$nFeature_RNA > 500 & scRNA$percent.mt < 20 & 
                                          scRNA$nCount_RNA > 500 & scRNA$percent.HB<0.1 & 
                                          scRNA$log10GenesPerUMI >0.78)]
scRNA=subset(scRNA, cell=cell)
saveRDS(scRNA,"scRNA_afterqc.RDS")

gc();scRNA <- SCTransform(scRNA,seed.use = 123456,variable.features.n = 1500,
                          conserve.memory = T);gc()
length(scRNA@assays$SCT@var.features)#

scRNA <- RunPCA(scRNA,seed.use = 123456)
scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident",
                    assay.use="SCT");gc()
####
ElbowPlot(scRNA,50,reduction="harmony")
DimHeatmap(scRNA,dims = 1:50,reduction="pca",nfeatures = 10,cell = 5000)
pct<-scRNA[["harmony"]]@stdev/sum(scRNA[["harmony"]]@stdev)*100
cumu<-cumsum(pct)
co1<-which(cumu >80 & pct<5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2 <- min(co1,co2)
co2#
pc.num=1:33

gc()
scRNA = RunTSNE(scRNA, reduction="harmony", dims=pc.num,
                seed.use=123456) %>% 
  RunUMAP(reduction="harmony", dims=pc.num,seed.use=123456L)
scRNA =  FindNeighbors(scRNA ,reduction="harmony", dims=pc.num) %>% 
  FindClusters(resolution=c(0.4,0.8,1.2,1.6),
               random.seed=123456)

clus_tree = clustree(scRNA)+theme(legend.position = "bottom")+
  scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low = "grey80",high = "red");clus_tree

scRNA = FindClusters(scRNA,resolution=0.8,
                     random.seed=123456)

DimPlot(scRNA,reduction = "tsne",group.by = "Cell_type")
VlnPlot(scRNA,features = c("EPCAM","CDH1","PTPRC","PECAM1","COL1A1","CD79A","CD3D","LAMP3","CD14","TPSAB1","CLEC4C"),pt.size = 0.1)
Idents(scRNA)="seurat_clusters"
scRNA = RenameIdents(scRNA,"0"="Epithelial cell","1"="Epithelial cell","2"="T cell",
                     "3"="T cell","4"="T cell","5"="B cell",
                     "6"="T cell","7"="Myeloid cell","8"="T cell",
                     "9"="Fibroblast cell","10"="Myeloid cell","11"="T cell",
                     "12"="Epithelial cell","13"="T cell","14"="Epithelial cell",
                     "15"="Endothelial cell","16"="Epithelial cell","17"="Fibroblast cell","18"="B cell","19"="Epithelial cell",
                     "20"="B cell","21"="Myeloid cell","22"="Epithelial cell","23"="Epithelial cell","24"="Fibroblast cell","25"="Myeloid cell",
                     "26"="Fibroblast cell","27"="Mast cell","28"="T cell","29"="Epithelial cell")

scRNA$cell_type_my = Idents(scRNA)

scRNA <- CellCycleScoring(scRNA,s.features=cc.genes$s.genes,
                          g2m.features=cc.genes$g2m.genes)

saveRDS(scRNA, file="scRNA_finnal.Rds")
diff = FindAllMarkers(scRNA,logfc.threshold = 0.5,test.use = "MAST",only.pos = T)
write.csv(diff,"diff.all.major.csv")


## epi ---------------------------------------------------------------------

diff = read.csv("diff.all.major.csv",row.names = 1)
scRNA = readRDS("scRNA_finnal.Rds")
Epithelial = SCTransform(Epithelial)

Epithelial <- RunPCA(Epithelial,seed.use = 123456)
Epithelial <- RunHarmony(Epithelial, group.by.vars="orig.ident",
                         assay.use="SCT");gc()
####
ElbowPlot(Epithelial,50,reduction="harmony")
DimHeatmap(Epithelial,dims = 1:50,reduction="pca",nfeatures = 10,cell = 5000)
pct<-Epithelial[["harmony"]]@stdev/sum(Epithelial[["harmony"]]@stdev)*100
cumu<-cumsum(pct)
co1<-which(cumu >80 & pct<5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2 <- min(co1,co2)
co2#最终使用的维数
pc.num=1:co2

gc()
Epithelial = RunTSNE(Epithelial, reduction="harmony", dims=pc.num,
                     seed.use=123456) %>% 
  RunUMAP(reduction="harmony", dims=pc.num,seed.use=123456L)
Epithelial =  FindNeighbors(Epithelial ,reduction="harmony", dims=pc.num) %>% 
  FindClusters(resolution=c(0.1,0.2,0.3,0.4,0.8,1.2,1.6),
               random.seed=123456)

clus_tree = clustree(Epithelial)+theme(legend.position = "bottom")+
  scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low = "grey80",high = "red");clus_tree

Epithelial = FindClusters(Epithelial,resolution=0.1,
                          random.seed=123456)
saveRDS(Epithelial,"../4single cell/Epithelial.RDS")


#Smart-seq scRNA-seq GSE81861 process ---------------------------------------------------------------------

if(T){
  gtf = rtracklayer::import("../software/genome_reference_sources_human/gencode.v32.primary_assembly.annotation.gtf")
  gtf_df <- as.data.frame(gtf)
  geneid_df <- dplyr::select(gtf_df,c(gene_id,gene_name,gene_type))
  geneid_df$gene_id2 = strsplit2(geneid_df$gene_id,"[.]")[,1]
  geneid_df2 = geneid_df[,c(2,3,4)]
  geneid_df2 = dplyr::distinct(geneid_df2)
  #tumor epithelial
  GSE81861=fread("4single cell/GSE81861_CRC_tumor_epithelial_cells_COUNT.csv.gz")
  GSE81861_info = fread("4single cell/GSE81861_info_tumor.csv")
  GSE81861_info = GSE81861_info[,-2]
  GSE81861_info$ENSG2 = strsplit2(GSE81861_info$ENSG,"[.]")[,1]
  GSE81861 = merge(GSE81861_info,GSE81861,by.x="gene",by.y="V1")
  GSE81861 = merge(geneid_df2,GSE81861,by.x="gene_id2",by.y="ENSG2")
  table(GSE81861$gene_type)
  GSE81861 = GSE81861[GSE81861$gene_type %in% c("lncRNA","protein_coding"),]
  saveRDS(GSE81861,"4single cell/GSE81861_tumor.RDS")
  
  GSE81861 <- readRDS("/media/dell/Mybio/CRC/4single cell/GSE81861_tumor.RDS")
  GSE81861 = avereps(GSE81861,GSE81861$gene_name)
  GSE81861 = as.data.frame(GSE81861)
  rownames(GSE81861)=GSE81861$gene_name
  GSE81861 = GSE81861[,-c(1,2,3,4,5,6)]
  count = apply(GSE81861, 1, as.numeric)
  count = as.data.frame(t(count))
  colnames(count)=strsplit2(colnames(GSE81861),"[_]")[,1]
  saveRDS(count,"4single cell/GSE81861_tumor_count.RDS")
  
  
  #normal epithelial
  GSE81861=fread("4single cell/GSE81861_CRC_NM_epithelial_cells_COUNT.csv.gz")
  GSE81861_info = fread("4single cell/GSE81861_info_normal.csv")
  GSE81861_info$ENSG2 = strsplit2(GSE81861_info$ENSG,"[.]")[,1]
  GSE81861 = merge(GSE81861_info,GSE81861,by.x="gene",by.y="V1")
  GSE81861 = merge(geneid_df2,GSE81861,by.x="gene_id2",by.y="ENSG2")
  table(GSE81861$gene_type)
  GSE81861 = GSE81861[GSE81861$gene_type %in% c("lncRNA","protein_coding"),]
  saveRDS(GSE81861,"4single cell/GSE81861_normal.RDS")
  
  GSE81861 <- readRDS("/media/dell/Mybio/CRC/4single cell/GSE81861_normal.RDS")
  GSE81861 = avereps(GSE81861,GSE81861$gene_name)
  GSE81861 = as.data.frame(GSE81861)
  rownames(GSE81861)=GSE81861$gene_name
  GSE81861 = GSE81861[,-c(1,2,3,4,5,6)]
  count = apply(GSE81861, 1, as.numeric)
  count = as.data.frame(t(count))
  colnames(count)=strsplit2(colnames(GSE81861),"[_]")[,1]
  saveRDS(count,"4single cell/GSE81861_normal_count.RDS")
}

GSE81861_normal_count <- readRDS("/media/dell/Mybio/CRC/4single cell/GSE81861_normal_count.RDS")
GSE81861_tumor_count <- readRDS("/media/dell/Mybio/CRC/4single cell/GSE81861_tumor_count.RDS")
scRNA_tumor = CreateSeuratObject(counts = GSE81861_tumor_count)
scRNA_tumor$Class = "tumor"
scRNA_noraml = CreateSeuratObject(counts = GSE81861_normal_count)
scRNA_noraml$Class = "normal"
scRNA=merge(scRNA_tumor,scRNA_noraml)

scRNA = SCTransform(scRNA,variable.features.n = 1500)

saveRDS(scRNA,"4single cell/GSE81861_normal_tumor_scRNA.RDS")



#Gene effect and RRA -----------------------------
geneeffect = read.csv("0essential data/geneeffect_all.csv",header = T,row.names = 1,check.names = F)
geneeffect= geneeffect[,-1]
# geneeffect_CERES = read.csv("0essential data/CRISPR_gene_effect_CERES.csv",header = T,row.names = 1)
# geneeffect_CERES = t(geneeffect_CERES);geneeffect_CERES[1:3,1:3]
# rownames(geneeffect_CERES) = unlist(strsplit2(rownames(geneeffect_CERES),split = "..",fixed = T))[,1]
# geneeffect_CERES[1:3,1:3]
# write.csv(geneeffect_CERES,"0essential data/CRISPR_geneeffect_CERES.csv")
geneeffect_CRISPR = read.csv("0essential data/CRISPR_geneeffect_CERES.csv",header = T,row.names = 1,check.names = F)
geneeffect_CRISPR_info = read.csv("0essential data/sample_info_CRISPR.csv")
geneeffect_CRISPR_info = geneeffect_CRISPR_info[which(geneeffect_CRISPR_info$sample_collection_site=="Colon"),]
geneeffect_CRISPR = geneeffect_CRISPR[,which(colnames(geneeffect_CRISPR) %in% geneeffect_CRISPR_info$DepMap_ID)]
geneeffect = merge(geneeffect,geneeffect_CRISPR,by="row.names")
write.csv(geneeffect,"0essential data/geneeffect_colon.csv")
##
geneeffect = read.csv("0essential data/geneeffect_colon_finnal.csv",header = T,row.names = 1,check.names = F)

geneeffect = as.matrix(geneeffect)
d = list()
count=c()
for (i in 1:ncol(geneeffect)){
  b=rank(geneeffect[,i],na.last = "keep")
  #c=rank(-b,na.last = "keep")
  d[[i]]=names(sort(b))
  count[i]=length(d[[i]])
}
#rank list
library(RobustRankAggreg)
rankmat = rankMatrix(d,N=count)
ranks = aggregateRanks(rmat = rankmat)
ranks$adj_pvalue <- apply(cbind(ranks$Score*max(count),1),1,min)
results <- subset(ranks,ranks$adj_pvalue<0.01)
write.csv(ranks,"0essential data/RRA_ranks.csv")
write.csv(results,"0essential data/RRA_result.csv")


## determine the optimized threshold -----------------------------------------------------------------------

if(T){
  ranks = read.csv("0essential data/RRA_result.csv",header = T,stringsAsFactors = F,row.names = 1)
  
  Hart_ESGs=read.csv("0essential data/other essential gene/hart common_essentials.csv",header = T,stringsAsFactors = F)
  Hart_non_ESGs=read.csv("0essential data/other essential gene/Hartnonessentials.csv",header = T,stringsAsFactors = F)
  Achilles_ESGs = read.csv("0essential data/other essential gene/Achilles_common_essentials_CERES.csv",header = T,stringsAsFactors = F)
  zhang_ESGs = read.csv("0essential data/other essential gene/pan-cancer ESGs.csv")
  zhang_ESlncRNA = read.table("0essential data/other essential gene/pan-cancer ESlncRNAs.txt",header = T)
  ###################get iESGs###############
  A=as.character(rownames(ranks))
  B=unique(c(as.character(Achilles_ESGs$hgnc_symbol),
             as.character(Hart_ESGs$hgnc_symbol),
             zhang_ESGs$X.Gene_symbol))

  d=c()
  e=c()
  for (i in 1:100) {
    a= A[1:round(length(A)*i/100)]
    c=intersect(a,B)
    RA=length(c)/round(length(A)*i/100)
    RB=length(c)/length(B)
    d[i]=(RA*RB)^(0.5)
    e[i] = length(c)
  }
  #plot(d)
  #plot(e)
  write.csv(cbind(d,e),"0essential data/overlap.csv")
  ESP= A[1:round(length(A)*(which(d==max(d)))/100)] #i = 37 d=0.790894筛选阈值
  ESP= A[1:round(length(A)*0.37)] #i = 37 d=0.790894筛选阈值
}


## filter ESG with smart-seq scRNA-seq -----------------------------------------------------------------------

GSE81861_normal_tumor_scRNA <- readRDS("/media/dell/Mybio/CRC/4single cell/GSE81861_normal_tumor_scRNA.RDS")
Idents(GSE81861_normal_tumor_scRNA)="Class"

diff_TN = FindAllMarkers(GSE81861_normal_tumor_scRNA,only.pos = T,logfc.threshold = 0.5)
diff_TN_all = FindAllMarkers(GSE81861_normal_tumor_scRNA,only.pos = T,logfc.threshold = 0)

write.csv(diff_TN,"4single cell/GSE81861_diff_epithelial_T_N.csv")
write.csv(diff_TN_all,"4single cell/GSE81861_diff_epithelial_T_N_all.csv")

diff_TN = read.csv("4single cell/GSE81861_diff_epithelial_T_N.csv",row.names = 1,stringsAsFactors = T)
diff_TN_all = read.csv("4single cell/GSE81861_diff_epithelial_T_N_all.csv",row.names = 1,stringsAsFactors = F)

exp = AverageExpression(GSE81861_normal_tumor_scRNA,group.by = "Class",
                        features = CaseMatch(rownames(GSE81861_normal_tumor_scRNA),ESP),
                        slot = "data")$SCT
exp = as.data.frame(log2(exp+1))
exp$gene = rownames(exp)

GSE81861_normal_tumor_scRNA_sub=subset(GSE81861_normal_tumor_scRNA,idents="tumor")
counts <- GetAssayData(object = GSE81861_normal_tumor_scRNA_sub, slot = "counts") %>% as.data.frame()
counts <- counts[ESP,]
counts = na.omit(counts)
nonzero <- counts > 0
keep_genes <- rowSums(nonzero) > 272*0.5#50%
filtered_counts <- counts[keep_genes,]
ESP2 = rownames(filtered_counts)

exp_sub = exp[which(exp$gene %in% ESP2),]
exp_sub = merge(exp_sub,diff_TN_all,by="gene",all.x=T)
exp_sub$cluster=as.character(exp_sub$cluster)
exp_sub$cluster[which(!exp_sub$cluster %in% c("tumor","normal"))] = "Consistently"
exp_sub$cluster[which(!exp_sub$p_val_adj<0.01)] = "Consistently"
exp_sub$cluster[which(exp_sub$cluster=="tumor")] = "Over-expression in tumor Eps"

# Consistently Over-expression in tumor Eps 
# 322                           74
grep("CTNNB1",ESP2)

write.csv(exp_sub,"0essential data/ESP_finnal.csv")

ESP = read.csv("0essential data/ESP_finnal.csv",header = T,row.names = 1,stringsAsFactors = F)[,1]



# Co-expression network construction --------------------------------------------------------------------


##TCGA -------------------------------------------------------------------------


COAD_rpkm =read.csv("2TCGA/COAD_rpkm.csv",header = T,row.names = 1,check.names = F)
READ_rpkm =read.csv("2TCGA/READ_rpkm.csv",header = T,row.names = 1,check.names = F)
rpkm = merge(COAD_rpkm,READ_rpkm,by = "row.names")
rpkm = column_to_rownames(rpkm,var = "Row.names")
rpkm = t(rpkm)
rpkm = rpkm[-grep("-11",rownames(rpkm)),]#去除癌旁
library(WGCNA)
library(multtest)
library(stringr)
enableWGCNAThreads(32)

#weighted
TOM = TOMsimilarityFromExpr(
  rpkm,
  corType = "pearson",
  power = 6,
  TOMDenom = "mean",
  verbose = 1)

geneid_tcga = data.frame(symbol=colnames(rpkm),num=1:length(colnames(rpkm)))
write.csv(geneid_tcga,"5network/geneid_tcga.csv")

probes = 1:length(colnames(rpkm))
dimnames(TOM) <- list(probes, probes)
saveRDS(TOM,"2TCGA/TCGA_CRC_T_network_TOM.RDS")

#unweighted
expdata.d<-cor(rpkm,
               method="pearson")
colnames(expdata.d)<-colnames(rpkm)
rownames(expdata.d)<-colnames(rpkm)

dim(expdata.d)#13007 13007
len<-ncol(rpkm)
expdata.p<-corPvalueFisher(expdata.d,len)#p value
expdata.q<-p.adjust(expdata.p,method="fdr",n=length(expdata.p))# adj pvalue
rm(expdata.p,len)
expdata.d[expdata.q>0.01]=0
saveRDS(expdata.d,"2TCGA/TCGA_CRC_T_network_expdata.d.RDS")

#Get node of network
TCGA_CRC_T_network_TOM <- readRDS("/media/dell/Mybio/CRC/2TCGA/TCGA_CRC_T_network_TOM.RDS")
TCGA_CRC_T_network_expdata.d <- readRDS("/media/dell/Mybio/CRC/2TCGA/TCGA_CRC_T_network_expdata.d.RDS")

cyt_TCGA_T_weight = exportNetworkToCytoscape(
  TCGA_CRC_T_network_TOM,
  weighted = TRUE,
  threshold = 0
)
cyt_TCGA_T_weight=cyt_TCGA_T_weight$edgeData[,-c(5,6)]
saveRDS(cyt_TCGA_T_weight,"5network/cyt_TCGA_T_weight.RDS")

cyt_TCGA_T_unweight = exportNetworkToCytoscape(
  TCGA_CRC_T_network_expdata.d,
  weighted = TRUE,
  threshold = 0.3
)
cyt_TCGA_T_unweight=cyt_TCGA_T_unweight$edgeData[,-c(5,6)]
saveRDS(cyt_TCGA_T_unweight,"5network/cyt_TCGA_T_unweight.RDS")


## smart-seq scRNA-seq -------------------------------------------------------------------------

Epithelial <- readRDS("/media/dell/Mybio/CRC/4single cell/GSE81861_normal_tumor_scRNA.RDS")
Idents(Epithelial)="Class"
Epithelial <- subset(Epithelial,idents = "tumor")
exprMat <- GetAssayData(Epithelial,assay = "SCT",slot ="counts") %>% as.matrix()
mydbDIR="/media/dell/Mybio/software/cisTarget_databases/"
mydbs <- c("hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
           "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
names(mydbs) <- c("500bp","10kb")
scenicOptions <- initializeScenic(org= "hgnc", 
                                  nCores = 32,
                                  dbDir = mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "hg38"
)
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat <- GetAssayData(Epithelial,assay = "SCT",slot ="data") %>% as.matrix()
exprMat_filtered <- exprMat[genesKept, ];rm(exprMat)
exprMat_filtered = t(exprMat_filtered)
exprMat_filtered[1:3,1:10]
saveRDS(exprMat_filtered,"4single cell/GSE81861_tumor_eps_exprMat_filtered.RDS")
###weighted 
exprMat_filtered <- readRDS("/media/dell/Mybio/CRC/4single cell/GSE81861_tumor_eps_exprMat_filtered.RDS")
enableWGCNAThreads(32)
##
TOM = TOMsimilarityFromExpr(
  exprMat_filtered,
  corType = "pearson",
  power = 6,
  TOMDenom = "mean",
  verbose = 1)
#for PPR
geneid_sc = data.frame(symbol=colnames(exprMat_filtered),num=1:length(colnames(exprMat_filtered)))
write.csv(geneid_sc,"5network/geneid_sc.csv")
probes = 1:length(colnames(exprMat_filtered))

dimnames(TOM) <- list(probes, probes)
saveRDS(TOM,"4single cell/TOM_GSE81861_tumor_eps.RDS")

#unweighted 
expdata.d<-cor(exprMat_filtered,
               method="pearson")
colnames(expdata.d)<-colnames(exprMat_filtered)
rownames(expdata.d)<-colnames(exprMat_filtered)

dim(expdata.d)#14603 14603
len<-ncol(exprMat_filtered)
expdata.p<-corPvalueFisher(expdata.d,len)#
expdata.q<-p.adjust(expdata.p,method="fdr",n=length(expdata.p))# adj pvalue
rm(expdata.p,len)
expdata.d[expdata.q>0.01]=0#
saveRDS(expdata.d,"4single cell/expdata.d_GSE81861_tumor_eps.RDS")

#Get node of network
expdata.d_GSE81861_tumor_eps <- readRDS("/media/dell/Mybio/CRC/4single cell/expdata.d_GSE81861_tumor_eps.RDS")
TOM_GSE81861_tumor_eps <- readRDS("/media/dell/Mybio/CRC/4single cell/TOM_GSE81861_tumor_eps.RDS")

cyt_single_T_weight = exportNetworkToCytoscape(
  TOM_GSE81861_tumor_eps,
  weighted = TRUE,
  threshold = 0
)
cyt_single_T_weight=cyt_single_T_weight$edgeData[,-c(5,6)]
saveRDS(cyt_single_T_weight,"5network/cyt_single_T_weight.RDS")

cyt_single_T_unweight = exportNetworkToCytoscape(
  expdata.d_GSE81861_tumor_eps,
  weighted = TRUE,
  threshold = 0.3
)
cyt_single_T_unweight=cyt_single_T_unweight$edgeData[,-c(5,6)]
saveRDS(cyt_single_T_unweight,"5network/cyt_single_T_unweight.RDS")

# a = subset(cyt_single_T_unweight,cyt_single_T_unweight$weight>0.3)
# b = geneid_df2[geneid_df2$gene_name %in% unique(c(a$fromNode,a$toNode)),]
# table(b$gene_type)
# intersect(b$gene_name,ESP)

if(T){
  gtf = rtracklayer::import("../software/genome_reference_sources_human/gencode.v32.primary_assembly.annotation.gtf")
  gtf_df <- as.data.frame(gtf)
  geneid_df <- dplyr::select(gtf_df,c(gene_id,gene_name,gene_type))
  geneid_df$gene_id2 = strsplit2(geneid_df$gene_id,"[.]")[,1]
  geneid_df2 = geneid_df[,c(2,3,4)]
  geneid_df2 = dplyr::distinct(geneid_df2)
}


# HT ----------------------------------------------------------------------
# Direction: A-B and B-A 
cyt_TCGA_T_unweight <- readRDS("/media/dell/Mybio/CRC/5network/cyt_TCGA_T_unweight.RDS")
cyt_single_T_unweight <- readRDS("/media/dell/Mybio/CRC/5network/cyt_single_T_unweight.RDS")

cyt_single_T_unweight= rbind(cyt_single_T_unweight[,1:3],
                             data.frame(fromNode=cyt_single_T_unweight$toNode,
                                        toNode=cyt_single_T_unweight$fromNode,
                                        weight=cyt_single_T_unweight$weight))
cyt_TCGA_T_unweight=rbind(cyt_TCGA_T_unweight[,1:3],
                          data.frame(fromNode=cyt_TCGA_T_unweight$toNode,
                                     toNode=cyt_TCGA_T_unweight$fromNode,
                                     weight=cyt_TCGA_T_unweight$weight))

Hart_ESGs=read.csv("0essential data/other essential gene/hart common_essentials.csv",header = T,stringsAsFactors = F)
Hart_non_ESGs=read.csv("0essential data/other essential gene/Hartnonessentials.csv",header = T,stringsAsFactors = F)
Achilles_ESGs = read.csv("0essential data/other essential gene/Achilles_common_essentials_CERES.csv",header = T,stringsAsFactors = F)
zhang_ESGs = read.csv("0essential data/other essential gene/pan-cancer ESGs.csv")
zhang_ESlncRNA = read.table("0essential data/other essential gene/pan-cancer ESlncRNAs.txt",header = T)
ESP = read.csv("0essential data/ESP_finnal.csv",header = T,row.names = 1,stringsAsFactors = F)[,1]
randomgenes = sample(unique(cyt_single_T_unweight$fromNode),length(ESP));length(unique(randomgenes))

network = list(cyt_single_T_unweight,cyt_TCGA_T_unweight)
essentialgene_all=ESP
# HT_results_10 = list()
# roc1=list()
# roc2=list()
# for (k in 1:10) {
#   essentialgene_test = ESP[((k-1)*69):(k*69)]
#   essentialgene_all = ESP[-which(ESP %in% essentialgene_test)]#train
  HT_results = list()

  for(j in 1:2){
    if(T){
      table_pcgs = table(network[[j]]$fromNode)
      query_pcg = names(table_pcgs)
      network_inuse = split(network[[j]][,1:2],network[[j]]$fromNode)
      essentialgene = essentialgene_all
      # essentialgene = sample(essentialgene_all,round(length(essentialgene_all)*0.7,0),replace = F)
      # essentialgene_test = setdiff(essentialgene_all,essentialgene)
      gc()
      if(T){
        Nt=length(query_pcg)
        Np = length(essentialgene)
        phy = c()
        for (i in 1:length(query_pcg)){
          Nn = table_pcgs[i]
          Nc = length(intersect(essentialgene,
                                network_inuse[[i]]$toNode))
          phy1= phyper(Nc-1,Np,(Nt-Np),Nn,lower.tail=F)
          phy = c(phy,phy1)
        }
        query_pcg=data.frame(gene=query_pcg,pvalue = phy)
        fdr=fdrtool(query_pcg$pvalue,statistic="pvalue",plot=F)
        query_pcg$fdr = fdr$qval
        HT_p=query_pcg
        colnames(HT_p)[1]="hgnc_symbol"
        rm(query_pcg)
        HT_p$log = -log10(HT_p$fdr)#-log10(p)
        HT_p = HT_p[order(HT_p$log,decreasing = T),]
        HT_p$cumu<-cumsum(HT_p$log/sum(HT_p$log)*100)
      }
    }
    HT_results[[j]]=HT_p
  }

saveRDS(HT_results,"5network/ESP_results.RDS")

# PPR ---------------------------------------------------------------------
write.table(cyt_TCGA_T_weight[,c(1,2,3)],"8lncRNA/rwr_tcga.txt",col.names = F,row.names = F,quote = F)
write.table(cyt_single_T_weight[,c(1,2,3)],"8lncRNA/rwr_tsingle.txt",col.names = F,row.names = F,quote = F)

write.table(na.omit(geneid_sc[match(ESP,geneid_sc$symbol),2]),"8lncRNA/ESP_sc.txt",col.names = F,row.names = F,quote = F)
geneid1 = read.csv("5network/geneid_tcga.csv",row.names = 1)
write.table(na.omit(geneid1[match(ESP,geneid1$symbol),2]),"8lncRNA/ESP_tcga.txt",col.names = F,row.names = F,quote = F)

# cd /media/dell/Mybio/CRC/8lncRNA
# nohup pyrwr --query-type ppr --graph-type undirected \
#   --input-path rwr_tsingle.txt \
#   --output-path rwr_single_score.csv \
#   --seeds ESP_sc.txt \
#   --c 0.5 &

# nohup pyrwr --query-type ppr --graph-type undirected \
# --input-path rwr_tcga.txt \
# --output-path rwr_tcga_score.csv \
# --seeds ESP_tcga.txt \
# --c 0.5 &


# ESLs ESP expression in  GSE81861----------------------------------------------------------------

GSE81861_normal_tumor_scRNA <- readRDS("/media/dell/Mybio/CRC/4single cell/GSE81861_normal_tumor_scRNA.RDS")
Idents(GSE81861_normal_tumor_scRNA)="Class"
ESLs = read.csv("8lncRNA/ESLs.csv")[,2]

exp = AverageExpression(GSE81861_normal_tumor_scRNA,group.by = "Class",
                        features = CaseMatch(rownames(GSE81861_normal_tumor_scRNA),ESLs),slot = "data")$SCT
exp = as.data.frame(log2(exp+1))
exp$gene = rownames(exp)

GSE81861_normal_tumor_scRNA_sub=subset(GSE81861_normal_tumor_scRNA,idents="tumor")
counts <- GetAssayData(object = GSE81861_normal_tumor_scRNA_sub, slot = "counts") %>% as.data.frame()
counts <- counts[ESLs,]
counts = na.omit(counts)
nonzero <- counts > 0
keep_genes <- rowSums(nonzero) > 272*0.3
filtered_counts <- counts[keep_genes,]
ESP2 = rownames(filtered_counts)

exp_sub = exp[which(exp$gene %in% ESP2),]
exp_sub = merge(exp_sub,diff_TN,by="gene",all.x=T)
exp_sub$cluster=as.character(exp_sub$cluster)
exp_sub$cluster[which(!exp_sub$cluster %in% c("tumor","normal"))] = "Consistently"
exp_sub$cluster[which(!exp_sub$p_val_adj<0.01)] = "Consistently"
exp_sub$cluster[which(exp_sub$cluster=="tumor")] = "Over-expression in tumor Eps"

table(exp_sub$cluster)
# Consistently           tumor 
# 607                  11
grep("PVT1",ESP2)

write.csv(exp_sub,"8lncRNA/ESLs_exp_sub.csv")
write.csv(exp_sub$gene,"8lncRNA/ESLs.csv")
intersect(exp_sub$gene,zhang_ESlncRNA$hgnc_symbol)

diff_TN = read.csv("4single cell/GSE81861_diff_epithelial_T_N.csv",row.names = 1,stringsAsFactors = T)
diff_TN_all = read.csv("4single cell/GSE81861_diff_epithelial_T_N_all.csv",row.names = 1,stringsAsFactors = F)

diff_TN = diff_TN[which(diff_TN$gene %in% ESLs),]

write.csv(diff_TN,"8lncRNA/diff_TN_ESL.csv")


## function of ESLs-------------------------------------------------------------------------
ESLs = read.csv("8lncRNA/ESLs.csv")[,2]

#lnc2cancer
data=read.csv("8lncRNA/ESLs_function.csv",row.names = 1)
pheatmap(t(data),cluster_rows = F,color = c("white","red"),angle_col = 45)

#No cluster
data = data[sort(rownames(data)),]
pheatmap(t(data),cluster_rows = F,cluster_cols = F,color = c("white","red"),angle_col = 45)

#lncRNA properties
character = data[,c('CRC.related','Co.exist.zhang_s.ESLs','Coding.Ability')]
pheatmap(t(character),cluster_rows = F,cluster_cols = T,color = c("white","#375093"),angle_col = 45)

#lncRNA clinical function mechainism
data = data[,c("Cell.Growth",'Survival','Apoptosis','Metastasis','EMT',
               'Circulating','Recurrence','Immune')]
pheatmap(t(data),cluster_rows = F,cluster_cols = T,color = c("white","#A13D3B"),angle_col = 45)


#based on co-expression network
cyt_single_T_weight <- readRDS("/media/dell/Mybio/CRC/5network/cyt_single_T_weight.RDS")
cyt_single_T_weight = rbind(cyt_single_T_weight[,1:3],
                             data.frame(fromNode=cyt_single_T_weight$toNode,
                                        toNode=cyt_single_T_weight$fromNode,
                                        weight=cyt_single_T_weight$weight))
network_inuse = split(cyt_single_T_weight[,1:3],cyt_single_T_weight$fromNode)
network_ESL = network_inuse[which(names(network_inuse) %in% ESLs)]

for (i in names(network_ESL)) {
  network_ESL[[i]]=network_ESL[[i]][order(network_ESL[[i]]$weight,decreasing = T),]
  network_ESL[[i]]=network_ESL[[i]][1:200,]
}
saveRDS(network_ESL,"8lncRNA/network_ESL.RDS")

gobp = list()
for (i in names(network_ESL)) {
  # list=network_ESL[[i]][,2]
  # list = bitr(list,fromType = "SYMBOL", toType = c("ENSEMBL","ENTREZID"),OrgDb = org.Hs.eg.db)
  go <- enrichGO(network_ESL[[i]][,2],keyType = "SYMBOL",'org.Hs.eg.db',
                   pAdjustMethod="BH",ont = "BP")
  gobp[[i]]=subset(go@result,go@result$p.adjust<0.05)
}
saveRDS(gobp,"8lncRNA/gpbp_ESLs.RDS")
data = data.frame()
for (i in names(gobp)) {
  if(nrow(gobp[[i]])!=0){
    data = rbind(data,data.frame(i,gobp[[i]]))
  }
}
data$GeneRatio = gsub("/","//",data$GeneRatio)
data$BgRatio = gsub("/","//",data$BgRatio)
write.csv(data,"8lncRNA/gobp_ESLs.csv")
a = as.data.frame.array(table(data$Description))

data2 = data[grep("telomere maintenance",data$Description),]
write.table(data2,"8lncRNA/data2.txt")


# ESLs ESP expression in different cell lineages (GSE132465) ---------------------------------------

Hart_ESGs=read.csv("0essential data/other essential gene/hart common_essentials.csv",header = T,stringsAsFactors = F)$hgnc_symbol
Hart_non_ESGs=read.csv("0essential data/other essential gene/Hartnonessentials.csv",header = T,stringsAsFactors = F)$hgnc_symbol
Achilles_ESGs = read.csv("0essential data/other essential gene/Achilles_common_essentials_CERES.csv",header = T,stringsAsFactors = F)$gene
zhang_ESGs = read.csv("0essential data/other essential gene/pan-cancer ESGs.csv")$X.Gene_symbol
zhang_ESlncRNA = read.table("0essential data/other essential gene/pan-cancer ESlncRNAs.txt",header = T)$hgnc_symbol
ESP = read.csv("0essential data/ESP_finnal.csv",header = T,row.names = 1,stringsAsFactors = F)[,1]
ESLs = read.csv("8lncRNA/ESLs.csv")[,2]

scRNA=readRDS("4single cell/scRNA_finnal.Rds")

gene = list(ESGs=c(ESP,ESLs)[c(ESP,ESLs) %in% rownames(scRNA)],
            ESP=c(ESP)[c(ESP) %in% rownames(scRNA)],
            ESLs=c(ESLs)[c(ESLs) %in% rownames(scRNA)],
            Hart_ESGs=c(Hart_ESGs)[c(Hart_ESGs) %in% rownames(scRNA)],
            Hart_non_ESGs=c(Hart_non_ESGs)[c(Hart_non_ESGs) %in% rownames(scRNA)],
            Achilles_ESGs=c(Achilles_ESGs)[c(Achilles_ESGs) %in% rownames(scRNA)],
            zhang_ESGs=c(zhang_ESGs)[c(zhang_ESGs) %in% rownames(scRNA)],
            zhang_ESlncRNA=c(zhang_ESlncRNA)[c(zhang_ESlncRNA) %in% rownames(scRNA)]
            )

Idents(scRNA)="Cell_type"
scRNA = AddModuleScore(scRNA,features=gene)
p=MySeuratWrappers::VlnPlot(scRNA,features = c("Cluster2","Cluster3","Cluster4",
                                             "Cluster5","Cluster6","Cluster7","Cluster8"),
                          pt.size = 0,split.by = "Class",direction = c("horizontal"),stacked = T)+
  stat_compare_means(label = "p.signif")+
  geom_boxplot(width = 0.2,position = position_dodge(0.9),col="#563624",outlier.shape = NA)+
  stat_summary(position = position_dodge(0.9),color="yellow",size=0.05);p

ggsave2("Fig/Fig.2f.pdf",p,width = 9,height = 5.5)

Idents(scRNA)="Cell_subtype"
scRNA = AddModuleScore(scRNA,features=gene)
p=MySeuratWrappers::VlnPlot(scRNA,features = c("Cluster2","Cluster3","Cluster4",
                                               "Cluster5","Cluster6","Cluster7","Cluster8"),
                            pt.size = 0,split.by = "Class",direction = c("horizontal"),stacked = T)+
  stat_compare_means(label = "p.signif")+
  geom_boxplot(width = 0.2,position = position_dodge(0.9),col="#563624",outlier.shape = NA)+
  stat_summary(position = position_dodge(0.9),color="yellow",size=0.05);p

ggsave2("Fig/Fig.2f1.pdf",p,width = 8,height = 8)



scRNA_sub = subset(scRNA,cells=which(scRNA$Class=="Tumor"))
diff= FindMarkers(scRNA_sub,ident.1 = "Epithelial cell",logfc.threshold = 0.5,only.pos = T)
diff$gene=rownames(diff)
write.csv(diff,"4single cell/scRNA_Epithelial_diff.csv")
diffESP =diff[diff$gene %in% c(ESP),]
diffESLs =diff[diff$gene %in% c(ESLs),]
diffESG =diff[diff$gene %in% c(ESLs,ESP),]
diffESG$type = ifelse(diffESG$gene %in% c(ESP),"ESPs","ESLs")
diffESG$time = diffESG$pct.1/diffESG$pct.2
diffESG$type2=ifelse(diffESG$time>3,"specifically","non-specifically")
write.csv(diffESG,"diffESG_specifically.csv")
p=MySeuratWrappers::VlnPlot(scRNA_sub,c("KRT8","MYC","PRELID3B","SCD","TUBA1C"),
                          pt.size = 0,stacked = T,direction = c("horizontal"));p


# TCGA ESP mutation -------------------------------------------------------

ESP = read.csv("0essential data/ESP_finnal.csv",header = T,row.names = 1,stringsAsFactors = F)[,1]
ESLs = read.csv("8lncRNA/ESLs.csv")[,2]
#COAD
library(maftools)
COAD_laml = fread("2TCGA/COAD_maf/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf")
COAD_laml$Tumor_Sample_Barcode2 = COAD_laml$Tumor_Sample_Barcode
COAD_laml$Tumor_Sample_Barcode = str_sub(COAD_laml$Tumor_Sample_Barcode,1,16)
clinicalData = fread("2TCGA/TCGA-COAD.survival.tsv")
colnames(clinicalData)[1]="Tumor_Sample_Barcode"
COAD_laml = read.maf(COAD_laml,
                     clinicalData = clinicalData)
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
COAD_all= getGeneSummary(COAD_laml)
COAD= getGeneSummary(COAD_laml)[which(getGeneSummary(COAD_laml)$Hugo_Symbol %in% ESP),]#ESLs没有
getFields(COAD_laml)# check value
getClinicalData(COAD_laml)
getSampleSummary(COAD_laml)
plotmafSummary(maf = COAD_laml, rmOutlier = TRUE,
               addStat = 'median', dashboard = TRUE,
               titvRaw=FALSE)
#oncoplot
oncoplot(maf = COAD_laml, genes = COAD$Hugo_Symbol[1:20],
         fontSize =0.6,colors =vc_cols,legendFontSize = 1)#clinicalFeatures = 'OS',sortByAnnotation = TRUE

oncoplot(maf = COAD_laml, genes = unique(c(COAD$Hugo_Symbol[1:20],COAD_all$Hugo_Symbol[1:4]),
         fontSize =0.6,colors =vc_cols,legendFontSize = 1))

#plotTiTv
laml.titv = titv(maf = COAD_laml, plot = T, useSyn = TRUE)
par(mar=c(10,10,10,10))
plotTiTv(res = laml.titv,plotType = "both",
         sampleOrder = NULL,
         color = NULL,
         showBarcodes = FALSE,
         textSize = 0.8,
         baseFontSize = 1,
         axisTextSize = c(1, 1),
         plotNotch = FALSE)
#Rainfall plots
rainfallPlot(maf = COAD_laml, detectChangePoints = TRUE, pointSize = 0.6)

#Somatic Interactions
par(mar=c(4,8,8,4))
Interact <- somaticInteractions(COAD_laml, genes = unique(c(COAD$Hugo_Symbol[1:20],COAD_all$Hugo_Symbol[1:4])),
                                pvalue = c(0.01),
                                returnAll = TRUE,
                                fontSize = 0.6,
                                showSigSymbols = TRUE,
                                sigSymbolsFontSize = 0.5,
                                showCounts = FALSE,
                                countStats = "sig",
                                countType = "all",
                                countsFontColor = "black",
                                colPal = "BrBG")
#lollipop plot
lollipopPlot(maf = COAD_laml, gene = c("CTNNB1"), showMutationRate = TRUE)#AACol = 'Protein_Change'
saveRDS(COAD,"2TCGA/COAD_maf/COAD.RDS")


#READ
library(maftools)
READ_laml = fread("2TCGA/READ_maf/TCGA.READ.mutect.faa5f62a-2731-4867-a264-0e85b7074e87.DR-10.0.somatic.maf")
READ_laml$Tumor_Sample_Barcode2 = READ_laml$Tumor_Sample_Barcode
READ_laml$Tumor_Sample_Barcode = str_sub(READ_laml$Tumor_Sample_Barcode,1,16)
clinicalData = fread("2TCGA/TCGA-READ.survival.tsv")
colnames(clinicalData)[1]="Tumor_Sample_Barcode"
READ_laml = read.maf(READ_laml,
                     clinicalData = clinicalData)
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
READ_all= getGeneSummary(READ_laml)
READ= getGeneSummary(READ_laml)[which(getGeneSummary(READ_laml)$Hugo_Symbol %in% ESP),]#ESLs没有
getFields(READ_laml)# check value
getClinicalData(READ_laml)
getSampleSummary(READ_laml)
plotmafSummary(maf = READ_laml, rmOutlier = TRUE,
               addStat = 'median', dashboard = TRUE,
               titvRaw=FALSE)
#oncoplot
oncoplot(maf = READ_laml, genes = READ$Hugo_Symbol[1:20],
         fontSize =0.6,colors =vc_cols,legendFontSize = 1)#clinicalFeatures = 'OS',sortByAnnotation = TRUE

oncoplot(maf = READ_laml, genes = unique(c(READ$Hugo_Symbol[1:20],READ_all$Hugo_Symbol[1:4]),
                                         fontSize =0.6,colors =vc_cols,legendFontSize = 1))

#plotTiTv
laml.titv = titv(maf = READ_laml, plot = T, useSyn = TRUE)
par(mar=c(10,10,10,10))
plotTiTv(res = laml.titv,plotType = "both",
         sampleOrder = NULL,
         color = NULL,
         showBarcodes = FALSE,
         textSize = 0.8,
         baseFontSize = 1,
         axisTextSize = c(1, 1),
         plotNotch = FALSE)
#Rainfall plots
rainfallPlot(maf = READ_laml, detectChangePoints = TRUE, pointSize = 0.6)

#Somatic Interactions
par(mar=c(4,8,8,4))
Interact <- somaticInteractions(READ_laml, genes = unique(c(READ$Hugo_Symbol[1:20],READ_all$Hugo_Symbol[1:4])),
                                pvalue = c(0.01),
                                returnAll = TRUE,
                                fontSize = 0.6,
                                showSigSymbols = TRUE,
                                sigSymbolsFontSize = 0.5,
                                showCounts = FALSE,
                                countStats = "sig",
                                countType = "all",
                                countsFontColor = "black",
                                colPal = "BrBG")
#W*H = 5*5
#lollipop plot
lollipopPlot(maf = READ_laml, gene = c("CTNNB1"), showMutationRate = TRUE)#AACol = 'Protein_Change'
saveRDS(READ,"2TCGA/READ_maf/READ.RDS")

#combined TCGA
all = rbind(COAD_laml,READ_laml)
CRC_laml = read.maf(all)
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
CRC_all= getGeneSummary(CRC_laml)[which(getGeneSummary(CRC_laml)$Hugo_Symbol %in% ESP),]

plotmafSummary(maf = CRC_laml, rmOutlier = TRUE,
               addStat = 'median', dashboard = TRUE,
               titvRaw=FALSE)
#oncoplot
oncoplot(maf = CRC_laml, genes = CRC_all$Hugo_Symbol[1:20],
         fontSize =0.6,colors =vc_cols,legendFontSize = 1)#clinicalFeatures = 'OS',sortByAnnotation = TRUE

#Somatic Interactions
par(mar=c(4,8,8,4))
Interact <- somaticInteractions(CRC_laml, genes = unique(CRC_all$Hugo_Symbol[1:20]),
                                pvalue = c(0.01),
                                returnAll = TRUE,
                                fontSize = 0.6,
                                showSigSymbols = TRUE,
                                sigSymbolsFontSize = 0.5,
                                showCounts = FALSE,
                                countStats = "sig",
                                countType = "all",
                                countsFontColor = "black",
                                colPal = "BrBG")

# PCAWG ESP ESL mutation -------------------------------------------------------

ESP = read.csv("0essential data/ESP_finnal.csv",header = T,row.names = 1,stringsAsFactors = F)[,1]
ESLs = read.csv("8lncRNA/ESLs.csv")[,2]
#COAD ESP
library(maftools)
CRC_laml = fread("12PCAWG/ESPssnp.maf")
CRC_laml = read.maf(CRC_laml)
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
CRC_all= getGeneSummary(CRC_laml)

plotmafSummary(maf = CRC_laml, rmOutlier = TRUE,
               addStat = 'median', dashboard = TRUE,
               titvRaw=FALSE)
#oncoplot
oncoplot(maf = CRC_laml, genes = CRC_all$Hugo_Symbol[330:348],
         fontSize =0.6,colors =vc_cols,legendFontSize = 1)#clinicalFeatures = 'OS',sortByAnnotation = TRUE

#Somatic Interactions
par(mar=c(4,8,8,4))
Interact <- somaticInteractions(CRC_laml, genes = unique(CRC_all$Hugo_Symbol[1:20]),
                                pvalue = c(0.01),
                                returnAll = TRUE,
                                fontSize = 0.6,
                                showSigSymbols = TRUE,
                                sigSymbolsFontSize = 0.5,
                                showCounts = FALSE,
                                countStats = "sig",
                                countType = "all",
                                countsFontColor = "black",
                                colPal = "BrBG")


#COAD ESLs
library(maftools)
CRC_laml = fread("12PCAWG/ESLcRNAsssnp.maf")
CRC_laml = read.maf(CRC_laml)
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
CRC_all= getGeneSummary(CRC_laml)

plotmafSummary(maf = CRC_laml, rmOutlier = TRUE,
               addStat = 'median', dashboard = TRUE,
               titvRaw=FALSE)#绘制整体的突变情况
#oncoplot
oncoplot(maf = CRC_laml, genes = CRC_all$Hugo_Symbol[1:20],
         fontSize =0.6,colors =vc_cols,legendFontSize = 1)#clinicalFeatures = 'OS',sortByAnnotation = TRUE



# Methylation -------------------------------------------------------------

ESP = read.csv("0essential data/ESP_finnal.csv",header = T,row.names = 1,stringsAsFactors = F)[,1]
ESLs = read.csv("8lncRNA/ESLs.csv")[,2]
#COAD
library(ChAMP)
library(impute)
methy = fread("2TCGA/TCGA-COAD.methylation450.tsv.gz",data.table = F,check.names = F)
rownames(methy) = methy[,1];methy=methy[,-1]
methy[1:4,1:4]

pd = fread("2TCGA/TCGA-COAD.GDC_phenotype.tsv.gz")
colnames(pd)[1]="Sample_Name"
colnames(pd)[116]="Sample_Group"
pd = as.data.frame(pd)
rownames(pd)=pd$Sample_Name
pd = pd[which(pd$Sample_Group %in% c("Primary Tumor","Solid Tissue Normal")),]
methy = methy[,which(colnames(methy) %in% pd$Sample_Name)]
pd=pd[colnames(methy),]
identical(rownames(pd),colnames(methy))
betaData =as.matrix(na.omit(methy))
# or full NA
# beta = as.matrix(methy)
# beta = impute.knn(beta)
# betaData = beta$data
# betaData = betaData+0.00001

myLoad = champ.filter(beta = betaData,pd=pd)
rm(betaData,beta)
dim(myLoad$beta)
champ.QC(beta = myLoad$beta, pheno = rownames(myLoad$pd))
library(doParallel)
detectCores()
myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=16)

saveRDS(myNorm,file = "10methy/COAD_myNorm.RDS")
saveRDS(myLoad,file = "10methy/COAD_myLoad.RDS")
#diff analysis
myLoad$pd$Sample_Group=ifelse(myLoad$pd$Sample_Group=="Primary Tumor","Tumor","Normal")
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Sample_Group)
head(myDMP[[1]])
df_DMP <- myDMP$Normal_to_Tumor
write.csv(df_DMP,"CRC/10methy/COAD_Normal_to_Tumor_DMP.csv")

n = which(as.character(probe.features$gene) %in% as.character(c(ESLs,ESP)))
ESG_prob = probe.features[as.numeric(n),];length(unique(ESG_prob$gene))
ESG_prob_beta = myNorm[which(rownames(myNorm) %in% rownames(ESG_prob)),]
ESG_prob = ESG_prob[which(rownames(ESG_prob) %in% rownames(ESG_prob_beta)),]
df_DMP = df_DMP[which(rownames(df_DMP) %in% rownames(ESG_prob)),]

write.csv(ESG_prob,"10methy/COAD_ESLs_ESP_methy_prob.csv")
write.csv(ESG_prob_beta,"10methy/COAD_ESLs_ESP_methy_beta.csv")
write.csv(df_DMP,"10methy/COAD_ESLs_ESP_Normal_to_Tumor_DMP.csv")

df_DMP = read.csv("10methy/COAD_ESLs_ESP_Normal_to_Tumor_DMP.csv",row.names = 1)

myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Sample_Group,
                   method="Bumphunter")

unique(df_DMP$gene) %in% unique(c(ESLs,ESP))

vc_cols = RColorBrewer::brewer.pal(n = 6, name = "Spectral")
df_DMP$genetype = ifelse(df_DMP$gene %in% ESLs,"lncRNA","PCG")
df_DMP$group = ifelse(df_DMP$logFC<0.2,ifelse(df_DMP$logFC>(-0.2),"not","hypo"),"hyper")
df_DMP$group2 = ifelse(df_DMP$logFC<0,"hypo","hyper")


df_DMP$label=""
label_num = match(c("cg19050555","cg26036626","cg19508622","cg10441379","cg15677364"),
                  rownames(df_DMP))
df_DMP$label[label_num]=rownames(df_DMP)[label_num]
# FBLIM1 MCM7 PVT1

p = ggplot(data=df_DMP,
       aes(x=logFC, y=-log10(adj.P.Val),color=feature)) +
  geom_point(alpha=0.6, size=1.5) +
  theme_set(theme_set(theme_bw(base_size=10)))+
  xlab("Mean difference methylation") + ylab("-log10(adjusted p-value)") +theme_bw()+
  scale_color_manual(values=vc_cols)+
  geom_vline(xintercept = 0.2,colour="red",size=0.1)+
  geom_vline(xintercept = -0.2,colour="red",size=0.1)+
  geom_text_repel(aes(label = label),max.overlaps = 300,
                size = 4,box.padding = unit(1, "lines"),
                point.padding = unit(1, "lines"),
                segment.color = "black",
                show.legend = FALSE);p
ggsave2("Fig/Fig.3e.pdf",p,height = 4,width = 5,device = "pdf")

a = as.data.frame.array(table(df_DMP$feature))
a$feature = rownames(a)
a$count = a$`table(df_DMP$feature)`
p = ggplot(data=a,
       aes(x=feature,y=count,fill=feature))+geom_bar(stat = "identity")+
  scale_fill_manual(values=vc_cols)+theme_classic2()
ggsave2("Fig/Fig.3f.pdf",p,height = 2,width = 5,device = "pdf")


# ESP were highly expressed, so we need to force hypomethylation. 
# hypermethylation, not likely happend on this genes.
df_DMP_orig=df_DMP
df_DMP = subset(df_DMP_orig,df_DMP$logFC<(-0.2))#

rpkm =read.csv("2TCGA/COAD_rpkm.csv",header = T,row.names = 1,check.names = F)
rpkm = rpkm[,which(colnames(rpkm) %in% colnames(ESG_prob_beta))]
beta = ESG_prob_beta[rownames(df_DMP),colnames(rpkm)]
identical(colnames(rpkm),colnames(beta))
group=rep("G",ncol(rpkm))
group[grep("-11",colnames(rpkm))]="Normal"
group[grep("-01",colnames(rpkm))]="Tumor"

c=data.frame()
for (i in rownames(beta)) {
  gene = as.character(ESG_prob[i,5])
  if(!is.na(rpkm[gene,])){
    a = data.frame(sample = names(beta[i,]),
                   sample2 = names(rpkm[gene,]),
                   prob_beta = as.numeric(beta[i,]),
                   gene_rpkm = as.numeric(rpkm[gene,]),
                   group=group)
    max = median(a$prob_beta)+3*mad(a$prob_beta)
    min = median(a$prob_beta)-3*mad(a$prob_beta)
    a = subset(a,a$prob_beta>min & a$prob_beta<max)
    b = cor.test(a$prob_beta,a$gene_rpkm,method = "pearson")
    c=rbind(c,data.frame(gene = gene,
                         prob=i,
                         cor=b[["estimate"]][["cor"]],
                         p=b$p.value))
  }
}
library(fdrtool)
c$adj.p=fdrtool(c$p,statistic = "pvalue",plot = F)$qval

for (i in c$prob) {
  i="cg23898497"
  gene = as.character(ESG_prob[i,5])
  a = data.frame(sample = names(beta[i,]),
                 sample2 = names(rpkm[gene,]),
                 prob_beta = as.numeric(beta[i,]),
                 gene_rpkm = as.numeric(rpkm[gene,]),
                 group=group)
  max = median(a$prob_beta)+3*mad(a$prob_beta)
  min = median(a$prob_beta)-3*mad(a$prob_beta)
  a = subset(a,a$prob_beta>min & a$prob_beta<max)
  p=ggplot(a) + geom_point(aes(x=prob_beta,y=gene_rpkm,col=group))+
    geom_smooth(aes(x=prob_beta,y=gene_rpkm),method = "lm",col="red")+
    xlab(paste0("Beta value of ",i))+ylab(paste0(gene," expression levels"))+
    theme_bw()+theme(axis.text=element_text(colour = 'black',size = 12),
                     axis.title=element_text(size = 14),
                     axis.text.x = element_text(size = 14),
                     axis.text.y = element_text(size = 14),
                     panel.border = element_blank(),
                     axis.line = element_line(colour = "black",size=0.15),
                     #legend.position = "none",
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank());p
  ggsave2(paste0("10methy/",i,".pdf"),p,height = 4,width = 5)
}
data = c[c$p<0.05&abs(c$cor)>0.3,]
df_DMP_data = df_DMP[which(rownames(df_DMP) %in% data$prob),]
df_DMP_data = merge(data,df_DMP_data,by.x= "prob",by.y="row.names")

write.csv(df_DMP_data,"10methy/df_DMP_data_COAD.csv")

#READ
ESP = read.csv("0essential data/ESP_finnal.csv",header = T,row.names = 1,stringsAsFactors = F)[,1]
ESLs = read.csv("8lncRNA/ESLs.csv")[,2]

library(ChAMP)
library(impute)
methy = fread("2TCGA/TCGA-READ.methylation450.tsv.gz",data.table = F,check.names = F)
rownames(methy) = methy[,1];methy=methy[,-1]
methy[1:4,1:4]

pd = fread("2TCGA/TCGA-READ.GDC_phenotype.tsv.gz")
colnames(pd)[1]="Sample_Name"
colnames(pd)[114]="Sample_Group"
pd = as.data.frame(pd)
rownames(pd)=pd$Sample_Name
pd = pd[which(pd$Sample_Group %in% c("Primary Tumor","Solid Tissue Normal")),]
methy = methy[,which(colnames(methy) %in% pd$Sample_Name)]
pd=pd[colnames(methy),]
identical(rownames(pd),colnames(methy))
betaData =as.matrix(na.omit(methy))
# or full NA
# beta = as.matrix(methy)
# beta = impute.knn(beta)
# betaData = beta$data
# betaData = betaData+0.00001
#过滤
myLoad = champ.filter(beta = betaData,pd=pd)
dim(myLoad$beta)
champ.QC(beta = myLoad$beta, pheno = rownames(myLoad$pd))
library(doParallel)
detectCores()
myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=16)

saveRDS(myNorm,file = "10methy/READ_myNorm.RDS")
saveRDS(myLoad,file = "10methy/READ_myLoad.RDS")
#diff analysis
myLoad$pd$Sample_Group=ifelse(myLoad$pd$Sample_Group=="Primary Tumor","Tumor","Normal")
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Sample_Group)
head(myDMP[[1]])
df_DMP <- myDMP$Tumor_to_Normal
write.csv(df_DMP,"10methy/READ_Normal_to_Tumor_DMP.csv")

n = which(as.character(probe.features$gene) %in% as.character(c(ESLs,ESP)))
ESG_prob = probe.features[as.numeric(n),];length(unique(ESG_prob$gene))
ESG_prob_beta = myNorm[which(rownames(myNorm) %in% rownames(ESG_prob)),]
ESG_prob = ESG_prob[which(rownames(ESG_prob) %in% rownames(ESG_prob_beta)),]
df_DMP = df_DMP[which(rownames(df_DMP) %in% rownames(ESG_prob)),]

df_DMP$genetype = ifelse(df_DMP$gene %in% ESLs,"lncRNA","PCG")
df_DMP$group = ifelse(df_DMP$logFC<0.2,ifelse(df_DMP$logFC>(-0.2),"not","hypo"),"hyper")
df_DMP$group2 = ifelse(df_DMP$logFC<0,"hypo","hyper")

write.csv(ESG_prob,"10methy/READ_ESLs_ESP_methy_prob.csv")
write.csv(ESG_prob_beta,"10methy/READ_ESLs_ESP_methy_beta.csv")
write.csv(df_DMP,"10methy/READ_ESLs_ESP_Normal_to_Tumor_DMP.csv")

df_DMP = read.csv("10methy/READ_ESLs_ESP_Normal_to_Tumor_DMP.csv",row.names = 1)

myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Sample_Group,
                   method="Bumphunter")

unique(df_DMP$gene) %in% unique(c(ESLs,ESP))

vc_cols = RColorBrewer::brewer.pal(n = 6, name = "Spectral")

df_DMP$label=""
label_num = match(c("cg09678212","cg08526705"),
                  rownames(df_DMP))
df_DMP$label[label_num]=rownames(df_DMP)[label_num]
# FBLIM1 MCM7 PVT1

p = ggplot(data=df_DMP,
           aes(x=logFC, y=-log10(adj.P.Val),color=feature)) +
  geom_point(alpha=0.6, size=1.5) +
  theme_set(theme_set(theme_bw(base_size=10)))+
  xlab("Mean difference methylation") + ylab("-log10(adjusted p-value)") +theme_bw()+
  scale_color_manual(values=vc_cols)+
  geom_vline(xintercept = 0.2,colour="red",size=0.1)+
  geom_vline(xintercept = -0.2,colour="red",size=0.1)+
  geom_text_repel(aes(label = label),max.overlaps = 300,
                  size = 4,box.padding = unit(1, "lines"),
                  point.padding = unit(1, "lines"),
                  segment.color = "black",
                  show.legend = FALSE);p

a = as.data.frame.array(table(df_DMP$feature))
a$feature = rownames(a)
a$count = a$`table(df_DMP$feature)`
p = ggplot(data=a,
           aes(x=feature,y=count,fill=feature))+geom_bar(stat = "identity")+
  scale_fill_manual(values=vc_cols)+theme_classic2()
ggsave2("Fig/supplementary Fig.3f.pdf",p,height = 2,width = 5,device = "pdf")


# ESP were highly expressed, so we need to force hypomethylation. 
# hypermethylation, not likely happend on this genes.
df_DMP_orig=df_DMP
df_DMP = subset(df_DMP_orig,df_DMP$logFC<(-0.2))#

rpkm =read.csv("2TCGA/READ_rpkm.csv",header = T,row.names = 1,check.names = F)
rpkm = rpkm[,which(colnames(rpkm) %in% colnames(ESG_prob_beta))]
beta = ESG_prob_beta[rownames(df_DMP),colnames(rpkm)]
identical(colnames(rpkm),colnames(beta))
group=rep("G",ncol(rpkm))
group[grep("-11",colnames(rpkm))]="Normal"
group[grep("-01",colnames(rpkm))]="Tumor"

c=data.frame()
for (i in rownames(beta)) {
  gene = as.character(ESG_prob[i,5])
  if(!is.na(rpkm[gene,])){
    a = data.frame(sample = names(beta[i,]),
                   sample2 = names(rpkm[gene,]),
                   prob_beta = as.numeric(beta[i,]),
                   gene_rpkm = as.numeric(rpkm[gene,]),
                   group=group)
    max = median(a$prob_beta)+3*mad(a$prob_beta)
    min = median(a$prob_beta)-3*mad(a$prob_beta)
    a = subset(a,a$prob_beta>min & a$prob_beta<max)
    b = cor.test(a$prob_beta,a$gene_rpkm,method = "pearson")
    c=rbind(c,data.frame(gene = gene,
                         prob=i,
                         cor=b[["estimate"]][["cor"]],
                         p=b$p.value))
  }
}
library(fdrtool)
c$adj.p=fdrtool(c$p,statistic = "pvalue",plot = F)$qval

for (i in c$prob) {
  # i="cg23898497"
  gene = as.character(ESG_prob[i,5])
  a = data.frame(sample = names(beta[i,]),
                 sample2 = names(rpkm[gene,]),
                 prob_beta = as.numeric(beta[i,]),
                 gene_rpkm = as.numeric(rpkm[gene,]),
                 group=group)
  max = median(a$prob_beta)+3*mad(a$prob_beta)
  min = median(a$prob_beta)-3*mad(a$prob_beta)
  a = subset(a,a$prob_beta>min & a$prob_beta<max)
  p=ggplot(a) + geom_point(aes(x=prob_beta,y=gene_rpkm,col=group))+
    geom_smooth(aes(x=prob_beta,y=gene_rpkm),method = "lm",col="red")+
    xlab(paste0("Beta value of ",i))+ylab(paste0(gene," expression levels"))+
    theme_bw()+theme(axis.text=element_text(colour = 'black',size = 12),
                     axis.title=element_text(size = 14),
                     axis.text.x = element_text(size = 14),
                     axis.text.y = element_text(size = 14),
                     panel.border = element_blank(),
                     axis.line = element_line(colour = "black",size=0.15),
                     #legend.position = "none",
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank());p
  ggsave2(paste0("10methy/",i,".pdf"),p,height = 4,width = 5)
}
data = c[c$p<0.05&abs(c$cor)>0.3,]
df_DMP_data = df_DMP[which(rownames(df_DMP) %in% data$prob),]
df_DMP_data = merge(data,df_DMP_data,by.x= "prob",by.y="row.names")

write.csv(df_DMP_data,"10methy/df_DMP_data_READ.csv")


# SCENIC and TF based on GSE132465 ----------------------------------------

ESP = read.csv("0essential data/ESP_finnal.csv",header = T,row.names = 1,stringsAsFactors = F)[,1]
ESLs = read.csv("8lncRNA/ESLs.csv")[,2]

# dir.create("7_scenic")
setwd("7_scenic")
aim_cell <- readRDS("/media/dell/Mybio/CRC/4single cell/GSE81861_normal_tumor_scRNA.RDS")
# Idents(aim_cell)="Class"
# aim_cell = subset(aim_cell,idents = "tumor")
exprMat <- GetAssayData(aim_cell,assay = "SCT",slot = "data") %>% as.matrix()


mydbDIR="/media/dell/Mybio/software/cisTarget_databases/"
mydbs <- c("hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
           "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")

names(mydbs) <- c("500bp","10kb")
scenicOptions <- initializeScenic(org= "hgnc", 
                                  nCores = 24,
                                  dbDir = mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "hg38"
)
scenicOptions@settings$verbose = TRUE
scenicOptions@settings$seed = 123456
saveRDS(scenicOptions,file = "output/scenicOptions.RDs")

### Co-expression network
genesKept <- geneFiltering(GetAssayData(aim_cell,
                                        assay = "SCT",slot = "counts") %>% as.matrix(),
                           scenicOptions)
interestingGenes <- c("SOX9", "SOX10", "PROM1")#check any gene that interesting missing?
interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_filtered <- exprMat[genesKept, ]
rm(exprMat);gc()

runCorrelation(exprMat_filtered, scenicOptions)

runGenie3(exprMat_filtered, scenicOptions,verbose=TRUE)

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)

detach("package:NMF", unload = TRUE)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
exprMat_all <-  GetAssayData(aim_cell,assay = "SCT",slot = "data") %>% as.matrix()
genesKept <- geneFiltering(exprMat_all, scenicOptions)
exprMat_all <- exprMat_all[genesKept, ]
gc()

saveRDS(aim_cell@meta.data, file="int/cellInfo.Rds")
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"

scenicOptions@settings$nCores=1
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
saveRDS(scenicOptions,file = "output/scenicOptions.RDs")


regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
regulonTargetsInfo = regulonTargetsInfo[which(regulonTargetsInfo$TF %in% ESP | regulonTargetsInfo$gene %in% c(ESLs,ESP)),]
a = as.data.frame.array(table(regulonTargetsInfo$TF))
a1 = as.data.frame.array(table(regulonTargetsInfo$gene))
write.csv(regulonTargetsInfo,"regulonTargetsInfo_ESP_ESLs.csv")

tableSubset <- regulonTargetsInfo[TF %in% c(ESP)]#
c = as.data.frame.array(table(tableSubset$TF,tableSubset$highConfAnnot))
c$TF=rownames(c)
c = c[order(c[,2],decreasing = T),]
c$TF = factor(c$TF,levels=rownames(c))
c1 = melt(c)
p = ggplot(c1,aes(x=TF,y=value,fill=variable))+
  geom_col()+theme_bw()+xlab("The number of relationships")+
  RotatedAxis();p
ggsave2("../../Fig/Fig.4a.pdf",p,height = 3,width = 4,device = "pdf")

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC = getAUC(regulonAUC)
regulonAUC<- data.frame(t(regulonAUC), check.names=F)
TF = strsplit2(colnames(regulonAUC),"[ ]|_")[,1]
regulonAUC=regulonAUC[,TF %in% regulonTargetsInfo$TF]
TF = unique(strsplit2(colnames(regulonAUC),"[ ]|_")[,1])

saveRDS(regulonAUC,"regulonAUC_finnal.RDS")

##DE regulon
aim_cell <- readRDS("/media/dell/Mybio/CRC/7_scenic/TandN/aim_cell.RDS")
regulonAUC <- readRDS("/media/dell/Mybio/CRC/7_scenic/TandN/regulonAUC_finnal.RDS")

cellInfo = aim_cell@meta.data
regulonAUC$group = cellInfo[rownames(regulonAUC),4]

diff2 = data.frame()
for (i in 1:(ncol(regulonAUC)-1)) {
  a = t.test(regulonAUC[,i]~group,regulonAUC)
  diff2 = rbind(diff2,
                data.frame(regulon = colnames(regulonAUC)[i],
                           p=a$p.value,
                           avgN=a[["estimate"]][["mean in group normal"]],
                           avgT=a[["estimate"]][["mean in group tumor"]]))
}
diff2$diff = diff2$avgT-diff2$avgN
diff2$group = ifelse(diff2$diff>0.01 & diff2$p<0.05,"activated in tumor Eps",
                     ifelse(diff2$diff<(-0.01) & diff2$p<0.05,"activated in normal Eps","Not"))
diff2$label=""
label_num = match(c("MYC_extended (5417g)","YY1 (2508g)","YY1_extended (4238g)","MYC (334g)","BCLAF1_extended (6775g)","BCLAF1 (3455g)",
                    "CDX1 (31g)","TCF7L2 (71g)"),
                  diff2$regulon)
diff2$label[label_num]=diff2$regulon[label_num]

p = ggplot(diff2,aes(x=diff,y=-log10(p),col=group))+geom_point()+
  scale_color_manual(values = c("#00BFC4","#F8766D","grey"))+
  geom_text_repel(aes(label = label),max.overlaps = 300,
                  size = 4,box.padding = unit(1, "lines"),
                  point.padding = unit(1, "lines"),
                  segment.color = "black",
                  show.legend = FALSE)+
  theme_bw()+coord_flip();p

ggsave2("../../Fig/Fig.4b.pdf",p,width = 6,height = 4.5)
write.csv(diff2,"diff2_regulon_AUC.csv")

diff2_TF = strsplit2(diff2$regulon[which(diff2$group=="activated in tumor Eps")],"[ ]|_")[,1]


library(RcisTarget)
data(motifAnnotations_hgnc)
TF= intersect(unique(motifAnnotations_hgnc$TF),c(ESP,ESLs))# 59个

# output/Step2_MotifEnrichment_preview.html in detail/subset:
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")


# output/Step2_regulonTargetsInfo.tsv in detail: 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")

aim_cell <- readRDS("/media/dell/Mybio/CRC/7_scenic/TandN/aim_cell.RDS")
for (i in c$TF) {
  x = grep(paste0(i),colnames(aim_cell@meta.data))
  y = colnames(aim_cell@meta.data)[x]
  for (k in y) {
    p = FeaturePlot(aim_cell, features=k, label=F, reduction = 'tsne')+
      theme(text = element_text(size=10),axis.text = element_text(size=10))+
      scale_color_distiller(palette = "Spectral")+coord_fixed();p
    ggsave2(paste0("tmp/Fig.4",k,".pdf"),p,height = 3.7,width = 4.7)
  }
}
p=DimPlot(aim_cell,reduction = 'tsne',label=T,group.by = "Class")
ggsave2("../../Fig/supplmentary Fig.4a.pdf",p,height = 3.7,width = 4.7)

tableSubset2 = regulonTargetsInfo[order(regulonTargetsInfo$CoexWeight,
                                        decreasing = T),]
tableSubset2 = tableSubset2[tableSubset2$gene %in% c("CTNNB1", "NDUFB9","SCNM1")]

write.csv(tableSubset2,"tableSubset_CRC_specific.csv")

#ESLs
# tableSubset <- regulonTargetsInfo[TF %in% ESLs]
tableSubset <- regulonTargetsInfo[gene %in% c(ESLs)]
a = as.data.frame.array(table(tableSubset$TF))
b = as.data.frame.array(table(tableSubset$gene,tableSubset$highConfAnnot))
b$gene = rownames(b)
b = b[order(b[,2],decreasing = T),]
level=unique(as.character(b$gene))
b =melt(b)
b$gene =factor(b$gene,levels = level)
p = ggplot(b,aes(x=gene,y=value,fill=variable))+
  geom_col()+theme_bw()+RotatedAxis();p
ggsave2("../../Fig/Fig.4d.pdf",p,height = 3.5,width =7)

tableSubset = tableSubset[highConfAnnot == "TRUE"]
write.csv(tableSubset,"tableSubset_ESLs.csv")


# visualization -----------------------------------------------------------

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")#int/3.4_regulonAUC.Rds
regulonAUC = getAUC(regulonAUC)#@assays@data@listData$AUC
regulonAUC <- data.frame(t(regulonAUC), check.names=F)
RegulonName_AUC <- gsub(' \\(','_',colnames(regulonAUC))
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(regulonAUC) <- RegulonName_AUC

aim_cell <- AddMetaData(aim_cell, regulonAUC)
aim_cell =RunPCA(aim_cell)
ElbowPlot(aim_cell)
aim_cell = RunTSNE(aim_cell,dims = 1:15)
aim_cell = RunUMAP(aim_cell,dims = 1:15)
saveRDS(aim_cell,"7_scenic/aim_cell.RDS")
for (i in c$TF) {
  i="EGR1"
  x = grep(paste0(i),colnames(aim_cell@meta.data))
  y = colnames(aim_cell@meta.data)[x]
  for (k in y) {
    p = FeaturePlot(aim_cell, features=k, label=F, reduction = 'tsne')+
      theme(text = element_text(size=10),axis.text = element_text(size=10))+
      scale_color_distiller(palette = "Spectral")+coord_fixed();p
    ggsave2(paste0("tmp/Fig.4",k,".png"),p,height = 3.7,width = 4.7)
  }
}

p = DimPlot(aim_cell,reduction = "tsne",group.by = "Class")
ggsave2("tmp/supplementary Fig.4a1.pdf",p,height = 3.7,width = 4.7)

p = FeaturePlot(aim_cell, features="EGR1", label=F, reduction = 'tsne')+
  theme(text = element_text(size=10),axis.text = element_text(size=10))+
  scale_color_distiller(palette = "Spectral")+coord_fixed();p
ggsave2(paste0("../Fig/Fig.4f.pdf"),p,height = 3.7,width = 4.7)


p1 = FeatureScatter(aim_cell, feature1 = "SPIB_776g", 
               feature2 = "MYC_334g",pt.size = 0.1)+NoLegend()+
  theme(text = element_text(size=20),axis.text = element_text(size=16));p1

pheatmap(regulonAUC,annotation_col = data.frame(row.names = rownames(aim_cell@meta.data),
                                                Class = aim_cell$Class))


BINmatrix <- loadInt(scenicOptions, "aucell_binary_full")#readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- gsub(' \\(','_',colnames(BINmatrix))
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN

# scRNAbin <- AddMetaData(scRNA, BINmatrix)
# saveRDS(scRNAbin, 'scRNAbin.rds')

BINmatrix = BINmatrix[rownames(cellInfo)[which(cellInfo$celltypesub=="SC-C2")],]
BINmatrix = BINmatrix[,which(colSums(BINmatrix)>(254*0.25))]
Heatmap(t(BINmatrix),show_column_names = F,col = c("#0099CC","#FF0033"))

#############integrating WCN-TEC, TF regulatory network, and CCN#############
ESP = read.csv("0essential data/ESP_finnal.csv",header = T,row.names = 1,stringsAsFactors = F)[,1]
ESLs = read.csv("8lncRNA/ESLs.csv")[,2]
scRNA = readRDS("4single cell/scRNA_finnal.Rds")

#LR from cellphonedb  V3
cpdb3=read.csv("11cellphonedb/cellphonedb-data-3.0.0/data/gene_input.csv")

LR3=unique(c(cpdb3$gene_name,cpdb3$hgnc_symbol))

LR3[LR3 %in% c(ESP,ESLs)]
#"COPA"  "RPS19" "TFRC"  "RACK1"


cyt_single_T_weight <- readRDS("5network/cyt_single_T_weight.RDS")

cyt_single_T_weight= rbind(cyt_single_T_weight[,1:3],
                           data.frame(fromNode=cyt_single_T_weight$toNode,
                                      toNode=cyt_single_T_weight$fromNode,
                                      weight=cyt_single_T_weight$weight))
network_inuse = split(cyt_single_T_weight[,1:3],cyt_single_T_weight$fromNode)
rm(cyt_single_T_weight)

network_ESG = network_inuse[which(names(network_inuse) %in% c(ESLs,ESP))]
saveRDS(network_ESG,"11LR_ESG/network_ESG.RDS")
rm(network_inuse)

if(T){#TOM weight >0.001 or top100 genes of each ESPs or ESLs.
  network_ESG <- readRDS("11LR_ESG/network_ESG.RDS")
  for (i in names(network_ESG)) {
    network_ESG[[i]]=network_ESG[[i]][order(network_ESG[[i]]$weight,decreasing = T),]
    a = subset(network_ESG[[i]],network_ESG[[i]]$weight>0.001)
    b = network_ESG[[i]][1:100,]
    c = distinct(rbind(a,b))
    network_ESG[[i]]=c
  }
}

network_ESG_ligand_receptor = data.frame()
for (i in names(network_ESG)) {
  network_ESG_ligand_receptor = rbind(network_ESG_ligand_receptor,
                                network_ESG[[i]])
}


regulon = read.csv("7_scenic/TandN/regulonTargetsInfo_ESP_ESLs.csv",row.names = 1)

regulon = regulon[which(regulon$TF %in% ESP | regulon$gene %in% c(ESLs,ESP)),]
regulon = regulon[,c(1,2,10)]
colnames(regulon)=colnames(network_ESG_ligand_receptor)

network_ESG_ligand_receptor = rbind(network_ESG_ligand_receptor,regulon)
network_ESG_ligand_receptor=network_ESG_ligand_receptor[,-3]
network_ESG_ligand_receptor=distinct(network_ESG_ligand_receptor)
saveRDS(network_ESG_ligand_receptor,"11LR_ESG/network_ESG_ligand_receptor_finnal.RDS")

network_ESG_ligand_receptor <- readRDS("11LR_ESG/network_ESG_ligand_receptor_finnal.RDS")

####cellphonedb V3####
scRNA = readRDS("4single cell/scRNA_finnal.Rds")
scRNA = subset(scRNA,cells = which(scRNA$Class=="Tumor"))

write.table(as.matrix(scRNA@assays$SCT@data), '11cellphonedb/cellphonedb_count.txt', sep='\t', quote=F)

meta_data1 <- cbind(rownames(scRNA@meta.data), scRNA@meta.data['cell_type_my'])
colnames(meta_data1)=c("barcode_sample",'cell_type')
write.table(meta_data1, '11cellphonedb/cellphonedb_meta1.txt', sep='\t', quote=F, row.names = F)

meta_data2 <- cbind(rownames(scRNA@meta.data), scRNA@meta.data['Cell_subtype'])
colnames(meta_data2)=c("barcode_sample",'cell_type')
write.table(meta_data2, '11cellphonedb/cellphonedb_meta2.txt', sep='\t', quote=F, row.names = F)


# nohup cellphonedb method statistical_analysis cellphonedb_meta2.txt cellphonedb_count.txt --counts-data=gene_name --threads=8 &
# nohup cellphonedb method statistical_analysis 11cellphonedb_meta2.txt cellphonedb_count.txt --counts-data=gene_name --threads=8 &

a= fread("11cellphonedb/out_meta1/significant_means.txt",data.table = F,stringsAsFactors = F)

d = data.frame()
for (i in 13:ncol(a)) {
  d = rbind(d,
            data.frame(source_target=colnames(a)[i],
                       a[,c(1:12)],
                       significant_means=a[,i]))
}
d$source = strsplit2(d$source_target,"[|]")[,1]
d$target = strsplit2(d$source_target,"[|]")[,2]

cellphonedb = d[,c("source","target","significant_means",
                   "gene_a","receptor_a","partner_a",
                   "gene_b","receptor_b","partner_b",
                   "secreted","is_integrin","source_target","id_cp_interaction","interacting_pair",
                   "annotation_strategy","rank")]

##significant_means>0 and associated with epi
cellphonedb = subset(cellphonedb,cellphonedb$significant_means>0)
cellphonedb = subset(cellphonedb,target=="Epithelial cell" | source=="Epithelial cell")
write.csv(cellphonedb,"11cellphonedb/out_tumoreps/cellphone_filted.csv")


#281pair
write.csv(cellphonedb,"11cellphonedb/cellphone_filted_Epi.csv")

#complex
data = unique(c(cellphonedb$partner_a,cellphonedb$partner_b))
data[grep('complex',data)]

L_R_pair = unique(c(cellphonedb$gene_a,cellphonedb$gene_b,
                    'FLT1','ITGA2','ITGB1','ITGA6','ITGAE','ITGB7','ITGA10','ITGA1','ITGA5','ITGA11','ITGA4','ITGA6','ITGB4',
                    'IFNGR1','IFNGR2'))

network_ESG_ligand_receptor_sub=subset(network_ESG_ligand_receptor,
                                       network_ESG_ligand_receptor$toNode %in% L_R_pair)
write.csv(network_ESG_ligand_receptor_sub,"network_ESG_ligand_receptor_sub.csv")
network_ESG_ligand_receptor_sub=read.csv("network_ESG_ligand_receptor_sub.csv",row.names = 1)


unique(network_ESG_ligand_receptor_sub$fromNode[network_ESG_ligand_receptor_sub$fromNode %in% ESP])
unique(network_ESG_ligand_receptor_sub$fromNode[network_ESG_ligand_receptor_sub$fromNode %in% ESLs])
unique(network_ESG_ligand_receptor_sub$fromNode[network_ESG_ligand_receptor_sub$fromNode %in% regulon$TF])


e = as.data.frame.array(table(network_ESG_ligand_receptor_sub$toNode))
e$group=rownames(e)
colnames(e)[1]="count"
e = e[order(e$count,decreasing = T),]

interestingGenes=c("MIF","CD47","CXADR","GPI","CDH1","RPS19")

cellphone_5=subset(cellphonedb,cellphonedb$gene_a %in% interestingGenes | cellphonedb$gene_b %in% interestingGenes)

cellphone_5$source_target2 = factor(cellphone_5$source_target,
                                       levels=c("Epithelial cell|B cell","Epithelial cell|Endothelial cell",
                                                "Epithelial cell|Epithelial cell","Epithelial cell|Fibroblast cell",
                                                "Epithelial cell|Mast cell","Epithelial cell|Myeloid cell",
                                                "Epithelial cell|T cell",
                                                "B cell|Epithelial cell","Endothelial cell|Epithelial cell",
                                                "Fibroblast cell|Epithelial cell",
                                                "Mast cell|Epithelial cell","Myeloid cell|Epithelial cell",
                                                "T cell|Epithelial cell"))


#subpopulation cellphonedb
a= fread("11cellphonedb/out_meta2/significant_means.txt",data.table = F,stringsAsFactors = F)

d = data.frame()
for (i in 13:ncol(a)) {
  d = rbind(d,
            data.frame(source_target=colnames(a)[i],
                       a[,c(1:12)],
                       significant_means=a[,i]))
}
d$source = strsplit2(d$source_target,"[|]")[,1]
d$target = strsplit2(d$source_target,"[|]")[,2]

cellphonedb = d[,c("source","target","significant_means",
                   "gene_a","receptor_a","partner_a",
                   "gene_b","receptor_b","partner_b",
                   "secreted","is_integrin","source_target","id_cp_interaction","interacting_pair",
                   "annotation_strategy","rank")]

#ignificant_means>0
cellphonedb = subset(cellphonedb,cellphonedb$significant_means>0)
cellphonedb = subset(cellphonedb,target %in% c("CMS1",'CMS2','CMS3','CMS4','Goblet cells','Intermediate','Mature Enterocytes type 1','Mature Enterocytes type 2','Stem-like/TA') |
                       source %in% c("CMS1",'CMS2','CMS3','CMS4','Goblet cells','Intermediate','Mature Enterocytes type 1','Mature Enterocytes type 2','Stem-like/TA'))

data = unique(c(cellphonedb$partner_a,cellphonedb$partner_b))
data[grep('complex',data)]

L_R_pair = unique(c(cellphonedb$gene_a,cellphonedb$gene_b,
                    'ITGB4','ITGA6',
                    'ITGB2','ITGAM','ITGA8','ITGA6',
                    'ITGA5','ITGAV','ITGB5','ITGAX',
                    'ITGA2','ITGA3','ITGA11','ITGA9',
                    'ITGA1','ITGB1',
                    'ITGAL','ITGB2','ITGA10','ITGAE','ITGB7',
                    'ITGA4','ITGB1',
                    'TGFBR2','TGFBR1',
                    'CD94','NKG2C',
                    'IL1R1','IFNGR1','IFNGR2',
                    'FLT1','LIFR','BMPR1A','BMPR2',
                    'LRP6','TREM2',
                    'FZD6','LRP6','FZD4',
                    'IL10RB','IL20RA'))

network_ESG_ligand_receptor_sub=subset(network_ESG_ligand_receptor,
                                       network_ESG_ligand_receptor$toNode %in% L_R_pair)
write.csv(network_ESG_ligand_receptor_sub,"network_ESG_ligand_receptor_sub2.csv")
network_ESG_ligand_receptor_sub=read.csv("network_ESG_ligand_receptor_sub.csv",row.names = 1)


unique(network_ESG_ligand_receptor_sub$fromNode[network_ESG_ligand_receptor_sub$fromNode %in% ESP])
unique(network_ESG_ligand_receptor_sub$fromNode[network_ESG_ligand_receptor_sub$fromNode %in% ESLs])
unique(network_ESG_ligand_receptor_sub$fromNode[network_ESG_ligand_receptor_sub$fromNode %in% regulon$TF])


e = as.data.frame.array(table(network_ESG_ligand_receptor_sub$toNode))
e$group=rownames(e)
colnames(e)[1]="count"
e = e[order(e$count,decreasing = T),]

interestingGenes=c("MIF","CD47","CXADR","GPI","CDH1","RPS19")
interestingGenes=c("MIF","CD47")

cellphone_5=subset(cellphonedb,cellphonedb$gene_a %in% interestingGenes | cellphonedb$gene_b %in% interestingGenes)

p=ggplot(cellphone_5,aes(x=source_target,y = interacting_pair))+
  geom_point(aes(fill=significant_means,stroke=0.3),size=5,shape=21)+
  scale_fill_distiller(palette = "Spectral")+
  theme_classic()+FontSize(y.text = 8,x.text = 4)+
  rotate_x_text(angle = 30);p
ggsave("Fig/Fig.5c1.pdf",p,width = 20,height = 4)


#####iLINCS######
#CD47,MIF related ESP ESL and TF
iLINCS = network_ESG_ligand_receptor_sub[which(network_ESG_ligand_receptor_sub$toNode %in% c("CD47","MIF")),1]

gene = read.csv("iLINCS_gene.csv",row.names = 1)
info = read.csv("iLINCS_info.csv",row.names = 1,stringsAsFactors = T)
gene = gene[,rownames(info)]
identical(rownames(info),colnames(gene))

# colnames(gene)=paste0(rownames(info)," (",info$Drug,")")
# rownames(info)=colnames(gene)
bk <- c(seq(-1.5,-0.1,by=0.1),seq(0,1.5,by=0.1))
gene2 = data.frame(gene = rownames(gene),
                   max = apply(gene, 1, max),
                   min=apply(gene, 1, min))
gene2$label=""
gene2$label=ifelse(gene2$max>3,gene2$gene,ifelse(gene2$min<(-3),gene2$gene,""))

ann_colors = list(
  DatasetID=c(GSE116436="#9932CC"),
  Drug= c("cisplatin"="#808000","dasatinib"="#FF00FF","doxorubicin"="#FA8072",
          "erlotinib"="#7B68EE","geldanamycin"="#9400D3","lapatinib"="#800080",
          "paclitaxel"="#A0522D","sunitinib"="#D2B48C","topotecan"="#D2691E","vorinostat"="#87CEEB"),
  Concentration = c("10nM"="#DC143C","100nM"="#0000FF","200nM"="#20B2AA",
                    "1000nM"="#FFA500","2000nM"="#9370DB","3000nM"="#98FB98",
                    "5000nM"="#F08080","10000nM"="#1E90FF","15000nM"="#7CFC00"),
  CellLine = c("COLO-205"="#40E0D0","HCC-2998"="#5F9EA0","HCT-15"="#FF1493","KM12"="#0000CD","SW-620"="#008B8B"),
  Time=c("2h"="#32CD32","6h"="#F0E68C","24h"="#FFFFE0")
)

pheatmap(gene,cluster_rows = T,cluster_cols = T,show_rownames = T,annotation_col = info,
         color = c(colorRampPalette(colors = c("#18499E","#FFFFFF"))(length(bk)/2),
                   colorRampPalette(colors = c("#FFFFFF","#FF003D"))(length(bk)/2)),
         breaks = bk,filename = "Fig/Fig.5e.pdf",width = 16,height = 10,annotation_colors =ann_colors,
         border=FALSE,na_col = "#000000",angle_col =45,fontsize_row = 2,fontsize_col = 10)
#colorRampPalette(colors = brewer.pal(12, "Spectral"))(length(bk))
