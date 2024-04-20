library(plyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(stringr)
library(Seurat)
library(HDF5Array)
library(rhdf5)
library(hdf5r)
library(DelayedArray)
library(DelayedMatrixStats)
require(gridExtra)
library(tidyr)
library("dendextend")
library(DESeq2)
source("smooth_DimPlot.R")
library(akima)
library(pheatmap)
#see http://www.cookbook-r.com/Graphs/ Plotting_means_and_error_bars_(ggplot2)
source("summarySE.R")
library(cowplot)
library(SeuratObject)
library(lme4)
library(ComplexHeatmap)
library(patchwork)
library(forcats)
library(DoubletFinder)


#read integrated object from cluster
pr.int <- readRDS("path/to/data/pr_seu_lung_new_combined_integrated.rds")
pr.int@meta.data = cbind(pr.int@meta.data, pr.int@reductions$umap@cell.embeddings)
pr.int@meta.data$timepoint <- gsub("pr_(.*D[0-9])_Z([0-9])_L","\\1",pr.int@meta.data$orig.ident)
pr.int@meta.data$hamster <- gsub("pr_(.*D[0-9])_Z([0-9])_L","Ha\\2",pr.int@meta.data$orig.ident)

DefaultAssay(pr.int) <- 'RNA'
SCoV2_rawcounts <- FetchData(pr.int, grep("SCoV2", pr.int@assays$RNA@counts@Dimnames[[1]], value="TRUE"), slot="counts")
pr.int@meta.data$SCoV2_sum <- rowSums(SCoV2_rawcounts)
SCoV2_rawcounts$SCoV2sum <- rowSums(SCoV2_rawcounts)
pr.int@meta.data$SCoV2_load <- SCoV2_rawcounts$SCoV2sum/pr.int@meta.data$nCount_RNA*100
DefaultAssay(pr.int) <- 'SCT'

ggplot()+
  geom_point(data=pr.int@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=timepoint), size=0.5, shape=16, alpha=0.5)+
  ggtitle("timepoint")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("pr_lung_new_integrated_timepoint.pdf", useDingbats=FALSE)
ggplot()+
  geom_point(data=pr.int@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=hamster), size=0.2, shape=16, alpha=0.5)+
  ggtitle("hamster")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("pr_lung_new_integrated_hamsters.pdf", useDingbats=FALSE)
UMAPPlot(pr.int, label=TRUE)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle("Clusters")
ggsave("pr_lung_new_integrated_clusters.pdf", useDingbats=FALSE)

ggplot()+geom_point(data=subset(pr.int@meta.data, SCoV2_load<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", size=0.5, shape = 16)+
  geom_point(data=subset(pr.int@meta.data, SCoV2_load>0), aes(x=UMAP_1, y=UMAP_2, colour=log10(SCoV2_load)), size=0.5, shape = 16)+
  scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle("log10 SCoV2 load")
ggsave(paste("pr_lung_new_integrated_log10SCoV2_load", "pdf", sep="."), useDingbats=FALSE)

unfilt_metadata <- read.table("path/to/data/pr_lung_new/pr_seu_lung_new_combined_250_metadata.txt")
pr.int@meta.data$cell_id <- rownames(pr.int@meta.data)
df <- dplyr::left_join(pr.int@meta.data, unfilt_metadata, by="cell_id")
pr.int@meta.data$DoubSing <- df$DoubSing
pr.int@meta.data$endocytotox <- df$endocytotox


#cell type annotation
#LOC101833790 = Ccr5

markers <- c("Ace2", "Ackr1", "Acta2", "Adgre1", "Arg1", "Aspn", "C1qb", "Camp", "Ccl21", "Ccr2", "Ccr5", "Cd3e", "Cd4", "Cd79b", "Cd8a", "Cldn5", "Col1a2", "Cx3cr1", "Cxcr2", "Dcn", "Ednrb", "Fcgr4", "Flt3", "Foxj1", "Gzma", "Il7r", "Irf8", "Lamp3", "Marco", "Mrc1", "Ms4a1", "Nkg7", "Pdpn", "Ppbp", "Plvap", "Retn", "Rtkn2", "S100a8", "Siglecf", "Tagln", "Tcf4", "Treml4")
markers <- c("S100a9", "Ceacam1", "Cd14")
for (gene in markers) {
  df <- FetchData(pr.int, gene)
  gene <- str_replace_all(gene, "-", "_")
  colnames(df) <- gene
  df = cbind(df, pr.int@reductions$umap@cell.embeddings)
  ggplot()+geom_point(data=subset(df, eval(parse(text = gene))<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", shape=16, size=0.5)+
    geom_point(data=subset(df, eval(parse(text = gene)) > 0), aes(x=UMAP_1, y=UMAP_2, colour=eval(parse(text = gene))), shape=16, size=0.5)+
    scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    coord_fixed(ratio=1)+
    ggtitle(gene)
  ggsave(paste("pr_lung_new_integrated", gene, "pdf", sep="."), useDingbats=FALSE)
}

avg <- AverageExpression(pr.int, assays=c("SCT"), features = markers, return.seurat = T, slot="data") 

DoHeatmap(avg, size=5, features=markers)+scale_fill_gradient2(low="midnightblue", mid= "white", high="red", na.value = "white")+NoLegend()
ggsave("pr_seu_lung_new_integrated_cluster_heatmap.pdf", width=12, height=7)

#Myofibroblasts, fibroblasts, smooth muscle
#LOC101844074 is Cox4i2 (smooth muscle marker)
musclemarkers <- c("Aspn", "Dcn", "Acta2", "Tagln", "Wif1", "Fgf18", "Col1a2", "Bsg", "Cox4i2", "Cnn1", "Myh11", "Actg2", "Gpc3", "Apoe", "Serpinf1", "Gpc3")
avg <- AverageExpression(pr.int, assays=c("SCT"), features = musclemarkers, return.seurat = T, slot="data") 
DoHeatmap(avg, size=5, features=musclemarkers)+scale_fill_gradient2(low="midnightblue", mid= "white", high="red", na.value = "white")+NoLegend()
ggsave("pr_lung_new_integrated_musclemarkers_heatmap.pdf")

Idents(pr.int) <- pr.int@meta.data$seurat_clusters

pr.int <- RenameIdents(pr.int, 
                       '0'='Endothelial',
                       '1'='Bcells',
                       '2'='Neutrophils',
                       '3'='Bcells',
                       '4'='Treml4+Macrophages',
                       '5'='MonocyticMacrophages',
                       '6'='Neutrophils',
                       '7'='Endothelial',
                       '8'='AT2',
                       '9'='TNKcells',
                       '10'='TNKcells',
                       '11'='AlveolarMacrophages',
                       '12'='mixed1',
                       '13'='Neutrophils',
                       '14'='InterstitialMacrophages',
                       '15'='TNKcells',
                       '16'='Myofibroblast',
                       '17'='SmoothMuscle',
                       '18'='TNKcells',
                       '19'='Endothelial',
                       '20'='TNKcells',
                       '21'='Endothelial',
                       '22'='Bcells',
                       '23'='Endothelial',
                       '24'='pDC',
                       '25'='MyeloidDendritic',
                       '26'='AT1',
                       '27'='MonocyticMacrophages',
                       '28'='InterstitialMacrophages',
                       '29'='TNKcells',
                       '30'='Ciliated',
                       '31'='Fibroblasts',
                       '32'='mixed2',
                       '33'='mixed3')
pr.int@meta.data$celltype <- Idents(pr.int)

saveRDS(pr.int, "path/to/data/pr_lung_new_combined_integrated_annotated.rds")
