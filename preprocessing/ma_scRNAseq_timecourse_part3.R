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


#read integrated object 
ma.int <- readRDS("/path/to/data/seu_lung_new_combined_integrated.rds")

#Add metadata columns
ma.int@meta.data = cbind(ma.int@meta.data, ma.int@reductions$umap@cell.embeddings)
ma.int@meta.data$timepoint <- gsub("ma_([de0-9]*)_lung_([0-9])","\\1",ma.int@meta.data$orig.ident)
ma.int@meta.data$hamster <- gsub("ma_([de0-9]*)_lung_([0-9])","Ha\\2",ma.int@meta.data$orig.ident)

DefaultAssay(ma.int) <- 'RNA'
SCoV2_rawcounts <- FetchData(ma.int, grep("SCoV2", ma.int@assays$RNA@counts@Dimnames[[1]], value="TRUE"), slot="counts")
ma.int@meta.data$SCoV2_sum <- rowSums(SCoV2_rawcounts)
SCoV2_rawcounts$SCoV2sum <- rowSums(SCoV2_rawcounts)
ma.int@meta.data$SCoV2_load <- SCoV2_rawcounts$SCoV2sum/ma.int@meta.data$nCount_RNA*100
DefaultAssay(ma.int) <- 'SCT'

ggplot()+
  geom_point(data=ma.int@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=timepoint), size=0.5, shape=16, alpha=0.5)+
  ggtitle("treatment")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("seu_lung_new_integrated_treatment.pdf", useDingbats=FALSE)
ggplot()+
  geom_point(data=ma.int@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=hamster), size=0.2, shape=16, alpha=0.5)+
  ggtitle("hamster")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("seu_lung_new_integrated_hamsters.pdf", useDingbats=FALSE)
UMAPPlot(ma.int, label=TRUE)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle("Clusters")
ggsave("seu_lung_new_integrated_clusters.pdf", useDingbats=FALSE)

ggplot()+geom_point(data=subset(ma.int@meta.data, SCoV2_load<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", size=0.5, shape = 16)+
  geom_point(data=subset(ma.int@meta.data, SCoV2_load>0), aes(x=UMAP_1, y=UMAP_2, colour=log10(SCoV2_load)), size=0.5, shape = 16)+
  scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle("log10 SCoV2 load")
ggsave(paste("seu_lung_new_integrated_log10SCoV2_load", "pdf", sep="."), useDingbats=FALSE)

unfilt_metadata <- read.table("/path/to/data/ma_new_lung/seu_lung_new_combined_250_metadata.txt")
ma.int@meta.data$cell_id <- rownames(ma.int@meta.data)
df <- dplyr::left_join(ma.int@meta.data, unfilt_metadata, by="cell_id")
ma.int@meta.data$DoubSing <- df$DoubSing
ma.int@meta.data$endocytotox <- df$endocytotox

#cell type annotation
#LOC101833790 = Ccr5

markers <- c("Ace2", "Ackr1", "Acta2", "Adgre1", "Arg1", "Aspn", "C1qb", "Camp", "Ccr2", "LOC101833790", "Cd3e", "Cd4", "Cd79b", "Cd8a", "Cldn5", "Col1a2", "COX1", "Cx3cr1", "Cxcr2", "Dcn", "Ednrb", "Fcgr4", "Flt3", "Foxj1", "Gzma", "Il7r", "Irf8", "Lamp3", "Marco", "Mrc1", "Ms4a1", "ND4", "Nkg7", "Pdpn", "Ppbp", "Plvap", "Retn", "Rtkn2", "S100a8", "Siglecf", "Tagln", "Tcf4", "Treml4")
for (gene in markers) {
  df <- FetchData(ma.int, gene)
  gene <- str_replace_all(gene, "-", "_")
  colnames(df) <- gene
  df = cbind(df, ma.int@reductions$umap@cell.embeddings)
  ggplot()+geom_point(data=subset(df, eval(parse(text = gene))<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", shape=16, size=0.5)+
    geom_point(data=subset(df, eval(parse(text = gene)) > 0), aes(x=UMAP_1, y=UMAP_2, colour=eval(parse(text = gene))), shape=16, size=0.5)+
    scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    coord_fixed(ratio=1)+
    ggtitle(gene)
  ggsave(paste("seu_lung_new_integrated", gene, "pdf", sep="."), useDingbats=FALSE)
}

avg <- AverageExpression(ma.int, assays=c("SCT"), features = markers, return.seurat = T, slot="data") 

DoHeatmap(avg, size=5, features=markers)+scale_fill_gradient2(low="midnightblue", mid= "white", high="red", na.value = "white")+NoLegend()
ggsave("seu_lung_new_integrated_cluster_heatmap.pdf", width=12, height=7)

#Myofibroblasts, fibroblasts, smooth muscle
#LOC101844074 is Cox4i2 (smooth muscle marker)
musclemarkers <- c("Aspn", "Dcn", "Acta2", "Tagln", "Wif1", "Fgf18", "Col1a2", "Bsg", "LOC101844074", "Cnn1", "Myh11", "Actg2", "Gpc3", "Apoe", "Serpinf1", "Gpc3")
avg <- AverageExpression(ma.int, assays=c("SCT"), features = musclemarkers, return.seurat = T, slot="data") 
DoHeatmap(avg, size=5, features=musclemarkers)+scale_fill_gradient2(low="midnightblue", mid= "white", high="red", na.value = "white")+NoLegend()
ggsave("seu_lung_new_integrated_musclemarkers_heatmap.pdf")

Idents(ma.int) <- ma.int@meta.data$seurat_clusters

#Annotate cells
ma.int <- RenameIdents(ma.int, 
                       '0'='Bcells',
                       '1'='Neutrophils',
                       '2'='MonocyticMacrophages',
                       '3'='AlveolarMacrophages',
                       '4'='TNKcells',
                       '5'='Endothelial',
                       '6'='AT2',
                       '7'='MonocyticMacrophages',
                       '8'='InterstitialMacrophages',
                       '9'='Endothelial',
                       '10'='TNKcells',
                       '11'='AT1',
                       '12'='TNKcells',
                       '13'='Treml4+Macrophages',
                       '14'='Endothelial',
                       '15'='TNKcells',
                       '16'='InterstitialMacrophages',
                       '17'='SmoothMuscle',
                       '18'='Myofibroblast',
                       '19'='mixed1',
                       '20'='Platelets',
                       '21'='mixed2',
                       '22'='MyeloidDendritic',
                       '23'='Endothelial',
                       '24'='Myofibroblast',
                       '25'='pDC',
                       '26'='mixed3',
                       '27'='AlveolarMacrophages',
                       '28'='Ciliated',
                       '29'='mixed4',
                       '30'='Fibroblasts',
                       '31'='TNKcells',
                       '32'='mixed5',
                       '33'='Bcells',
                       '34'='mixed6',
                       '35'='Bcells',
                       '36'='mixed7',
                       '37'='Bcells',
                       '38'='Bcells')
ma.int@meta.data$celltype <- Idents(ma.int)

ma.int[["percent.mt"]] <- PercentageFeatureSet(object = ma.int, features=c("ATP6", "COX1", "COX2", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"))
VlnPlot(ma.int, features = "percent.mt", pt.size = 0, y.max = 20)

saveRDS(ma.int, "/path/to/data/seu_lung_new_combined_integrated_annotated.rds")
