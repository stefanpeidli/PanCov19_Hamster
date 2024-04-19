library(ggplot2)
library(stringr)
library(Seurat)
library(HDF5Array)
library(DelayedArray)
library(DelayedMatrixStats)
require(gridExtra)
library(dplyr)
library(DoubletFinder)

#Create Seurat objects from individual cellranger output folders
#Use very low depth per cell thresholds when reading in data to preserve cells with low RNA content such as neutrophils

data_dir <- '/path/to/data/scRNAseq_Uninfectedanimal_1_lung/scRNAseq_Uninfectedanimal_1_lung/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d0_lung_1 = CreateSeuratObject(counts = data, project="ma_d0_lung_1", min.cells=5, min.features=250)

data_dir <- '/path/to/data/scRNAseq_Uninfectedanimal_2_lung/scRNAseq_Uninfectedanimal_2_lung/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d0_lung_2 = CreateSeuratObject(counts = data, project="ma_d0_lung_2", min.cells=5, min.features=250)

data_dir <- '/path/to/data/scRNAseq_Uninfectedanimal_3_lung/scRNAseq_Uninfectedanimal_3_lung/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d0_lung_3 = CreateSeuratObject(counts = data, project="ma_d0_lung_3", min.cells=5, min.features=250)

data_dir <- '/path/to/data/scRNAseq_2dpianimal_1_lung/scRNAseq_2dpianimal_1_lung/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d2_lung_1 = CreateSeuratObject(counts = data, project="ma_d2_lung_1", min.cells=5, min.features=250)

data_dir <- '/path/to/data/scRNAseq_2dpianimal_2_lung/scRNAseq_2dpianimal_2_lung/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d2_lung_2 = CreateSeuratObject(counts = data, project="ma_d2_lung_2", min.cells=5, min.features=250)

data_dir <- '/path/to/data/scRNAseq_2dpianimal_3_lung/scRNAseq_2dpianimal_3_lung/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d2_lung_3 = CreateSeuratObject(counts = data, project="ma_d2_lung_3", min.cells=5, min.features=250)

data_dir <- '/path/to/data/scRNAseq_3dpianimal_1_lung/scRNAseq_3dpianimal_1_lung/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d3_lung_1 = CreateSeuratObject(counts = data, project="ma_d3_lung_1", min.cells=5, min.features=250)

data_dir <- '/path/to/data//scRNAseq_3dpianimal_2_lung/scRNAseq_3dpianimal_2_lung/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d3_lung_2 = CreateSeuratObject(counts = data, project="ma_d3_lung_2", min.cells=5, min.features=250)

data_dir <- '/path/to/data/scRNAseq_3dpianimal_3_lung/scRNAseq_3dpianimal_3_lung/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d3_lung_3 = CreateSeuratObject(counts = data, project="ma_d3_lung_3", min.cells=5, min.features=250)

data_dir <- '/path/to/data/scRNAseq_5dpianimal_1_lung/scRNAseq_5dpianimal_1_lung/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d5_lung_1 = CreateSeuratObject(counts = data, project="ma_d5_lung_1", min.cells=5, min.features=250)

data_dir <- '/path/to/data/scRNAseq_5dpianimal_2_lung/scRNAseq_5dpianimal_2_lung/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d5_lung_2 = CreateSeuratObject(counts = data, project="ma_d5_lung_2", min.cells=5, min.features=250)

data_dir <- '/path/to/data/scRNAseq_5dpianimal_3_lung/scRNAseq_5dpianimal_3_lung/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d5_lung_3 = CreateSeuratObject(counts = data, project="ma_d5_lung_3", min.cells=5, min.features=250)

data_dir <- '/path/to/data/scRNAseq_14dpianimal_1_lung/scRNAseq_14dpianimal_1_lung/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_e14_lung_1 = CreateSeuratObject(counts = data, project="ma_e14_lung_1", min.cells=5, min.features=250)

data_dir <- '/path/to/data/scRNAseq_14dpianimal_2_lung/scRNAseq_14dpianimal_2_lung/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_e14_lung_2 = CreateSeuratObject(counts = data, project="ma_e14_lung_2", min.cells=5, min.features=250)

data_dir <- '/path/to/data/scRNAseq_14dpianimal_3_lung/scRNAseq_14dpianimal_3_lung/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_e14_lung_3 = CreateSeuratObject(counts = data, project="ma_e14_lung_3", min.cells=5, min.features=250)


#Merge objects in object list

obj_list <- list(ma_d0_lung_1, ma_d0_lung_2, ma_d0_lung_3, ma_d2_lung_1, ma_d2_lung_2, ma_d2_lung_3, ma_d3_lung_1, ma_d3_lung_2, ma_d3_lung_3, ma_d5_lung_1, ma_d5_lung_2, ma_d5_lung_3, ma_e14_lung_1, ma_e14_lung_2, ma_e14_lung_3)
names(obj_list) <- c("ma_d0_lung_1", "ma_d0_lung_2", "ma_d0_lung_3", "ma_d2_lung_1", "ma_d2_lung_2", "ma_d2_lung_3", "ma_d3_lung_1", "ma_d3_lung_2", "ma_d3_lung_3", "ma_d5_lung_1", "ma_d5_lung_2", "ma_d5_lung_3", "ma_e14_lung_1", "ma_e14_lung_2", "ma_e14_lung_3") 

#Iterate through the object lists and write results to a text file  

expr <- list()
for (the_file in names(obj_list)) {
  seu <- obj_list[[the_file]]
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu)
  seu <- RunUMAP(seu, dims = 1:18)
  seu <- FindNeighbors(seu, dims = 1:18)
  seu <- FindClusters(seu, resolution = 0.9)
  pdf(paste0(the_file, "_clusters.pdf"))
  DimPlot(seu, reduction = "umap", label=TRUE)
  dev.off()
  
  seu@meta.data = cbind(seu@meta.data, seu@reductions$umap@cell.embeddings)
  sweep.res.list_lung <- paramSweep_v3(seu, PCs = 1:18, sct = FALSE)
  sweep.stats_lung <- summarizeSweep(sweep.res.list_lung, GT = FALSE)
  bcmvn_lung <- find.pK(sweep.stats_lung)
  homotypic.prop <- modelHomotypic(seu@meta.data$seurat_clusters)

  # Assuming 5% doublet formation rate
  nExp_poi <- round(0.05*nrow(seu@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  seu <- doubletFinder_v3(seu, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  colnames(seu@meta.data) <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "RNA_snn_res.0.9", "seurat_clusters", "UMAP_1", "UMAP_2", "pAnn", "DoubSing")
  ggplot()+geom_point(data=seu@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=DoubSing), size=0.25)+theme_bw()
  ggsave(paste0(the_file, "_DoubSing.pdf"))
  expr[[the_file]] <- seu@meta.data
}
DoubSing <- do.call(rbind,expr)
DoubSing$cell_id <- rownames(DoubSing)
DoubSing <- DoubSing %>% select(c("cell_id", "DoubSing")) %>% mutate (cell_id = gsub("\\.", "\\_", cell_id))

write.table(DoubSing, "./DoubSing.txt", sep="\t", quote=FALSE)


#Merge into single Seurat object, and write it to disk for subsequent processing
hamster_all <- merge(ma_d0_lung_1, y = c(ma_d0_lung_2, ma_d0_lung_3, ma_d2_lung_1, ma_d2_lung_2, ma_d2_lung_3, ma_d3_lung_1, ma_d3_lung_2, ma_d3_lung_3, ma_d5_lung_1, ma_d5_lung_2, ma_d5_lung_3, ma_e14_lung_1, ma_e14_lung_2, ma_e14_lung_3), add.cell.ids = c("ma_d0_lung_1", "ma_d0_lung_2", "ma_d0_lung_3", "ma_d2_lung_1", "ma_d2_lung_2", "ma_d2_lung_3", "ma_d3_lung_1", "ma_d3_lung_2", "ma_d3_lung_3", "ma_d5_lung_1", "ma_d5_lung_2", "ma_d5_lung_3", "ma_e14_lung_1", "ma_e14_lung_2", "ma_e14_lung_3"), project = "ma_lung")

saveRDS(hamster_all, "./seu_lung_new_combined_250.rds")


#Print number of cells per sample
hamster_all@meta.data %>% group_by(orig.ident) %>% tally()


