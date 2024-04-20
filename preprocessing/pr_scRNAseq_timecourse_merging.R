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

data_dir <- "./pr-D0-Z1-L/pr-D0-Z1-L/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_D0_Z1_L = CreateSeuratObject(counts = data, project="pr_D0_Z1_L", min.cells=5, min.features=250)

data_dir <- "./pr-D0-Z2-L/pr-D0-Z2-L/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_D0_Z2_L = CreateSeuratObject(counts = data, project="pr_D0_Z2_L", min.cells=5, min.features=250)

data_dir <- "./pr-D0-Z3-L/pr-D0-Z3-L/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_D0_Z3_L = CreateSeuratObject(counts = data, project="pr_D0_Z3_L", min.cells=5, min.features=250)

data_dir <- "./pr-hd-D2-Z1-L/pr-hd-D2-Z1-L/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_hd_D2_Z1_L = CreateSeuratObject(counts = data, project="pr_hd_D2_Z1_L", min.cells=5, min.features=250)

data_dir <- "./pr-hd-D2-Z2-L/pr-hd-D2-Z2-L/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_hd_D2_Z2_L = CreateSeuratObject(counts = data, project="pr_hd_D2_Z2_L", min.cells=5, min.features=250)

data_dir <- "./pr-hd-D2-Z3-L/pr-hd-D2-Z3-L/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_hd_D2_Z3_L = CreateSeuratObject(counts = data, project="pr_hd_D2_Z3_L", min.cells=5, min.features=250)

data_dir <- "./pr-hd-D3-Z1-L/pr-hd-D3-Z1-L/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_hd_D3_Z1_L = CreateSeuratObject(counts = data, project="pr_hd_D3_Z1_L", min.cells=5, min.features=250)

data_dir <- "./pr-hd-D3-Z2-L/pr-hd-D3-Z2-L/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_hd_D3_Z2_L = CreateSeuratObject(counts = data, project="pr_hd_D3_Z2_L", min.cells=5, min.features=250)

data_dir <- "./pr-hd-D3-Z3-L/pr-hd-D3-Z3-L/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_hd_D3_Z3_L = CreateSeuratObject(counts = data, project="pr_hd_D3_Z3_L", min.cells=5, min.features=250)

data_dir <- "./pr-ld-D2-Z1-L/pr-ld-D2-Z1-L/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_ld_D2_Z1_L = CreateSeuratObject(counts = data, project="pr_ld_D2_Z1_L", min.cells=5, min.features=250)

data_dir <- "./pr-ld-D2-Z2-L/pr-ld-D2-Z2-L/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_ld_D2_Z2_L = CreateSeuratObject(counts = data, project="pr_ld_D2_Z2_L", min.cells=5, min.features=250)

data_dir <- "./pr-ld-D2-Z3-L/pr-ld-D2-Z3-L/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_ld_D2_Z3_L = CreateSeuratObject(counts = data, project="pr_ld_D2_Z3_L", min.cells=5, min.features=250)

data_dir <- "./pr-ld-D3-Z1-L/pr-ld-D3-Z1-L/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_ld_D3_Z1_L = CreateSeuratObject(counts = data, project="pr_ld_D3_Z1_L", min.cells=5, min.features=250)

data_dir <- "./pr-ld-D3-Z2-L/pr-ld-D3-Z2-L/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_ld_D3_Z2_L = CreateSeuratObject(counts = data, project="pr_ld_D3_Z2_L", min.cells=5, min.features=250)

data_dir <- "./pr-ld-D3-Z3-L/pr-ld-D3-Z3-L/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_ld_D3_Z3_L = CreateSeuratObject(counts = data, project="pr_ld_D3_Z3_L", min.cells=5, min.features=250)

#Merge objects in object list

obj_list <- list(pr_D0_Z1_L, pr_D0_Z2_L, pr_D0_Z3_L, pr_hd_D2_Z1_L, pr_hd_D2_Z2_L, pr_hd_D2_Z3_L, pr_hd_D3_Z1_L, pr_hd_D3_Z2_L, pr_hd_D3_Z3_L, pr_ld_D2_Z1_L, pr_ld_D2_Z2_L, pr_ld_D2_Z3_L, pr_ld_D3_Z1_L, pr_ld_D3_Z2_L, pr_ld_D3_Z3_L)

names(obj_list) <- c("pr_D0_Z1_L", "pr_D0_Z2_L", "pr_D0_Z3_L", "pr_hd_D2_Z1_L", "pr_hd_D2_Z2_L", "pr_hd_D2_Z3_L", "pr_hd_D3_Z1_L", "pr_hd_D3_Z2_L", "pr_hd_D3_Z3_L", "pr_ld_D2_Z1_L", "pr_ld_D2_Z2_L", "pr_ld_D2_Z3_L", "pr_ld_D3_Z1_L", "pr_ld_D3_Z2_L", "pr_ld_D3_Z3_L") 


#Iterate through the object lists and write DoubletFinder results to a text file  

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

write.table(DoubSing, "./DoubSing_pr_lung.txt", sep="\t", quote=FALSE)

#Merge into single Seurat object, and write it to disk for subsequent processing
hamster_all <- merge(pr_D0_Z1_L, y = c(pr_D0_Z2_L, pr_D0_Z3_L, pr_hd_D2_Z1_L, pr_hd_D2_Z2_L, pr_hd_D2_Z3_L, pr_hd_D3_Z1_L, pr_hd_D3_Z2_L, pr_hd_D3_Z3_L, pr_ld_D2_Z1_L, pr_ld_D2_Z2_L, pr_ld_D2_Z3_L, pr_ld_D3_Z1_L, pr_ld_D3_Z2_L, pr_ld_D3_Z3_L), add.cell.ids = c("pr_D0_Z1_L", "pr_D0_Z2_L", "pr_D0_Z3_L", "pr_hd_D2_Z1_L", "pr_hd_D2_Z2_L", "pr_hd_D2_Z3_L", "pr_hd_D3_Z1_L", "pr_hd_D3_Z2_L", "pr_hd_D3_Z3_L", "pr_ld_D2_Z1_L", "pr_ld_D2_Z2_L", "pr_ld_D2_Z3_L", "pr_ld_D3_Z1_L", "pr_ld_D3_Z2_L", "pr_ld_D3_Z3_L"), project = "ma_lung")

saveRDS(hamster_all, "./pr_seu_lung_new_combined_250.rds")

hamster_all@meta.data %>% group_by(orig.ident) %>% tally()

