library(ggplot2)
library(stringr)
library(Seurat)
library(HDF5Array)
library(DelayedArray)
library(DelayedMatrixStats)
require(gridExtra)
library(dplyr)
library(DoubletFinder)

data_dir <- "./pr-D0-Z1-B/pr-D0-Z1-B/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_D0_Z1_B = CreateSeuratObject(counts = data, project="pr_D0_Z1_B", min.cells=5, min.features=250)

data_dir <- "./pr-D0-Z2-B/pr-D0-Z2-B/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_D0_Z2_B = CreateSeuratObject(counts = data, project="pr_D0_Z2_B", min.cells=5, min.features=250)

data_dir <- "./pr-D0-Z3-B/pr-D0-Z3-B/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_D0_Z3_B = CreateSeuratObject(counts = data, project="pr_D0_Z3_B", min.cells=5, min.features=250)

data_dir <- "./pr-hd-D2-Z1-B/pr-hd-D2-Z1-B/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_hd_D2_Z1_B = CreateSeuratObject(counts = data, project="pr_hd_D2_Z1_B", min.cells=5, min.features=250)

data_dir <- "./pr-hd-D2-Z2-B/pr-hd-D2-Z2-B/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_hd_D2_Z2_B = CreateSeuratObject(counts = data, project="pr_hd_D2_Z2_B", min.cells=5, min.features=250)

data_dir <- "./pr-hd-D2-Z3-B/pr-hd-D2-Z3-B/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_hd_D2_Z3_B = CreateSeuratObject(counts = data, project="pr_hd_D2_Z3_B", min.cells=5, min.features=250)

data_dir <- "./pr-hd-D3-Z1-B/pr-hd-D3-Z1-B/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_hd_D3_Z1_B = CreateSeuratObject(counts = data, project="pr_hd_D3_Z1_B", min.cells=5, min.features=250)

data_dir <- "./pr-hd-D3-Z2-B/pr-hd-D3-Z2-B/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_hd_D3_Z2_B = CreateSeuratObject(counts = data, project="pr_hd_D3_Z2_B", min.cells=5, min.features=250)

data_dir <- "./pr-hd-D3-Z3-B/pr-hd-D3-Z3-B/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_hd_D3_Z3_B = CreateSeuratObject(counts = data, project="pr_hd_D3_Z3_B", min.cells=5, min.features=250)

data_dir <- "./pr-ld-D2-Z1-B/pr-ld-D2-Z1-B/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_ld_D2_Z1_B = CreateSeuratObject(counts = data, project="pr_ld_D2_Z1_B", min.cells=5, min.features=250)

data_dir <- "./pr-ld-D2-Z2-B/pr-ld-D2-Z2-B/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_ld_D2_Z2_B = CreateSeuratObject(counts = data, project="pr_ld_D2_Z2_B", min.cells=5, min.features=250)

data_dir <- "./pr-ld-D2-Z3-B/pr-ld-D2-Z3-B/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_ld_D2_Z3_B = CreateSeuratObject(counts = data, project="pr_ld_D2_Z3_B", min.cells=5, min.features=250)

data_dir <- "./pr-ld-D3-Z1-B/pr-ld-D3-Z1-B/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_ld_D3_Z1_B = CreateSeuratObject(counts = data, project="pr_ld_D3_Z1_B", min.cells=5, min.features=250)

data_dir <- "./pr-ld-D3-Z2-B/pr-ld-D3-Z2-B/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_ld_D3_Z2_B = CreateSeuratObject(counts = data, project="pr_ld_D3_Z2_B", min.cells=5, min.features=250)

data_dir <- "./pr-ld-D3-Z3-B/pr-ld-D3-Z3-B/outs/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir, gene.column = 2)
pr_ld_D3_Z3_B = CreateSeuratObject(counts = data, project="pr_ld_D3_Z3_B", min.cells=5, min.features=250)



hamster_all <- merge(pr_D0_Z1_B, y = c(pr_D0_Z2_B, pr_D0_Z3_B, pr_hd_D2_Z1_B, pr_hd_D2_Z2_B, pr_hd_D2_Z3_B, pr_hd_D3_Z1_B, pr_hd_D3_Z2_B, pr_hd_D3_Z3_B, pr_ld_D2_Z1_B, pr_ld_D2_Z2_B, pr_ld_D2_Z3_B, pr_ld_D3_Z1_B, pr_ld_D3_Z2_B, pr_ld_D3_Z3_B), add.cell.ids = c("pr_D0_Z1_B", "pr_D0_Z2_B", "pr_D0_Z3_B", "pr_hd_D2_Z1_B", "pr_hd_D2_Z2_B", "pr_hd_D2_Z3_B", "pr_hd_D3_Z1_B", "pr_hd_D3_Z2_B", "pr_hd_D3_Z3_B", "pr_ld_D2_Z1_B", "pr_ld_D2_Z2_B", "pr_ld_D2_Z3_B", "pr_ld_D3_Z1_B", "pr_ld_D3_Z2_B", "pr_ld_D3_Z3_B"), project = "ma_blood")


hamster_all@meta.data %>% group_by(orig.ident) %>% tally()

cells_to_keep <- read.table("./pr_blood_cells_to_keep.txt", sep="\t") 
hamster_all <- subset(hamster_all, cells=cells_to_keep$cell_id)

hamster_all@meta.data %>% group_by(orig.ident) %>% tally()


DefaultAssay(hamster_all) <- "RNA"


## integrate by hamster to remove batch effects

hamster.list <- SplitObject(hamster_all, split.by = "orig.ident")
for (i in 1:length(hamster.list)){
  hamster.list[[i]] <- SCTransform(hamster.list[[i]], verbose =T)
}

hamster.features <- SelectIntegrationFeatures(object.list = hamster.list, nfeatures = 3000)
hamster.list <- PrepSCTIntegration(object.list = hamster.list, anchor.features = hamster.features, 
                                   verbose = T)

hamster.anchors <- FindIntegrationAnchors(object.list = hamster.list, normalization.method = "SCT", 
                                          anchor.features = hamster.features, verbose = T)
hamster.integrated <- IntegrateData(anchorset = hamster.anchors, normalization.method = "SCT", 
                                    verbose = T)

## run dimensional reductions
#   PCA
hamster.integrated<- RunPCA(hamster.integrated, verbose = FALSE)
#   UMAP
hamster.integrated<- RunUMAP(hamster.integrated, dims = 1:30, verbose = FALSE)

hamster.integrated <- FindNeighbors(hamster.integrated, dims = 1:30)
hamster.integrated <- FindClusters(hamster.integrated, resolution = 0.5)


saveRDS(hamster.integrated, "./pr_seu_blood_new_combined_integrated.rds")


