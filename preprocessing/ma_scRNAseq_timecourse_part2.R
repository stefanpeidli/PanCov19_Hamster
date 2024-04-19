library(ggplot2)
library(stringr)
library(Seurat)
library(HDF5Array)
library(DelayedArray)
library(DelayedMatrixStats)
require(gridExtra)
library(dplyr)


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


hamster_all <- merge(ma_d0_lung_1, y = c(ma_d0_lung_2, ma_d0_lung_3, ma_d2_lung_1, ma_d2_lung_2, ma_d2_lung_3, ma_d3_lung_1, ma_d3_lung_2, ma_d3_lung_3, ma_d5_lung_1, ma_d5_lung_2, ma_d5_lung_3, ma_e14_lung_1, ma_e14_lung_2, ma_e14_lung_3), add.cell.ids = c("ma_d0_lung_1", "ma_d0_lung_2", "ma_d0_lung_3", "ma_d2_lung_1", "ma_d2_lung_2", "ma_d2_lung_3", "ma_d3_lung_1", "ma_d3_lung_2", "ma_d3_lung_3", "ma_d5_lung_1", "ma_d5_lung_2", "ma_d5_lung_3", "ma_e14_lung_1", "ma_e14_lung_2", "ma_e14_lung_3"), project = "ma_lung")

#Print number of cells per sample before filtering
hamster_all@meta.data %>% group_by(orig.ident) %>% tally()

cells_to_keep <- read.table("./ma_lung_cells_to_keep.txt") 
hamster_all <- subset(hamster_all, cells=cells_to_keep$cell_id)

#Print number of cells per sample before filtering
hamster_all@meta.data %>% group_by(orig.ident) %>% tally()


DefaultAssay(hamster_all) <- "RNA"


#Integrate by hamster to remove batch effects

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

hamster.integrated<- RunPCA(hamster.integrated, verbose = FALSE)
hamster.integrated<- RunUMAP(hamster.integrated, dims = 1:30, verbose = FALSE)

hamster.integrated <- FindNeighbors(hamster.integrated, dims = 1:30)
hamster.integrated <- FindClusters(hamster.integrated, resolution = 0.5)


saveRDS(hamster.integrated, "./seu_lung_new_combined_integrated.rds")

