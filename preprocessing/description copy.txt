Description Syrian/Golden hamster (Mesocricetus auratus) and Roborovski dwarf hamster (Phodopus roborovskii) lung scRNA-seq object processing

Step 1 – cellranger

For Syrian hamsters, cellranger (10x Genomics) was run using the raw sequencing data (fastq files) from Nouailles et al. 2021 (https://doi.org/10.1038/s41467-021-25030-7) and the MesAurJ_rev3_ext_S2 reference using standard settings for Syrian hamster data. The reference is described in Nouailles et al. 2023 (https://doi.org/10.1038/s41564-023-01352-8), NCBI GEO entry GSE200596. For dwarf hamsters, the fastq files are available through the NCBI GEO entry accompanying this paper (GSE241133). The reference phorob_curated_rev3_S2.gtf is describe in Andreotti et al. 2022 (https://doi.org/10.1093/gbe/evac100), figshare https://figshare.com/articles/dataset/Phodopus_roborovskii_assembly/16695457 .

Step 2 – merging of cellranger output and doublet detection

In script ma_scRNAseq_timecourse_merging.R/pr_scRNAseq_timecourse_merging.R, the cellranger output folders are converted to Seurat objects and merged. A very low depth per cell thresholds is uwed when reading in the data in order to preserve cells with low RNA content such as neutrophils. Doublets are identified using the DoubletFinder package, and written to a text file.

Step 3 - initial annotation of cell types and cell type-specific thresholdings

In script ma_scRNAseq_timecourse_part1.R/pr_scRNAseq_timecourse_part1.R, there is a cell type annotation of the cells with the low threshold. Subsequently, cell type-specific tresholds are calculated and written to a text file ma_lung_cells_to_keep.txt. This means that e.g. for neutrophils the threshold is much lower as for macrophages.

Step 4 – filter Seurat objects and integrate

In script ma_scRNAseq_timecourse_part2.R/pr_scRNAseq_timecourse_part2.R using the cells_to_keep text file from the previous step, cells below the threshold are filtered out (cellranger input is read in again). Subsequently, samples are integrated using the standard sctransform procedure. Note: this step needs a lot of memory.

Step 5 – annotate cell types on integrated object again

In script ma_scRNAseq_timecourse_part3.R/pr_scRNAseq_timecourse_part3.R, the integrated object is read in and annotated again, and the annotated object saved as seu_lung_new_combined_integrated_annotated.rds