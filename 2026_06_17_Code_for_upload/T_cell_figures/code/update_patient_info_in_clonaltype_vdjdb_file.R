library(Seurat)
library(SeuratDisk)
library(anndata)
library(reticulate)
library(zellkonverter)
library(hdf5r)
library(dplyr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(readxl)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in the single cell object which Anne sent me
setwd('/Users/henryhampton/OneDrive - UNSW/2026_06_17_Code_for_upload/T_cell_figures/code/')
data <- readRDS("../data/clonaltype_vdjdb_annotation.rds")


# Add the patient_alias and remove the 'patient' column
alias = readxl::read_excel('../data/patient_alias_2025.xlsx', sheet = 'Sheet1')
alias_dict <- setNames(alias$alias, alias$PID)
data@meta.data['alias'] = alias_dict[data@meta.data$patient]
data@meta.data <- data@meta.data %>% dplyr::select(-patient)

parts <- strsplit(data@meta.data$sample, "_")
data@meta.data$sample = sapply(parts, function(x) {
  paste(c(alias_dict[[x[1]]], x[-1]), collapse = "_")
})

parts <- strsplit(rownames(data@meta.data), "_")
rownames(data@meta.data) = sapply(parts, function(x) {
  paste(c(alias_dict[[x[1]]], x[-1]), collapse = "_")
})

parts <- strsplit((data@assays$RNA@counts@Dimnames[[2]]), "_")
data@assays$RNA@counts@Dimnames[[2]] = sapply(parts, function(x) {
  paste(c(alias_dict[[x[1]]], x[-1]), collapse = "_")
})

parts <- strsplit((data@assays$RNA@data@Dimnames[[2]]), "_")
data@assays$RNA@data@Dimnames[[2]] = sapply(parts, function(x) {
  paste(c(alias_dict[[x[1]]], x[-1]), collapse = "_")
})


parts <- strsplit(names(data@active.ident), "_")
names(data@active.ident) = sapply(parts, function(x) {
     paste(c(alias_dict[[x[1]]], x[-1]), collapse = "_")
})

parts <- strsplit(rownames(data@reductions$XPCA_@cell.embeddings), "_")
rownames(data@reductions$XPCA_@cell.embeddings) = sapply(parts, function(x) {
  paste(c(alias_dict[[x[1]]], x[-1]), collapse = "_")
})

parts <- strsplit(rownames(data@reductions$XPCAHARMONY_@cell.embeddings), "_")
rownames(data@reductions$XPCAHARMONY_@cell.embeddings) = sapply(parts, function(x) {
  paste(c(alias_dict[[x[1]]], x[-1]), collapse = "_")
})

parts <- strsplit(rownames(data@reductions$XUMAP_@cell.embeddings), "_")
rownames(data@reductions$XUMAP_@cell.embeddings) = sapply(parts, function(x) {
  paste(c(alias_dict[[x[1]]], x[-1]), collapse = "_")
})

annotations = read_xlsx('../data/Copy of annotations_for_chan.xlsx' )

anno_dict = setNames(annotations$Chans_annotation, annotations$Annes_Annotation)

anno_dict
data@meta.data$cell_type1 = anno_dict[data@meta.data$celltype_0527]


new_meta <- data@meta.data %>% dplyr::select(-c("leiden_res0_5", "leiden_res0_6",          
                                                "leiden_res0_7", "leiden_res0_8", "leiden_res0_9",          
                                                "leiden_res0_10", "leiden_res0_12", 
                                                'category', 'celltype', 'barcode', 
                                                'celltype_group', "celltype_0425", "celltype_merge","celltype_0527", 
                                                "BEST", 'celltype_0516', "total_counts_ribo", "pct_counts_ribo", "BEST",
                                                "specific_outcome_C12D29", "specific_outcome_C6D28" ))

print(dim(new_meta))            # confirm columns actually dropped here
print(colnames(new_meta)) 
data@meta.data <- new_meta


saveRDS(data, "../data/clonaltype_vdjdb_annotation_patient_updated.RDS")
