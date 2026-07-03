# Examine ERV expression in T cells from patients with MDS
library(Seurat)
library(dplyr)
library(ggplot2)
library(gghalves)
library(ggpattern)
library(svglite)
library(tidyr)
library(readxl)

# Set the working directory
setwd('/Users/henryhampton/OneDrive - UNSW/2026_06_17_Code_for_upload/T_cell_figures/code')

# Read in the TE annotation file which Anne sent me
df <- read.csv('../data/TE_family.csv',  header = TRUE)

# Read in the Seurat object whcih Anne sent me
data<-readRDS("../data/combined_scTE_all.RDS")

# Add the patient_alias and remove the 'patient' column
alias = read_excel('../data/patient_alias_2025.xlsx')
alias_dict <- setNames(alias$alias, alias$PID)
data@meta.data['alias'] = alias_dict[data@meta.data$patient]
data@meta.data <- data@meta.data %>% dplyr::select(-patient)
rm(alias)

parts <- strsplit(data@meta.data$sample, "_")
data@meta.data$sample = sapply(parts, function(x) {
  paste(c(alias_dict[[x[1]]], x[-1]), collapse = "_")
})

# update the clustering

annotations = read_xlsx('../data/Copy of annotations_for_chan.xlsx' )

anno_dict = setNames(annotations$Chans_annotation, annotations$Annes_Annotation)

anno_dict
data@meta.data$cell_type1 = anno_dict[data@meta.data$celltype_0527]


new_meta <- data@meta.data %>% dplyr::select(-c("leiden_res0_5", "leiden_res0_6",          
                                                "leiden_res0_7", "leiden_res0_8", "leiden_res0_9",          
                                                "leiden_res0_10", "leiden_res0_12", 
                                                '_index', 'category', 'celltype', 'barcode', 
                                                'celltype_group', "celltype_0425", "celltype_merge","celltype_0527", 
                                                "BEST", 'celltype_0516', "total_counts_ribo", "pct_counts_ribo", "BEST",
                                                "specific_outcome_C12D29", "specific_outcome_C6D28" ))

print(dim(new_meta))            # confirm columns actually dropped here
print(colnames(new_meta)) 
data@meta.data <- new_meta

print(dim(data@meta.data))

                                            

saveRDS(data, "../data/combined_scTE_all_patient_updated.RDS")
getwd()

