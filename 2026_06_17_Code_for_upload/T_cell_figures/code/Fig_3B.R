library(scRepertoire)
library(Seurat)
library(scCustomize)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(dplyr)
library(patchwork)
library(svglite)
library(SeuratDisk)
library(zellkonverter)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd('/Users/henryhampton/OneDrive - UNSW/2026_06_17_Code_for_upload/T_cell_figures/code')

# Read in data
sce <- readH5AD("../data/t_cell_obj_v2.h5ad")
seurat_obj <- as.Seurat(sce, counts = "counts", data = "logcounts")

# 
# # Read in Chan's updated annotations
annotations = read_excel('../data/Copy of annotations_for_chan.xlsx')
annotations$Chans_annotation <- gsub('Memory1', 'Memory_1', annotations$Chans_annotation, fixed = TRUE)
annotations$Chans_annotation <- gsub('Memory2', 'Memory_2', annotations$Chans_annotation, fixed = TRUE)
annotations_dict_rev = setNames(annotations$Chans_annotation, annotations$Annes_Annotation)


# Grab the recent color dictionary
color = read_excel('../data/t_cell_color_scheme.xlsx')  
colors <- setNames(color$color, annotations_dict_rev[color$cell_type])

# Extract UMAP coordinates and metadata
umap_data <- as.data.frame(Embeddings(seurat_obj, reduction = "XUMAP_"))
umap_data$celltype <- seurat_obj@meta.data$annotation #annotations_dict_rev[seurat_obj@meta.data$annotation]

# Create the plot WITHOUT legend
p <- ggplot(umap_data, aes(x = Xumap_1, y = Xumap_2, color = celltype)) +
  geom_point(size = 0.2, alpha = 1, stroke = 0) +
  scale_color_manual(values = colors) +
  theme_classic() +
  labs(title = "UMAP by Cell Type",
       x = "UMAP 1",
       y = "UMAP 2",
       color = "Cell Type") +
  theme(legend.position = "none",  # Remove legend
        aspect.ratio = 1)          # Make it square

# Save the plot
ggsave("../results/umap_plot_v2.pdf", plot = p,
       width = 5, height = 5, units = "cm", dpi = 600)

p
