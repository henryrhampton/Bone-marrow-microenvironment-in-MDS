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

# Grab the UMAP data
plot_data <- as.data.frame(Embeddings(seurat_obj, reduction = "XUMAP_"))
plot_data$cloneType <- seurat_obj@meta.data$cloneType

color_names = c('Hyperexpanded (0.1 < X <= 1)' = '#ebe82e',
                'Large (0.01 < X <= 0.1)' = '#f69542', 
                'Medium (0.001 < X <= 0.01)'  =  '#ca4879', 
                'Small (1e-04 < X <= 0.001)' ='#793092')



p <- ggplot(plot_data, aes(x = Xumap_1, y = Xumap_2, color = cloneType)) +
  geom_point(size = 0.2, alpha = 1, stroke = 0) +
  #scale_color_viridis_c(option = "magma", trans = "log10") +  # log scale often helps since clone sizes are usually skewed
  scale_color_manual(values = color_names, na.value = "grey50") +
  theme_classic(base_size = 7) +
  labs(title = "UMAP by Clonotype Frequency",
       x = "UMAP 1",
       y = "UMAP 2",
       color = "Clonotype\nFrequency") +
  theme(aspect.ratio = 1, legend.position = "none")
p



# Save the plot
ggsave("../results/umap_plot_clone_size.pdf", plot = p,
       width = 5, height = 5, units = "cm", dpi = 600)

