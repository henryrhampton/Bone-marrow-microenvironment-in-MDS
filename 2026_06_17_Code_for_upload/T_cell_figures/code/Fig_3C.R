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
data <- as.Seurat(sce, counts = "counts", data = "logcounts")
rm(sce)


# Read in Chan's updated annotations
annotations = read_excel('../data/Copy of annotations_for_chan.xlsx')
annotations_dict_rev = setNames(annotations$Chans_annotation, annotations$Annes_Annotation)


# Grab the recent color dictionary
color = read_excel('../data/t_cell_color_scheme.xlsx')  
colors <- setNames(color$color, annotations_dict_rev[color$cell_type])


marker_genes = c('CD4', 
                 'CCR7', 'LEF1',
                 'SELL',
                 #'PRKCA',  #'ARID1B', 
                 #'GATA3', #'ANXA1',
                 'FOXP3', 
                 'CD8A',  "GZMK", 'KLRB1', 'TRAV1-2',
                 #'LTB',  
                 'MTRNR2L8', 'IL7R', 
                 'GZMH', 'ITGB1', 'CMC1', 
                 
                 #'S1PR1',
                 'KLF2', 'CD69',
                 'GZMB', 'ZNF683', 'THEMIS',
                 'FCGR3A',  'KLRF1',   
                 'TRDV2')


celltype_order <- c("CD4_Naïve", "CD4 TCM", "CD4 TEM", "Treg", "CD4_Cytotoxic",
                    "MAIT", "GZMK+ Effector_Memory", "GZMH+ Effector_Memory_1",
                    "GZMH+ Effector_Memory_2", "Trm", "TEMRA", "TEMRA_NK-like", "gdT")

p <- p + scale_y_discrete(limits = rev(celltype_order))

p <- DotPlot(
  data,
  features = marker_genes,
  group.by = "annotation",
  cols = c("white", "black"),
  scale = TRUE,
  dot.scale = 6
) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7)) +
  labs(x = NULL, y = NULL)

p$layers[[1]]$aes_params$shape <- 21
p$layers[[1]]$aes_params$stroke <- 0.25
p$layers[[1]]$mapping$fill <- p$layers[[1]]$mapping$colour
p$layers[[1]]$mapping$colour <- NULL

p <- p + scale_fill_gradient(low = "white", high = "black")
p <- p + scale_y_discrete(limits = rev(celltype_order), labels = annotations_dict_rev)
p

# Specify the celltypes to color the plot by
celltypes <- levels(factor(annotations_dict_rev))

colorbar_df <- data.frame(
  celltype = factor(celltype_order, levels = rev(celltype_order)),
  x = 1
)


colorbar <- ggplot(colorbar_df, aes(x = x, y = celltype, fill = celltype)) +
  geom_tile() +
  scale_fill_manual(values = colors) +
  theme_void() +
  theme(legend.position = "none")

colorbar + p + plot_layout(widths = c(0.05, 1))
combined <- colorbar + p + plot_layout(widths = c(0.05, 1))

ggsave(filename = '../results/dotplot_v3.pdf', plot = p, dpi = 1200, width = 7.25, height = 3.5)
ggsave("../results/dotplot_v3_combined.pdf", plot = combined, width = 7.25, height = 3.5, dpi = 1200)
#ggsave(filename = 'results/dotplot_v2.pdf', plot = p, dpi = 1200, width = 9.25, height = 4.5)

