library("scRepertoire")
library("Seurat")
library("scCustomize")
library("ggplot2")
library("readxl")

#
results = read.csv('../data/percent_non_unique.csv')


# Grab the recent color dictionary
color = readxl::read_excel('../data/t_cell_color_scheme.xlsx')  
colors <- setNames(color$color, annotations_dict_rev[color$cell_type]) # match colors to Chan's annotations

# Create a vector to specify the plotting order
celltype_order = c("CD4_Naïve", "CD4 TCM", "CD4 TEM", "Treg", "CD4_Cytotoxic", "MAIT",
                   "GZMK+ Effector_Memory","GZMH+ Effector_Memory", "Activated/exhausted Memory", "Trm",
                   "TEMRA",  "TEMRA_NK-like")

# Add a number to each name in celltype_order
numbered_labels <- setNames(
  paste0(seq_along(celltype_order), ".", celltype_order),
  celltype_order
)

# Make the graph
ggplot() +
  geom_col(data = filename_means, 
           aes(x = filename, y = percent_non_unique, fill = filename), 
           alpha = 1, color = "black",
           linewidth = 0.5, width = 0.8) +
  geom_point(data = results, 
             aes(x = filename, y = percent_non_unique), 
             color = "black",
             fill = "white", size = 1.0, shape = 21, stroke = 0.5, 
             position = position_jitter(width = 0.2, height = 0)) +
  scale_fill_manual(values = colors) +
  scale_x_discrete(limits = celltype_order, labels = numbered_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  theme_minimal() +
  theme(
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    text = element_text(size = 7),
    plot.title = element_text(size = 7),
    axis.title = element_text(size = 7),
    axis.title.y = element_text(margin = margin(r = 0)),
    axis.text.y = element_text(margin = margin(r = 0)),
    legend.text = element_text(size = 7), 
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.position = "none") +
  labs(x = "", 
       y = "% non-unique clones",
       title = "Mean Percent Non-Unique (bars) with Individual Values (dots) by Cell Type")


ggsave("../results/percent_non_unique.pdf", width = 6.5, height = 6, units = "cm", dpi = 300)
