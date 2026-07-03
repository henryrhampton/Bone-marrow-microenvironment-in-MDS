library(limma)
library(reshape2)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)
library(biomaRt)
library(readxl)
library(edgeR)
library(Seurat)
library(tidyr)
library(SeuratDisk)
library(zellkonverter)
# 
setwd('/Users/henryhampton/OneDrive - UNSW/2026_06_17_Code_for_upload/T_cell_figures/code')


# Read in data
sce <- readH5AD("../data/t_cell_obj_v2.h5ad")
data <- as.Seurat(sce, counts = "counts", data = "logcounts")


rownames(data)[grep("^IFN", rownames(data))]


# Get counts matrix
counts <- GetAssayData(data, layer = "counts")

# Number of cells per cluster
table(data$leiden)

data@meta.data$IFN_group <- NA

data@meta.data[data$sample%in%c('P12_C1D1','P03_C1D1','P09_C1D8','P01_C1D1',
                                'P11_C1D1','P17_C1D1',
                                'P18_C1D1'),]$IFN_group<-"pre_R"
data@meta.data[data$sample%in%c('P12_C6D8','P12_C7D1','P12_C7D22','P03_C7D22',
                                'P03_C12D29','P09_C12D29','P01_C7D1',
                                'P01_C12D29','P17_C7D1','P18_C7D1','P11_C6D8'),]$IFN_group<-"post_R"
data@meta.data[data$sample%in%c('P09_C1D8','P02_C1D1',
                                'P18_C1D1'),]$IFN_group<-"pre_NR"
data@meta.data[data$sample%in%c('P09_C7D1','P02_C7D1','P02_Progression',
                                'P17_Progression','P18_C7D22'),]$IFN_group<-"post_NR"




celltypes_to_plot <- c("GZMK+ Effector_Memory", "GZMH+ Effector_Memory_1","GZMH+ Effector_Memory_2",
                       "TEMRA", "TEMRA_NK-like")

data@meta.data$group_order <- paste(data@meta.data$annotation, data@meta.data$IFN_group, sep = "_")

data@meta.data$timepoint_group <- ifelse(data@meta.data$IFN_group %in% c('pre_NR', 'pre_R'), 'Pre', 'Post')

ifn_genes = c("IFNG", "IFNL1")

data <- subset(data, cells = colnames(data)[!is.na(data$IFN_group)])


# 
# Build combined data across all cell types
all_long <- lapply(celltypes_to_plot, function(ct) {
  sub_obj  <- subset(data, subset = annotation == ct)
  expr_mat <- FetchData(sub_obj, vars = c(ifn_genes, "IFN_group", "timepoint_group"))

  ct_data <- expr_mat %>%
    filter(timepoint_group %in% c("Pre", "Post"),
           grepl("_R$|_NR$", IFN_group))

  if (nrow(ct_data) == 0) return(NULL)

  ct_data %>%
    mutate(outcome  = ifelse(grepl("_R$", IFN_group), "Responder", "Non-Responder"),
           celltype = ct) %>%
    pivot_longer(cols = all_of(ifn_genes), names_to = "gene", values_to = "expression")
}) %>% bind_rows()

all_long$gene            <- factor(all_long$gene,     levels = ifn_genes)
all_long$timepoint_group <- factor(all_long$timepoint_group, levels = c("Pre", "Post"))
all_long$outcome         <- factor(all_long$outcome,  levels = c("Responder", "Non-Responder"))
all_long$celltype        <- factor(all_long$celltype, levels = celltypes_to_plot)

# Build combined stats
all_stats <- lapply(celltypes_to_plot, function(ct) {
  ct_data <- all_long %>% filter(celltype == ct)
  lapply(ifn_genes, function(g) {
    lapply(c("Responder", "Non-Responder"), function(out) {
      pattern <- ifelse(out == "Responder", "_R$", "_NR$")
      pre  <- ct_data[grepl(pattern, ct_data$IFN_group) & ct_data$timepoint_group == "Pre"  & ct_data$gene == g, "expression", drop = TRUE]
      post <- ct_data[grepl(pattern, ct_data$IFN_group) & ct_data$timepoint_group == "Post" & ct_data$gene == g, "expression", drop = TRUE]
      if (length(pre) > 1 & length(post) > 1) {
        res <- wilcox.test(pre, post)
        data.frame(celltype = ct, gene = factor(g, levels = ifn_genes),
                   outcome = out, p_value = res$p.value,
                   mean_pre = mean(pre, na.rm = TRUE), mean_post = mean(post, na.rm = TRUE))
      }
    }) %>% bind_rows()
  }) %>% bind_rows()
}) %>% bind_rows() %>%
  mutate(
    p_adj      = p.adjust(p_value, method = "BH"),
    signif     = case_when(p_adj < 0.001 ~ "*\n*\n*", p_adj < 0.01 ~ "*\n*", p_adj < 0.05 ~ "*", TRUE ~ ""),
    star_color = ifelse(p_adj < 0.05, ifelse(mean_post > mean_pre, "Post", "Pre"), "Not significant")
  )

all_stats$celltype <- factor(all_stats$celltype, levels = celltypes_to_plot)

all_long$facet_label  <- paste(all_long$celltype,  all_long$gene,  sep = "\n")
all_stats$facet_label <- paste(all_stats$celltype, all_stats$gene, sep = "\n")

# Set factor order so celltypes are grouped together
facet_levels <- paste(rep(celltypes_to_plot, each = length(ifn_genes)),
                      rep(ifn_genes, times = length(celltypes_to_plot)), sep = "\n")
all_long$facet_label  <- factor(all_long$facet_label,  levels = facet_levels)
all_stats$facet_label <- factor(all_stats$facet_label, levels = facet_levels)



library(ggpattern)
bg_df <- tibble(
  outcome = factor(c("Responder", "Non-Responder"),
                   levels = c("Responder", "Non-Responder")),
  xmin    = c(0.8, 1.8) - 0.2,   # Pre centres: 0.8 (R), 1.8 (NR)
  xmax    = c(0.8, 1.8) + 0.2
)

# Create a for loop to plot the expression of IFNL1 and IFNG
for (gene in ifn_genes) {
  long_sub  <- all_long  %>% filter(gene == !!gene)
  stats_sub <- all_stats %>% filter(gene == !!gene)
  
  # Reset factor levels to control panel order
  long_sub$celltype  <- factor(long_sub$celltype,  levels = celltypes_to_plot)
  stats_sub$celltype <- factor(stats_sub$celltype, levels = celltypes_to_plot)
  
  
  
  p <- ggplot() +
    # 1. Hatched background behind Pre columns
    geom_rect(
      data        = bg_df,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill        = "#e5e5e5",   # pale grey — adjust to taste
      colour      = NA
    ) +
    # 2. Violins — both timepoints together as before
    geom_violin(
      data      = long_sub,
      aes(x = outcome, y = expression, fill = outcome, alpha = timepoint_group),
      position  = position_dodge(width = 0.8),
      width     = 0.6,
      scale     = "width",
      trim      = TRUE,
      linewidth = 0.2,
      colour    = "black"
    ) +
    # 3. Jitter
    geom_jitter(
      data     = long_sub,
      aes(x = outcome, y = expression,
          group = interaction(outcome, timepoint_group)),
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
      size     = 0.3,
      alpha    = 1.0,
      shape    = 21,
      fill     = NA,
      stroke   = 0.2
    ) +
    # 4. Significance stars
    geom_text(
      data        = stats_sub,
      aes(x = outcome, y = Inf, label = signif, color = star_color),
      inherit.aes = FALSE,
      vjust       = 1.5,
      size        = 4
    ) +
    facet_wrap(~ celltype, nrow = 1) +
    scale_fill_manual(values = c("Responder" = "#0e5485", "Non-Responder" = "#f79120")) +
    scale_alpha_manual(values = c("Pre" = 0.4, "Post" = 1.0)) +
    scale_color_manual(values = c("Pre" = "#f7912080", "Post" = "#0e5485",
                                  "Not significant" = "black")) +
    scale_x_discrete(labels = c("Responder" = "R", "Non-Responder" = "NR")) +
    labs(title = gene, y = "Normalised Expression", x = NULL) +
    theme_classic(base_size = 7) +
    theme(
      axis.text.x      = element_text(angle = 45, vjust = 1, hjust = 1),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold"),
      legend.position  = "none"
    )
  out = '../results/'
  ggsave(paste0(out, gene, "_expression.pdf"), plot = p, 
         width = 8, height = 3.8, units = "cm")
}
p
