# Examine ERV expression in T cells from patients with MDS
library(Seurat)
library(dplyr)
library(ggplot2)
library(gghalves)
library(ggpattern)
library(svglite)
library(tidyr)

# Set the working directory

setwd('/Users/henryhampton/OneDrive - UNSW/2026_06_17_Code_for_upload/T_cell_figures/code')

# Read in the TE annotation file which Anne sent me
df <- read.csv('../data/TE_family.csv',  header = TRUE)

# Read in the Seurat object whcih Anne sent me
obj<-readRDS("../data/combined_scTE_all_patient_updated.RDS")

# Convert the seurat object to a Seurat v5 object
current_idents <- Idents(obj)
names(current_idents) <- colnames(obj)
Idents(obj) <- current_idents
obj[["RNA"]] <- as(object = obj[["RNA"]], Class = "Assay5")

if(!"TE_group" %in% colnames(obj@meta.data)) {
  obj$TE_group <- NA
}

obj@meta.data[obj$sample%in%c('P12_C1D1','P03_C1D1','P01_C1D1',
                              'P11_C1D1','P17_C1D1',
                              'P18_C1D1'),]$TE_group<-"pre_R"
obj@meta.data[obj$sample%in%c('P12_C6D8','P12_C7D1','P12_C7D22','P03_C7D22',
                              'P03_C12D29','P09_C12D29','P01_C7D1',
                              'P01_C12D29','P17_C7D1','P18_C7D1','P11_C6D8'),]$TE_group<-"post_R"
obj@meta.data[obj$sample%in%c('P09_C1D8','P02_C1D1'),]$TE_group<-"pre_NR"
obj@meta.data[obj$sample%in%c('P09_C7D1','P02_C7D1','P02_Progression',
                              'P17_Progression','P18_C7D22'),]$TE_group<-"post_NR"


sine <- df[df$V2 == "SINE", "V1"]
line <- df[df$V2 == "LINE", "V1"]
ltr <- df[df$V2 == "LTR", "V1"]

te_gene_sets <- list(
  sine_score = sine,
  line_score = line,
  ltr_score = ltr
)

# Add all scores in one call
obj <- AddModuleScore(
  object = obj,
  features = te_gene_sets,
  name = names(te_gene_sets)
)


celltypes_to_plot <- c("GZMH+ Effector_Memory1", "TEMRA_NK-like" ,"GZMK+ Effector_Memory" ,     
                   "GZMH+ Effector_Memory2", "TEMRA")

# Filter metadata
plot_data <- obj@meta.data[obj@meta.data$cell_type1 %in% celltypes_to_plot, ]
plot_data <- plot_data[!is.na(plot_data$TE_group), ]

plot_data$group_order <- paste(plot_data$cell_type1, plot_data$TE_group, sep = "_")
plot_data$group_order <- gsub('_pre', '', plot_data$group_order)
plot_data$group_order <- gsub('_post', '', plot_data$group_order)
plot_data$group_order <- factor(plot_data$group_order, 
                                levels = c("GZMH+ Effector_Memory1_NR", "GZMH+ Effector_Memory1_R", 
                                           "TEMRA_NK-like_NR", "TEMRA_NK-like_R",
                                           "GZMK+ Effector_Memory_NR", "GZMK+ Effector_Memory_R" ,     
                                           "GZMH+ Effector_Memory2_NR", "GZMH+ Effector_Memory2_R", 
                                           "TEMRA_NR", "TEMRA_R"))

plot_data$timepoint_group <- ifelse(plot_data$TE_group %in% c('pre_NR', 'pre_R'), 'Pre', 'Post')


# Define the scores to analyze uasing a statistical test
scores_to_analyze <- c("sine_score1", "line_score2",  "ltr_score3")

# Create a list to store all results
all_ttest_results <- list()

getwd()


library(ggpattern)


score_labels <- c(
  "sine_score1" = "SINE",
  "line_score2" = "LINE"
)

for (ct in celltypes_to_plot) {
  
  cat("\n========== Plotting", ct, "==========\n")
  ct_data <- plot_data[plot_data$cell_type1 == ct, ]
  
  # Wilcoxon tests: SINE and LINE only
  ttest_df <- lapply(c("sine_score1", "line_score2"), function(s) {
    lapply(c("Responder", "Non-Responder"), function(out) {
      grp  <- ifelse(out == "Responder", "_R$", "_NR$")
      pre  <- ct_data[grepl(grp, ct_data$TE_group) & ct_data$timepoint_group == "Pre",  s]
      post <- ct_data[grepl(grp, ct_data$TE_group) & ct_data$timepoint_group == "Post", s]
      if (length(pre) > 1 & length(post) > 1) {
        res <- wilcox.test(pre, post)
        data.frame(
          TE_score  = factor(s, levels = c("sine_score1", "line_score2")),
          outcome   = factor(out, levels = c("Responder", "Non-Responder")),
          p_value   = res$p.value,
          mean_pre  = mean(pre,  na.rm = TRUE),
          mean_post = mean(post, na.rm = TRUE)
        )
      }
    }) %>% bind_rows()
  }) %>% bind_rows() %>%
    mutate(
      p_adj      = p.adjust(p_value, method = "BH"),
      signif     = case_when(
        p_adj < 0.001 ~ "*\n*\n*",
        p_adj < 0.01  ~ "*\n*",
        p_adj < 0.05  ~ "*",
        TRUE          ~ ""
      ),
      star_color = ifelse(p_adj < 0.05,
                          ifelse(mean_post > mean_pre, "Post", "Pre"),
                          "Not significant"),
      # Stars centred over outcome group integer positions (1 and 3)
      x_numeric  = ifelse(outcome == "Responder", 1, 3)
    )
  
  # Outcome groups spaced at x_int = 1 (Responder) and 3 (Non-Responder)
  # This gives a clear gap and prevents any rect/box overlap between groups.
  # Within each group: Pre = x_int - 0.4, Post = x_int + 0.4
  #   R Pre = 0.6,  R Post = 1.4
  #   NR Pre = 2.6, NR Post = 3.4
  ct_long <- ct_data %>%
    mutate(
      outcome         = factor(ifelse(grepl("_R$", TE_group), "Responder", "Non-Responder"),
                               levels = c("Responder", "Non-Responder")),
      timepoint_group = factor(timepoint_group, levels = c("Pre", "Post"))
    ) %>%
    pivot_longer(cols = all_of(c("sine_score1", "line_score2")),
                 names_to = "TE_score", values_to = "score_value") %>%
    mutate(
      TE_score  = factor(TE_score, levels = c("sine_score1", "line_score2")),
      x_int     = ifelse(outcome == "Responder", 1, 3),
      x_numeric = x_int + ifelse(timepoint_group == "Pre", -0.4, 0.4)
    )
  
  ct_long_pre  <- ct_long %>% filter(timepoint_group == "Pre")
  ct_long_post <- ct_long %>% filter(timepoint_group == "Post")
  
  # Whisker-based y limits
  whisker_sl <- ct_long %>%
    group_by(TE_score, outcome, timepoint_group) %>%
    summarise(
      q1    = quantile(score_value, 0.25, na.rm = TRUE),
      q3    = quantile(score_value, 0.75, na.rm = TRUE),
      iqr   = IQR(score_value, na.rm = TRUE),
      lower = min(score_value[score_value >= q1 - 1.5 * iqr], na.rm = TRUE),
      upper = max(score_value[score_value <= q3 + 1.5 * iqr], na.rm = TRUE),
      .groups = "drop"
    )
  y_min_sl    <- min(whisker_sl$lower, na.rm = TRUE)
  y_max_sl    <- max(whisker_sl$upper, na.rm = TRUE)
  y_buffer_sl <- (y_max_sl - y_min_sl) * 0.03
  
  # Background rects: span full width of each outcome group (Pre + Post)
  # R:  x_int=1 → Pre centre 0.6, Post centre 1.4 → xmin=0.3, xmax=1.7
  # NR: x_int=3 → Pre centre 2.6, Post centre 3.4 → xmin=2.3, xmax=3.7
  # No overlap between groups (gap from 1.7 to 2.3)
  bg_df <- tidyr::expand_grid(
    tibble(x_int = c(1, 3)),
    TE_score = factor(c("sine_score1", "line_score2"), levels = c("sine_score1", "line_score2"))
  ) %>%
    mutate(
      xmin = (x_int - 0.4) - 0.4,   # = x_int - 0.7
      xmax = (x_int - 0.4) + 0.4    # = x_int - 0.1
    )
  
  p <- ggplot() +
    # 1. Crosshatch background — white fill so pattern is visible on white canvas
    geom_rect( data = bg_df, aes(xmin = xmin, xmax = xmax, 
                                 ymin = -Inf, ymax = Inf), 
               inherit.aes = FALSE, 
               fill = "#e5e5e5", colour = NA ) +
    # 2. Pre boxes: semi-transparent so crosshatch shows through
    geom_boxplot(
      data          = ct_long_pre,
      aes(x = x_numeric, y = score_value, fill = outcome, group = outcome),
      width         = 0.6,
      alpha         = 0.5,
      outlier.shape = NA,
      colour        = "grey30",
      linewidth     = 0.3
    ) +
    # 3. Post boxes: fully opaque, covers crosshatch behind them
    geom_boxplot(
      data          = ct_long_post,
      aes(x = x_numeric, y = score_value, fill = outcome, group = outcome),
      width         = 0.6,
      alpha         = 1.0,
      outlier.shape = NA,
      colour        = "grey30",
      linewidth     = 0.3
    ) +
    # 4. Significance stars centred between each Pre/Post pair
    geom_text(
      data        = ttest_df,
      aes(x = x_numeric, y = Inf, label = signif, color = star_color),
      inherit.aes = FALSE,
      vjust       = 1.5,
      size        = 4
    ) +
    facet_wrap(~ TE_score, scales = "fixed",
               labeller = labeller(TE_score = score_labels)) +
    scale_fill_manual(values  = c("Responder" = "#0e5485", "Non-Responder" = "#f79120")) +
    scale_color_manual(values = c("Pre" = "#f7912080", "Post" = "#0e5485", "Not significant" = "black")) +
    # Breaks at group centres 1 and 3
    scale_x_continuous(breaks = c(1, 3), labels = c("R", "NR"),
                       limits = c(0.1, 4.1)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    coord_cartesian(ylim = c(y_min_sl - y_buffer_sl, y_max_sl + y_buffer_sl)) +
    labs(title = ct, y = "TE Score", x = NULL) +
    theme_classic(base_size = 7) +
    theme(
      axis.text.x      = element_text(angle = 45, vjust = 1, hjust = 1),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold"),
      legend.position  = "none"
    )
  
  print(p)
  
  ggsave(paste0('../results/', ct, "_TE_scores_boxplot_SINE_LINE.svg"), p,
         width = 1.75, height = 1.5, device = svglite)
}

library(ggpattern)

for (ct in celltypes_to_plot) {
  
  cat("\n========== Plotting", ct, "==========\n")
  ct_data <- plot_data[plot_data$cell_type1 == ct, ]
  
  # Wilcoxon tests: SINE and LINE only
  ttest_df <- lapply(c("ltr_score3"), function(s) {
    lapply(c("Responder", "Non-Responder"), function(out) {
      grp  <- ifelse(out == "Responder", "_R$", "_NR$")
      pre  <- ct_data[grepl(grp, ct_data$TE_group) & ct_data$timepoint_group == "Pre",  s]
      post <- ct_data[grepl(grp, ct_data$TE_group) & ct_data$timepoint_group == "Post", s]
      if (length(pre) > 1 & length(post) > 1) {
        res <- wilcox.test(pre, post)
        data.frame(
          TE_score  = factor(s, levels = c("ltr_score3")),
          outcome   = factor(out, levels = c("Responder", "Non-Responder")),
          p_value   = res$p.value,
          mean_pre  = mean(pre,  na.rm = TRUE),
          mean_post = mean(post, na.rm = TRUE)
        )
      }
    }) %>% bind_rows()
  }) %>% bind_rows() %>%
    mutate(
      p_adj      = p.adjust(p_value, method = "BH"),
      signif     = case_when(
        p_adj < 0.001 ~ "*\n*\n*",
        p_adj < 0.01  ~ "*\n*",
        p_adj < 0.05  ~ "*",
        TRUE          ~ ""
      ),
      star_color = ifelse(p_adj < 0.05,
                          ifelse(mean_post > mean_pre, "Post", "Pre"),
                          "Not significant"),
      # Stars centred over outcome group integer positions (1 and 3)
      x_numeric  = ifelse(outcome == "Responder", 1, 3)
    )
  
  # Outcome groups spaced at x_int = 1 (Responder) and 3 (Non-Responder)
  # This gives a clear gap and prevents any rect/box overlap between groups.
  # Within each group: Pre = x_int - 0.4, Post = x_int + 0.4
  #   R Pre = 0.6,  R Post = 1.4
  #   NR Pre = 2.6, NR Post = 3.4
  ct_long <- ct_data %>%
    mutate(
      outcome         = factor(ifelse(grepl("_R$", TE_group), "Responder", "Non-Responder"),
                               levels = c("Responder", "Non-Responder")),
      timepoint_group = factor(timepoint_group, levels = c("Pre", "Post"))
    ) %>%
    pivot_longer(cols = all_of(c("ltr_score3")),
                 names_to = "TE_score", values_to = "score_value") %>%
    mutate(
      TE_score  = factor(TE_score, levels = c("ltr_score3")),
      x_int     = ifelse(outcome == "Responder", 1, 3),
      x_numeric = x_int + ifelse(timepoint_group == "Pre", -0.4, 0.4)
    )
  
  ct_long_pre  <- ct_long %>% filter(timepoint_group == "Pre")
  ct_long_post <- ct_long %>% filter(timepoint_group == "Post")
  
  # Whisker-based y limits
  whisker_sl <- ct_long %>%
    group_by(TE_score, outcome, timepoint_group) %>%
    summarise(
      q1    = quantile(score_value, 0.25, na.rm = TRUE),
      q3    = quantile(score_value, 0.75, na.rm = TRUE),
      iqr   = IQR(score_value, na.rm = TRUE),
      lower = min(score_value[score_value >= q1 - 1.5 * iqr], na.rm = TRUE),
      upper = max(score_value[score_value <= q3 + 1.5 * iqr], na.rm = TRUE),
      .groups = "drop"
    )
  y_min_sl    <- min(whisker_sl$lower, na.rm = TRUE)
  y_max_sl    <- max(whisker_sl$upper, na.rm = TRUE)
  y_buffer_sl <- (y_max_sl - y_min_sl) * 0.03
  
  # Background rects: span full width of each outcome group (Pre + Post)
  # R:  x_int=1 → Pre centre 0.6, Post centre 1.4 → xmin=0.3, xmax=1.7
  # NR: x_int=3 → Pre centre 2.6, Post centre 3.4 → xmin=2.3, xmax=3.7
  # No overlap between groups (gap from 1.7 to 2.3)
  bg_df <- tidyr::expand_grid(
    tibble(x_int = c(1, 3)),
    TE_score = factor(c("ltr_score3"), levels = c("ltr_score3"))
  ) %>%
    mutate(
      xmin = (x_int - 0.4) - 0.4,   # = x_int - 0.7
      xmax = (x_int - 0.4) + 0.4    # = x_int - 0.1
    )
  
  p <- ggplot() +
    # 1. Crosshatch background — white fill so pattern is visible on white canvas
    geom_rect(
      data        = bg_df,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill        = "grey90",
      colour      = NA,
      alpha       = 0.25
    ) +
    
    
    # 2. Pre boxes: semi-transparent so crosshatch shows through
    geom_boxplot(
      data          = ct_long_pre,
      aes(x = x_numeric, y = score_value, fill = outcome, group = outcome),
      width         = 0.6,
      alpha         = 0.5,
      outlier.shape = NA,
      colour        = "grey30",
      linewidth     = 0.3
    ) +
    # 3. Post boxes: fully opaque, covers crosshatch behind them
    geom_boxplot(
      data          = ct_long_post,
      aes(x = x_numeric, y = score_value, fill = outcome, group = outcome),
      width         = 0.6,
      alpha         = 1.0,
      outlier.shape = NA,
      colour        = "grey30",
      linewidth     = 0.3
    ) +
    # 4. Significance stars centred between each Pre/Post pair
    geom_text(
      data        = ttest_df,
      aes(x = x_numeric, y = Inf, label = signif, color = star_color),
      inherit.aes = FALSE,
      vjust       = 1.5,
      size        = 4
    ) +
    facet_wrap(~ TE_score, scales = "fixed",
               labeller = labeller(TE_score = score_labels)) +
    scale_fill_manual(values  = c("Responder" = "#0e5485", "Non-Responder" = "#f79120")) +
    scale_color_manual(values = c("Pre" = "#f7912080", "Post" = "#0e5485", "Not significant" = "black")) +
    # Breaks at group centres 1 and 3
    scale_x_continuous(breaks = c(1, 3), labels = c("R", "NR"),
                       limits = c(0.1, 4.1)) +
    scale_y_continuous(
      expand = expansion(mult = c(0.02, 0.02)),
      labels = scales::label_number(accuracy = 0.1)
    ) +
    coord_cartesian(ylim = c(y_min_sl - y_buffer_sl, y_max_sl + y_buffer_sl)) +
    labs(title = ct, y = "TE Score", x = NULL) +
    theme_classic(base_size = 7) +
    theme(
      axis.text.x      = element_text(angle = 45, vjust = 1, hjust = 1),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold"),
      legend.position  = "none"
    )
  
  print(p)
  
  ggsave(paste0('../results/', ct, "_TE_scores_boxplot_LTR.svg"), p,
         width = 1.75, height = 1.5, device = svglite)
}


