library(dplyr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(readxl)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
setwd('/Users/henryhampton/Library/CloudStorage/OneDrive-UNSW/2026_06_17_Code_for_upload/T_cell_figures/code')

# Read in files
tcr = read.csv('../data/human_reactive_tcrs_summary_table.csv')

# Read in a sheet to specify the color of ea
color = read_excel('../data/patient_color_scheme.xlsx')
color_dict = setNames(color$color, color$patient)

# 
# combined_data$source <- factor(combined_data$source, 
#                                levels = c("vdj_db", "ie_db"))

# Create dot plot with jitter and white fill with colored edges
p = ggplot(combined_data, aes(x = factor(source, levels = c("vdj_db", "ie_db")), 
                              y = percent_non_na, color = patient)) +
  geom_jitter(size = 1, shape = 21, fill = "white", stroke = 1,
              width = 0.1, height = 0) +  # width controls x-jitter, height controls y-jitter
  scale_color_manual(values = color_dict) +
  theme_classic() +
  labs(color = "Patient")

ggsave(filename = '../results/percent_self_reactive.pdf',
       plot = p,
       scale = 1,
       width = 7,
       height = 5,
       units = c( "cm"),
       dpi = 600)
p

mean(combined_data$percent_non_na[combined_data$source == "ie_db"])

mean(combined_data$percent_non_na[combined_data$source == "vdj_db"])
