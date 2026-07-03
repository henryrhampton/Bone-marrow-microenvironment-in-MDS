library(Spectre)
library(ggplot2)
library(CytoNorm)
library(flowCore)
library(data.table)
library(scales)
library(colorRamps)
library(ggthemes)
library(RColorBrewer)
library(ggpointdensity)
Spectre::package.load()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd('/home/henry/Cytof_Analysis/Myeloid_Panel/second_batch_alignment')
# Read in data
aligned_cells = fread('aligned_myeloid_panel_clustered.csv', sep = ',')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd('/home/henry/Cytof_Analysis/Myeloid_Panel/second_batch_alignment')
# Read in data
aligned_cells = fread('aligned_myeloid_panel.csv', sep = ',')

# Delete superfluous columns these are batch control and columns that have been transformed but not aligned
aligned_cells = aligned_cells[ , (24:53) := NULL ]
aligned_cells = aligned_cells[ , (36:62) := NULL ]
aligned_cells <- aligned_cells[aligned_cells[["Batch_Control"]] != 'Batch_Control',]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compare using all markers versus cellular markers
# Make graphs using cellular markers
setwd("/home/henry/Cytof_Analysis/Myeloid_Panel/second_batch_alignment/comparision_cellular_cluster_markers/cellular_markers")
# Note I deleted Ki67, PD1, CD141 and CD47 from clustering
visualisation.cols = names(aligned_cells)[c(36:62)]
aligned.cellular.cols = names(aligned_cells)[c(36:39, 41:44, 46:58, 61:62)]
aligned_sub <- do.subsample(aligned_cells, targets= c(rep(2000, times = 118)), divide.by = 'FileName')
aligned_sub <- run.flowsom(dat = aligned_sub, use.cols = aligned.cellular.cols, meta.k = 20) 
aligned_sub <- run.umap(aligned_sub, aligned.cellular.cols)

make.colour.plot(aligned_sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", col.type = 'factor', add.label = TRUE)
make.multi.plot(aligned_sub, "UMAP_X", "UMAP_Y", visualisation.cols, dot.size = 0.25, colours = 'viridis')

cell_table <- table(aligned_cells$FileName, aligned_cells$FlowSOM_metacluster)

# Write to CSV
write.csv(cell_table, "myeloid_cell_counts_by_cluster.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd('/home/henry/Cytof_Analysis/Myeloid_Panel/second_batch_alignment/HSPC_analysis')
# Isolate  HSPCs
CD34 = subset(aligned_cells, FlowSOM_metacluster %in% c(1,2))

# remove sample due to a low number of CD34+ cells
CD34 = CD34[FileName != 'Control_12_Myeloid_Panel',] 
CD34 = CD34[FileName != 'P02_C1_D1_Myeloid_Panel',]
CD34 = CD34[FileName != 'Control_6_HSA1385_Myeloid_Panel',]
CD34 = CD34[FileName != 'P06_C7_D1_Myeloid_Panel',]
CD34 = CD34[ , (64:65) := NULL ]

aligned.cellular.cols = names(CD34)[c(36:39, 42:43, 46, 48, 51, 54,55, 57, 58, 61:62)]
# CD117,CD172ab_SIRPab, CD34, CD38, CD123, CD68, CD33, CD16, CD64, CD11b, LILRB4_FITC, CD45RA, CD40, HLA-DR, CD56 were used to cluster cells
visualization.cols = names(CD34)[c(36:62)]

CD34 = run.flowsom(dat = CD34, use.cols = aligned.cellular.cols, meta.k = 20)
CD34_sub <- do.subsample(CD34, targets= c(rep(121, times = 115)), divide.by = 'FileName')
CD34_sub <- run.umap(CD34_sub, aligned.cellular.cols)

make.colour.plot(CD34_sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", 
		col.type = 'factor', dot.size = 0.1, add.label = TRUE, title = 'CD34_FlowSOM_Metaclusters')
make.multi.plot(CD34_sub, "UMAP_X", "UMAP_Y", visualization.cols, 
		dot.size = 0.25, colours = 'viridis')	


cell_table <- table(CD34$FileName, CD34$FlowSOM_metacluster)

# Write to CSV
write.csv(cell_table, "CD34_cell_counts_by_cluster.csv")
