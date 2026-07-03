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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd('/home/henry/Cytof_Analysis/T_Cell_Panel/aligned_data')
#~~~~~~~~~~
# Read in data and delete superfluous columns
aligned = fread('aligned_t_cell_panel.csv', sep = ',')
aligned = aligned[ , (22:52) := NULL ]
aligned = aligned[ , (33:63) := NULL ]
aligned = aligned[aligned[["Batch_Control"]] != 'TRUE',]


#~~~~~~~
# Run FlowSOM so that I can select populations of interest
aligned.cellular.cols = names(aligned)[c(33:36, 38:40, 42:49, 51:63)]
aligned = run.flowsom(dat = aligned, use.cols = aligned.cellular.cols, meta.k = 40) 
aligned_sub = do.subsample(aligned, targets= c(rep(2500, times = 118)), divide.by = 'FileName')
aligned_sub = run.umap(aligned_sub, aligned.cellular.cols)

cell_table <- table(aligned$patient_id, aligned$FlowSOM_metacluster)

# Write to CSV
write.csv(cell_table, "cell_counts_by_cluster.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd('/home/henry/Cytof_Analysis/T_Cell_Panel/aligned_data/clustering_cell_markers/CD8_T_cells') 
CD8_T_cells = aligned[FlowSOM_metacluster %in% c(12, 18, 19, 25, 26, 27, 30, 32, 36, 38)]
CD8_T_cells = CD8_T_cells[, c("FlowSOM_cluster","FlowSOM_metacluster"):=NULL] 
aligned.cellular.cols = names(CD8_T_cells)[c(35, 38:40, 42:44, 46:47, 49:52, 54:55, 57, 59, 63)]
# markers used in clustering
# CD183_CXCR3, CD3, CD39, CD314, CD194_CCR4, CD197_CCR7, CD127, CD57, Perforin, CD45RA, CD69, CD4, CD8, KLRG1, CD159c, CD45RO, Granzyme_B, CD56
# markers not used in clustering
#CD33, Ki67, CD279_PD1, CD11c, CD19, CD69, CD16, Tim3, CD38, CD14, HLA-DR, CD8a, CD25

visualization.cols = names(aligned)[c(33:63)]
CD8_T_cells = run.flowsom(dat = CD8_T_cells, use.cols = aligned.cellular.cols, meta.k = 40)
# Delete clusters of contaminating CD4 positive and CD8 negative cells
CD8_T_cells = CD8_T_cells[!(FlowSOM_metacluster %in% c(15, 32, 34))]
CD8_T_cells = CD8_T_cells[, c("FlowSOM_cluster","FlowSOM_metacluster"):=NULL] 
CD8_T_cells = run.flowsom(dat = CD8_T_cells, use.cols = aligned.cellular.cols, meta.k = 40)
CD8_T_cells_sub = do.subsample(CD8_T_cells, targets= c(rep(500, times = 118)), divide.by = 'FileName')
CD8_T_cells_sub = run.umap(CD8_T_cells_sub, aligned.cellular.cols)
make.colour.plot(CD8_T_cells_sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", title = 'FlowSOM_metacluster', dot.size =0.25, col.type = 'factor')
make.multi.plot(CD8_T_cells_sub, "UMAP_X", "UMAP_Y", visualization.cols, dot.size = 0.15, colours = 'viridis')

CD8_T_cells_sub_1_10 = CD8_T_cells_sub[FlowSOM_metacluster %in% c(1:10)]
make.colour.plot(CD8_T_cells_sub_1_10, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", title = 'FlowSOM_metaclusters 1 to 10', dot.size =0.25, col.type = 'factor')

CD8_T_cells_sub_11_20 = CD8_T_cells_sub[FlowSOM_metacluster %in% c(11:20)]
make.colour.plot(CD8_T_cells_sub_11_20, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", title = 'FlowSOM_metaclusters 11 to 20', dot.size =0.25, col.type = 'factor')

CD8_T_cells_sub_21_30 = CD8_T_cells_sub[FlowSOM_metacluster %in% c(21:30)]
make.colour.plot(CD8_T_cells_sub_21_30, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", title = 'FlowSOM_metaclusters 21 to 30', dot.size =0.25, col.type = 'factor')

CD8_T_cells_sub_31_40 = CD8_T_cells_sub[FlowSOM_metacluster %in% c(31:40)]
make.colour.plot(CD8_T_cells_sub_31_40, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", title = 'FlowSOM_metaclusters 31 to 40', dot.size =0.25, col.type = 'factor')



cell_table <- table(CD8_T_cells$patient_id, CD8_T_cells$FlowSOM_metacluster)

# Write to CSV
write.csv(cell_table, "cd8_t_cell_counts_by_cluster.csv")
