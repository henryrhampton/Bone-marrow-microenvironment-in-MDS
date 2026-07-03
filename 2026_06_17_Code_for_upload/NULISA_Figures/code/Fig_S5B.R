library(ggplot2)
library(edgeR) #create DElist to calculate FPKM
library(Signac);library(readxl)
library(RColorBrewer)
library(reticulate)
library(dplyr);library(ggpubr)
library(stringr)
library("VennDiagram")
library(glmnet)
library(ComplexHeatmap)
library(circlize)
library(ggthemes)
library(ggrepel)
###############################################################
########## NULISA data for MDS patients from Sydeny
########## Author Bofei; Sept 2024
##############################################################
setwd('OneDrive - UNSW/2026_06_17_Code_for_upload/NULISA_Figures/code/')
npq_mds=read_excel("../data/NULISA_for_upload.xlsx", sheet = "NPQ") %>% as.data.frame()
npq_mds=npq_mds[,1:82]
npq_mds = t(npq_mds)
colnames(npq_mds) = npq_mds[1,]
npq_mds = npq_mds[-1,]
# Read in the cytokine pairs and create a named vector
cytokine_pairs = read_excel('../data/2026_02_25_receptor_cytokine_combinations.xlsx',
                            sheet = 'receptor_cytokine_combo')


named_vec <- setNames(cytokine_pairs$Receptor, cytokine_pairs$Cytokine)

rn <- rownames(npq_mds)
npq_mds <- apply(npq_mds, 2, as.numeric)
rownames(npq_mds) <- rn

npq_mds <- as.data.frame(npq_mds)

for (i in seq_along(named_vec)) {
  ligand   <- names(named_vec)[i]
  receptor <- named_vec[[i]]
  new_col  <- paste0(ligand, "_", receptor)
  
  npq_mds[, new_col] <- npq_mds[, ligand] - npq_mds[, receptor]
}

npq_mds = t(npq_mds)


clinical=read.csv("../data/unsw_australia_plasma_protein_concentration_patient_updated.csv") %>% as.data.frame()
clinical=clinical[1:46,];colnames(clinical)[1]="sample"
clinical[,1]=paste0(clinical$PID,"_",clinical$Cycle_day)
##There are bone marrow and pb samples, need to analyze separately
index=unlist(lapply(colnames(npq_mds),function(x){endsWith(x,"_PB")}))
mds_pb=npq_mds[,index];mds_bm=npq_mds[,!index]

mds_pb=as.data.frame(t(mds_pb));mds_bm=as.data.frame(t(mds_bm))


###add patient ID and cycle to the column
mds_pb$PID=unlist(lapply(rownames(mds_pb),function(x){substr(x,1,3)}))
mds_bm$PID=unlist(lapply(rownames(mds_bm),function(x){substr(x,1,3)}))
mds_pb$cycle=unlist(lapply(rownames(mds_pb),function(x){substr(x,5,9)}))
mds_bm$cycle=unlist(lapply(rownames(mds_bm),function(x){ifelse(nchar(x)>14,substr(x,5,11),substr(x,5,9))}))
colnames(mds_pb)[c(99,102,122,174)]=c("HLA_DRA","IFNA1_IFNA13","IL17A_IL17F","LTA_LTB")
colnames(mds_bm)[c(99,102,122,174)]=c("HLA_DRA","IFNA1_IFNA13","IL17A_IL17F","LTA_LTB")

my_dict <- setNames(clinical$outcome_6, clinical$PID)
mds_pb['outcome_6'] = my_dict[mds_pb$PID]
mds_bm['outcome_6'] = my_dict[mds_bm$PID]
mds_pb=mds_pb[mds_pb$cycle != "C12_D9",]


folder_path <- '../data/differential_expression_tables/'
list.files(folder_path, pattern = "*.csv")



# Specify the files to make a heatmap
files_info <- data.frame(
  filename = c("blood_c1d1_r_v_nr_cytokine_receptor_pairs.csv", "blood_c7d1_r_v_nr_cytokine_receptor_pairs.csv", 
               "bone_marrow_c1d1_r_v_nr_cytokine_receptor_pairs.csv", "bone_marrow_c7d1_r_v_nr_cytokine_receptor_pairs.csv"),
  data_source = c("mds_pb", "mds_pb", "mds_bm", "mds_bm"),
  cycle = c("C1_D1", "C7_D1", "C1_D1", "C7_D1"),
  title = c("Blood C1D1", "Blood C7D1", "Bone Marrow C1D1", "Bone Marrow C7D1")
)

files_info

# Loop through each file
for(i in 1:nrow(files_info)) {
  cat("\n=== Creating heatmap for:", files_info$title[i], "===\n")
  
  # Load differential expression data
  file_path <- file.path(folder_path, files_info$filename[i])
  deg_data <- read.csv(file_path)
  deg_data <- deg_data[deg_data$adj.P.Val < '0.05', ]
  deg_data <- deg_data[grepl('_', deg_data$Gene), ]
  deg_data <- deg_data[deg_data$Gene != 'IFNA1_IFNA13', ]
  deg_data <- deg_data[order(deg_data$logFC, decreasing = FALSE), ]
  
  # Select data source
  if(files_info$data_source[i] == "mds_pb") {
    mds_data <- mds_pb
  } else {
    mds_data <- mds_bm
  }
  
  # Filter for cycle and create heatmap
  cycle_filter <- mds_data$cycle == files_info$cycle[i]
  response_data <- mds_data$outcome_6[cycle_filter]
  responder_first <- order(response_data)
  
  # Create matrix and scaling
  mat <- mds_data[cycle_filter, deg_data$Gene]
  mat_ordered <- mat[responder_first, ]
  mat_scaled <- t(scale(mat_ordered))
  
  # Apply aliases to columns
  sample_names_ordered <- rownames(mat_ordered)
  # pids_ordered <- sapply(sample_names_ordered, extract_pid)
  # aliases_ordered <- pid_to_alias[pids_ordered]
  # aliases_ordered[is.na(aliases_ordered)] <- sample_names_ordered[is.na(aliases_ordered)]
  # colnames(mat_scaled) <- aliases_ordered
  
  # Create annotation and heatmap
  column_ha <- HeatmapAnnotation(
    Response = response_data[responder_first],
    col = list(Response = c("non-responder_2" = '#fe9003', "responder_1" = '#115284')),
    height = unit(2.5, "mm"), show_legend = FALSE
  )
  
  p <- Heatmap(
    as.matrix(mat_scaled), 
    name = "Gene Z-score",
    col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
    cluster_columns = TRUE, 
    cluster_rows = TRUE, 
    column_names_rot = 90, 
    top_annotation = column_ha,
    show_column_names = TRUE,
    column_title = files_info$title[i],
    show_heatmap_legend = FALSE
    
  )
  plot(p)
  #draw(p, annotation_legend_side = 'right', heatmap_legend_side = "right")
  #Create output filename
  output_filename <- paste0("heatmap_", gsub("\\.csv$", "", files_info$filename[i]), "_cytokine_receptor.pdf")
  output_path <- file.path("../graphs/nulisa_heatmaps/", output_filename)
  
  # Create directory if it doesn't exist
  if(!dir.exists("heatmaps")) dir.create("heatmaps")
  
  # Save ComplexHeatmap
  min_height <- 4  # Minimum height in inches
  height_per_row <- 0.125  # Reduce from 0.25 to 0.125
  calculated_height <- max(min_height, height_per_row * nrow(deg_data))
  
  pdf(output_path, width = 5, height = calculated_height)
  draw(p, annotation_legend_side = 'right', heatmap_legend_side = "right",
       padding = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) # Use draw() for ComplexHeatmap objects
  dev.off()
  
  cat("Saved heatmap to:", output_path, "\n")
}



