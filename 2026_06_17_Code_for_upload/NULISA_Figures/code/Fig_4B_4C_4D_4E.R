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

#~~~~~~~~
setwd("/Users/henryhampton/Library/CloudStorage/OneDrive-UNSW/2026_06_17_Code_for_upload/NULISA_Figures/code/")

npq_mds=read_excel("../data/NULISA_for_upload.xlsx", sheet = "NPQ") %>% as.data.frame()
npq_mds=npq_mds[,1:82]

# Read in the clinical data associated with the NULISA measurements
clinical=read.csv("../data/unsw_australia_plasma_protein_concentration_patient_updated.csv") %>% as.data.frame()
clinical=clinical[1:46,];colnames(clinical)[1]="sample"
clinical[,1]=paste0(clinical$PID,"_",clinical$Cycle_day)


##There are bone marrow and pb samples, these need to be analyzed separately
index=unlist(lapply(colnames(npq_mds),function(x){endsWith(x,"_PB")}))
mds_pb=npq_mds[,index];mds_bm=npq_mds[,!index]
rownames(mds_pb)=npq_mds$targetName
rownames(mds_bm)=mds_bm$targetName;mds_bm=mds_bm[,2:47]
mds_pb=as.data.frame(t(mds_pb));mds_bm=as.data.frame(t(mds_bm))

# ## Add patient ID and cycle to the column
extract_pid <- function(x) sub("_C\\d+_D\\d+.*", "", x) # Grab the patient id and remove the timepoint
mds_pb$PID <- extract_pid(rownames(mds_pb))
mds_bm$PID <- extract_pid(rownames(mds_bm))

mds_pb$cycle=unlist(lapply(rownames(mds_pb),function(x){substr(x,5,9)}))
mds_bm$cycle=unlist(lapply(rownames(mds_bm),function(x){ifelse(nchar(x)>5,substr(x,5,nchar(x)-4),substr(x,5,15))}))

colnames(mds_pb)[c(99,102,122,174)]=c("HLA_DRA","IFNA1_IFNA13","IL17A_IL17F","LTA_LTB")
colnames(mds_bm)[c(99,102,122,174)]=c("HLA_DRA","IFNA1_IFNA13","IL17A_IL17F","LTA_LTB")
mds_bm$sample = rownames(mds_bm)
mds_bm$cycle

# Add clinical info
my_dict <- setNames(clinical$outcome_6, clinical$PID)

mds_pb['outcome_6'] = my_dict[mds_pb$PID]
mds_bm['outcome_6'] = my_dict[mds_bm$PID]

clinical$log2_protein <- log2(clinical$Protein_concentration..mg.ml.)

folder_path <- '../data/differential_expression_tables/'

list.files(folder_path, pattern = "*.csv")

# Specify the files to make a heatmap
files_info <- data.frame(
  filename = c("blood_c7d1_r_v_nr.csv"),
  data_source = c("mds_pb"),
  cycle = c("C7_D1"),
  title = c("Blood C7D1")
)

library(readxl)

# Read the molecule annotation once (outside the loop)
molecule_annotation <- read_excel("../data/molecule_color_annotation.xlsx")

# Specify the files to make a heatmap
files_info <- data.frame(
  filename = c("blood_c1d1_r_v_nr.csv", "blood_c7d1_r_v_nr.csv", 
               "bm_protein_protein_corrected_c1d1_r_v_nr.csv", "bm_protein_protein_corrected_c7d1_r_v_nr.csv"),
  data_source = c("mds_pb", "mds_pb", "mds_bm", "mds_bm"),
  cycle = c("C1_D1", "C7_D1", "C1_D1", "C7_D1"),
  title = c("Blood C1D1", "Blood C7D1", "Bone Marrow C1D1", "Bone Marrow C7D1")
)



# Loop through each file
for(i in 1:nrow(files_info)) {
  cat("\n=== Creating heatmap for:", files_info$title[i], "===\n")
  
  # Load differential expression data
  file_path <- file.path(folder_path, files_info$filename[i])
  deg_data <- read.csv(file_path)
  deg_data <- deg_data[deg_data$P.Value < 0.05, ]
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
  rownames(mat_ordered) <- substr(rownames(mat_ordered), 1, 3)
  mat_scaled <- scale(mat_ordered)  # REMOVE t() here - don't transpose
  
  # Print range before capping
  cat("Before capping:", range(mat_scaled, na.rm = TRUE), "\n")
  
  # Cap values at -3 to 3
  mat_scaled[mat_scaled < -3] <- -3
  mat_scaled[mat_scaled > 3] <- 3
  
  # Print range after capping
  cat("After capping:", range(mat_scaled, na.rm = TRUE), "\n")
  
  # Apply aliases to rows (instead of columns)
  #sample_names_ordered <- rownames(mat_ordered)
  #pids_ordered <- sapply(sample_names_ordered, extract_pid)
  # aliases_ordered <- pid_to_alias[pids_ordered]
  # aliases_ordered[is.na(aliases_ordered)] <- sample_names_ordered[is.na(aliases_ordered)]
  # rownames(mat_scaled) <- aliases_ordered  # Change from colnames to rownames
  
  
  color_breaks <- c(-2, 0, 2)
  
  
  left_ha <- rowAnnotation(
    Response = response_data[responder_first],
    col = list(Response = c("non-responder_2" = '#fe9003', "responder_1" = '#115284')),
    width = unit(2.5, "mm"), 
    show_legend = TRUE
  )
  
  
  # ===== CREATE ROW ANNOTATION =====
  # Get all molecules from the heatmap (these are now the column names of mat_scaled)
  all_molecules <- colnames(mat_scaled)  # Changed from rownames to colnames
  
  # Create annotation vectors for each group
  group1 <- rep(NA, length(all_molecules))
  group2 <- rep(NA, length(all_molecules))
  group3 <- rep(NA, length(all_molecules))
  
  # Loop through each column and assign group labels
  for(j in 1:nrow(molecule_annotation)) {
    mol1 <- molecule_annotation[[1]][j]
    mol2 <- molecule_annotation[[2]][j]
    mol3 <- molecule_annotation[[3]][j]
    
    # Check and assign for column 1
    if(!is.na(mol1) && mol1 %in% all_molecules) {
      idx <- which(all_molecules == mol1)
      group1[idx] <- "Blood counts"
    }
    
    # Check and assign for column 2
    if(!is.na(mol2) && mol2 %in% all_molecules) {
      idx <- which(all_molecules == mol2)
      group2[idx] <- "Stromal targeting"
    }
    
    # Check and assign for column 3
    if(!is.na(mol3) && mol3 %in% all_molecules) {
      idx <- which(all_molecules == mol3)
      group3[idx] <- "Stromal produced"
    }
  }
  
  # Create colors for each group
  col_list <- list(
    Group1 = c("Blood counts" = "#66C2A5"),
    Group2 = c("Stromal targeting" = "#FC8D62"),
    Group3 = c("Stromal produced" = "#8DA0CB")
  )
  
  # Create a custom annotation function that ONLY draws dots
  anno_dot <- function(group_vector, dot_color) {
    anno_simple(
      x = ifelse(is.na(group_vector), NA, 1),
      col = c("1" = "white"),
      pch = ifelse(is.na(group_vector), NA, 19),
      pt_size = unit(3, "mm"),
      pt_gp = gpar(col = ifelse(is.na(group_vector), NA, dot_color)),
      simple_anno_size = unit(0.5, "mm")
    )
  }
  
  
  # Create TOP annotation using HeatmapAnnotation
  top_ha <- HeatmapAnnotation(
    Group1 = anno_dot(group1, "#66C2A5"),
    Group2 = anno_dot(group2, "#FC8D62"),
    Group3 = anno_dot(group3, "#8DA0CB"),
    annotation_name_side = "left",
    height = unit(1, "cm")
  )
  
  
  # Create the heatmap with annotation inside
  # Then use this in your Heatmap call:
  p <- Heatmap(
    as.matrix(mat_scaled), 
    name = "Gene Z-score",
    col = colorRamp2(color_breaks, c("#2166AC", "white", "#B2182B")),
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    column_names_rot = 90,
    left_annotation = left_ha,    # Changed variable name
    top_annotation = top_ha,      # Changed variable name
    show_column_names = TRUE,
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8),
    column_title = files_info$title[i],
    heatmap_width = unit(10, "cm"),
    column_names_max_height = unit(6, "cm"),
    heatmap_legend_param = list(
      at = c(-2, -1, 0, 1, 2),
      labels = c("-2", "-1", "0", "1", "2")
    )
  )
  
  # Create output filename
  output_filename <- paste0("heatmap_annotated_", gsub("\\.csv$", "", files_info$filename[i]), ".pdf")
  output_path <- file.path("../graphs/nulisa_heatmaps", output_filename)
  
  # Create directory if it doesn't exist
  output_dir <- "../graphs/nulisa_heatmaps"
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Calculate height
  min_height <- 4
  height_per_row <- 0.125
  calculated_height <- max(min_height, height_per_row * nrow(deg_data))
  
  # Save ComplexHeatmap
  pdf(output_path, width = 10, height = calculated_height)
  draw(p, annotation_legend_side = 'right', heatmap_legend_side = "right", 
       padding = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  dev.off()
  
  cat("Saved heatmap to:", output_path, "\n")
}  
  