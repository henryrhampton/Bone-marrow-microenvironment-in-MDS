library(Seurat)
library(SeuratDisk)
library(anndata)
library(reticulate)
library(zellkonverter)
library(hdf5r)
library(dplyr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in the single cell object which Anne sent me
setwd('/Users/henryhampton/OneDrive - UNSW/2026_06_17_Code_for_upload/T_cell_figures/code/')
sce <- readRDS("../data/clonaltype_vdjdb_annotation_patient_updated.RDS")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in a csv file detailing human antigens
vdj <- read.csv("../data/vdjdb_paired_tcr.csv", stringsAsFactors = FALSE, row.names = 1)

# Filter the VDJdb object for human antigens
vdj <- vdj[vdj$antigen.species == 'HomoSapiens', ]

tcr_cols <- 'CTaa'

# Split the combined TCR into alpha and beta chains
sce@meta.data[[tcr_cols]] <- as.character(sce@meta.data[[tcr_cols]])
tcr_combined <- as.character(sce@meta.data[[tcr_cols]])
tcr_split <- strsplit(tcr_combined, "_")

# Extract the first TCR sequence from the alpha chain (a small number of cells have dual TCRs)
sce@meta.data$TRA_CDR3 <- sapply(tcr_split, function(x) {
  if (length(x) >= 2 && !is.na(x[1])) x[1] else NA
})

# Extract the first TCR sequence from the beta chain (a small number of cells have dual TCRs)
sce@meta.data$TRB_CDR3 <- sapply(tcr_split, function(x) {
  if (length(x) >= 2 && !is.na(x[2])) x[2] else NA
})

# Update CTaa to be just the BETA chain
sce@meta.data$CTaa <- sce@meta.data$TRB_CDR3

# Identify rows to annotate (HomoSapiens only)
homo_sapiens_rows <- sce@meta.data$antigen == "HomoSapiens" & !is.na(sce@meta.data$antigen)
cat("Cells with antigen == 'HomoSapiens':", sum(homo_sapiens_rows), "\n")

# === BETA CHAIN MATCHING ===
human_self_beta <- vdj[!is.na(vdj$TRB_junction_aa), ]
beta_antigen_lookup <- setNames(human_self_beta$antigen.epitope, human_self_beta$TRB_junction_aa)
beta_gene_lookup <- setNames(human_self_beta$antigen.gene, human_self_beta$TRB_junction_aa)

# === ALPHA CHAIN MATCHING ===
human_self_alpha <- vdj[!is.na(vdj$TRA_junction_aa), ]
alpha_antigen_lookup <- setNames(human_self_alpha$antigen.epitope, human_self_alpha$TRA_junction_aa)
alpha_gene_lookup <- setNames(human_self_alpha$antigen.gene, human_self_alpha$TRA_junction_aa)

# Initialize columns
sce@meta.data$human_self_antigen <- NA
sce@meta.data$human_self_gene <- NA
sce@meta.data$match_strategy <- NA

# Match TCRs using BETA chain first
sce@meta.data$human_self_antigen[homo_sapiens_rows] <- 
  beta_antigen_lookup[sce@meta.data$CTaa[homo_sapiens_rows]]
sce@meta.data$human_self_gene[homo_sapiens_rows] <- 
  beta_gene_lookup[sce@meta.data$CTaa[homo_sapiens_rows]]
sce@meta.data$match_strategy[homo_sapiens_rows & !is.na(sce@meta.data$human_self_antigen)] <- "beta"

# For cells that didn't match with beta, try ALPHA chain
unmatched_homo <- homo_sapiens_rows & is.na(sce@meta.data$human_self_antigen)

sce@meta.data$human_self_antigen[unmatched_homo] <- 
  alpha_antigen_lookup[sce@meta.data$TRA_CDR3[unmatched_homo]]
sce@meta.data$human_self_gene[unmatched_homo] <- 
  alpha_gene_lookup[sce@meta.data$TRA_CDR3[unmatched_homo]]
sce@meta.data$match_strategy[unmatched_homo & !is.na(sce@meta.data$human_self_antigen)] <- "alpha"

# === HANDLE DUAL TCR CELLS ===
# Identify cells with multiple chains (semicolon-separated)
has_dual_tcr <- homo_sapiens_rows & 
  (grepl(";", sce@meta.data$TRA_CDR3) | grepl(";", sce@meta.data$TRB_CDR3))

if (sum(has_dual_tcr, na.rm = TRUE) > 0) {
  cat("\nFound", sum(has_dual_tcr), "cells with dual TCRs\n")
  
  # For each dual TCR cell, try to match individual chains
  dual_tcr_cells <- which(has_dual_tcr)
  
  for (idx in dual_tcr_cells) {
    # Skip if already matched
    if (!is.na(sce@meta.data$human_self_antigen[idx])) next
    
    # Split the chains by semicolon
    tra_chains <- strsplit(sce@meta.data$TRA_CDR3[idx], ";")[[1]]
    trb_chains <- strsplit(sce@meta.data$TRB_CDR3[idx], ";")[[1]]
    
    # Try matching each TRB chain
    for (trb in trb_chains) {
      if (trb %in% names(beta_antigen_lookup)) {
        sce@meta.data$human_self_antigen[idx] <- beta_antigen_lookup[trb]
        sce@meta.data$human_self_gene[idx] <- beta_gene_lookup[trb]
        sce@meta.data$match_strategy[idx] <- "dual_tcr_beta"
        break
      }
    }
    
    # If still not matched, try TRA chains
    if (is.na(sce@meta.data$human_self_antigen[idx])) {
      for (tra in tra_chains) {
        if (tra %in% names(alpha_antigen_lookup)) {
          sce@meta.data$human_self_antigen[idx] <- alpha_antigen_lookup[tra]
          sce@meta.data$human_self_gene[idx] <- alpha_gene_lookup[tra]
          sce@meta.data$match_strategy[idx] <- "dual_tcr_alpha"
          break
        }
      }
    }
  }
}

# Flag dual TCR cells
sce@meta.data$is_dual_tcr <- grepl(";", sce@meta.data$TRA_CDR3) | 
  grepl(";", sce@meta.data$TRB_CDR3)

# === RESULTS ===
cat("\n=== MATCHING RESULTS ===\n")
cat("Total HomoSapiens cells:", sum(homo_sapiens_rows), "\n")
cat("Matched to VDJdb:", sum(!is.na(sce@meta.data$human_self_antigen[homo_sapiens_rows])), "\n\n")

cat("Matching strategy breakdown:\n")
print(table(sce@meta.data$match_strategy[homo_sapiens_rows], useNA = "ifany"))

cat("\nSelf-antigens recognized:\n")
print(sort(table(sce@meta.data$human_self_gene[homo_sapiens_rows]), decreasing = TRUE))
print(table(sce@meta.data$human_self_gene[homo_sapiens_rows], sce@meta.data$alias[homo_sapiens_rows]))

# Add validation flag
sce@meta.data$vdjdb_validated <- homo_sapiens_rows & 
  !is.na(sce@meta.data$human_self_antigen)

cat("\nValidation rate:", 
    round(100 * sum(!is.na(sce@meta.data$human_self_antigen[homo_sapiens_rows])) / 
            sum(homo_sapiens_rows), 1), "%\n")

# Save results
validated_cells <- sce@meta.data[sce@meta.data$vdjdb_validated == TRUE, 
                                 c("TRA_CDR3", "TRB_CDR3", 
                                   "human_self_antigen", "human_self_gene", 
                                   "match_strategy", "is_dual_tcr", "alias", 'timepoint')]

write.csv(validated_cells, file = '../results/human_reactive_tcrs_from_vdjdb.csv')

# Create a summary table with counts and percentages
summary_table <- sce@meta.data %>%
  group_by(alias) %>%
  summarise(
    total_cells = n(),
    non_na_count = sum(!is.na(human_self_gene)),
    percent_non_na = round((sum(!is.na(human_self_gene)) / n()) * 100, 2)
  ) %>%
  arrange(desc(percent_non_na))

write.csv(summary_table, file = '../results/human_reactive_tcrs_summary_table_vdjdb.csv')

# Print the table
print(summary_table)


n_antigens <- length(unique(na.omit(sce@meta.data$human_self_gene)))
#seurat_colors <- scales::hue_pal()(n_antigens)
unique_genes <- unique(na.omit(sce@meta.data$human_self_gene))

# Generate colors for your genes
color_palette <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique_genes))
names(color_palette) <- unique_genes

# Add grey for NA
all_colors <- c(color_palette, "NA" = "grey90")

# Create display column
sce@meta.data$human_gene_display <- ifelse(
  is.na(sce@meta.data$human_self_gene),
  "NA",
  as.character(sce@meta.data$human_self_gene)  
)

# Plot with VDJdb-matched cells on top
DimPlot(sce, 
        reduction = "XUMAP_",
        group.by = "human_gene_display",
        cols = all_colors,
        pt.size = 0.1, 
        order = unique_genes) +  # Plot colored cells on top
  labs(title = "Human Self-Reactive T Cells") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.line = element_line(color = "black"))  # Keep axis lines

ggsave("../results/umap_self_reactive.pdf", width = 4.923, height = 4) # I spent ages trying to obtain a square pl

######
library(readr)

receptor_data <- read.csv("../data/receptor_full_v3/tcr_full_v3.csv", 
                          stringsAsFactors = FALSE)

# Set the first row as column names
col_names <- as.character(receptor_data[1, ])
receptor_data <- receptor_data[-1, ]  # Remove the header row
colnames(receptor_data) <- col_names

# Get original names
orig_names <- colnames(receptor_data)

# Fix duplicate names by adding Chain 1 and Chain 2 prefixes
new_names <- orig_names
new_names[13:41] <- paste("Chain1", orig_names[13:41], sep = "_")
new_names[42:70] <- paste("Chain2", orig_names[42:70], sep = "_")

# Apply new names
colnames(receptor_data) <- new_names
# Check the structure
dim(receptor_data)
colnames(receptor_data)

# Step 1: Extract TCR data from Seurat object
tcr_data <- sce@meta.data[, c("BARCODE", "TRA_CDR3", "TRB_CDR3", "is_dual_tcr")]
# Add the full cell barcode (rownames) as a new column
tcr_data$full_Barcode <- rownames(tcr_data)

# Filter for cells that have at least one TCR chain
tcr_data <- tcr_data[!is.na(tcr_data$TRA_CDR3) | !is.na(tcr_data$TRB_CDR3), ]


# Step 2: Prepare IEDB receptor data (human only)
human_tcr <- receptor_data[receptor_data$`Chain1_Organism IRI` == "http://purl.obolibrary.org/obo/NCBITaxon_9606", ]

human_tcr_tra <- human_tcr
human_tcr_tra$IEDB_TRA_CDR3 <- human_tcr_tra$`Chain1_CDR3 Curated`  # Save it first

matches_tra <- merge(
  tcr_data,
  human_tcr_tra,
  by.x = "TRA_CDR3",
  by.y = "Chain1_CDR3 Curated",
  all.x = FALSE
)

cat("TRA CDR3 matches:", nrow(matches_tra), "\n")

# Step 4: Match TRB (beta chain)
# For TRB matches
human_tcr_trb <- human_tcr
human_tcr_trb$IEDB_TRB_CDR3 <- human_tcr_trb$`Chain2_CDR3 Curated`  # Save it first

matches_trb <- merge(
  tcr_data,
  human_tcr_trb,
  by.x = "TRB_CDR3",
  by.y = "Chain2_CDR3 Curated",
  all.x = FALSE
)


cat("TRB CDR3 matches:", nrow(matches_trb), "\n\n")

# Now combine
matches_tra$Match_Type <- "TRA"
matches_trb$Match_Type <- "TRB"

# Step 5: Display results
if(nrow(matches_tra) > 0) {
  cat("=== ALPHA CHAIN MATCHES ===\n")
  print(matches_tra[, c("BARCODE", "TRA_CDR3", "Name", "Source Molecule", "MHC Allele Names")])
}

if(nrow(matches_trb) > 0) {
  cat("\n=== BETA CHAIN MATCHES ===\n")
  print(matches_trb[, c("BARCODE", "TRB_CDR3", "Name", "Source Molecule", "MHC Allele Names")])
}

# Step 6: Combine and save results
# Define the columns we want to keep
key_columns <- c("BARCODE", "TRA_CDR3", "TRB_CDR3", "is_dual_tcr",
                 "Reference Name", "Name", "Source Molecule", "Source Organism",
                 "MHC Allele Names", "Chain1_Type", "Chain2_Type",
                 "Chain1_CDR3 Curated", "Chain2_CDR3 Curated")

# Add match type to each
if(nrow(matches_tra) > 0) {
  matches_tra$Match_Type <- "TRA"
}

if(nrow(matches_trb) > 0) {
  matches_trb$Match_Type <- "TRB"
}

# Combine with bind_rows (handles column mismatches)
all_matches <- bind_rows(matches_tra, matches_trb)


# Summary by source molecule
source_table <- sort(table(all_matches$`Source Molecule`), decreasing = TRUE)

# Show first few key columns
cat("\nFirst 10 matches:\n")
print(all_matches[1:min(10, nrow(all_matches)), 
                  c("BARCODE", "Match_Type", "TRA_CDR3", "TRB_CDR3", 
                    "Name", "Source Molecule", "Source Organism")])

all_matches = all_matches[all_matches$`Source Organism` == 'Homo sapiens (human)', ]

# Some TCRs have multiple characterised antigens. Therefore, pick one of the antigens
# top_epitope_per_cell <- all_matches %>%
#   group_by(full_Barcode) %>%
#   slice(1) %>%
#   select(full_Barcode, 'Source Molecule')

top_epitope_per_cell <- all_matches %>%
  group_by(full_Barcode) %>%
  dplyr::slice(1) %>%
  dplyr::select(full_Barcode, `Source Molecule`)


# Add to Seurat
sce$IEDB_Top_Source <- NA
sce@meta.data[top_epitope_per_cell$full_Barcode, "IEDB_Top_Source"] <- top_epitope_per_cell$`Source Molecule`

# Read in csv file detailing the name
gene_name = read_csv('../data/protein_to_gene_mapping.csv')
protein_to_gene_vector <- setNames(gene_name$Gene_Symbol, 
                                   gene_name$Protein_Name)

# Replace protein names with gene symbols in your Seurat object
sce@meta.data$IEDB_Top_Gene <- protein_to_gene_vector[sce@meta.data$IEDB_Top_Source]

# Check the result
table(sce@meta.data$IEDB_Top_Gene, useNA = "ifany")
table(sce@meta.data$IEDB_Top_Gene, sce@meta.data$alias, useNA = "ifany")

x = table(sce@meta.data$IEDB_Top_Gene, sce@meta.data$alias, useNA = "ifany")
write.csv(x, file = '../results/human_reactive_tcrs_from_iedb.csv')

rm(x)
unique_genes <- unique(na.omit(sce@meta.data$IEDB_Top_Gene))

# Generate colors for your genes
color_palette <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique_genes))
names(color_palette) <- unique_genes

# Add grey for NA
all_colors <- c(color_palette, "NA" = "grey90")

# Create display column
sce@meta.data$IEDB_Gene_Display <- ifelse(
  is.na(sce@meta.data$IEDB_Top_Gene),
  "NA",
  as.character(sce@meta.data$IEDB_Top_Gene)
)

# Plot with IEDB-matched cells on top
DimPlot(sce, 
        reduction = "XUMAP_",
        group.by = "IEDB_Gene_Display",
        cols = all_colors,
        pt.size = 0.1, 
        order = unique_genes) +  # Plot colored cells on top
  labs(title = "Human Self-Reactive T Cells") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.line = element_line(color = "black"))


ggsave("../results/iedbumap_self_reactive.pdf", width = 7.536, height = 4) # I spent ages trying to obtain a square pl

# Create a summary table with counts and percentages
summary_table <- sce@meta.data %>%
  group_by(patient) %>%
  summarise(
    total_cells = n(),
    non_na_count = sum(!is.na(IEDB_Top_Gene)),
    percent_non_na = round((sum(!is.na(IEDB_Top_Gene)) / n()) * 100, 2)
  ) %>%
  arrange(desc(percent_non_na))

# Print the table
print(summary_table)

write.csv(summary_table, file = '../results/human_reactive_tcrs_summary_table_from_iedb.csv')

# Read in human antigen colors
antigen_color = readxl::read_excel('../data/t_cell_human_antigen_color.xlsx', sheet = 'Sheet1')
antigen_color_dict <- setNames(antigen_color[['Color']], antigen_color[['Antigen']])


DimPlot(sce, 
           reduction = "XUMAP_",
           group.by = "antigen",
           pt.size = 0.1, 
           order = unique_genes) +
  scale_colour_manual(
    values = antigen_color_dict,
    na.value = "grey90"
  ) +
  labs(title = "Human Self-Reactive T Cells") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.line = element_line(color = "black"))
