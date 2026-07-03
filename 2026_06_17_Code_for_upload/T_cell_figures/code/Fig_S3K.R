library(limma)
library(reshape2)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)
library(biomaRt)
library(readxl)
library(edgeR)
# 
setwd('/Users/henryhampton/OneDrive - UNSW/2026_06_17_Code_for_upload/T_cell_figures/code')

# Read in the data which was downloaded
df = read.table('../data/GSE274999_oralAZA-STAR-counts.txt', header = TRUE)

count <- 0
for (i in colnames(df)) {
  if (grepl('JT', i)) {
    count <- count + 1
  }
}

# List of interferon gene symbols
ifn_symbols <- c("IFNA1", "IFNA2", "IFNA4", "IFNA5", "IFNA6", "IFNA7", 
                 "IFNA8", "IFNA10", "IFNA13", "IFNA14", "IFNA16", "IFNA17", 
                 "IFNA21", "IFNB1", "IFNE", "IFNK", "IFNW1", "IFNG", 
                 "IFNL1", "IFNL2", "IFNL3", "IFNL4", "TGFB1")

# Get the correct ENSG IDs
correct_ensg <- mapIds(org.Hs.eg.db, 
                       keys = ifn_symbols,
                       column = "ENSEMBL",
                       keytype = "SYMBOL",
                       multiVals = "first")

all_mappings <- select(org.Hs.eg.db, 
                       keys = ifn_symbols,
                       columns = c("ENSEMBL", "SYMBOL"),
                       keytype = "SYMBOL")

# Create corrected named vector
interferon_genes_correct <- setNames(correct_ensg, ifn_symbols)
print(interferon_genes_correct)


interferon_genes_correct <- c(
  # Type I Interferons - IFN-α family
  "IFNA1" = "ENSG00000197919",
  "IFNA2" = "ENSG00000188379",
  "IFNA4" = "ENSG00000236637",
  "IFNA5" = "ENSG00000147873",
  "IFNA6" = "ENSG00000120235",
  "IFNA7" = "ENSG00000214042",
  "IFNA8" = "ENSG00000120242",
  "IFNA10" = "ENSG00000186803",
  "IFNA13" = "ENSG00000233816",
  "IFNA14" = "ENSG00000228083",
  "IFNA16" = "ENSG00000147885",
  "IFNA17" = "ENSG00000234829",
  "IFNA21" = "ENSG00000137080",
  
  # Type I Interferons - other
  "IFNB1" = "ENSG00000171855",
  "IFNE" = "ENSG00000184995",
  "IFNK" = "ENSG00000147896",
  "IFNW1" = "ENSG00000177047",
  
  # Type II Interferon
  "IFNG" = "ENSG00000111537",
  
  # Type III Interferons - both isoforms
  "IFNL1_iso1" = "ENSG00000182393",
  "IFNL1_iso2" = "ENSG00000291872",
  "IFNL2_iso1" = "ENSG00000183709",
  "IFNL2_iso2" = "ENSG00000291875",
  "IFNL3_iso1" = "ENSG00000197110",
  "IFNL3_iso2" = "ENSG00000291876",
  "IFNL4_iso1" = "ENSG00000272395",
  "IFNL4_iso2" = "ENSG00000292250", 
  
  # Positive Control Cytokine
  "TGFB1" = "ENSG00000105329")

# Find all IFN genes
pattern <- paste0("^(", paste(interferon_genes_correct, collapse = "|"), ")")
ifn_idx <- grep(pattern, df$Geneid)

# For plotting, you might want to collapse isoforms back to the main gene name
ensg_clean <- sub("\\.\\d+$", "", df$Geneid[ifn_idx])
gene_symbols_corrected <- names(interferon_genes_correct)[match(ensg_clean, interferon_genes_correct)]
gene_symbols_plotting <- sub("_iso\\d+$", "", gene_symbols_corrected)


# Convert to long format in one go
ifn_long <- df[ifn_idx, -1] %>%
  setNames(c(colnames(.)[-1])) %>%  # if Geneid is first column
  `rownames<-`(gene_symbols_corrected) %>%
  as.matrix() %>%
  melt() %>%
  filter(grepl("JT", Var2)) %>%
  mutate(value = as.numeric(value))

# Then apply log transformation to allow for easy plotting
ifn_long$value_log <- log10(ifn_long$value + 1)

# Remove and outlier sample which expresses high levels of all IFNs
ifn_long = ifn_long[ifn_long$Var2 != 'JT184', ]

ifn_long$Var1 <- gsub('_iso1', '', ifn_long$Var1)

# Specify exact order you want
gene_order <- c("IFNA1", "IFNA2", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", 
                "IFNA10", "IFNA13", "IFNA14", "IFNA16", "IFNA17", "IFNA21",
                "IFNB1",  "IFNE", "IFNK", "IFNW1", "IFNG",
                "IFNL1", "IFNL2", "IFNL3", "IFNL4", 'TGFB1')

ifn_long$Var1 <- factor(ifn_long$Var1, levels = gene_order)

# Make a graph which shows IFN expression across all CD34 cells
ggplot(ifn_long, aes(x = Var1, y = value_log, fill = Var1)) +
  geom_violin(scale = "width", alpha = 0.7, linewidth = 0.2) +
  geom_jitter(width = 0.2, size = 0.75, shape = 21, fill = "white", color = "black", 
              stroke = 0.5, alpha = 1.0) +
  theme_classic() +  # Already has no gridlines and axis lines
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(title = "Expression of IFNs",
       x = "Interferon Gene", 
       y = "Expression Level (log10 + 1)")

ggsave("../results/ifn_expression_across_all_samples_violin.pdf", width = 4, height = 2)
