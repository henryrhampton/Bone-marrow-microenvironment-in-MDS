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
setwd('~/OneDrive - UNSW/2026_06_17_Code_for_upload/NULISA_Figures/code/')
npq_mds=read_excel("../data/NULISA_for_upload.xlsx", sheet = "NPQ") %>% as.data.frame()
npq_mds=npq_mds[,1:82]
clinical=read.csv("../data/unsw_australia_plasma_protein_concentration_patient_updated.csv") %>% as.data.frame()
clinical=clinical[1:46,];colnames(clinical)[1]="sample"
clinical[,1]=paste0(clinical$PID,"_",clinical$Cycle_day)

##There are bone marrow and pb samples, need to analyze separately
index=unlist(lapply(colnames(npq_mds),function(x){endsWith(x,"_PB")}))
mds_pb=npq_mds[,index];mds_bm=npq_mds[,!index]
rownames(mds_pb)=npq_mds$targetName
rownames(mds_bm)=mds_bm$targetName;mds_bm=mds_bm[,2:47]
mds_pb=as.data.frame(t(mds_pb));mds_bm=as.data.frame(t(mds_bm))
###add patient ID and cycle to the column
mds_pb$PID=unlist(lapply(rownames(mds_pb),function(x){substr(x,1,3)}))
mds_bm$PID=unlist(lapply(rownames(mds_bm),function(x){substr(x,1,3)}))
mds_pb$cycle=unlist(lapply(rownames(mds_pb),function(x){substr(x,5,9)}))
mds_bm$cycle=unlist(lapply(rownames(mds_bm),function(x){ifelse(nchar(x)>14,substr(x,5,nchar(x)-4),substr(x,5,9))}))
colnames(mds_pb)[c(99,102,122,174)]=c("HLA_DRA","IFNA1_IFNA13","IL17A_IL17F","LTA_LTB")
colnames(mds_bm)[c(99,102,122,174)]=c("HLA_DRA","IFNA1_IFNA13","IL17A_IL17F","LTA_LTB")

my_dict <- setNames(clinical$outcome_6, clinical$PID)
mds_pb['outcome_6'] = my_dict[mds_pb$PID]
mds_bm['outcome_6'] = my_dict[mds_bm$PID]
mds_pb=mds_pb[mds_pb$cycle != "C12_D",]


######limma to run differential expression on blood samples
tp = c("C1_D1", "C7_D1")

for(cycle in tp) {
  cat("\n=== Processing cycle:", cycle, "===\n")
  
  # Filter data for current cycle
  temp <- mds_pb[mds_pb$cycle == cycle, ]
  
  # Clean outcome labels
  temp$outcome_6[temp$outcome_6 == "non-responder_2"] <- "Nonresponder"
  temp$outcome_6[temp$outcome_6 == "responder_1"] <- "Responder"
  temp$outcome_6 <- factor(temp$outcome_6, levels = c("Responder", "Nonresponder"))
  
  # Prepare expression and metadata
  exp <- temp[, 1:247]
  meta <- data.frame(contrast = factor(temp$outcome_6))
  exp <- as.data.frame(t(exp))
  colnames(exp) <- rownames(meta)
  
  # Differential expression analysis
  design <- model.matrix(~ 0 + contrast, data = meta)
  fit <- lmFit(exp, design)
  contrast <- makeContrasts(NR_CR = contrastNonresponder - contrastResponder, levels = design)
  fits <- contrasts.fit(fit, contrast)
  ebFit <- eBayes(fits)
  
  # Get results
  limma.res <- topTable(ebFit, coef = "NR_CR", adjust.method = 'fdr', number = Inf)
  limma.res <- limma.res %>% 
    filter(!is.na(adj.P.Val)) %>% 
    mutate(logP = -log10(adj.P.Val)) %>%
    mutate(tag = "NR -vs- CR") %>%
    mutate(Gene = rownames(limma.res))
  
  # Create output filename based on cycle
  output_filename <- paste0('blood_', tolower(gsub("_", "", cycle)), '_r_v_nr.csv')
  output_path <- file.path('../data/differential_expression_tables/', 
                           output_filename)
  
  # Save results
  write.csv(limma.res, file = output_path, row.names = FALSE)
  
  cat("Saved results to:", output_path, "\n")
  cat("Number of genes:", nrow(limma.res), "\n")
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clinical$Protein_ave = mean(clinical$Protein_concentration..mg.ml.)
clinical$rel_protein = (clinical$Protein_concentration..mg.ml.)/(clinical$Protein_ave)
clinical$NPQ = log2(clinical$rel_protein)

mds_bm2 = mds_bm
clinical=clinical[order(match(clinical$sample,mds_bm2$sample)),]
for (i in 1:dim(clinical)[1]){
  mds_bm2[i,1:247]=mds_bm2[i,1:247] - clinical$NPQ[i]
}



tp = c("C1_D1", "C7_D1")

for(cycle in tp) {
  cat("\n=== Processing cycle:", cycle, "===\n")
  
  # Use the cycle variable instead of hardcoded "C1_D1"
  temp <- mds_bm2[mds_bm2$cycle == cycle, ]
  
  temp$outcome_6[temp$outcome_6 == "non-responder_2"] <- "Nonresponder"
  temp$outcome_6[temp$outcome_6 == "responder_1"] <- "Responder"
  temp$outcome_6 <- factor(temp$outcome_6, levels = c("Responder", "Nonresponder"))
  
  exp <- temp[, 1:247]
  meta <- data.frame(contrast = factor(temp$outcome_6))
  exp <- as.data.frame(t(exp))
  colnames(exp) <- rownames(meta)
  
  design <- model.matrix(~ 0 + contrast, data = meta)
  fit <- lmFit(exp, design)
  contrast <- makeContrasts(NR_CR = contrastNonresponder - contrastResponder, levels = design)
  fits <- contrasts.fit(fit, contrast)
  ebFit <- eBayes(fits)
  
  ## result
  limma.res <- topTable(ebFit, coef = "NR_CR", adjust.method = 'fdr', number = Inf)
  limma.res <- limma.res %>% 
    filter(!is.na(adj.P.Val)) %>% 
    mutate(logP = -log10(adj.P.Val)) %>%
    mutate(tag = "NR -vs- CR") %>%
    mutate(Gene = rownames(limma.res))
  
  # Create dynamic filename based on cycle
  cycle_clean <- tolower(gsub("_", "", cycle))  # "C1_D1" -> "c1d1", "C7_D1" -> "c7d1"
  output_filename <- paste0('bm_protein_protein_scaled_', cycle_clean, '_r_v_nr.csv')
  output_path <- file.path('../data/differential_expression_tables', 
                           output_filename)
  
  write.csv(limma.res, file = output_path, row.names = FALSE)
  cat("Saved results to:", output_path, "\n")
  cat("Number of genes processed:", nrow(limma.res), "\n")
}
