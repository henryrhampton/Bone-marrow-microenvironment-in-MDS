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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd('/home/henry/Cytof_Analysis/T_Cell_Panel/')
dat.list <- read.files(file.type = ".csv", do.embed.file.names = TRUE)

#Change column names so that I can merge the two sets of experiments (i.e some markers were give a single or a double underscore)
for (df in seq_along(dat.list)) {
  colnames(dat.list[[df]]) = gsub(x = colnames(dat.list[[df]]), pattern = '156Gd_CD279__PD1', replacement = '156Gd_CD279_PD1')
}
for (df in seq_along(dat.list)) {
  colnames(dat.list[[df]]) = gsub(x = colnames(dat.list[[df]]), pattern = '149Sm_cd366_tim3', replacement = '149Sm_Tim3')
}
for (df in seq_along(dat.list)) {
  colnames(dat.list[[df]]) = gsub(x = colnames(dat.list[[df]]), pattern = '149Sm_CD366', replacement = '149Sm_Tim3')
}
for (df in seq_along(dat.list)) {
  colnames(dat.list[[df]]) = gsub(x = colnames(dat.list[[df]]), pattern = '161Dy_NKp46', replacement = '161Dy_CD335_Nkp46')
}

for (df in seq_along(dat.list)) {
  colnames(dat.list[[df]]) = gsub(x = colnames(dat.list[[df]]), pattern = '^113In$', replacement = '113In_CD45')
}
for (df in seq_along(dat.list)) {
  colnames(dat.list[[df]]) = gsub(x = colnames(dat.list[[df]]), pattern = '^104Pd$', replacement = '104Pd_CD45')
}
for (df in seq_along(dat.list)) {
  colnames(dat.list[[df]]) = gsub(x = colnames(dat.list[[df]]), pattern = '^106Pd$', replacement = '106Pd_CD45')
}
for (df in seq_along(dat.list)) {
  colnames(dat.list[[df]]) = gsub(x = colnames(dat.list[[df]]), pattern = '^108Pd$', replacement = '108Pd_CD45')
}
for (df in seq_along(dat.list)) {
  colnames(dat.list[[df]]) = gsub(x = colnames(dat.list[[df]]), pattern = '^110Pd$', replacement = '110Pd_CD45')
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delete superfluous columns to save memory
for (df in dat.list){
  df = df[, c('190BCKG'):=NULL]
}
for (df in dat.list){
  df = df[, c('140Ce'):=NULL]
}
for (df in dat.list){
  df = df[, c('Offset'):=NULL]
}
for (df in dat.list){
  df = df[, c('Residual'):=NULL]
}
for (df in dat.list){
  df = df[, c('190Os'):=NULL]
}
for (df in dat.list){
  df = df[, c('SampleID'):=NULL]
}
for (df in dat.list){
  df = df[, c('192Pt'):=NULL]
}
for (df in dat.list){
  df = df[, c('141Pr_CD196_CCR6'):=NULL]
}
for (df in dat.list){
  df = df[, c('155Gd_CD196_CCR6'):=NULL]
}
for (df in dat.list){
  df = df[, c('155Gd'):=NULL]
}
for (df in dat.list){
  df = df[, c('141Pr'):=NULL]
}
# Merge metadata with data
cell.dat <- do.merge.files(dat.list)

cell.dat <- cell.dat[cell.dat[["FileName"]] != 'PBMC_Control_2021_02_16_Myeloid_Panel',]
cell.dat <- cell.dat[cell.dat[["FileName"]] != 'PBMC_Control_2021_02_24_Myeloid_Panel',]
cell.dat <- cell.dat[cell.dat[["FileName"]] != 'PBMC_Control_2021_03_02_Myeloid_Panel',]
head(cell.dat)
dim(cell.dat)
any(is.na(cell.dat))



#Choose columns to transform
as.matrix(names(cell.dat))
transform.cols.nums <- c(3:5, 8:12, 14, 16:18, 20, 23:31, 38:40, 43, 45:49)
transform.cols.nums <- names(cell.dat)[transform.cols.nums] 


# Transform data
cell.dat = do.asinh(cell.dat, use.cols = transform.cols.nums, cofactor = 5, append.cf = FALSE, reduce.noise = FALSE)
cell.dat <- cell.dat [, !transform.cols.nums, with = FALSE]

meta.dat <- fread("t_cell_panel_patient_information.txt")
cell.dat <- do.add.cols(cell.dat, "FileName", meta.dat, "File_Name", rmv.ext = TRUE)

ref.dat <- cell.dat[cell.dat[["Reference"]] == TRUE,]

batch.col <- "Batch"
method <- '95p'
crs.append <- '_coarseAlign'
align.cols = names(cell.dat)[c(22:52)]

cell.dat <- run.align(ref.dat = ref.dat,
                            target.dat = cell.dat,
                            batch.col = batch.col,
                            align.cols = align.cols,
                            method = method,
                            append.name = crs.append)


ref.dat <- cell.dat[cell.dat[["Reference"]] == TRUE,]

cellular.cols = names(cell.dat)[c(64:94)]
cluster.cols = names(cell.dat)[c(64:94)]
batch.col <- "Batch"
align.model <- prep.cytonorm(dat = ref.dat,
                                 cellular.cols = cellular.cols,
                                 cluster.cols = cluster.cols,
                                 batch.col = batch.col, 
                                 xdim = 14,
                                 ydim = 14,
                                 meta.k = 8)

cytonorm.goal <- 'mean'
cytonorm.nQ <- 101
align.cols = names(cell.dat)[c(64:94)]
align.model <- train.cytonorm(model = align.model,
                                  align.cols = align.cols,
                                  cytonorm.goal = cytonorm.goal,
                                  cytonorm.nQ = cytonorm.nQ)
fine.append <- '_fineAlign'
cell.dat <- run.cytonorm(dat = cell.dat,
                             model = align.model,
                             batch.col = batch.col,
                             append.name = fine.append,
                             )


cell.dat = cell.dat[ , (66:92) := NULL ]

setwd('/home/henry/Cytof_Analysis/T_Cell_Panel/aligned_data')
fwrite(cell.dat, file = 'aligned_t_cell_panel_quantile.csv')
