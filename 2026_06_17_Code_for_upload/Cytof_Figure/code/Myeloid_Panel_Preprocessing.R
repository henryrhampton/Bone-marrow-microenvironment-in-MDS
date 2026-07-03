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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PrimaryDirectory = setwd('/Users/henryhampton/Documents/2021_Cytof_Data/')
## Set input directory
InputDirectory <- getwd()
### Set metadata directory
setwd("metadata")
getwd()
MetadataDirectory <- getwd()
setwd(PrimaryDirectory)

### Create output directory

dir.create("output-align")
setwd("output-align")
getwd()
OutputDirectory <- getwd()
setwd(PrimaryDirectory)

setwd('/Users/henryhampton/Documents/2021_Cytof_Data/Myeloid_Panel/')
dat.list <- read.files(file.type = ".csv", do.embed.file.names = TRUE)

#Change column names so that I can merge the two sets of experiments (i.e some markers were give a single or a double underscore)
for (df in seq_along(dat.list)) {
  colnames(dat.list[[df]]) = gsub(x = colnames(dat.list[[df]]), pattern = '147Sm_Cy5', replacement = '147Sm_NKG2DFc_Cy5')
}
for (df in seq_along(dat.list)) {
  colnames(dat.list[[df]]) = gsub(x = colnames(dat.list[[df]]), pattern = '164Dy_CD172ab__SIRPab_', replacement = '164Dy_CD172ab_SIRPab')
}
for (df in seq_along(dat.list)) {
  colnames(dat.list[[df]]) = gsub(x = colnames(dat.list[[df]]), pattern = '161Dy_CD274__PDL1_', replacement = '161Dy_CD274_PDL1')
}
for (df in seq_along(dat.list)) {
  colnames(dat.list[[df]]) = gsub(x = colnames(dat.list[[df]]), pattern = '161Dy_CD274__PDL1', replacement = '161Dy_CD274_PDL1')
}
for (df in seq_along(dat.list)) {
  colnames(dat.list[[df]]) = gsub(x = colnames(dat.list[[df]]), pattern = '175Lu_PE', replacement = '175Lu_STAT1_PE')
}
for (df in seq_along(dat.list)) {
  colnames(dat.list[[df]]) = gsub(x = colnames(dat.list[[df]]), pattern = '156Gd_CD279__PD1_', replacement = '156Gd_CD279_PD1')
}
for (df in seq_along(dat.list)) {
  colnames(dat.list[[df]]) = gsub(x = colnames(dat.list[[df]]), pattern = '156Gd_CD279__PD1', replacement = '156Gd_CD279_PD1')
}
for (df in seq_along(dat.list)) {
  colnames(dat.list[[df]]) = gsub(x = colnames(dat.list[[df]]), pattern = '173Yb_CD47_biotin', replacement = '173Yb_CD47_Biotin')
}
for (df in seq_along(dat.list)) {
  colnames(dat.list[[df]]) = gsub(x = colnames(dat.list[[df]]), pattern = '173Yb_biotin', replacement = '173Yb_CD47_Biotin')
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
# Delete superfluous columns
for (df in dat.list){
  df = df[, c('190BCKG'):=NULL]
}
for (df in dat.list){
  df = df[, c('140Ce'):=NULL]
}
for (df in dat.list){
  df = df[, c('139La'):=NULL]
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

# Merge metadata with data
cell.dat <- do.merge.files(dat.list)


#~~~~~~~~~~~~~~~~~~~~~~~~~~
#Choose columns to transform 
as.matrix(names(cell.dat))
transform.cols.nums <- c(5:12, 14:15, 17:19, 21, 24:31, 40:42, 45, 47:51)
transform.cols.nums <- names(cell.dat)[transform.cols.nums] 

# Transform data
cell.dat = do.asinh(cell.dat, use.cols = transform.cols.nums, cofactor = 5, append.cf = FALSE, reduce.noise = FALSE)
cell.dat <- cell.dat [, !transform.cols.nums, with = FALSE]

meta.dat <- fread("myeloid_panel_patient_information.txt")
cell.dat <- do.add.cols(cell.dat, "FileName", meta.dat, "File_Name", rmv.ext = TRUE)

ref.dat <- cell.dat[cell.dat[["Reference"]] == TRUE,]

batch.col <- "Batch"
method <- '95p'
crs.append <- '_coarseAlign'
align.cols = names(cell.dat)[c(22:53)]

cell.dat <- run.align(ref.dat = ref.dat,
                            target.dat = cell.dat,
                            batch.col = batch.col,
                            align.cols = align.cols,
                            method = method,
                            append.name = crs.append)


ref.dat <- cell.dat[cell.dat[["Reference"]] == TRUE,]

cellular.cols = names(cell.dat)[c(64:95)]
cluster.cols = names(cell.dat)[c(64:95)]
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
align.cols = names(cell.dat)[c(64:95)]
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


cell.dat = cell.dat[ , (64:95) := NULL ]

setwd('/home/henry/Cytof_Analysis/Myeloid_Panel/aligned_data')
fwrite(cell.dat, file = 'aligned_myeloid_panel.csv')


