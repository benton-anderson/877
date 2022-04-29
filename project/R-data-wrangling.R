# LIBD vignette https://bioconductor.org/packages/release/data/experiment/vignettes/spatialLIBD/inst/doc/spatialLIBD.html

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Dino")

# devtools::install_github('satijalab/seurat-data')
set.seed(1)

library(spatialLIBD)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(sctransform)
library(Dino)
library(Matrix)
library(here)

options(max.print=100)

# Get the sce data from spatialLIBD
ehub <- ExperimentHub::ExperimentHub()
sce <- fetch_data(type = "sce", eh = ehub)
# spe <- sce_to_spe(sce)  # Convert to a SpatialExperiment object
scecounts <- assays(sce)$counts
# dim(mycounts)
# colnames(mycounts)
# colData(sce)$sample_name %in% c(151507, 151669, 151673)
# x <- grep("151507|151669|151673", colData(sce)$sample_name)
# length(colData(sce)$sample_name)
# table(colData(sce)$sample_name %in% c(151507, 151669, 151673))
# assays(sce)$counts[, colData(sce)$sample_name]
# as.matrix(assays(sce)$counts[, colData(sce)$sample_name])

# write.csv(rowData(sce), 'rowdata.csv')
# write.csv(colData(sce), 'coldata.csv')

### Get the column indices for the 6 columns
sample_names <- c(151507, 151508, 151669, 151670, 151673, 151674)
coldata <- colData(sce)
# scecounts <- scecounts[, coldata$sample_name %in% sample_names]
# dim(scecounts)
rowsums <- rowSums(scecounts) > 10000
table(rowsums)
scecounts <- scecounts[rowsums, ]
dim(scecounts)

## SCT NORMALIZE INDIVIDUAL MATRICES
for (sample in sample_names) {
  print(sample)
  data_subset <- scecounts[, coldata$sample_name == sample]
  seuratobject <- CreateSeuratObject(data_subset)
  sct_norm_matrix <- as.matrix(GetAssayData(SCTransform(seuratobject)))
  sct_filepath <- paste(here(), "/data/SCT_norms/", as.character(sample), "_SCT_norm_sparse.csv", sep = "")
  print(sct_filepath)
  write.csv(sct_norm_matrix, file = sct_filepath)
  # break
}

## DINO NORMALIZE INDIVIDUAL MATRICES 
for (sample in sample_names) {
  print(sample)
  dino_norm_matrix <- Dino(as.matrix(scecounts[, coldata$sample_name == sample]), nCores=4)
  dino_filepath <- paste(here(), "/data/Dino_norms/", as.character(sample), 
                        "_Dino_norm_sparse.csv", sep = "")
  print(dino_filepath)
  write.csv(dino_norm_matrix, file = dino_filepath)
  # break
}

test <- Dino(as.matrix(scecounts[, coldata$sample_name == 151507]))

# $meta.data
# $assays$SCT
# $assays$RNA
# $meta.data has interesting before&after stats going from the original "RNA" assay to "SCT" normalized data
attributes(sct)$meta.data

attributes(sct)$assays$RNA[1:10, 1:10]
attributes(sct)$assays$SCT[1:10, 1:10]
sct_norm_matrix <- as.matrix(GetAssayData(sct))
dim(sct_norm_matrix)
sct_norm_matrix[1:5, 1:5]

sct_filepath <- paste(here(), "/project/data/SCT_norm_sparse.mtx", sep = "")
Matrix:writeMM(sct_norm_matrix, file = sct_filepath)  # sparse, doesn't work because writeMM doesn't recognize a matrix object
write.csv(sct_norm_matrix, file = sct_filepath)

attributes(attributes(sct)$assays$SCT)


dinot <- Dino(f, ncores=0)  # ncores = 0 uses all cores 

dim(dinot)

write.csv(dinot, file = paste(here(), "/project/data/dino_norm.csv", sep = ""))               # Write dense matrix, about 5 mins and 720 MB
Matrix::writeMM(dinot, file = paste(here(), "/project/data/dino_norm_sparse.mtx", sep = ""))  # Write sparse matrix, 5 mins and 1.4 GB !!!
# Sparse matrix is bad because the data, after Dino, is no longer sparse

dinot[1:20, 1:20]






















# Get brain object from Seurat 
# InstallData("stxBrain")
# brain <- LoadData("stxBrain", type = "anterior1")
# typeof(brain)
# data("pbmcSmall")
# brain_sctnorm <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
# brain_sctnorm
# plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
# plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
# wrap_plots(plot1, plot2)











