# load necessary libraries
library(Seurat)
library(SeuratObject)
library(tidyr)
library(ggplot2)
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
install.packages("BiocManager")
BiocManager::install("decontX")
BiocManager::install("SingleCellExperiment")
library(decontX)
library(SingleCellExperiment)

# load raw data into counts matrix format
dirs <- list.dirs('C:/Users/tanta/Downloads/GSE220797_Processed_V2/', 
                  recursive = F, full.names = F)
for (x in dirs) {
  cts <- ReadMtx(
    mtx = paste0('C:/Users/tanta/Downloads/GSE220797_Processed_V2/', x,
                       '/matrix.mtx.gz'), 
    features = paste0('C:/Users/tanta/Downloads/GSE220797_Processed_V2/',
                             x,'/features.tsv.gz'), 
    cells = paste0('C:/Users/tanta/Downloads/GSE220797_Processed_V2/',
                         x,'/barcodes.tsv.gz'))
  name <- basename(x)
  
  assign(name, CreateSeuratObject(counts = cts))
}

# create list of object names
object_names <- c(
  "patient1_2hr", "patient2_2hr", "patient3_2hr", "patient4_2hr",
  "patient5_2hr", "patient6_2hr", "patient1_cit", "patient2_cit",
  "patient3_cit", "patient4_cit", "patient5_cit", "patient6_cit"
)

# save objects as RDS files
for (name in object_names) {
  obj <- get(name)
  saveRDS(obj, paste0(name, ".rds"))  
}

# read RDS files
for (name in object_names) {
  obj <- readRDS(file = paste0(name, ".rds"))
  assign(name, obj, envir = .GlobalEnv) 
}

# Define an QC function that filters based on 1) high mDNA content 2) small
# library size and 3) outliers 
qc_filter <- function(obj, mito_pattern = "^MT-", mads_thresh = 4, 
                      hard_thresh = 50, min_features = 200, 
                      filt_intercept = 100, filt_slope = 0.055) {
  obj <- PercentageFeatureSet(obj, pattern = mito_pattern, 
                              col.name = "pct_counts_Mito")
  mito_thresh <- median(obj$pct_counts_Mito) + mad(obj$pct_counts_Mito) * 
    mads_thresh
  drop_mito <- obj$pct_counts_Mito > mito_thresh | obj$pct_counts_Mito > 
    hard_thresh
  obj <- obj[, !drop_mito]
  drop_Lib <- obj$nFeature_RNA < min_features
  obj <- obj[, !drop_Lib]
  to_inspect <- obj$nFeature_RNA < (obj$nCount_RNA * filt_slope + filt_intercept)
  obj <- obj[, !to_inspect]
  return(obj)
}

# standard pre-processing workflow
for (name in object_names) {
  
  # print the sample we are on
  print(name)
  
  # retrieve object from name
  obj <- get(name)
  
  # standard pre-processing workflow
  obj_filtered <- qc_filter(obj)
  obj_filtered <- NormalizeData(obj_filtered)
  obj_filtered <- FindVariableFeatures(obj_filtered)
  obj_filtered <- ScaleData(obj_filtered)
  obj_filtered <- RunPCA(obj_filtered, nfeatures.print = 10)

  # Find significant PCs
  stdv <- obj_filtered[["pca"]]@stdev
  sum.stdv <- sum(obj_filtered[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  min.pc
  
  # finish pre-processing
  obj_filtered <- RunUMAP(obj_filtered, dims = 1:min.pc)
  obj_filtered <- FindNeighbors(object = obj_filtered, dims = 1:min.pc)              
  obj_filtered <- FindClusters(object = obj_filtered, resolution = 0.1)
  
  # assign filtered object back to original name
  assign(name, obj_filtered, envir = .GlobalEnv)
}

# run a for loop to identify doublets in each object
for (name in object_names) {
  
  # print the sample we are on
  print(name)
  
  # retrieve object from name
  obj <- get(name)
  
  # pK identification (no ground-truth)
  sweep.list <- paramSweep(obj, PCs = 1:min.pc)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  # homotypic doublet proportion estimate
  annotations <- obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  
  # adjust doublet rate based on total cell count in each sample
  if (nrow(obj@meta.data) < 5000) {
    doublet_rate <- 0.025  
  } else if (nrow(obj@meta.data) < 10000) {
    doublet_rate <- 0.05  
  } else {
    doublet_rate <- 0.075  
  }
  nExp.poi <- round(doublet_rate * nrow(obj@meta.data))
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  obj <- doubletFinder(seu = obj, PCs = 1:min.pc, pK = optimal.pk, 
                       nExp = nExp.poi.adj)
  
  # assign filtered object back to original name
  assign(name, obj, envir = .GlobalEnv)
}  

# remove doublets from each object
for (name in object_names) {
  
  # print the sample we are on
  print(name)
  
  # retrieve object from name
  obj <- get(name)
  
  # change column header "DF.classifications_0.25_0.26_799" to "doublet_finder"
  metadata <- obj@meta.data
  colnames(metadata)[8] <- "doublet_finder"
  obj@meta.data <- metadata 
  
  # subset and save
  obj_singlets <- subset(obj, doublet_finder == "Singlet")
  obj <- obj_singlets
  remove(obj_singlets)
  
  # assign filtered object back to original name
  assign(name, obj, envir = .GlobalEnv)
}

# run decontx to remove ambient RNA contamination
for (name in object_names) {
  
  # print the sample we are on
  print(name)
  
  # retrieve object from name
  obj <- get(name)
  
  counts <- GetAssayData(object = obj, assay = "RNA", layer = "counts")
  sce <- SingleCellExperiment(list(counts = counts))
  sce <- decontX(sce)
  obj[["decontXcounts"]] <- CreateAssayObject(counts = decontXcounts(sce))
  
  # assign filtered object back to original name
  assign(name, obj, envir = .GlobalEnv)
}

# change default assay from "RNA" to "decontXcounts"
DefaultAssay(patient1_2hr) <- "decontXcounts"

# prepare for integration
for (name in object_names) {
  
  # print the sample we are on
  print(name)
  
  # retrieve object from name
  obj <- get(name)
  
  # standard pre-processing workflow
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  
  # assign filtered object back to original name
  assign(name, obj, envir = .GlobalEnv)
}

# create an object list for integration
obj_list <- mget(object_names)

# select integration features
features <- SelectIntegrationFeatures(object.list = obj_list)

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj_list, 
                                  anchor.features = features)

# integrate data
seurat_integrated <- IntegrateData(anchorset = anchors)

# scale data, run PCA and UMAP and visualize integrated data
seurat_integrated <- ScaleData(seurat_integrated)
seurat_integrated <- RunPCA(seurat_integrated)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:50)
