library(Seurat)
library(zellkonverter)

BASE_PATH <- '../../../data/'
data <- readRDS(file = paste(BASE_PATH, 'pancreas/Pancreas_sce_quantifications_combined.rds', sep=''))

writeH5AD(data, file=paste(BASE_PATH, 'pancreas/Pancreas_adata_quantifications_combined.h5ad', sep=''), X_name='velocyto_spliced')
