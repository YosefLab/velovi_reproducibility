library(Seurat)
library(zellkonverter)

BASE_PATH <- '../../../data/'
data <- readRDS(file = paste(BASE_PATH, 'pfc/PFC_sce_quantifications_combined.rds', sep=''))

writeH5AD(data, file=paste(BASE_PATH, 'pfc/PFC_adata_quantifications_combined.h5ad', sep=''), X_name='velocyto_spliced')
