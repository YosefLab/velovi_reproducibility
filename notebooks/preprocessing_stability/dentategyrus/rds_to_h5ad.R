library(Seurat)
library(zellkonverter)

BASE_PATH <- '../../../data/'
data <- readRDS(file = paste(BASE_PATH, 'dentategyrus/Dentate_gyrus_sce_quantifications_combined.rds', sep=''))
writeH5AD(data, file=paste(BASE_PATH, 'dentategyrus/Dentate_gyrus_adata_quantifications_combined.h5ad', sep=''), X_name='kallisto_bustools_prepref_isocollapse_exclude_spliced')