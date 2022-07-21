# Preprocessing analysis for RNA velocity in mouse brain

To reproduce results, the RDS file `OldBrain_sce_quantifications_combined.rds` provided by [Sonesen *et al.*](https://doi.org/10.1371/journal.pcbi.1008585) is needed. The scripts and notebooks need to be run in the following order:

1. `rds_to_h5ad.R`
2. `adata_generation.ipynb`
3. `velocyto_vi.ipynb`
4. Every notebook specific to one preprocessing algorithm and RNA velocity inference approach: `*_vi.ipynb`, `*_em_ss.ipynb`
5. `correlation_analysis.ipynb`
