Data source: https://www.10xgenomics.com/resources/datasets/10k-human-pbmcs-3-ht-v3-1-chromium-x-3-1-high

kb index source: https://www.kallistobus.tools/tutorials/kb_velocity/python/kb_velocity/#download-a-pre-built-human-rna-velocity-index

```
curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/6.1.0/10k_PBMC_3p_nextgem_Chromium_X/10k_PBMC_3p_nextgem_Chromium_X_fastqs.tar

tar -xvf 10k_PBMC_3p_nextgem_Chromium_X_fastqs.tar

kb ref -d linnarsson -i index.idx -g t2g.txt -c1 spliced_t2c.txt -c2 unspliced_t2c.txt

kb count --h5ad -i index.idx -g t2g.txt -x 10xv3 -o 10k_PBMC_3p_nextgem_Chromium_X \
-c1 spliced_t2c.txt -c2 unspliced_t2c.txt --workflow lamanno --filter bustools -t 2 \
10k_PBMC_3p_nextgem_Chromium_X_fastqs/10k_PBMC_3p_nextgem_Chromium_X_S1_L001_R1_001.fastq.gz \
10k_PBMC_3p_nextgem_Chromium_X_fastqs/10k_PBMC_3p_nextgem_Chromium_X_S1_L001_R2_001.fastq.gz \
10k_PBMC_3p_nextgem_Chromium_X_fastqs/10k_PBMC_3p_nextgem_Chromium_X_S1_L002_R1_001.fastq.gz \
10k_PBMC_3p_nextgem_Chromium_X_fastqs/10k_PBMC_3p_nextgem_Chromium_X_S1_L002_R2_001.fastq.gz \
10k_PBMC_3p_nextgem_Chromium_X_fastqs/10k_PBMC_3p_nextgem_Chromium_X_S1_L003_R1_001.fastq.gz \
10k_PBMC_3p_nextgem_Chromium_X_fastqs/10k_PBMC_3p_nextgem_Chromium_X_S1_L003_R2_001.fastq.gz \
10k_PBMC_3p_nextgem_Chromium_X_fastqs/10k_PBMC_3p_nextgem_Chromium_X_S1_L004_R1_001.fastq.gz \
10k_PBMC_3p_nextgem_Chromium_X_fastqs/10k_PBMC_3p_nextgem_Chromium_X_S1_L004_R2_001.fastq.gz
```
