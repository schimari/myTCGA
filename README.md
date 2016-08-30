# myTCGA
R scripts for TCGA data retrieval and analysis

r_retrieve.r                        - retrieves TCGA data via RTCGA

r_analyse.survival.r                - calculates z-score and KM survival curve based on RNA-seq and clinical data (up and down regulation)

r_analyse.survival_logfc.r          - calculates logFC and KM survival curve based on RNA-seq and clinical data (up and down regulation)

r_TCGA_limma.r                      - uses a linear model via limma to calculate topTreat and topTable output from RNA-seq data

r_TCGA_limma_tumor_alive_vs_dead.r  - same as above but only data from tumor patients. Contrast is between alive or dead
