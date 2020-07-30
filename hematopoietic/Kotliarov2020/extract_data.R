#script for converting the filtered raw counts to a 10X object to work with in python

#install pacakges
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Seurat")
#BiocManager::install("DropletUtils")

#load libraries --------------------
library(Seurat)
library(DropletUtils)

#set working directory -------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#get session info ------------------
sessionInfo()

# R version 3.6.3 (2020-02-29)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Catalina 10.15.3
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods  
# [9] base     
# 
# other attached packages:
#   [1] DropletUtils_1.6.1          SingleCellExperiment_1.8.0  SummarizedExperiment_1.16.1
# [4] DelayedArray_0.12.3         BiocParallel_1.20.1         matrixStats_0.56.0         
# [7] Biobase_2.46.0              GenomicRanges_1.38.0        GenomeInfoDb_1.22.1        
# [10] IRanges_2.20.2              S4Vectors_0.24.4            BiocGenerics_0.32.0        
# [13] Seurat_3.1.5               
# 
# loaded via a namespace (and not attached):
#   [1] TH.data_1.0-10         Rtsne_0.15             colorspace_1.4-1      
# [4] ellipsis_0.3.0         ggridges_0.5.2         XVector_0.26.0        
# [7] rstudioapi_0.11        leiden_0.3.3           listenv_0.8.0         
# [10] npsurv_0.4-0.1         ggrepel_0.8.2          fansi_0.4.1           
# [13] mvtnorm_1.1-0          R.methodsS3_1.8.0      codetools_0.2-16      
# [16] splines_3.6.3          mnormt_1.5-7           lsei_1.2-0            
# [19] TFisher_0.2.0          jsonlite_1.6.1         ica_1.0-2             
# [22] cluster_2.1.0          R.oo_1.23.0            png_0.1-7             
# [25] uwot_0.1.8             HDF5Array_1.14.4       sctransform_0.2.1     
# [28] BiocManager_1.30.10    compiler_3.6.3         httr_1.4.1            
# [31] dqrng_0.2.1            assertthat_0.2.1       Matrix_1.2-18         
# [34] lazyeval_0.2.2         limma_3.42.2           cli_2.0.2             
# [37] htmltools_0.4.0        tools_3.6.3            rsvd_1.0.3            
# [40] igraph_1.2.5           gtable_0.3.0           glue_1.4.0            
# [43] GenomeInfoDbData_1.2.2 RANN_2.6.1             reshape2_1.4.4        
# [46] dplyr_0.8.5            Rcpp_1.0.4             vctrs_0.2.4           
# [49] multtest_2.42.0        gdata_2.18.0           ape_5.3               
# [52] nlme_3.1-147           gbRd_0.4-11            lmtest_0.9-37         
# [55] stringr_1.4.0          globals_0.12.5         lifecycle_0.2.0       
# [58] irlba_2.3.3            gtools_3.8.2           future_1.17.0         
# [61] edgeR_3.28.1           zlibbioc_1.32.0        MASS_7.3-51.6         
# [64] zoo_1.8-8              scales_1.1.0           sandwich_2.5-1        
# [67] rhdf5_2.30.1           RColorBrewer_1.1-2     yaml_2.2.1            
# [70] reticulate_1.15        pbapply_1.4-2          gridExtra_2.3         
# [73] ggplot2_3.3.0          stringi_1.4.6          mutoss_0.1-12         
# [76] plotrix_3.7-8          caTools_1.18.0         bibtex_0.4.2.2        
# [79] Rdpack_0.11-1          rlang_0.4.6            pkgconfig_2.0.3       
# [82] bitops_1.0-6           lattice_0.20-41        Rhdf5lib_1.8.0        
# [85] ROCR_1.0-11            purrr_0.3.4            patchwork_1.0.0       
# [88] htmlwidgets_1.5.1      cowplot_1.0.0          tidyselect_1.0.0      
# [91] RcppAnnoy_0.0.16       plyr_1.8.6             magrittr_1.5          
# [94] R6_2.4.1               gplots_3.0.3           multcomp_1.4-13       
# [97] pillar_1.4.4           sn_1.6-1               fitdistrplus_1.0-14   
# [100] survival_3.1-12        RCurl_1.98-1.2         tibble_3.0.1          
# [103] future.apply_1.5.0     tsne_0.1-3             crayon_1.3.4          
# [106] KernSmooth_2.23-17     plotly_4.9.2.1         locfit_1.5-9.4        
# [109] grid_3.6.3             data.table_1.12.8      metap_1.3             
# [112] digest_0.6.25          tidyr_1.0.2            numDeriv_2016.8-1.1   
# [115] R.utils_2.9.2          munsell_0.5.0          viridisLite_0.3.0 

#read R Object ---------------------
obj <- readRDS('R Objects/H1_day0_demultilexed_singlets.RDS')

#write to mtx file -----------------

#write out gene_expression data
#cant use the "raw" data because we do not have all of the metadata for each of the cells
write10xCounts(path = 'extracted_data/gene_expression', obj@data)

#write out citeseq data
write10xCounts(path = 'extracted_data/CITEseq', obj@assay$CITE@raw.data)

#write out metadata ---------------

write.table(obj@meta.data,file = "extracted_data/gene_expression/metadata.tsv", row.names = FALSE, sep = '\t')
write.table(obj@meta.data,file = "extracted_data/CITEseq/metadata.tsv", row.names = FALSE, sep = '\t')
