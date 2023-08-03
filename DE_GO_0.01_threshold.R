# AMD server 
load("DE/ED_contrast_gene_n122.Rdata")
library(org.Hs.eg.db)
library(clusterProfiler)
library(cowplot)
library(writexl)

packageVersion("clusterProfiler")
# [1] ‘4.6.2’
packageVersion("org.Hs.eg.db")
# [1] ‘3.16.0’

geneUniverse = as.character(sigGeneOV$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

id=which(sigGeneCNT$`p_ED-Control`<0.01)
sigem=as.character(sigGeneCNT$EntrezID[id])
sigem=sigem[!is.na(sigem)]
length(sigem)
# [1] 951

goBP=enrichGO(gene          = sigem,
               universe      = geneUniverse,
               OrgDb         = org.Hs.eg.db,
               ont           = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.1,
               qvalueCutoff  = 0.5,
               readable      = TRUE)
goMF=enrichGO(gene          = sigem,
              universe      = geneUniverse,
              OrgDb         = org.Hs.eg.db,
              ont           = "MF",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.1,
              qvalueCutoff  = 0.5,
              readable      = TRUE)
# no sig
goCC=enrichGO(gene          = sigem,
              universe      = geneUniverse,
              OrgDb         = org.Hs.eg.db,
              ont           = "CC",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.1,
              qvalueCutoff  = 0.5,
              readable      = TRUE) 
# no sig

# plot
dbp=dotplot(goBP,showCategory = 20,title="ED-CT only p<0.01 Biological Process GO")

cowplot::plot_grid(dbp,dmf,dcc,ncol=3)

write_xlsx(list("BP"=goBP@result,"MF"=goMF@result,"CC"=goCC@result), 
				'./20230803_01_ED-CT.xlsx')

print(sessionInfo())
[1] rappdirs_0.3.3           R.methodsS3_1.8.2        tidyr_1.3.0             
  [4] bumphunter_1.40.0        bit64_4.0.5              DelayedArray_0.24.0     
  [7] R.utils_2.12.2           data.table_1.14.8        rpart_4.1.19            
 [10] KEGGREST_1.38.0          RCurl_1.98-1.12          GEOquery_2.66.0         
 [13] derfinder_1.32.0         generics_0.1.3           GenomicFeatures_1.50.4  
 [16] callr_3.7.3              cowplot_1.1.1            RSQLite_2.3.1           
 [19] shadowtext_0.1.2         bit_4.0.5                tzdb_0.4.0              
 [22] enrichplot_1.18.3        webshot_0.5.5            xml2_1.3.5              
 [25] httpuv_1.6.11            viridis_0.6.3            gargle_1.5.2            
 [28] xfun_0.39                hms_1.1.3                evaluate_0.21           
 [31] promises_1.2.0.1         fansi_1.0.4              restfulr_0.0.15         
 [34] progress_1.2.2           dbplyr_2.3.3             igraph_1.5.0            
 [37] htmlwidgets_1.6.2        googledrive_2.1.1        purrr_1.0.1             
 [40] ellipsis_0.3.2           dplyr_1.1.2              backports_1.4.1         
 [43] annotate_1.76.0          biomaRt_2.54.1           vctrs_0.6.3             
 [46] remotes_2.4.2.1          cachem_1.0.8             withr_2.5.0             
 [49] ggforce_0.4.1            HDO.db_0.99.1            BSgenome_1.66.3         
 [52] checkmate_2.2.0          treeio_1.22.0            GenomicAlignments_1.34.1
 [55] prettyunits_1.1.1        svglite_2.1.1            cluster_2.1.4           
 [58] lazyeval_0.2.2           ape_5.7-1                segmented_1.6-4         
 [61] crayon_1.5.2             labeling_0.4.2           pkgconfig_2.0.3         
 [64] tweenr_2.0.2             pkgload_1.3.2.1          nnet_7.3-18             
 [67] rlang_1.1.1              lifecycle_1.0.3          miniUI_0.1.1.1          
 [70] downloader_0.4           filelock_1.0.2           BiocFileCache_2.6.1     
 [73] cellranger_1.1.0         polyclip_1.10-4          rngtools_1.5.2          
 [76] aplot_0.1.10             base64enc_0.1-3          processx_3.8.2          
 [79] png_0.1-8                viridisLite_0.4.2        rjson_0.2.21            
 [82] bitops_1.0-7             gson_0.1.0               R.oo_1.25.0             
 [85] Biostrings_2.66.0        blob_1.2.4               doRNG_1.8.6             
 [88] stringr_1.5.0            qvalue_2.30.0            readr_2.1.4             
 [91] gridGraphics_0.5-1       scales_1.2.1             memoise_2.0.1           
 [94] magrittr_2.0.3           plyr_1.8.8               derfinderHelper_1.32.0  
 [97] zlibbioc_1.44.0          scatterpie_0.1.8         compiler_4.2.2          
[100] BiocIO_1.8.0             RColorBrewer_1.1-3       Rsamtools_2.14.0        
[103] cli_3.6.1                XVector_0.38.0           urlchecker_1.0.1        
[106] patchwork_1.1.2          ps_1.7.5                 htmlTable_2.4.1         
[109] Formula_1.2-5            MASS_7.3-60              tidyselect_1.2.0        
[112] stringi_1.7.12           highr_0.10               yaml_2.3.7              
[115] GOSemSim_2.24.0          locfit_1.5-9.8           ggrepel_0.9.3           
[118] grid_4.2.2               VariantAnnotation_1.44.1 fastmatch_1.1-3         
[121] tools_4.2.2              parallel_4.2.2           rstudioapi_0.15.0       
[124] foreach_1.5.2            foreign_0.8-84           gridExtra_2.3           
[127] farver_2.1.1             ggraph_2.1.0             digest_0.6.33           
[130] shiny_1.7.4.1            Rcpp_1.0.11              later_1.3.1             
[133] httr_1.4.6               colorspace_2.1-0         rvest_1.0.3             
[136] XML_3.99-0.14            fs_1.6.3                 splines_4.2.2           
[139] yulab.utils_0.0.6        tidytree_0.4.4           graphlayouts_1.0.0      
[142] ggplotify_0.1.1          sessioninfo_1.2.2        systemfonts_1.0.4       
[145] xtable_1.8-4             ggtree_3.6.2             jsonlite_1.8.7          
[148] GenomicFiles_1.34.0      tidygraph_1.2.3          ggfun_0.1.1             
[151] R6_2.5.1                 profvis_0.3.8            Hmisc_5.0-1             
[154] pillar_1.9.0             htmltools_0.5.5          mime_0.12               
[157] glue_1.6.2               fastmap_1.1.1            codetools_0.2-18        
[160] fgsea_1.24.0             pkgbuild_1.4.2           utf8_1.2.3              
[163] lattice_0.21-8           tibble_3.2.1             curl_5.0.1              
[166] rentrez_1.2.3            GO.db_3.16.0             rmarkdown_2.23          
[169] munsell_0.5.0            GenomeInfoDbData_1.2.9   iterators_1.0.14        
[172] reshape2_1.4.4           gtable_0.3.3            
