# cytoChain
This is a first prototype of cytoChain, a Shiny web app to handle high-dimensional cytometry datasets. A beta release is currently deployed at the `shinyapps.io` server and it is available at the following [link](//abbatidanilo.shinyapps.io/cytoChain)

In case you are interested in cytoChain for your cytometric analisys, drop a mail with a description of your goals, experiments and aims to abbati.danilo@hsr.it. 

cytoChain has been deployed by Danilo Abbati at the Experimental Hematology department of the San Raffaele Hospital - Milan (Italy). 

## Get the package and run the application
You should install **shiny** and all the other libraries listed below in the sessionInfo() section and then run:

```
shiny::runGitHub("abbatidanilo/cytoChain", username = abbatidanilo, ref = "master", subdir = NULL, port = NULL, launch.browser = getOption("shiny.launch.browser", interactive()))
```

## Operation 
Load your compensated samples, namely your *.FCS* (standard 3.0) files with the *File Input* control slot in the `Pre-clustering workflow` or in the `Metadata and assays` main side-panel and go throught step by step, all the tabs you need. It is possible to set all the relevant parameters offered in the various side-panels and you may also edit some meta-data related to your samples in the `Metadata and assays` main side-panel. The various analysis are offered within the `Clustering workflow` main side-panel and the related sub-choices.

[Manual available](//github.com/abbatidanilo/cytoChain/tree/main/www). The application will be keep on updating with new features

## Citation
If you use cytoChain, please use the following citation: "Flow cytometry data mining by cytoChain identifies determinants of exhaustion and stemness in TCR-engineered T cells" by **Francesco Manfredi**, **Danilo Abbati** et al. (doi: 10.1002/eji.202049103). This paper was published in the European Journal of Immunology in August 2021.

## sessionInfo()
sessionInfo()
R version 4.3.1 (2023-06-16 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default


locale:
[1] LC_COLLATE=English_United Kingdom.utf8  LC_CTYPE=English_United Kingdom.utf8    LC_MONETARY=English_United Kingdom.utf8 LC_NUMERIC=C                           
[5] LC_TIME=English_United Kingdom.utf8    

time zone: Europe/Rome
tzcode source: internal

attached base packages:
[1] grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] flowVS_1.34.0               umap_0.2.10.0               Rtsne_0.16                  ConsensusClusterPlus_1.66.0 FastPG_0.0.8               
 [6] Rphenograph_0.99.1          DDoutlier_0.1.0             FlowSOM_2.10.0              kohonen_3.0.12              mclust_6.0.0               
[11] factoextra_1.0.7            matrixStats_1.0.0           MASS_7.3-60                 spade_1.10.4                Rclusterpp_0.2.6           
[16] Rcpp_1.0.11                 igraph_1.5.1                DepecheR_1.18.0             flowStats_4.14.0            flowTrans_1.54.0           
[21] flowClust_3.40.0            flowViz_1.66.0              lattice_0.22-5              flowAI_1.32.0               flowWorkspace_4.14.0       
[26] flowCore_2.14.0             audio_0.1-11                fs_1.6.3                    ggpubr_0.6.0                streamgraph_0.9.0          
[31] gridExtra_2.3               png_0.1-8                   cluster_2.1.4               pheatmap_1.0.12             ggrepel_0.9.4              
[36] rhandsontable_0.3.8         reshape2_1.4.4              RColorBrewer_1.1-3          plotly_4.10.3               lubridate_1.9.3            
[41] forcats_1.0.0               stringr_1.5.0               dplyr_1.1.3                 purrr_1.0.2                 readr_2.1.4                
[46] tidyr_1.3.0                 tibble_3.2.1                ggplot2_3.4.4               tidyverse_2.0.0             DT_0.30                    
[51] shinyjs_2.1.0               shinyBS_0.61.1              shinyWidgets_0.8.0          shinycssloaders_1.0.0       shinydashboard_0.7.2       
[56] shinythemes_1.2.0           shiny_1.7.5.1              

loaded via a namespace (and not attached):
  [1] splines_4.3.1       later_1.3.1         bitops_1.0-7        polyclip_1.10-6     graph_1.80.0        xts_0.13.1          XML_3.99-0.15       deSolve_1.38       
  [9] lifecycle_1.0.4     rstatix_0.7.2       crosstalk_1.2.0     backports_1.4.1     magrittr_2.0.3      sass_0.4.7          rmarkdown_2.25      yaml_2.3.7         
 [17] jquerylib_0.1.4     httpuv_1.6.12       collapse_2.0.3      askpass_1.2.0       reticulate_1.34.0   abind_1.4-5         zlibbioc_1.48.0     BiocGenerics_0.48.1
 [25] RCurl_1.98-1.13     pracma_2.4.2        tweenr_2.0.2        S4Vectors_0.40.1    gmodels_2.18.1.1    gdata_3.0.0         moments_0.14.1      ellipse_0.5.0      
 [33] RSpectra_0.16-1     codetools_0.2-19    ggforce_0.4.1       tidyselect_1.2.0    farver_2.1.1        viridis_0.6.4       beanplot_1.3.1      stats4_4.3.1       
 [41] jsonlite_1.8.7      ks_1.14.1           ellipsis_0.3.2      iterators_1.0.14    foreach_1.5.2       ggnewscale_0.4.9    tools_4.3.1         snow_0.4-4         
 [49] rARPACK_0.11-0      glue_1.6.2          mnormt_2.1.1        xfun_0.41           mixOmics_6.26.0     withr_2.5.2         BiocManager_1.30.22 fastmap_1.1.1      
 [57] latticeExtra_0.6-30 fansi_1.0.5         openssl_2.1.1       caTools_1.18.2      digest_0.6.33       timechange_0.2.0    R6_2.5.1            mime_0.12          
 [65] colorspace_2.1-0    gtools_3.9.4        jpeg_0.1-10         utf8_1.2.4          generics_0.1.3      hexbin_1.28.3       data.table_1.14.8   corpcor_1.6.10     
 [73] FNN_1.1.3.2         robustbase_0.99-0   httr_1.4.7          htmlwidgets_1.6.2   IDPmisc_1.1.20      pkgconfig_2.0.3     gtable_0.3.4        changepoint_2.2.4  
 [81] RProtoBufLib_2.14.0 pcaPP_2.0-3         htmltools_0.5.7     carData_3.0-5       clue_0.3-65         scales_1.2.1        Biobase_2.62.0      doSNOW_1.0.20      
 [89] colorRamps_2.3.1    knitr_1.45          rstudioapi_0.15.0   tzdb_0.4.0          cachem_1.0.8        zoo_1.8-12          KernSmooth_2.23-22  fda_6.1.4          
 [97] pillar_1.9.0        vctrs_0.6.4         RANN_2.6.1          gplots_3.1.3        promises_1.2.1      car_3.1-2           cytolib_2.14.0      xtable_1.8-4       
[105] Rgraphviz_2.46.0    evaluate_0.23       mvtnorm_1.2-3       cli_3.6.1           compiler_4.3.1      rlang_1.1.2         ncdfFlow_2.48.0     rrcov_1.7-4        
[113] ggsignif_0.6.4      interp_1.1-4        fds_1.8             plyr_1.8.9          stringi_1.7.12      rainbow_3.7         BiocParallel_1.36.0 viridisLite_0.4.2  
[121] deldir_1.0-9        hdrcde_3.4          munsell_0.5.0       lazyeval_0.2.2      Matrix_1.6-1.1      hms_1.1.3           fontawesome_0.5.2   memoise_2.0.1      
[129] broom_1.0.5         RcppParallel_5.1.7  bslib_0.5.1         DEoptimR_1.1-3     
