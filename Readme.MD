# cytoChain

cytoChain has been deployed by Danilo Abbati at the Experimental Hematology department of the San Raffaele Hospital - Milan (Italy).

In case you are interested in cytoChain for your cytometric analisys, drop a mail with a description of your goals, experiments and aims to abbati.danilo@hsr.it.

Please notice: this is still a prototype of cytoChain, a Shiny web app to handle high-dimensional cytometry datasets. A beta release is currently deployed at the `shinyapps.io` server and it is available at the following [link](//abbatidanilo.shinyapps.io/cytoChain). Notice though that this current public server implementation does not support heavy experiments. A more stable solution could be achieved deploying locally, downloading these source R files.

## Get the package and run the application

You should install **shiny** and all the other libraries listed below in the sessionInfo() section and then run:

```
shiny::runGitHub("abbatidanilo/cytoChain", username = abbatidanilo, ref = "master", subdir = NULL, port = NULL, launch.browser = getOption("shiny.launch.browser", interactive()))
```

## Operation

Load your compensated samples, namely your _.FCS_ (standard 3.0) files with the _File Input_ control slot in the `Pre-clustering workflow` or in the `Metadata and assays` main side-panel and go throught step by step, all the tabs you need. It is possible to set all the relevant parameters offered in the various side-panels and you may also edit some meta-data related to your samples in the `Metadata and assays` main side-panel. The various analysis are offered within the `Clustering workflow` main side-panel and the related sub-choices.

[Manual available](//github.com/abbatidanilo/cytoChain/tree/main/www). The application will be keep on updating with new features

## Citation

If you use cytoChain, please use the following citation: "Flow cytometry data mining by cytoChain identifies determinants of exhaustion and stemness in TCR-engineered T cells" by **Francesco Manfredi**, **Danilo Abbati** et al. (doi: 10.1002/eji.202049103). This paper was published in the European Journal of Immunology in August 2021.

## sessionInfo()

> sessioninfo::session_info()
> ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
> setting value
> version R version 4.4.1 (2024-06-14 ucrt)
> os Windows 10 x64 (build 19045)
> system x86_64, mingw32
> ui RStudio
> language (EN)
> collate English_United Kingdom.utf8
> ctype English_United Kingdom.utf8
> tz Europe/Rome
> date 2025-04-15
> rstudio 2024.12.1+563 Kousa Dogwood (desktop)
> pandoc 3.2 @ C:/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)

─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
! package _ version date (UTC) lib source
abind 1.4-8 2024-09-12 [1] CRAN (R 4.4.1)
askpass 1.2.1 2024-10-04 [1] CRAN (R 4.4.1)
audio _ 0.1-11 2023-08-18 [1] CRAN (R 4.4.0)
backports 1.5.0 2024-05-23 [1] CRAN (R 4.4.0)
beanplot 1.3.1 2022-04-09 [1] CRAN (R 4.4.1)
Biobase 2.64.0 2024-05-01 [1] Bioconductor 3.19 (R 4.4.0)
BiocGenerics 0.50.0 2024-05-01 [1] Bioconductor 3.19 (R 4.4.0)
BiocManager 1.30.25 2024-08-28 [1] CRAN (R 4.4.1)
BiocParallel 1.38.0 2024-05-01 [1] Bioconductor 3.19 (R 4.4.0)
bitops 1.0-9 2024-10-03 [1] CRAN (R 4.4.1)
broom 1.0.7 2024-09-26 [1] CRAN (R 4.4.1)
bslib 0.8.0 2024-07-29 [1] CRAN (R 4.4.1)
cachem 1.1.0 2024-05-16 [1] CRAN (R 4.4.1)
car 3.1-3 2024-09-27 [1] CRAN (R 4.4.1)
carData 3.0-5 2022-01-06 [1] CRAN (R 4.4.1)
caTools 1.18.3 2024-09-04 [1] CRAN (R 4.4.1)
changepoint 2.2.4 2022-11-03 [1] CRAN (R 4.4.1)
cli 3.6.3 2024-06-21 [1] CRAN (R 4.4.1)
clue 0.3-65 2023-09-23 [1] CRAN (R 4.4.1)
cluster _ 2.1.6 2023-12-01 [1] CRAN (R 4.4.1)
ClusterR 1.3.3 2024-06-18 [1] CRAN (R 4.4.1)
codetools 0.2-20 2024-03-31 [1] CRAN (R 4.4.1)
collapse 2.0.16 2024-08-21 [1] CRAN (R 4.4.1)
colorRamps 2.3.4 2024-03-07 [1] CRAN (R 4.4.0)
colorspace 2.1-1 2024-07-26 [1] CRAN (R 4.4.1)
ConsensusClusterPlus _ 1.68.0 2024-05-01 [1] Bioconductor 3.19 (R 4.4.0)
corpcor 1.6.10 2021-09-16 [1] CRAN (R 4.4.0)
crosstalk 1.2.1 2023-11-23 [1] CRAN (R 4.4.1)
cytolib 2.16.0 2024-05-01 [1] Bioconductor 3.19 (R 4.4.0)
data.table 1.16.2 2024-10-10 [1] CRAN (R 4.4.1)
DDoutlier _ 0.1.0 2018-05-30 [1] CRAN (R 4.4.1)
deldir 2.0-4 2024-02-28 [1] CRAN (R 4.4.0)
DEoptimR 1.1-3 2023-10-07 [1] CRAN (R 4.4.0)
DepecheR _ 1.20.0 2024-05-01 [1] Bioconductor 3.19 (R 4.4.0)
deSolve 1.40 2023-11-27 [1] CRAN (R 4.4.1)
devtools _ 2.4.5 2022-10-11 [1] CRAN (R 4.4.1)
digest 0.6.37 2024-08-19 [1] CRAN (R 4.4.1)
doSNOW 1.0.20 2022-02-04 [1] CRAN (R 4.4.1)
dplyr _ 1.1.4 2023-11-17 [1] CRAN (R 4.4.1)
DT _ 0.33 2024-04-04 [1] CRAN (R 4.4.1)
ellipse 0.5.0 2023-07-20 [1] CRAN (R 4.4.1)
ellipsis 0.3.2 2021-04-29 [1] CRAN (R 4.4.1)
evaluate 1.0.1 2024-10-10 [1] CRAN (R 4.4.1)
factoextra _ 1.0.7 2020-04-01 [1] CRAN (R 4.4.1)
fansi 1.0.6 2023-12-08 [1] CRAN (R 4.4.1)
farver 2.1.2 2024-05-13 [1] CRAN (R 4.4.1)
fastmap 1.2.0 2024-05-15 [1] CRAN (R 4.4.1)
FastPG _ 0.0.8 2024-10-23 [1] Github (sararselitsky/FastPG@44c9282)
fda 6.2.0 2024-09-17 [1] CRAN (R 4.4.1)
fds 1.8 2018-10-31 [1] CRAN (R 4.4.1)
flowAI _ 1.34.0 2024-05-01 [1] Bioconductor 3.19 (R 4.4.0)
flowClust _ 3.42.0 2024-05-01 [1] Bioconductor 3.19 (R 4.4.0)
flowCore _ 2.16.0 2024-05-01 [1] Bioconductor 3.19 (R 4.4.0)
FlowSOM _ 2.12.0 2024-05-01 [1] Bioconductor 3.19 (R 4.4.0)
flowStats _ 4.16.0 2024-05-01 [1] Bioconductor 3.19 (R 4.4.0)
flowTrans _ 1.56.0 2024-05-01 [1] Bioconductor 3.19 (R 4.4.0)
flowViz _ 1.68.0 2024-05-01 [1] Bioconductor 3.19 (R 4.4.0)
flowVS _ 1.36.0 2024-05-01 [1] Bioconductor 3.19 (R 4.4.0)
flowWorkspace _ 4.16.0 2024-05-01 [1] Bioconductor 3.19 (R 4.4.0)
FNN 1.1.4.1 2024-09-22 [1] CRAN (R 4.4.1)
fontawesome 0.5.2 2023-08-19 [1] CRAN (R 4.4.1)
forcats _ 1.0.0 2023-01-29 [1] CRAN (R 4.4.1)
foreach 1.5.2 2022-02-02 [1] CRAN (R 4.4.1)
Formula 1.2-5 2023-02-24 [1] CRAN (R 4.4.0)
fs _ 1.6.4 2024-04-25 [1] CRAN (R 4.4.1)
gdata 3.0.1 2024-10-22 [1] CRAN (R 4.4.1)
generics 0.1.3 2022-07-05 [1] CRAN (R 4.4.1)
ggforce 0.4.2 2024-02-19 [1] CRAN (R 4.4.1)
ggnewscale 0.5.0 2024-07-19 [1] CRAN (R 4.4.1)
ggplot2 _ 3.5.1 2024-04-23 [1] CRAN (R 4.4.1)
ggpubr _ 0.6.0 2023-02-10 [1] CRAN (R 4.4.1)
ggrepel _ 0.9.6 2024-09-07 [1] CRAN (R 4.4.1)
ggsignif 0.6.4 2022-10-13 [1] CRAN (R 4.4.1)
glue 1.8.0 2024-09-30 [1] CRAN (R 4.4.1)
gmodels 2.19.1 2024-03-06 [1] CRAN (R 4.4.1)
gmp 0.7-5 2024-08-23 [1] CRAN (R 4.4.1)
gplots 3.2.0 2024-10-05 [1] CRAN (R 4.4.1)
graph 1.82.0 2024-05-01 [1] Bioconductor 3.19 (R 4.4.0)
gridExtra _ 2.3 2017-09-09 [1] CRAN (R 4.4.1)
gtable 0.3.5 2024-04-22 [1] CRAN (R 4.4.1)
gtools 3.9.5 2023-11-20 [1] CRAN (R 4.4.1)
hdrcde 3.4 2021-01-18 [1] CRAN (R 4.4.1)
hexbin 1.28.4 2024-09-04 [1] CRAN (R 4.4.1)
hms 1.1.3 2023-03-21 [1] CRAN (R 4.4.1)
htmltools 0.5.8.1 2024-04-04 [1] CRAN (R 4.4.1)
htmlwidgets 1.6.4 2023-12-06 [1] CRAN (R 4.4.1)
httpuv 1.6.15 2024-03-26 [1] CRAN (R 4.4.1)
httr 1.4.7 2023-08-15 [1] CRAN (R 4.4.1)
IDPmisc 1.1.21 2024-02-08 [1] CRAN (R 4.4.1)
igraph _ 2.1.1 2024-10-19 [1] CRAN (R 4.4.1)
interp 1.1-6 2024-01-26 [1] CRAN (R 4.4.1)
iterators 1.0.14 2022-02-05 [1] CRAN (R 4.4.1)
jpeg 0.1-10 2022-11-29 [1] CRAN (R 4.4.0)
jquerylib 0.1.4 2021-04-26 [1] CRAN (R 4.4.1)
jsonlite 1.8.9 2024-09-20 [1] CRAN (R 4.4.1)
KernSmooth 2.23-24 2024-05-17 [1] CRAN (R 4.4.1)
knitr _ 1.48 2024-07-07 [1] CRAN (R 4.4.1)
kohonen _ 3.0.12 2023-06-09 [1] CRAN (R 4.4.1)
ks 1.14.3 2024-09-20 [1] CRAN (R 4.4.1)
later 1.3.2 2023-12-06 [1] CRAN (R 4.4.1)
lattice _ 0.22-6 2024-03-20 [1] CRAN (R 4.4.1)
latticeExtra 0.6-30 2022-07-04 [1] CRAN (R 4.4.1)
lazyeval 0.2.2 2019-03-15 [1] CRAN (R 4.4.1)
lifecycle 1.0.4 2023-11-07 [1] CRAN (R 4.4.1)
lubridate _ 1.9.3 2023-09-27 [1] CRAN (R 4.4.1)
magrittr 2.0.3 2022-03-30 [1] CRAN (R 4.4.1)
MASS _ 7.3-61 2024-06-13 [1] CRAN (R 4.4.1)
Matrix 1.7-1 2024-10-18 [1] CRAN (R 4.4.1)
matrixStats _ 1.4.1 2024-09-08 [1] CRAN (R 4.4.1)
mclust _ 6.1.1 2024-04-29 [1] CRAN (R 4.4.1)
memoise 2.0.1 2021-11-26 [1] CRAN (R 4.4.1)
mime 0.12 2021-09-28 [1] CRAN (R 4.4.0)
miniUI 0.1.1.1 2018-05-18 [1] CRAN (R 4.4.1)
mixOmics 6.28.0 2024-05-01 [1] Bioconductor 3.19 (R 4.4.0)
mnormt 2.1.1 2022-09-26 [1] CRAN (R 4.4.0)
moments 0.14.1 2022-05-02 [1] CRAN (R 4.4.0)
munsell 0.5.1 2024-04-01 [1] CRAN (R 4.4.1)
mvtnorm 1.3-1 2024-09-03 [1] CRAN (R 4.4.1)
ncdfFlow 2.50.0 2024-05-01 [1] Bioconductor 3.19 (R 4.4.0)
openssl 2.2.2 2024-09-20 [1] CRAN (R 4.4.1)
pcaPP 2.0-5 2024-08-19 [1] CRAN (R 4.4.1)
pheatmap _ 1.0.12 2019-01-04 [1] CRAN (R 4.4.1)
pillar 1.9.0 2023-03-22 [1] CRAN (R 4.4.1)
pkgbuild 1.4.4 2024-03-17 [1] CRAN (R 4.4.1)
pkgconfig 2.0.3 2019-09-22 [1] CRAN (R 4.4.1)
pkgload 1.4.0 2024-06-28 [1] CRAN (R 4.4.1)
plotly _ 4.10.4 2024-01-13 [1] CRAN (R 4.4.1)
plyr 1.8.9 2023-10-02 [1] CRAN (R 4.4.1)
png _ 0.1-8 2022-11-29 [1] CRAN (R 4.4.0)
polyclip 1.10-7 2024-07-23 [1] CRAN (R 4.4.1)
pracma 2.4.4 2023-11-10 [1] CRAN (R 4.4.1)
profvis 0.4.0 2024-09-20 [1] CRAN (R 4.4.1)
promises 1.3.0 2024-04-05 [1] CRAN (R 4.4.1)
purrr _ 1.0.2 2023-08-10 [1] CRAN (R 4.4.1)
R6 2.5.1 2021-08-19 [1] CRAN (R 4.4.1)
rainbow 3.8 2024-01-23 [1] CRAN (R 4.4.1)
RANN 2.6.2 2024-08-25 [1] CRAN (R 4.4.1)
rARPACK 0.11-0 2016-03-10 [1] CRAN (R 4.4.1)
Rclusterpp _ 0.2.6 2023-11-07 [1] Github (nolanlab/Rclusterpp@a073806)
RColorBrewer _ 1.1-3 2022-04-03 [1] CRAN (R 4.4.0)
Rcpp _ 1.0.13 2024-07-17 [1] CRAN (R 4.4.1)
D RcppParallel 5.1.9 2024-08-19 [1] CRAN (R 4.4.1)
RCurl 1.98-1.16 2024-07-11 [1] CRAN (R 4.4.1)
readr _ 2.1.5 2024-01-10 [1] CRAN (R 4.4.1)
remotes 2.5.0 2024-03-17 [1] CRAN (R 4.4.1)
reshape2 _ 1.4.4 2020-04-09 [1] CRAN (R 4.4.1)
reticulate 1.39.0 2024-09-05 [1] CRAN (R 4.4.1)
Rgraphviz 2.48.0 2024-05-01 [1] Bioconductor 3.19 (R 4.4.0)
rhandsontable _ 0.3.8 2021-05-27 [1] CRAN (R 4.4.1)
rlang 1.1.4 2024-06-04 [1] CRAN (R 4.4.1)
rmarkdown 2.28 2024-08-17 [1] CRAN (R 4.4.1)
robustbase 0.99-4-1 2024-09-27 [1] CRAN (R 4.4.1)
Rphenograph _ 0.99.1 2024-10-23 [1] Github (JinmiaoChenLab/Rphenograph@0298487)
RProtoBufLib 2.16.0 2024-05-01 [1] Bioconductor 3.19 (R 4.4.0)
rrcov 1.7-6 2024-08-19 [1] CRAN (R 4.4.1)
RSpectra 0.16-2 2024-07-18 [1] CRAN (R 4.4.1)
rstatix _ 0.7.2 2023-02-01 [1] CRAN (R 4.4.1)
rstudioapi 0.17.1 2024-10-22 [1] CRAN (R 4.4.1)
Rtsne _ 0.17 2023-12-07 [1] CRAN (R 4.4.1)
S4Vectors 0.42.1 2024-07-03 [1] Bioconductor 3.19 (R 4.4.1)
sass 0.4.9 2024-03-15 [1] CRAN (R 4.4.1)
scales 1.3.0 2023-11-28 [1] CRAN (R 4.4.1)
sessioninfo 1.2.2 2021-12-06 [1] CRAN (R 4.4.1)
shiny _ 1.9.1 2024-08-01 [1] CRAN (R 4.4.1)
shinyBS _ 0.61.1 2022-04-17 [1] CRAN (R 4.4.1)
shinycssloaders _ 1.1.0 2024-07-30 [1] CRAN (R 4.4.1)
shinydashboard _ 0.7.2 2021-09-30 [1] CRAN (R 4.4.1)
shinyjs _ 2.1.0 2021-12-23 [1] CRAN (R 4.4.1)
shinythemes _ 1.2.0 2021-01-25 [1] CRAN (R 4.4.1)
shinyWidgets _ 0.8.7 2024-09-23 [1] CRAN (R 4.4.1)
snow 0.4-4 2021-10-27 [1] CRAN (R 4.4.0)
spade _ 1.10.4 2021-06-10 [1] Github (nolanlab/spade@356955f)
streamgraph _ 0.9.0 2024-10-23 [1] Github (hrbrmstr/streamgraph@76f7173)
stringi 1.8.4 2024-05-06 [1] CRAN (R 4.4.0)
stringr _ 1.5.1 2023-11-14 [1] CRAN (R 4.4.1)
tibble _ 3.2.1 2023-03-20 [1] CRAN (R 4.4.1)
tidyr _ 1.3.1 2024-01-24 [1] CRAN (R 4.4.1)
tidyselect 1.2.1 2024-03-11 [1] CRAN (R 4.4.1)
tidyverse _ 2.0.0 2023-02-22 [1] CRAN (R 4.4.1)
timechange 0.3.0 2024-01-18 [1] CRAN (R 4.4.1)
tweenr 2.0.3 2024-02-26 [1] CRAN (R 4.4.1)
tzdb 0.4.0 2023-05-12 [1] CRAN (R 4.4.1)
umap _ 0.2.10.0 2023-02-01 [1] CRAN (R 4.4.1)
urlchecker 1.0.1 2021-11-30 [1] CRAN (R 4.4.1)
usethis _ 3.0.0 2024-07-29 [1] CRAN (R 4.4.1)
utf8 1.2.4 2023-10-22 [1] CRAN (R 4.4.1)
vctrs 0.6.5 2023-12-01 [1] CRAN (R 4.4.1)
viridis 0.6.5 2024-01-29 [1] CRAN (R 4.4.1)
viridisLite 0.4.2 2023-05-02 [1] CRAN (R 4.4.1)
withr 3.0.1 2024-07-31 [1] CRAN (R 4.4.1)
xfun 0.48 2024-10-03 [1] CRAN (R 4.4.1)
XML 3.99-0.17 2024-06-25 [1] CRAN (R 4.4.1)
xtable 1.8-4 2019-04-21 [1] CRAN (R 4.4.1)
xts 0.14.1 2024-10-15 [1] CRAN (R 4.4.1)
yaml 2.3.10 2024-07-26 [1] CRAN (R 4.4.1)
zlibbioc 1.50.0 2024-05-01 [1] Bioconductor 3.19 (R 4.4.0)
zoo 1.8-12 2023-04-13 [1] CRAN (R 4.4.1)

[1] C:/R/R-4.4.1/library
