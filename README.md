# cytoChain
This first prototype of cytoChain, a Shiny web app to handle high-dimensional cytofluorometry datasets, is currently deployed at the shinyapps.io servers and it is available at the following link:  https://abbatidanilo.shinyapps.io/cytoChain/
In case you are interested in cytoChain for your cytometric analisys drop a mail with a description of your goals, expriments and aims to abbati.danilo@hsr.it. 
cytoChain has been deployed by Danilo Abbati at the Experimental Hematology department of the San Raffaele Hospital - Milan (Italy). 

## Get the package
Install this package using the devtools library.

```
devtools::install_github("[abbatidanilo/cytoChain]")
```

## Operation 
Load your compensated samples, namely your *.FCS* (standard 3.0) files with the *File Input* control slot in the `Pre-clustering workflow` or in the `Metadata and assays` main side-panel and go throught step by step, all the tabs you need. It is possible to set all the relevant parameters offered in the various side-panels and you may also edit some meta-data related to your samples in the `Metadata and assays` main side-panel. The various analysis are offered within the `Clustering workflow` main side-panel and the related sub-choices.

[Manual available](//github.com/abbatidanilo/cytoChain/tree/main/www). The application will be keep on updating with new features

## Citation
If you use cytoChain, please use the following citation: "Flow cytometry data mining by cytoChain identifies determinants of exhaustion and stemness in TCR-engineered T cells" by **Francesco Manfredi**, **Danilo Abbati** et al. (doi: 10.1002/eji.202049103). This paper was published in the European Journal of Immunology in August 2021.

