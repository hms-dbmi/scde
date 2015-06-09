# SCDE

[![Build Status](https://travis-ci.org/hms-dbmi/scde.svg?branch=master)](https://travis-ci.org/hms-dbmi/scde)

The `scde` package implements a set of statistical methods for analyzing single-cell RNA-seq data. `scde` fits individual error models for single-cell RNA-seq measurements. These models can then be used for assessment of differential expression between groups of cells, as well as other types of analysis. The `scde` package also contains the `pagoda` framework which applies pathway and gene set overdispersion analysis to identify aspects of transcriptional heterogeneity among single cells. 
  
The overall approach to the differential expression analysis is detailed in the following publication:  
["Bayesian approach to single-cell differential expression analysis" (Kharchenko PV, Silberstein L, Scadden DT, Nature Methods, doi:10.1038/nmeth.2967)](http://www.nature.com/nmeth/journal/v11/n7/abs/nmeth.2967.html)

The overall approach to pathawys and gene set overdispersion analysis is detailed in the following publication:
"Characterizing transcriptional heterogeneity through pathway and gene set overdispersion analysis" (Fan J, Salathia N, Liu R, Kaeser G, Yung Y, Herman J, Kaper F, Fan JB, Zhang K, Chun J, and Kharchenko PV) COMING SOON!

# Sample analyses and images

Single cell error models

Differential expression analysis

Pathway and gene set overdispersion analysis with interactive GUI

# Installation 

```
require(devtools)
devtools::install_github("hms-dbmi/scde")
```

# Tutorials

Please refer to the following vignettes to help you get started with using `scde`:
- [Single cell differential expression analysis](vignettes/diffexp.md)
- [Characterizing single cell transcriptional heterogeneity using pathway and gene set overdispersion analysis](vignettes/pagoda.md)

Additional tutorials are available on the [Kharchenko Lab website](http://pklab.med.harvard.edu/scde/index.html). 
