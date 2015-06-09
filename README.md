# SCDE

[![Build Status](https://travis-ci.org/hms-dbmi/scde.svg?branch=master)](https://travis-ci.org/hms-dbmi/scde)

The `scde` package implements a set of statistical methods for analyzing single-cell RNA-seq data. `scde` fits individual error models for single-cell RNA-seq measurements. These models can then be used for assessment of differential expression between groups of cells, as well as other types of analysis. The `scde` package also contains the `pagoda` framework which applies pathway and gene set overdispersion analysis to identify aspects of transcriptional heterogeneity among single cells. 
  
The overall approach to the differential expression analysis is detailed in the following publication:  
["Bayesian approach to single-cell differential expression analysis" (Kharchenko PV, Silberstein L, Scadden DT, Nature Methods, doi:10.1038/nmeth.2967)](http://www.nature.com/nmeth/journal/v11/n7/abs/nmeth.2967.html)

The overall approach to pathawys and gene set overdispersion analysis is detailed in the following publication:
"Characterizing transcriptional heterogeneity through pathway and gene set overdispersion analysis" (Fan J, Salathia N, Liu R, Kaeser G, Yung Y, Herman J, Kaper F, Fan JB, Zhang K, Chun J, and Kharchenko PV) COMING SOON!

# Sample analyses and images

## Single cell error modeling
<table>
  <tr>
    <td width=600px>
      <img src="vignettes/figures/pagoda-cell.model.fits-0.png" width="600px">
    </td>
    <td>
      <code>scde</code> fits individual error models for single cells using counts derived from single-cell RNA-seq data to estimate drop-out and amplification biases on gene expression magnitude.
    </td>
  </tr>
</table>

## Differential expression analysis
<table>
  <tr>
    <td width=300px>
      <img src="vignettes/figures/scde-diffexp3-1.png" width="300">
    </td>
    <td width=300px>
      <pre>
                   lb      mle        ub       ce        Z       cZ
Dppa5a        8.075160 9.965929 11.541570 8.075160 7.160813 5.968921
Pou5f1        5.357179 7.208557  9.178109 5.357179 7.160333 5.968921
Gm13242       5.672307 7.681250  9.768974 5.672307 7.159987 5.968921
Tdh           5.829872 8.075160 10.281057 5.829872 7.159599 5.968921
Ift46         5.435961 7.366121  9.217500 5.435961 7.150271 5.968921
4930509G22Rik 5.435961 7.484295  9.808365 5.435961 7.115804 5.957784 </pre>
      <br>
      <code>scde</code> compares groups of single cells and tests for differential expression, taking into account variability in the single cell RNA-seq data due to drop-out and amplification biases in order to identify more robustly differentially expressed genes. 
    </td>
  </tr>
</table>

## Pathway and gene set overdispersion analysis with GUI
<table>
  <tr>
    <td width=600px>
      <img src="vignettes/figures/pagoda-Screen_Shot_2015-06-07_at_4.53.46_PM.png" width="600"> 
    </td>
    <td>
      <code>scde</code> contains <code>pagoda</code> routines that characterize aspects of transcriptional heterogeneity in populations of single cells using pre-defined gene sets as well as 'de novo' gene sets derived from the data. Significant aspects are used to cluster cells into subpopulations. A graphical user interface can be deployed to interactively explore results. 
    <td>
  </tr>
</table>
    
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
