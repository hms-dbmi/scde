Creating custom pathway annotations or gene sets
================================================

In this vignette, we show you how to create and use your own custom pathway annotations or gene sets with pagoda.

GO annotations
==============

``` r
# Use the org.Hs.eg.db package for GO annotations
library(org.Hs.eg.db)
# Translate gene names to ids
ids <- unlist(lapply(mget(rownames(cd), org.Hs.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
# Reverse map
rids <- names(ids)
names(rids) <- ids
# Convert ids per GO category to gene names
go.env <- eapply(org.Hs.egGO2ALLEGS, function(x) as.character(na.omit(rids[x])))
go.env <- clean.gos(go.env) # Remove GOs with too few or too many genes
go.env <- list2env(go.env)  # Convert to an environment

# Test
class(go.env)
```

    ## [1] "environment"

``` r
head(ls(go.env)) # Look at gene set names
```

    ## [1] "GO:0000002" "GO:0000003" "GO:0000012" "GO:0000014" "GO:0000018"
    ## [6] "GO:0000022"

``` r
head(get(ls(go.env)[1], go.env)) # Look at one gene set
```

    ## [1] "SLC25A4" "DNA2"    "TYMP"    "LIG3"    "LIG3"    "MEF2A"

BioMart
=======

Alternatively, we can use Ensembl's BioMart service to get the GO annotations.

``` r
library(biomaRt)
library(GO.db)

# Initialize the connection to the Ensembl BioMart Service
# Available datasets can be listed with 
# listDatasets(useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org"))
# Use mmusculus_gene_ensembl for mouse
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")

# Constructs a dataframe with two columns: hgnc_symbol and go_id
# If rownames are Ensembl IDs, use ensembl_gene_id as filter value
go <- getBM(attributes = c("hgnc_symbol", "go_id"), filters = "hgnc_symbol", values = rownames(cd), mart = ensembl)

# Use the GO.db library to add a column with the GO-term to the dataframe
go$term <- Term(go$go_id)

# Create a named list of character vectors out of the df
s = split(go$hgnc_symbol, paste(go$go_id,go$term))

# Saves the list as a R environment
go.env <- list2env(s)

# Test
class(go.env)
```

    ## [1] "environment"

``` r
head(ls(go.env)) # Look at gene set names
```

    ## [1] " NA"                                                 
    ## [2] "GO:0000002 mitochondrial genome maintenance"         
    ## [3] "GO:0000003 reproduction"                             
    ## [4] "GO:0000009 alpha-1,6-mannosyltransferase activity"   
    ## [5] "GO:0000010 trans-hexaprenyltranstransferase activity"
    ## [6] "GO:0000012 single strand break repair"

``` r
head(get(ls(go.env)[1], go.env)) # Look at one gene set
```

    ## [1] "HLA-DOB" "CLDN6"   "YPEL5"   "DYNLRB1" "LRRC41"  "RSPH3"

From GMT
========

The GMT file format is a tab delimited file format that describes gene sets. GMT files for Broad's MSigDB and other gene sets can be downloaded from the [Broad Website](http://www.broadinstitute.org/gsea/downloads.jsp).

``` r
## read in Broad gmt format
library(GSA)
filename <- 'https://raw.githubusercontent.com/JEFworks/genesets/master/msigdb.v5.0.symbols.gmt'
gs <- GSA.read.gmt(filename)

## number of gene sets
n <- length(gs$geneset.names)

## create environment
env <- new.env(parent=globalenv())
invisible(lapply(1:n,function(i) {
  genes <- as.character(unlist(gs$genesets[i]))
  name <- as.character(gs$geneset.names[i])
  assign(name, genes, envir = env)
}))

go.env <- env

# Test
class(go.env)
```

    ## [1] "environment"

``` r
head(ls(go.env)) # Look at gene set names
```

    ## [1] "3_5_CYCLIC_NUCLEOTIDE_PHOSPHODIESTERASE_ACTIVITY"
    ## [2] "3_5_EXONUCLEASE_ACTIVITY"                        
    ## [3] "AAACCAC,MIR-140"                                 
    ## [4] "AAAGACA,MIR-511"                                 
    ## [5] "AAAGGAT,MIR-501"                                 
    ## [6] "AAAGGGA,MIR-204,MIR-211"

``` r
head(get(ls(go.env)[1], go.env)) # Look at one gene set
```

    ## [1] "PDE3B"  "PDE4D"  "PDE3A"  "PDE10A" "PDE4C"  "PDE7B"
