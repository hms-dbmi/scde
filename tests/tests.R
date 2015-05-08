# tests for travis.ci
library(scde)

######
# Basic diff exp and batch correction tests
######

# load example dataset
data(es.mef.small)
# factor determining cell types
sg <- factor(gsub("(MEF|ESC).*", "\\1", colnames(es.mef.small)), levels=c("ESC", "MEF"))
# the group factor should be named accordingly
names(sg) <- colnames(es.mef.small)
table(sg)

# clean up the dataset
cd <- es.mef.small
# omit genes that are never detected
cd <- cd[rowSums(cd)>0, ]
# omit cells with very poor coverage
cd <- cd[, colSums(cd)>1e4]

# calculate models
# takes too long to run on travis...
# o.ifm <- scde.error.models(counts=cd, groups=sg, n.cores=1, threshold.segmentation=T, save.crossfit.plots=F, save.model.plots=F, verbose=1)
# devtools::use_data(o.ifm)  # save for later since this step takes a long time
data(o.ifm)

# filter out cells that don't show positive correlation with
# the expected expression magnitudes (very poor fits)
valid.cells <- o.ifm$corr.a > 0
table(valid.cells)
o.ifm <- o.ifm[valid.cells, ]

# estimate gene expression prior
o.prior <- scde.expression.prior(models=o.ifm, counts=cd, length.out=400, show.plot=F)

# define two groups of cells
groups <- factor(gsub("(MEF|ESC).*", "\\1", rownames(o.ifm)), levels=c("ESC", "MEF"))
names(groups) <- row.names(o.ifm)
# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups=groups, n.randomizations=100, n.cores=1, verbose=1)
# top upregulated genes (tail would show top downregulated ones)

scde.test.gene.expression.difference("Tdh", models=o.ifm, counts=cd, prior=o.prior)

batch <- as.factor(ifelse(rbinom(nrow(o.ifm), 1, 0.5)==1, "batch1", "batch2"))
# check the interaction between batches and cell types (shouldn't be any)
table(groups, batch)
# test the Tdh gene again
scde.test.gene.expression.difference("Tdh", models=o.ifm, counts=cd, prior=o.prior, batch=batch)

# test for all of the genes
ediff.batch <- scde.expression.difference(o.ifm, cd, o.prior, groups=groups, batch=batch, n.randomizations=100, n.cores=1, return.posteriors=T, verbose=1)
