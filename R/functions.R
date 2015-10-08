##' Single-cell Differential Expression (with Pathway And Gene set Overdispersion Analysis)
##'
##' The scde package implements a set of statistical methods for analyzing single-cell RNA-seq data.
##' scde fits individual error models for single-cell RNA-seq measurements. These models can then be used for
##' assessment of differential expression between groups of cells, as well as other types of analysis.
##' The scde package also contains the pagoda framework which applies pathway and gene set overdispersion analysis
##' to identify and characterize putative cell subpopulations based on transcriptional signatures.
##' See vignette("diffexp") for a brief tutorial on differential expression analysis.
##' See vignette("pagoda") for a brief tutorial on pathway and gene set overdispersion analysis to identify and characterize cell subpopulations.
##' More extensive tutorials are available at \url{http://pklab.med.harvard.edu/scde/index.html}.
##'  (test)
##' @name scde
##' @docType package
##' @author Peter Kharchenko \email{Peter_Kharchenko@@hms.harvard.edu}
##' @author Jean Fan \email{jeanfan@@fas.harvard.edu}
NULL

################################# Sample data

##' Sample data
##'
##' A subset of Saiful et al. 2011 dataset containing first 20 ES and 20 MEF cells.
##'
##' @name es.mef.small
##' @docType data
##' @references \url{http://www.ncbi.nlm.nih.gov/pubmed/21543516}
##' @export
NULL

##' Sample data
##'
##' Single cell data from Pollen et al. 2014 dataset.
##'
##' @name pollen
##' @docType data
##' @references \url{www.ncbi.nlm.nih.gov/pubmed/25086649}
##' @export
NULL

##' Sample error model
##'
##' SCDE error model generated from a subset of Saiful et al. 2011 dataset containing first 20 ES and 20 MEF cells.
##'
##' @name o.ifm
##' @docType data
##' @references \url{http://www.ncbi.nlm.nih.gov/pubmed/21543516}
##' @export
NULL

##' Sample error model
##'
##' SCDE error model generated from the Pollen et al. 2014 dataset.
##'
##' @name knn
##' @docType data
##' @references \url{www.ncbi.nlm.nih.gov/pubmed/25086649}
##' @export
NULL

# Internal model data
#
# Numerically-derived correction for NB->chi squared approximation stored as an local regression model
#
# @name scde.edff

################################# SCDE Methods

##' Fit single-cell error/regression models
##'
##' Fit error models given a set of single-cell data (counts) and an optional grouping factor (groups). The cells (within each group) are first cross-compared to determine a subset of genes showing consistent expression. The set of genes is then used to fit a mixture model (Poisson-NB mixture, with expression-dependent concomitant).
##'
##' Note: the default implementation has been changed to use linear-scale fit with expression-dependent NB size (overdispersion) fit. This represents an interative improvement on the originally published model. Use linear.fit=F to revert back to the original fitting procedure.
##'
##' @param counts read count matrix. The rows correspond to genes (should be named), columns correspond to individual cells. The matrix should contain integer counts
##' @param groups an optional factor describing grouping of different cells. If provided, the cross-fits and the expected expression magnitudes will be determined separately within each group. The factor should have the same length as ncol(counts).
##' @param min.nonfailed minimal number of non-failed observations required for a gene to be used in the final model fitting
##' @param threshold.segmentation use a fast threshold-based segmentation during cross-fit (default: TRUE)
##' @param min.count.threshold the number of reads to use to guess which genes may have "failed" to be detected in a given measurement during cross-cell comparison (default: 4)
##' @param zero.count.threshold threshold to guess the initial value (failed/non-failed) during error model fitting procedure (defaults to the min.count.threshold value)
##' @param zero.lambda the rate of the Poisson (failure) component (default: 0.1)
##' @param save.crossfit.plots whether png files showing cross-fit segmentations should be written out (default: FALSE)
##' @param save.model.plots whether pdf files showing model fits should be written out (default = TRUE)
##' @param n.cores number of cores to use
##' @param min.size.entries minimum number of genes to use when determining expected expression magnitude during model fitting
##' @param max.pairs maximum number of cross-fit comparisons that should be performed per group (default: 5000)
##' @param min.pairs.per.cell minimum number of pairs that each cell should be cross-compared with
##' @param verbose 1 for increased output
##' @param linear.fit Boolean of whether to use a linear fit in the regression (default: TRUE).
##' @param local.theta.fit Boolean of whether to fit the overdispersion parameter theta, ie. the negative binomial size parameter, based on local regression (default: set to be equal to the linear.fit parameter)
##' @param theta.fit.range Range of valid values for the overdispersion parameter theta, ie. the negative binomial size parameter (default: c(1e-2, 1e2))
##'
##' @return a model matrix, with rows corresponding to different cells, and columns representing different parameters of the determined models
##'
##' @useDynLib scde
##'
##' @examples
##' \donttest{
##' data(es.mef.small)
##' cd <- es.mef.small
##' cd <- cd[rowSums(cd) > 0, ]
##' cd <- cd[, colSums(cd) > 1e4]
##' sg <- factor(gsub("(MEF|ESC).*", "\\1", colnames(cd)), levels = c("ESC", "MEF"))
##' names(sg) <- colnames(cd)
##' o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = 10, threshold.segmentation = TRUE)
##' }
##'
##' @export
scde.error.models <- function(counts, groups = NULL, min.nonfailed = 3, threshold.segmentation = TRUE, min.count.threshold = 4, zero.count.threshold = min.count.threshold, zero.lambda = 0.1, save.crossfit.plots = FALSE, save.model.plots = TRUE, n.cores = 12, min.size.entries = 2e3, max.pairs = 5000, min.pairs.per.cell = 10, verbose = 0, linear.fit = TRUE, local.theta.fit = linear.fit, theta.fit.range = c(1e-2, 1e2)) {
    # default same group
    if(is.null(groups)) {
        groups <- as.factor(rep("cell", ncol(counts)))
    }
    # check for integer counts
    if(any(!unlist(lapply(counts,is.integer)))) {
      stop("Some of the supplied counts are not integer values (or stored as non-integer types). Aborting!\nThe method is designed to work on read counts - do not pass normalized read counts (e.g. FPKM values). If matrix contains read counts, but they are stored as numeric values, use counts<-apply(counts,2,function(x) {storage.mode(x) <- 'integer'; x}) to recast.");
    }

    # crossfit
    if(verbose) {
        cat("cross-fitting cells.\n")
    }
    cfm <- calculate.crossfit.models(counts, groups, n.cores = n.cores, threshold.segmentation = threshold.segmentation, min.count.threshold = min.count.threshold, zero.lambda = zero.lambda, max.pairs = max.pairs, save.plots = save.crossfit.plots, min.pairs.per.cell = min.pairs.per.cell, verbose = verbose)
    # error model for each cell
    if(verbose) {
        cat("building individual error models.\n")
    }
    ifm <- calculate.individual.models(counts, groups, cfm, min.nonfailed = min.nonfailed, zero.count.threshold = zero.count.threshold, n.cores = n.cores, save.plots = save.model.plots, linear.fit = linear.fit, return.compressed.models = TRUE, verbose = verbose, min.size.entries = min.size.entries, local.theta.fit = local.theta.fit, theta.fit.range = theta.fit.range)
    rm(cfm)
    gc()
    return(ifm)
}


##' Estimate prior distribution for gene expression magnitudes
##'
##' Use existing count data to determine a prior distribution of genes in the dataset
##'
##' @param models models determined by \code{\link{scde.error.models}}
##' @param counts count matrix
##' @param length.out number of points (resolution) of the expression magnitude grid (default: 400). Note: larger numbers will linearly increase memory/CPU demands.
##' @param show.plot show the estimate posterior
##' @param pseudo.count pseudo-count value to use (default 1)
##' @param bw smoothing bandwidth to use in estimating the prior (default: 0.1)
##' @param max.quantile determine the maximum expression magnitude based on a quantile (default : 0.999)
##' @param max.value alternatively, specify the exact maximum expression magnitude value
##'
##' @return a structure describing expression magnitude grid ($x, on log10 scale) and prior ($y)
##'
##' @examples
##' data(es.mef.small)
##' cd <- es.mef.small
##' cd <- cd[rowSums(cd) > 0, ]
##' cd <- cd[, colSums(cd) > 1e4]
##' data(o.ifm)  # Load precomputed model. Use ?scde.error.models to see how o.ifm was generated
##' o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
##'
##' @export
scde.expression.prior <- function(models, counts, length.out = 400, show.plot = FALSE, pseudo.count = 1, bw = 0.1, max.quantile = 1-1e-3, max.value = NULL) {
    fpkm <- scde.expression.magnitude(models, counts)
    fail <- scde.failure.probability(models, counts = counts)
    fpkm <- log10(exp(as.matrix(fpkm))+1)
    wts <- as.numeric(as.matrix(1-fail[, colnames(fpkm)]))
    wts <- wts/sum(wts)

    # fit density on a mirror image
    if(is.null(max.value)) {
        x <- as.numeric(fpkm)
        max.value <- as.numeric(quantile(x[x<Inf], p = max.quantile))
    }
    md <- density(c(-1*as.numeric(fpkm), as.numeric(fpkm)), bw = bw, weights = c(wts/2, wts/2), n = 2*length.out+1, from = -1*max.value, to = max.value)

    gep <- data.frame(x = md$x[-seq_len(length.out)], y = md$y[-seq_len(length.out)])
    gep$y[is.na(gep$y)] <- 0
    gep$y <- gep$y+pseudo.count/nrow(fpkm) # pseudo-count
    gep$y <- gep$y/sum(gep$y)
    if(show.plot) {
        par(mfrow = c(1, 1), mar = c(3.5, 3.5, 3.5, 0.5), mgp = c(2.0, 0.65, 0), cex = 0.9)
        plot(gep$x, gep$y, col = 4, panel.first = abline(h = 0, lty = 2), type = 'l', xlab = "log10( signal+1 )", ylab = "probability density", main = "signal prior")
    }
    gep$lp <- log(gep$y)

    # grid weighting (for normalization)
    gep$grid.weight <- diff(10^c(gep$x[1], gep$x+c(diff(gep$x)/2, 0))-1)

    return(gep)
    plot(x)
}


##' Test for expression differences between two sets of cells
##'
##' Use the individual cell error models to test for differential expression between two groups of cells.
##'
##' @param models models determined by \code{\link{scde.error.models}}
##' @param counts read count matrix
##' @param prior gene expression prior as determined by \code{\link{scde.expression.prior}}
##' @param groups a factor determining the two groups of cells being compared. The factor entries should correspond to the rows of the model matrix. The factor should have two levels. NAs are allowed (cells will be omitted from comparison).
##' @param batch a factor (corresponding to rows of the model matrix) specifying batch assignment of each cell, to perform batch correction
##' @param n.randomizations number of bootstrap randomizations to be performed
##' @param n.cores number of cores to utilize
##' @param batch.models (optional) separate models for the batch data (if generated using batch-specific group argument). Normally the same models are used.
##' @param return.posteriors whether joint posterior matrices should be returned
##' @param verbose integer verbose level (1 for verbose)
##'
##' @return \subsection{default}{
##' a data frame with the following fields:
##' \itemize{
##' \item{lb, mle, ub} {lower bound, maximum likelihood estimate, and upper bound of the 95% confidence interval for the expression fold change on log2 scale.}
##' \item{ce} { conservative estimate of expression-fold change (equals to the min(abs(c(lb, ub))), or 0 if the CI crosses the 0}
##' \item{Z} { uncorrected Z-score of expression difference}
##' \item{cZ} {expression difference Z-score corrected for multiple hypothesis testing using Holm procedure}
##' }
##'  If batch correction has been performed (\code{batch} has been supplied), analogous data frames are returned in slots \code{$batch.adjusted} for batch-corrected results, and \code{$batch.effect} for the differences explained by batch effects alone.
##' }}
##' \subsection{return.posteriors = TRUE}{
##' A list is returned, with the default results data frame given in the \code{$results} slot.
##' \code{difference.posterior} returns a matrix of estimated expression difference posteriors (rows - genes, columns correspond to different magnitudes of fold-change - log2 values are given in the column names)
##' \code{joint.posteriors} a list of two joint posterior matrices (rows - genes, columns correspond to the expression levels, given by prior$x grid)
##' }
##'
##' @examples
##' \donttest{
##' data(es.mef.small)
##' cd <- es.mef.small
##' cd <- cd[rowSums(cd) > 0, ]
##' cd <- cd[, colSums(cd) > 1e4]
##' sg <- factor(gsub("(MEF|ESC).*", "\\1", colnames(cd)), levels = c("ESC", "MEF"))
##' names(sg) <- colnames(cd)
##' o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = 10, threshold.segmentation = TRUE)
##' o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
##' # make sure groups corresponds to the models (o.ifm)
##' groups <- factor(gsub("(MEF|ESC).*", "\\1", rownames(o.ifm)), levels = c("ESC", "MEF"))
##' names(groups) <- row.names(o.ifm)
##' ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups, n.randomizations = 100, n.cores = n.cores, verbose = 1)
##' }
##'
##' @export
scde.expression.difference <- function(models, counts, prior, groups = NULL, batch = NULL, n.randomizations = 150, n.cores = 10, batch.models = models, return.posteriors = FALSE, verbose = 0) {
    if(!all(rownames(models) %in% colnames(counts))) {
        stop("ERROR: provided count data does not cover all of the cells specified in the model matrix")
    }

    ci <- match(rownames(models), colnames(counts))
    counts <- as.matrix(counts[, ci])

    if(is.null(groups)) { # recover groups from models
        groups <- as.factor(attr(models, "groups"))
        if(is.null(groups)) stop("ERROR: groups factor is not provided, and models structure is lacking groups attribute")
        names(groups) <- rownames(models)
    }
    if(length(levels(groups)) != 2) {
        stop(paste("ERROR: wrong number of levels in the grouping factor (", paste(levels(groups), collapse = " "), "), but must be two.", sep = ""))
    }

    correct.batch <- FALSE
    if(!is.null(batch)) {
        if(length(levels(batch)) > 1) {
            correct.batch <- TRUE
        } else {
            if(verbose) {
                cat("WARNING: only one batch level detected. Nothing to correct for.")
            }
        }
    }

    # batch control
    if(correct.batch) {
        batch <- as.factor(batch)
        # check batch-group interactions
        bgti <- table(groups, batch)
        bgti.ft <- fisher.test(bgti)
        if(verbose) {
            cat("controlling for batch effects. interaction:\n")
            print(bgti)
        }
        #if(any(bgti == 0)) {
        #  cat("ERROR: cannot control for batch effect, as some batches are found only in one group:\n")
        #  print(bgti)
        #}
        if(bgti.ft$p.value < 1e-3) {
            cat("WARNING: strong interaction between groups and batches! Correction may be ineffective:\n")
            print(bgti.ft)
        }

        # calculate batch posterior
        if(verbose) {
            cat("calculating batch posteriors\n")
        }
        batch.jpl <- tapply(seq_len(nrow(models)), groups, function(ii) {
            scde.posteriors(models = batch.models, counts = counts, prior = prior, batch = batch, composition = table(batch[ii]), n.cores = n.cores, n.randomizations = n.randomizations, return.individual.posteriors = FALSE)
        })
        if(verbose) {
            cat("calculating batch differences\n")
        }
        batch.bdiffp <- calculate.ratio.posterior(batch.jpl[[1]], batch.jpl[[2]], prior, n.cores = n.cores)
        batch.bdiffp.rep <- quick.distribution.summary(batch.bdiffp)
    } else {
        if(verbose) {
            cat("comparing groups:\n")
            print(table(as.character(groups)))
        }
    }


    # fit joint posteriors for each group
    jpl <- tapply(seq_len(nrow(models)), groups, function(ii) {
        scde.posteriors(models = models[ii, , drop = FALSE], counts = counts[, ii, drop = FALSE], prior = prior, n.cores = n.cores, n.randomizations = n.randomizations)
    })
    if(verbose) {
        cat("calculating difference posterior\n")
    }
    # calculate difference posterior
    bdiffp <- calculate.ratio.posterior(jpl[[1]], jpl[[2]], prior, n.cores = n.cores)

    if(verbose) {
        cat("summarizing differences\n")
    }
    bdiffp.rep <- quick.distribution.summary(bdiffp)

    if(correct.batch) {
        if(verbose) {
            cat("adjusting for batch effects\n")
        }
        # adjust for batch effects
        a.bdiffp <- calculate.ratio.posterior(bdiffp, batch.bdiffp, prior = data.frame(x = as.numeric(colnames(bdiffp)), y = rep(1/ncol(bdiffp), ncol(bdiffp))), skip.prior.adjustment = TRUE, n.cores = n.cores)
        a.bdiffp.rep <- quick.distribution.summary(a.bdiffp)

        # return with batch correction info
        if(return.posteriors) {
            return(list(batch.adjusted = a.bdiffp.rep, results = bdiffp.rep, batch.effect = batch.bdiffp.rep, difference.posterior = bdiffp, batch.adjusted.difference.posterior = a.bdiffp, joint.posteriors = jpl))
        } else {
            return(list(batch.adjusted = a.bdiffp.rep, results = bdiffp.rep, batch.effect = batch.bdiffp.rep))
        }
    } else {
        # no batch correction return
        if(return.posteriors) {
            return(list(results = bdiffp.rep, difference.posterior = bdiffp, joint.posteriors = jpl))
        } else {
            return(bdiffp.rep)
        }
    }
}


##' View differential expression results in a browser
##'
##' Launches a browser app that shows the differential expression results, allowing to sort, filter, etc.
##' The arguments generally correspond to the \code{scde.expression.difference()} call, except that the results of that call are also passed here. Requires \code{Rook} and \code{rjson} packages to be installed.
##'
##' @param results result object returned by \code{scde.expression.difference()}. Note to browse group posterior levels, use \code{return.posteriors = TRUE} in the \code{scde.expression.difference()} call.
##' @param models model matrix
##' @param counts count matrix
##' @param prior prior
##' @param groups group information
##' @param batch batch information
##' @param geneLookupURL The URL that will be used to construct links to view more information on gene names. By default (if can't guess the organism) the links will forward to ENSEMBL site search, using \code{geneLookupURL = "http://useast.ensembl.org/Multi/Search/Results?q = {0}"}. The "{0}" in the end will be substituted with the gene name. For instance, to link to GeneCards, use \code{"http://www.genecards.org/cgi-bin/carddisp.pl?gene = {0}"}.
##' @param server optional previously returned instance of the server, if want to reuse it.
##' @param name app name (needs to be altered only if adding more than one app to the server using \code{server} parameter)
##' @param port Interactive browser port
##'
##' @return server instance, on which $stop() function can be called to kill the process.
##'
##' @examples
##' \donttest{
##' data(es.mef.small)
##' cd <- es.mef.small
##' cd <- cd[rowSums(cd) > 0, ]
##' cd <- cd[, colSums(cd) > 1e4]
##' sg <- factor(gsub("(MEF|ESC).*", "\\1", colnames(cd)), levels = c("ESC", "MEF"))
##' names(sg) <- colnames(cd)
##' o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = 10, threshold.segmentation = TRUE)
##' o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
##' # make sure groups corresponds to the models (o.ifm)
##' groups <- factor(gsub("(MEF|ESC).*", "\\1", rownames(o.ifm)), levels = c("ESC", "MEF"))
##' names(groups) <- row.names(o.ifm)
##' ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups, n.randomizations = 100, n.cores = 10, verbose = 1)
##' scde.browse.diffexp(ediff, o.ifm, cd, o.prior, groups = groups, geneLookupURL="http://www.informatics.jax.org/searchtool/Search.do?query={0}")  # creates browser
##' }
##'
##' @export
scde.browse.diffexp <- function(results, models, counts, prior, groups = NULL, batch = NULL, geneLookupURL = NULL, server = NULL, name = "scde", port = NULL) {
    #require(Rook)
    #require(rjson)
    if(is.null(server)) { server <- get.scde.server(port) }
    sa <- ViewDiff$new(results, models, counts, prior, groups = groups, batch = batch, geneLookupURL = geneLookupURL)
    server$add(app = sa, name = name)
    browseURL(paste(server$full_url(name), "index.html", sep = "/"))
    return(server)
}


##' View PAGODA application
##'
##' Installs a given pagoda app (or any other rook app) into a server, optionally
##' making a call to show it in the browser.
##'
##' @param app pagoda app (output of make.pagoda.app()) or another rook app
##' @param name URL path name for this app
##' @param browse whether a call should be made for browser to show the app
##' @param port optional port on which the server should be initiated
##' @param ip IP on which the server should listen (typically localhost)
##' @param server an (optional) Rook server instance (defaults to ___scde.server)
##'
##' @examples
##' \donttest{
##' app <- make.pagoda.app(tamr2, tam, varinfo, go.env, pwpca, clpca, col.cols=col.cols, cell.clustering=hc, title="NPCs")
##' # show app in the browser (port 1468)
##' show.app(app, "pollen", browse = TRUE, port=1468)
##' }
##'
##' @return Rook server instance
##'
##' @export
show.app <- function(app, name, browse = TRUE, port = NULL, ip = '127.0.0.1', server = NULL) {
    if(is.null(server)) { server <- get.scde.server(port) }
    server$add(app = app, name = name)
    if(browse) {
        browseURL(paste(server$full_url(name), "index.html", sep = "/"))
    }
    return(server)
}
# get SCDE server from saved session
get.scde.server <- function(port = NULL, ip = '127.0.0.1') {
    if(exists("___scde.server", envir = globalenv())) {
        server <- get("___scde.server", envir = globalenv())
    } else {
        require(Rook)
        server <- Rhttpd$new()
        assign("___scde.server", server, envir = globalenv())
        server$start(listen = ip, port = port)
    }
    return(server)
}


# calculate individual and joint posterior information
# models - all or a subset of models belonging to a particular group
#
##' Calculate joint expression magnitude posteriors across a set of cells
##'
##' Calculates expression magnitude posteriors for the individual cells, and then uses bootstrap resampling to calculate a joint expression posterior for all the specified cells. Alternatively during batch-effect correction procedure, the joint posterior can be calculated for a random composition of cells of different groups (see \code{batch} and \code{composition} parameters).
##'
##' @param models models models determined by \code{\link{scde.error.models}}
##' @param counts read count matrix
##' @param prior gene expression prior as determined by \code{\link{scde.expression.prior}}
##' @param n.randomizations number of bootstrap iterations to perform
##' @param batch a factor describing which batch group each cell (i.e. each row of \code{models} matrix) belongs to
##' @param composition a vector describing the batch composition of a group to be sampled
##' @param return.individual.posteriors whether expression posteriors of each cell should be returned
##' @param return.individual.posterior.modes whether modes of expression posteriors of each cell should be returned
##' @param ensemble.posterior Boolean of whether to calculate the ensemble posterior (sum of individual posteriors) instead of a joint (product) posterior. (default: FALSE)
##' @param n.cores number of cores to utilize
##'
##' @return \subsection{default}{ a posterior probability matrix, with rows corresponding to genes, and columns to expression levels (as defined by \code{prior$x})
##' }
##' \subsection{return.individual.posterior.modes}{ a list is returned, with the \code{$jp} slot giving the joint posterior matrix, as described above. The \code{$modes} slot gives a matrix of individual expression posterior mode values on log scale (rows - genes, columns -cells)}
##' \subsection{return.individual.posteriors}{ a list is returned, with the \code{$post} slot giving a list of individual posterior matrices, in a form analogous to the joint posterior matrix, but reported on log scale }
##'
##' @examples
##' data(es.mef.small)
##' cd <- es.mef.small
##' cd <- cd[rowSums(cd) > 0, ]
##' cd <- cd[, colSums(cd) > 1e4]
##' data(o.ifm)  # Load precomputed model. Use ?scde.error.models to see how o.ifm was generated
##' o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
##' # calculate joint posteriors
##' jp <- scde.posteriors(o.ifm, cd, o.prior, n.cores = 1)
##'
##' @export
scde.posteriors <- function(models, counts, prior, n.randomizations = 100, batch = NULL, composition = NULL, return.individual.posteriors = FALSE, return.individual.posterior.modes = FALSE, ensemble.posterior = FALSE, n.cores = 20) {
    if(!all(rownames(models) %in% colnames(counts))) { stop("ERROR: provided count data does not cover all of the cells specified in the model matrix") }
    if(!is.null(batch)) { # calculating batch-sampled posteriors instead of evenly sampled ones
        if(is.null(composition)) { stop("ERROR: group composition must be provided if the batch argument is passed") }
        batchil <- tapply(c(1:nrow(models))-1, batch, I)
    }
    # order counts according to the cells
    ci <- match(rownames(models), colnames(counts))
    counts <- as.matrix(counts[, ci, drop = FALSE])
    marginals <- 10^prior$x - 1
    marginals[marginals<0] <- 0
    marginals <- log(marginals)

    min.slope <- 1e-10
    if(any(models$corr.a<min.slope)) {
        cat("WARNING: the following cells have negatively-correlated or 0-slope fits: ", paste(rownames(models)[models$corr.a<min.slope], collapse = " "), ". Setting slopes to 1e-10.\n")
        models$corr.a[models$corr.a<min.slope] <- min.slope
    }

    postflag <- 0
    if(return.individual.posteriors) {
        postflag <- 2
        if(return.individual.posterior.modes) {
            postflag <- 3
        }
    } else if(return.individual.posterior.modes) {
        postflag <- 1
    }

    ensembleflag <- ifelse(ensemble.posterior, 1, 0)

    localthetaflag <- "corr.ltheta.b" %in% colnames(models)
    squarelogitconc <- "conc.a2" %in% colnames(models)

    # prepare matrix models
    mn <- c("conc.b", "conc.a", "fail.r", "corr.b", "corr.a", "corr.theta", "corr.ltheta.b", "corr.ltheta.t", "corr.ltheta.m", "corr.ltheta.s", "corr.ltheta.r", "conc.a2")
    mc <- match(c(mn), colnames(models))
    mm <- matrix(NA, nrow(models), length(mn))
    mm[, which(!is.na(mc))] <- as.matrix(models[, mc[!is.na(mc)], drop = FALSE])

    chunk <- function(x, n) split(x, sort(rank(x) %% n.cores))
    if(n.cores > 1 && nrow(counts) > n.cores) { # split by genes
        xl <- papply(chunk(seq_len(nrow(counts)), n.cores), function(ii) {
            ucl <- lapply(seq_len(ncol(counts)), function(i) as.vector(unique(counts[ii, i, drop = FALSE])))
            uci <- do.call(cbind, lapply(seq_len(ncol(counts)), function(i) match(counts[ii, i, drop = FALSE], ucl[[i]])-1))
            #x <- logBootPosterior(models, ucl, uci, marginals, n.randomizations, 1, postflag)
            if(!is.null(batch)) {
                x <- .Call("logBootBatchPosterior", mm, ucl, uci, marginals, batchil, composition, n.randomizations, ii[1], postflag, localthetaflag, squarelogitconc, PACKAGE = "scde")
            } else {
                x <- .Call("logBootPosterior", mm, ucl, uci, marginals, n.randomizations, ii[1], postflag, localthetaflag, squarelogitconc, ensembleflag, PACKAGE = "scde")
            }
        }, n.cores = n.cores)
        if(postflag == 0) {
            x <- do.call(rbind, xl)
        } else if(postflag == 1) {
            x <- list(jp = do.call(rbind, lapply(xl, function(d) d$jp)), modes = do.call(rbind, lapply(xl, function(d) d$modes)))
        } else if(postflag == 2) {
            x <- list(jp = do.call(rbind, lapply(xl, function(d) d$jp)), post = lapply(seq_along(xl[[1]]$post), function(pi) { do.call(rbind, lapply(xl, function(d) d$post[[pi]])) }))
        } else if(postflag == 3) {
            x <- list(jp = do.call(rbind, lapply(xl, function(d) d$jp)), modes = do.call(rbind, lapply(xl, function(d) d$modes)), post = lapply(seq_along(xl[[1]]$post), function(pi) { do.call(rbind, lapply(xl, function(d) d$post[[pi]])) }))
        }
        rm(xl)
        gc()
    } else {
        # unique count lists with matching indices
        ucl <- lapply(seq_len(ncol(counts)), function(i) as.vector(unique(counts[, i, drop = FALSE])))
        uci <- do.call(cbind, lapply(seq_len(ncol(counts)), function(i) match(counts[, i, drop = FALSE], ucl[[i]])-1))
        #x <- logBootPosterior(models, ucl, uci, marginals, n.randomizations, 1, postflag)
        if(!is.null(batch)) {
            x <- .Call("logBootBatchPosterior", mm, ucl, uci, marginals, batchil, composition, n.randomizations, 1, postflag, localthetaflag, squarelogitconc, PACKAGE = "scde")
        } else {
            x <- .Call("logBootPosterior", mm, ucl, uci, marginals, n.randomizations, 1, postflag, localthetaflag, squarelogitconc, ensembleflag, PACKAGE = "scde")
        }
    }
    if(postflag == 0) {
        rownames(x) <- rownames(counts)
        colnames(x) <- as.character(exp(marginals))
    } else if(postflag == 1) {
        rownames(x$jp) <- rownames(counts)
        colnames(x$jp) <- as.character(exp(marginals))
        rownames(x$modes) <- rownames(counts)
        colnames(x$modes) <- rownames(models)
    } else if(postflag == 2) {
        rownames(x$jp) <- rownames(counts)
        colnames(x$jp) <- as.character(exp(marginals))
        names(x$post) <- rownames(models)
        x$post <- lapply(x$post, function(d) {
            rownames(d) <- rownames(counts)
            colnames(d) <- as.character(exp(marginals))
            return(d)
        })
    } else if(postflag == 3) {
        rownames(x$jp) <- rownames(counts)
        colnames(x$jp) <- as.character(exp(marginals))
        rownames(x$modes) <- rownames(counts)
        colnames(x$modes) <- rownames(models)
        names(x$post) <- rownames(models)
        x$post <- lapply(x$post, function(d) {
            rownames(d) <- rownames(counts)
            colnames(d) <- as.character(exp(marginals))
            return(d)
        })
    }
    return(x)
}


# get estimates of expression magnitude for a given set of models
# models - entire model matrix, or a subset of cells (i.e. select rows) of the model matrix for which the estimates should be obtained
# counts - count data that covers the desired set of genes (rows) and all specified cells (columns)
# return - a matrix of log(FPM) estimates with genes as rows and cells  as columns (in the model matrix order).
##' Return scaled expression magnitude estimates
##'
##' Return point estimates of expression magnitudes of each gene across a set of cells, based on the regression slopes determined during the model fitting procedure.
##'
##' @param models models determined by \code{\link{scde.error.models}}
##' @param counts count matrix
##'
##' @return a matrix of expression magnitudes on a log scale (rows - genes, columns - cells)
##'
##' @examples
##' data(es.mef.small)
##' cd <- es.mef.small
##' cd <- cd[rowSums(cd) > 0, ]
##' cd <- cd[, colSums(cd) > 1e4]
##' data(o.ifm)  # Load precomputed model. Use ?scde.error.models to see how o.ifm was generated
##' # get expression magnitude estimates
##' lfpm <- scde.expression.magnitude(o.ifm, cd)
##'
##' @export
scde.expression.magnitude <- function(models, counts) {
    if(!all(rownames(models) %in% colnames(counts))) { stop("ERROR: provided count data does not cover all of the cells specified in the model matrix") }
    t((t(log(counts[, rownames(models), drop = FALSE]))-models$corr.b)/models$corr.a)
}


# calculate drop-out probability given either count data or magnitudes (log(FPM))
# magnitudes can either be a per-cell matrix or a single vector of values which will be evaluated for each cell
# returns a probability of a drop out event for every gene (rows) for every cell (columns)
##' Calculate drop-out probabilities given a set of counts or expression magnitudes
##'
##' Returns estimated drop-out probability for each cell (row of \code{models} matrix), given either an expression magnitude
##' @param models models determined by \code{\link{scde.error.models}}
##' @param magnitudes a vector (\code{length(counts) == nrows(models)}) or a matrix (columns correspond to cells) of expression magnitudes, given on a log scale
##' @param counts a vector (\code{length(counts) == nrows(models)}) or a matrix (columns correspond to cells) of read counts from which the expression magnitude should be estimated
##'
##' @return a vector or a matrix of drop-out probabilities
##'
##' @examples
##' data(es.mef.small)
##' cd <- es.mef.small
##' cd <- cd[rowSums(cd) > 0, ]
##' cd <- cd[, colSums(cd) > 1e4]
##' data(o.ifm)  # Load precomputed model. Use ?scde.error.models to see how o.ifm was generated
##' o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
##' # calculate probability of observing a drop out at a given set of magnitudes in different cells
##' mags <- c(1.0, 1.5, 2.0)
##' p <- scde.failure.probability(o.ifm, magnitudes = mags)
##' # calculate probability of observing the dropout at a magnitude corresponding to the
##' # number of reads actually observed in each cell
##' self.p <- scde.failure.probability(o.ifm, counts = cd)
##'
##' @export
scde.failure.probability <- function(models, magnitudes = NULL, counts = NULL) {
    if(is.null(magnitudes)) {
        if(!is.null(counts)) {
            magnitudes <- scde.expression.magnitude(models, counts)
        } else {
            stop("ERROR: either magnitudes or counts should be provided")
        }
    }
    if(is.matrix(magnitudes)) { # a different vector for every cell
        if(!all(rownames(models) %in% colnames(magnitudes))) { stop("ERROR: provided magnitude data does not cover all of the cells specified in the model matrix") }
        if("conc.a2" %in% names(models)) {
            x <- t(1/(exp(t(magnitudes)*models$conc.a +t(magnitudes^2)*models$conc.a2 + models$conc.b)+1))
        } else {
            x <- t(1/(exp(t(magnitudes)*models$conc.a + models$conc.b)+1))
        }
    } else { # a common vector of magnitudes for all cells
        if("conc.a2" %in% names(models)) {
            x <- t(1/(exp((models$conc.a %*% t(magnitudes)) + (models$conc.a2 %*% t(magnitudes^2)) + models$conc.b)+1))
        } else {
            x <- t(1/(exp((models$conc.a %*% t(magnitudes)) + models$conc.b)+1))
        }
    }
    x[is.nan(x)] <- 0
    colnames(x) <- rownames(models)
    x
}


##' Test differential expression and plot posteriors for a particular gene
##'
##' The function performs differential expression test and optionally plots posteriors for a specified gene.
##'
##' @param gene name of the gene to be tested
##' @param models models
##' @param counts read count matrix (must contain the row corresponding to the specified gene)
##' @param prior expression magnitude prior
##' @param groups a two-level factor specifying between which cells (rows of the models matrix) the comparison should be made
##' @param batch optional multi-level factor assigning the cells (rows of the model matrix) to different batches that should be controlled for (e.g. two or more biological replicates). The expression difference estimate will then take into account the likely difference between the two groups that is explained solely by their difference in batch composition. Not all batch configuration may be corrected this way.
##' @param batch.models optional set of models for batch comparison (typically the same as models, but can be more extensive, or recalculated within each batch)
##' @param n.randomizations number of bootstrap/sampling iterations that should be performed
##' @param show.plots whether the plots should be shown
##' @param return.details whether the posterior should be returned
##' @param verbose set to T for some status output
##' @param ratio.range optionally specifies the range of the log2 expression ratio plot
##' @param show.individual.posteriors whether the individual cell expression posteriors should be plotted
##' @param n.cores number of cores to use (default = 1)
##'
##' @return by default returns MLE of log2 expression difference, 95% CI (upper, lower bound), and a Z-score testing for expression difference. If return.details = TRUE, a list is returned containing the above structure, as well as the expression fold difference posterior itself.
##'
##' @examples
##' data(es.mef.small)
##' cd <- es.mef.small
##' cd <- cd[rowSums(cd) > 0, ]
##' cd <- cd[, colSums(cd) > 1e4]
##' data(o.ifm)  # Load precomputed model. Use ?scde.error.models to see how o.ifm was generated
##' o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
##' scde.test.gene.expression.difference("Tdh", models = o.ifm, counts = cd, prior = o.prior)
##'
##' @export
scde.test.gene.expression.difference <- function(gene, models, counts, prior, groups = NULL, batch = NULL, batch.models = models, n.randomizations = 1e3, show.plots = TRUE, return.details = FALSE, verbose = FALSE, ratio.range = NULL, show.individual.posteriors = TRUE, n.cores = 1) {
    if(!gene %in% rownames(counts)) {
        stop("ERROR: specified gene (", gene, ") is not found in the count data")
    }

    ci <- match(rownames(models), colnames(counts))
    counts <- as.matrix(counts[gene, ci, drop = FALSE])


    if(is.null(groups)) { # recover groups from models
        groups <- as.factor(attr(models, "groups"))
        if(is.null(groups)) stop("ERROR: groups factor is not provided, and models structure is lacking groups attribute")
        names(groups) <- rownames(models)
    }
    if(length(levels(groups)) != 2) {
        stop(paste("ERROR: wrong number of levels in the grouping factor (", paste(levels(groups), collapse = " "), "), but must be two.", sep = ""))
    }

    if(verbose) {
        cat("comparing gene ", gene, " between groups:\n")
        print(table(as.character(groups)))
    }

    # calculate joint posteriors
    jpl <- tapply(seq_len(nrow(models)), groups, function(ii) {
        scde.posteriors(models = models[ii, , drop = FALSE], counts = counts[, ii, drop = FALSE], prior = prior, n.cores = n.cores, n.randomizations = n.randomizations, return.individual.posteriors = TRUE)
    })

    bdiffp <- calculate.ratio.posterior(jpl[[1]]$jp, jpl[[2]]$jp, prior, n.cores = n.cores)

    bdiffp.rep <- quick.distribution.summary(bdiffp)

    nam1 <- levels(groups)[1]
    nam2 <- levels(groups)[2]

    # batch control
    correct.batch <- !is.null(batch) && length(levels(batch)) > 1
    if(correct.batch) {
        batch <- as.factor(batch)
        # check batch-group interactions
        bgti <- table(groups, batch)
        bgti.ft <- fisher.test(bgti)
        if(verbose) {
            cat("controlling for batch effects. interaction:\n")
        }
        if(any(bgti == 0)) {
            cat("ERROR: cannot control for batch effect, as some batches are found only in one group:\n")
            print(bgti)
        }
        if(bgti.ft$p.value<1e-3) {
            cat("WARNING: strong interaction between groups and batches! Correction may be ineffective:\n")
            print(bgti)
            print(bgti.ft)
        }
        # calculate batch posterior
        batch.jpl <- tapply(seq_len(nrow(models)), groups, function(ii) {
            scde.posteriors(models = batch.models, counts = counts, prior = prior, batch = batch, composition = table(batch[ii]), n.cores = n.cores, n.randomizations = n.randomizations, return.individual.posteriors = FALSE)
        })
        batch.bdiffp <- calculate.ratio.posterior(batch.jpl[[1]], batch.jpl[[2]], prior, n.cores = n.cores)
        a.bdiffp <- calculate.ratio.posterior(bdiffp, batch.bdiffp, prior = data.frame(x = as.numeric(colnames(bdiffp)), y = rep(1/ncol(bdiffp), ncol(bdiffp))), skip.prior.adjustment = TRUE)
        a.bdiffp.rep <- quick.distribution.summary(a.bdiffp)
    }


    if(show.plots) {
        # show each posterior
        layout(matrix(c(1:3), 3, 1, byrow = TRUE), heights = c(2, 1, 2), widths = c(1), FALSE)
        par(mar = c(2.5, 3.5, 2.5, 3.5), mgp = c(1.5, 0.65, 0), cex = 0.9)
        #par(mar = c(2.5, 3.5, 0.5, 3.5), mgp = c(1.5, 0.65, 0), cex = 0.9)

        pp <- exp(do.call(rbind, lapply(jpl[[1]]$post, as.numeric)))
        cols <- rainbow(nrow(pp), s = 0.8)
        plot(c(), c(), xlim = range(prior$x), ylim = range(c(0, pp)), xlab = "expression level", ylab = "individual posterior", main = nam1)
        if(show.individual.posteriors) {
            lapply(seq_len(nrow(pp)), function(i) lines(prior$x, pp[i, ], col = rgb(1, 0.5, 0, alpha = 0.25)))
        }
        #legend(x = ifelse(which.max(na.omit(pjpc)) > length(pjpc)/2, "topleft", "topright"), bty = "n", col = cols, legend = rownames(pp), lty = rep(1, nrow(pp)))
        if(correct.batch) {
            par(new = TRUE)
            plot(prior$x, batch.jpl[[1]][1, ], axes = FALSE, ylab = "", xlab = "", type = 'l', col = 8, lty = 1, lwd = 2)
        }
        pjpc <- jpl[[1]]$jp
        par(new = TRUE)
        jpr <- range(c(0, na.omit(pjpc)))
        plot(prior$x, pjpc, axes = FALSE, ylab = "", xlab = "", ylim = jpr, type = 'l', col = 1, lty = 1, lwd = 2)
        axis(4, pretty(jpr, 5), col = 1)
        mtext("joint posterior", side = 4, outer = FALSE, line = 2)


        # ratio plot
        if(is.null(ratio.range)) { ratio.range <- range(as.numeric(colnames(bdiffp))/log10(2)) }

        par(mar = c(2.5, 3.5, 0.5, 3.5), mgp = c(1.5, 0.65, 0), cex = 0.9)
        rv <- as.numeric(colnames(bdiffp))/log10(2)
        rp <- as.numeric(bdiffp[1, ])
        plot(rv, rp, xlab = "log2 expression ratio", ylab = "ratio posterior", type = 'l', lwd = ifelse(correct.batch, 1, 2), main = "", axes = FALSE, xlim = ratio.range, ylim = c(0, max(bdiffp)))
        axis(1, pretty(ratio.range, 5), col = 1)
        abline(v = 0, lty = 2, col = 8)
        if(correct.batch) { # with batch correction
            # show batch difference
            par(new = TRUE)
            plot(as.numeric(colnames(batch.bdiffp))/log10(2), as.numeric(batch.bdiffp[1, ]), xlab = "", ylab = "", type = 'l', lwd = 1, main = "", axes = FALSE, xlim = ratio.range, col = 8, ylim = c(0, max(batch.bdiffp)))
            # fill out the a.bdiffp confidence interval
            par(new = TRUE)
            rv <- as.numeric(colnames(a.bdiffp))/log10(2)
            rp <- as.numeric(a.bdiffp[1, ])
            plot(rv, rp, xlab = "", ylab = "", type = 'l', lwd = 2, main = "", axes = FALSE, xlim = ratio.range, col = 2, ylim = c(0, max(rp)))
            axis(2, pretty(c(0, max(a.bdiffp)), 2), col = 1)
            r.lb <- which.min(abs(rv-a.bdiffp.rep$lb))
            r.ub <- which.min(abs(rv-a.bdiffp.rep$ub))
            polygon(c(rv[r.lb], rv[r.lb:r.ub], rv[r.ub]), y = c(-10, rp[r.lb:r.ub], -10), col = rgb(1, 0, 0, alpha = 0.2), border = NA)
            abline(v = a.bdiffp.rep$mle, col = 2, lty = 2)
            abline(v = c(rv[r.ub], rv[r.lb]), col = 2, lty = 3)

            legend(x = ifelse(a.bdiffp.rep$mle > 0, "topleft", "topright"), legend = c(paste("MLE: ", round(a.bdiffp.rep$mle, 2), sep = ""), paste("95% CI: ", round(a.bdiffp.rep$lb, 2), " : ", round(a.bdiffp.rep$ub, 2), sep = ""), paste("Z = ", round(a.bdiffp.rep$Z, 2), sep = ""), paste("cZ = ", round(a.bdiffp.rep$cZ, 2), sep = "")), bty = "n")

        } else {  # without batch correction
            # fill out the bdiffp confidence interval
            axis(2, pretty(c(0, max(bdiffp)), 2), col = 1)

            r.lb <- which.min(abs(rv-bdiffp.rep$lb))
            r.ub <- which.min(abs(rv-bdiffp.rep$ub))
            polygon(c(rv[r.lb], rv[r.lb:r.ub], rv[r.ub]), y = c(-10, rp[r.lb:r.ub], -10), col = rgb(1, 0, 0, alpha = 0.2), border = NA)
            abline(v = bdiffp.rep$mle, col = 2, lty = 2)
            abline(v = c(rv[r.ub], rv[r.lb]), col = 2, lty = 3)

            legend(x = ifelse(bdiffp.rep$mle > 0, "topleft", "topright"), legend = c(paste("MLE: ", round(bdiffp.rep$mle, 2), sep = ""), paste("95% CI: ", round(bdiffp.rep$lb, 2), " : ", round(bdiffp.rep$ub, 2), sep = ""), paste("Z = ", round(bdiffp.rep$Z, 2), sep = ""), paste("aZ = ", round(bdiffp.rep$cZ, 2), sep = "")), bty = "n")
        }

        # distal plot
        par(mar = c(2.5, 3.5, 2.5, 3.5), mgp = c(1.5, 0.65, 0), cex = 0.9)
        #par(mar = c(2.5, 3.5, 0.5, 3.5), mgp = c(1.5, 0.65, 0), cex = 0.9)
        dp <- exp(do.call(rbind, lapply(jpl[[2]]$post, as.numeric)))
        cols <- rainbow(nrow(dp), s = 0.8)
        plot(c(), c(), xlim = range(prior$x), ylim = range(c(0, dp)), xlab = "expression level", ylab = "individual posterior", main = nam2)
        if(show.individual.posteriors) {
            lapply(seq_len(nrow(dp)), function(i) lines(prior$x, dp[i, ], col = rgb(0, 0.5, 1, alpha = 0.25)))
        }
        if(correct.batch) {
            par(new = TRUE)
            plot(prior$x, batch.jpl[[2]][1, ], axes = FALSE, ylab = "", xlab = "", type = 'l', col = 8, lty = 1, lwd = 2)
        }
        djpc <- jpl[[2]]$jp
        #legend(x = ifelse(which.max(na.omit(djpc)) > length(djpc)/2, "topleft", "topright"), bty = "n", col = cols, legend = rownames(dp), lty = rep(1, nrow(dp)))
        par(new = TRUE)
        jpr <- range(c(0, na.omit(djpc)))
        plot(prior$x, djpc, axes = FALSE, ylab = "", xlab = "", ylim = jpr, type = 'l', col = 1, lty = 1, lwd = 2)
        axis(4, pretty(jpr, 5), col = 1)
        mtext("joint posterior", side = 4, outer = FALSE, line = 2)
    }

    if(return.details) {
        if(correct.batch) { # with batch correction
            return(list(results = a.bdiffp.rep, difference.posterior = a.bdiffp, results.nobatchcorrection = bdiffp.rep))
        } else {
            return(list(results = bdiffp.rep, difference.posterior = bdiffp, posteriors = jpl))
        }
    } else {
        if(correct.batch) { # with batch correction
            return(a.bdiffp.rep)
        } else {
            return(bdiffp.rep)
        }
    }
}


# fit models to external (bulk) reference
##' Fit scde models relative to provided set of expression magnitudes
##'
##' If group-average expression magnitudes are available (e.g. from bulk measurement), this method can be used
##' to fit individual cell error models relative to that reference
##'
##' @param counts count matrix
##' @param reference a vector of expression magnitudes (read counts) corresponding to the rows of the count matrix
##' @param min.fpm minimum reference fpm of genes that will be used to fit the models (defaults to 1). Note: fpm is calculated from the reference count vector as reference/sum(reference)*1e6
##' @param n.cores number of cores to use
##' @param zero.count.threshold read count to use as an initial guess for the zero threshold
##' @param nrep number independent of mixture fit iterations to try (default = 1)
##' @param save.plots whether to write out a pdf file showing the model fits
##' @param plot.filename model fit pdf filename
##' @param verbose verbose level
##'
##' @return matrix of scde models
##'
##' @examples
##' \donttest{
##' data(es.mef.small)
##' cd <- es.mef.small
##' cd <- cd[rowSums(cd) > 0, ]
##' cd <- cd[, colSums(cd) > 1e4]
##' data(o.ifm)  # Load precomputed model. Use ?scde.error.models to see how o.ifm was generated
##' o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
##' # calculate joint posteriors across all cells
##' jp <- scde.posteriors(models = o.ifm, cd, o.prior, n.cores = 10, return.individual.posterior.modes = TRUE, n.randomizations = 100)
##' # use expected expression magnitude for each gene
##' av.mag <- as.numeric(jp$jp %*% as.numeric(colnames(jp$jp)))
##' # translate into counts
##' av.mag.counts <- as.integer(round(av.mag))
##' # now, fit alternative models using av.mag as a reference (normally this would correspond to bulk RNA expression magnitude)
##' ref.models <- scde.fit.models.to.reference(cd, av.mag.counts, n.cores = 1)
##' }
##'
##' @export
scde.fit.models.to.reference <- function(counts, reference, n.cores = 10, zero.count.threshold = 1, nrep = 1, save.plots = FALSE, plot.filename = "reference.model.fits.pdf", verbose = 0, min.fpm = 1) {
    return.compressed.models <- TRUE
    verbose <- 1
    ids <- colnames(counts)
    ml <- papply(seq_along(ids), function(i) {
        df <- data.frame(count = counts[, ids[i]], fpm = reference/sum(reference)*1e6)
        df <- df[df$fpm > min.fpm, ]
        m1 <- fit.nb2.mixture.model(df, nrep = nrep, verbose = verbose, zero.count.threshold = zero.count.threshold)
        if(return.compressed.models) {
            v <- get.compressed.v1.model(m1)
            cl <- clusters(m1)
            rm(m1)
            gc()
            return(list(model = v, clusters = cl))
        } else {
            return(m1)
        }
    }, n.cores = n.cores)
    names(ml) <- ids

    # check if there were errors in the multithreaded portion
    lapply(seq_along(ml), function(i) {
        if(class(ml[[i]]) == "try-error") {
            message("ERROR encountered in building a model for cell ", ids[i], ":")
            message(ml[[i]])
            tryCatch(stop(paste("ERROR encountered in building a model for cell ", ids[i])), error = function(e) stop(e))
        }
    })

    if(save.plots) {
        # model fits
        #CairoPNG(file = paste(group, "model.fits.png", sep = "."), width = 1024, height = 300*length(ids))
        pdf(file = plot.filename, width = 13, height = 4)
        #l <- layout(matrix(seq(1, 4*length(ids)), nrow = length(ids), byrow = TRUE), rep(c(1, 1, 1, 0.5), length(ids)), rep(1, 4*length(ids)), FALSE)
        l <- layout(matrix(seq(1, 4), nrow = 1, byrow = TRUE), rep(c(1, 1, 1, 0.5), 1), rep(1, 4), FALSE)
        par(mar = c(3.5, 3.5, 3.5, 0.5), mgp = c(2.0, 0.65, 0), cex = 0.9)
        invisible(lapply(seq_along(ids), function(i) {
            df <- data.frame(count = counts[, ids[i]], fpm = reference/sum(reference)*1e6)
            df <- df[df$fpm > min.fpm, ]
            plot.nb2.mixture.fit(ml[[i]], df, en = ids[i], do.par = FALSE, compressed.models = return.compressed.models)
        }))
        dev.off()
    }

    if(return.compressed.models) {
        # make a joint model matrix
        jmm <- data.frame(do.call(rbind, lapply(ml, function(m) m$model)))
        rownames(jmm) <- names(ml)
        jmm
        return(jmm)
    } else {
        return(ml)
    }
}


##' Determine principal components of a matrix using per-observation/per-variable weights
##'
##' Implements a weighted PCA
##'
##' @param mat matrix of variables (columns) and observations (rows)
##' @param matw  corresponding weights
##' @param npcs number of principal components to extract
##' @param nstarts number of random starts to use
##' @param smooth smoothing span
##' @param em.tol desired EM algorithm tolerance
##' @param em.maxiter maximum number of EM iterations
##' @param seed random seed
##' @param center whether mat should be centered (weighted centering)
##' @param n.shuffles optional number of per-observation randomizations that should be performed in addition to the main calculations to determine the lambda1 (PC1 eigenvalue) magnitude under such randomizations (returned in $randvar)
##'
##' @return a list containing eigenvector matrix ($rotation), projections ($scores), variance (weighted) explained by each component ($var), total (weighted) variance of the dataset ($totalvar)
##'
##' @examples
##' set.seed(0)
##' mat <- matrix( c(rnorm(5*10,mean=0,sd=1), rnorm(5*10,mean=5,sd=1)), 10, 10)  # random matrix
##' base.pca <- bwpca(mat)  # non-weighted pca, equal weights set automatically
##' matw <- matrix( c(rnorm(5*10,mean=0,sd=1), rnorm(5*10,mean=5,sd=1)), 10, 10)  # random weight matrix
##' matw <- abs(matw)/max(matw)
##' base.pca.weighted <- bwpca(mat, matw)  # weighted pca
##'
##' @export
bwpca <- function(mat, matw = NULL, npcs = 2, nstarts = 1, smooth = 0, em.tol = 1e-6, em.maxiter = 25, seed = 1, center = TRUE, n.shuffles = 0) {
    if(smooth<4) { smooth <- 0 }
    if(any(is.nan(matw))) {
      stop("bwpca: weight matrix contains NaN values")
    }
    if(any(is.nan(mat))) {
      stop("bwpca: value matrix contains NaN values")
    }
    if(is.null(matw)) {
        matw <- matrix(1, nrow(mat), ncol(mat))
        nstarts <- 1
    }
    if(center) { mat <- t(t(mat)-colSums(mat*matw)/colSums(matw)) }

    res <- .Call("baileyWPCA", mat, matw, npcs, nstarts, smooth, em.tol, em.maxiter, seed, n.shuffles, PACKAGE = "scde")
    #res <- bailey.wpca(mat, matw, npcs, nstarts, smooth, em.tol, em.maxiter, seed)
    rownames(res$rotation) <- colnames(mat)
    rownames(res$scores) <- rownames(mat)
    colnames(res$rotation) <- paste("PC", seq(1:ncol(res$rotation)), sep = "")
    res$sd <- t(sqrt(res$var))
    res
}


##' Winsorize matrix
##'
##' Sets the ncol(mat)*trim top outliers in each row to the next lowest value same for the lowest outliers
##'
##' @param mat matrix
##' @param trim fraction of outliers (on each side) that should be Winsorized, or (if the value is  >= 1) the number of outliers to be trimmed on each side
##'
##' @return Winsorized matrix
##'
##' @examples
##' set.seed(0)
##' mat <- matrix( c(rnorm(5*10,mean=0,sd=1), rnorm(5*10,mean=5,sd=1)), 10, 10)  # random matrix
##' mat[1,1] <- 1000  # make outlier
##' range(mat)  # look at range of values
##' win.mat <- winsorize.matrix(mat, 0.1)
##' range(win.mat)  # note outliers removed
##'
##' @export
winsorize.matrix <- function(mat, trim) {
    if(trim  >  0.5) { trim <- trim/ncol(mat)  }
    wm <- .Call("winsorizeMatrix", mat, trim, PACKAGE = "scde")
    rownames(wm) <- rownames(mat)
    colnames(wm) <- colnames(mat)
    return(wm)
}


############################ PAGODA functions


##' Build error models for heterogeneous cell populations, based on K-nearest neighbor cells.
##'
##' Builds cell-specific error models assuming that there are multiple subpopulations present
##' among the measured cells. The models for each cell are based on average expression estimates
##' obtained from K closest cells within a given group (if groups = NULL, then within the entire
##' set of measured cells). The method implements fitting of both the original log-fit models
##' (when linear.fit = FALSE), or newer linear-fit models (linear.fit = TRUE, default) with locally
##' fit overdispersion coefficient (local.theta.fit = TRUE, default).
##'
##' @param counts count matrix (integer matrix, rows- genes, columns- cells)
##' @param groups optional groups partitioning known subpopulations
##' @param cor.method correlation measure to be used in determining k nearest cells
##' @param k number of nearest neighbor cells to use during fitting. If k is set sufficiently high, all of the cells within a given group will be used.
##' @param min.nonfailed minimum number of non-failed measurements (within the k nearest neighbor cells) required for a gene to be taken into account during error fitting procedure
##' @param min.size.entries minimum number of genes to use for model fitting
##' @param min.count.threshold minimum number of reads required for a measurement to be considered non-failed
##' @param save.model.plots whether model plots should be saved (file names are (group).models.pdf, or cell.models.pdf if no group was supplied)
##' @param max.model.plots maximum number of models to save plots for (saves time when there are too many cells)
##' @param n.cores number of cores to use through the calculations
##' @param min.fpm optional parameter to restrict model fitting to genes with group-average expression magnitude above a given value
##' @param verbose level of verbosity
##' @param fpm.estimate.trim trim fraction to be used in estimating group-average gene expression magnitude for model fitting (0.5 would be median, 0 would turn off trimming)
##' @param linear.fit whether newer linear model fit with zero intercept should be used (T), or the log-fit model published originally (F)
##' @param local.theta.fit whether local theta fitting should be used (only available for the linear fit models)
##' @param theta.fit.range allowed range of the theta values
##' @param alpha.weight.power 1/theta weight power used in fitting theta dependency on the expression magnitude
##'
##' @return a data frame with parameters of the fit error models (rows- cells, columns- fitted parameters)
##'
##' @examples
##' \donttest{
##' data(pollen)
##' cd <- pollen
##' cd <- cd[,colSums(cd>0)>1.8e3]
##' cd <- cd[rowSums(cd)>10,]
##' cd <- cd[rowSums(cd>0)>5,]
##' knn <- knn.error.models(cd, k=ncol(cd)/4, n.cores=10, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
##' }
##'
##' @export
knn.error.models <- function(counts, groups = NULL, k = round(ncol(counts)/2), min.nonfailed = 5, min.count.threshold = 1, save.model.plots = TRUE, max.model.plots = 50, n.cores = parallel::detectCores(), min.size.entries = 2e3, min.fpm = 0, cor.method = "pearson", verbose = 0, fpm.estimate.trim = 0.25, linear.fit = TRUE, local.theta.fit = linear.fit, theta.fit.range = c(1e-2, 1e2), alpha.weight.power = 1/2) {
    threshold.prior = 1-1e-6

    # check for integer counts
    if(any(!unlist(lapply(counts,is.integer)))) {
      stop("Some of the supplied counts are not integer values (or stored as non-integer types). Aborting!\nThe method is designed to work on read counts - do not pass normalized read counts (e.g. FPKM values). If matrix contains read counts, but they are stored as numeric values, use counts<-apply(counts,2,function(x) {storage.mode(x) <- 'integer'; x}) to recast.");
    }

    # TODO:
    #  - implement check for k >= n.cells (to avoid correlation calculations)
    #  - implement error reporting/handling for failed cell fits

    if(is.null(groups)) {
        groups <- as.factor(rep("cell", ncol(counts)))
    }
    names(groups) <- colnames(counts)

    if(k >  ncol(counts)-1) {
        message("the value of k (", k, ") is too large, setting to ", (ncol(counts)-1))
        k <- ncol(counts)-1
    }

    ls <- estimate.library.sizes(counts, NULL, groups, min.size.entries, verbose = verbose, return.details = TRUE, vil = counts >= min.count.threshold)
    ca <- counts
    ca[ca<min.count.threshold] <- NA # a version of counts with all "drop-out" components set to NA
    mll <- tapply(colnames(counts), groups, function(ids) {
        # use Spearman rank correlation on pairwise complete observations to establish distance relationships between cells
        group <- as.character(groups[ids[1]])

        if(verbose > 0) {
            cat(group, ": calculating cell-cell similarities ...")
        }

        #if(n.cores > 1) { allowWGCNAThreads(n.cores) } else { disableWGCNAThreads() }
        #celld <- WGCNA::cor(log10(matrix(as.numeric(as.matrix(ca)), nrow = nrow(ca), ncol = ncol(ca))+1), method = cor.method, use = "p", nThreads = n.cores)
        if(is.element("WGCNA", installed.packages()[, 1])) {
            celld <- WGCNA::cor(sqrt(matrix(as.numeric(as.matrix(ca[, ids])), nrow = nrow(ca), ncol = length(ids))), method = cor.method, use = "p", nThreads = n.cores)
        } else {
            celld <- stats::cor(sqrt(matrix(as.numeric(as.matrix(ca[, ids])), nrow = nrow(ca), ncol = length(ids))), method = cor.method, use = "p")
        }
        rownames(celld) <- colnames(celld) <- ids

        if(verbose > 0) {
            cat(" done\n")
        }

        # TODO: correct for batch effect in cell-cell similarity matrix
        if(FALSE) {
            # number batches 10^(seq(0, n)) compute matrix of id sums, NA the diagonal,
            bid <- 10^(as.integer(batch)-1)
            bm <- matrix(bid, byrow = TRUE, nrow = length(bid), ncol = length(bid))+bid
            diag(bm) <- NA

            # use tapply to calculate means shifts per combination reconstruct shift vector, matrix, subtract
            # select the upper triangle, tapply to it to correct celld vector directly
        }

        if(verbose)  message(paste("fitting", group, "models:"))

        ml <- papply(seq_along(ids), function(i) { try({
            if(verbose)  message(paste(group, '.', i, " : ", ids[i], sep = ""))
            # determine k closest cells
            oc <- ids[-i][order(celld[ids[i], -i, drop = FALSE], decreasing = TRUE)[1:min(k, length(ids)-1)]]
            #set.seed(i)   oc <- sample(ids[-i], k)
            # determine a subset of genes that show up sufficiently often
            #fpm <- rowMeans(t(t(counts[, oc, drop = FALSE])/(ls$ls[oc])))
            fpm <- apply(t(ca[, oc, drop = FALSE])/(ls$ls[oc]), 2, mean, trim = fpm.estimate.trim, na.rm = TRUE)
            # rank genes by the number of non-zero occurrences, take top genes
            vi <- which(rowSums(counts[, oc] > min.count.threshold)  >=  min(ncol(oc)-1, min.nonfailed) & fpm > min.fpm)
            if(length(vi)<40)  message("WARNING: only ", length(vi), " valid genes were found to fit ", ids[i], " model")
            df <- data.frame(count = counts[vi, ids[i]], fpm = fpm[vi])

            # determine failed-component posteriors for each gene
            #fp <- ifelse(df$count <=  min.count.threshold, threshold.prior, 1-threshold.prior)
            fp <- ifelse(df$count <=  min.count.threshold & df$fpm  >=  median(df$fpm[df$count <=  min.count.threshold]), threshold.prior, 1-threshold.prior)
            cp <- cbind(fp, 1-fp)

            if(linear.fit) {
                # use a linear fit (nb2gth)
                m1 <- fit.nb2gth.mixture.model(df, prior = cp, nrep = 1, verbose = verbose, zero.count.threshold = min.count.threshold, full.theta.range = theta.fit.range, theta.fit.range = theta.fit.range, use.constant.theta.fit = !local.theta.fit, alpha.weight.power = alpha.weight.power)

            }  else {
                # mixture fit (the originally published method)
                m1 <- fit.nb2.mixture.model(df, prior = cp, nrep = 1, verbose = verbose, zero.count.threshold = min.count.threshold)
            }
            v <- get.compressed.v1.model(m1)
            cl <- clusters(m1)
            m1<-list(model = v, clusters = cl)
            #plot.nb2.mixture.fit(m1, df, en = ids[i], do.par = FALSE, compressed.models = TRUE)
            return(m1)
            #})
        })}, n.cores = n.cores)
        vic <- which(unlist(lapply(seq_along(ml), function(i) {
            if(class(ml[[i]]) == "try-error") {
                message("ERROR encountered in building a model for cell ", ids[i], " - skipping the cell. Error:")
                message(ml[[i]])
                #tryCatch(stop(paste("ERROR encountered in building a model for cell ", ids[i])), error = function(e) stop(e))
                return(FALSE);
            }
            return(TRUE);
        })))
        ml <- ml[vic]; names(ml) <- ids[vic];

        if(length(vic)<length(ids)) {
          message("ERROR fitting of ", (length(ids)-length(vic)), " out of ", length(ids), " cells resulted in errors reporting remaining ", length(vic), " cells")
        }
        if(length(vic)<length(ids)) {
                # model fits
                if(verbose)  message("plotting ", group, " model fits... ")
                tryCatch( {
                    pdf(file = paste(group, "model.fits.pdf", sep = "."), width = ifelse(local.theta.fit, 13, 15), height = 4)
                    l <- layout(matrix(seq(1, 4), nrow = 1, byrow = TRUE), rep(c(1, 1, 1, ifelse(local.theta.fit, 1, 0.5)), 1), rep(1, 4), FALSE)
                    par(mar = c(3.5, 3.5, 3.5, 0.5), mgp = c(2.0, 0.65, 0), cex = 0.9)
                    invisible(lapply(vic[1:min(max.model.plots, length(vic))], function(i) {
                        oc <- ids[-i][order(celld[ids[i], -i, drop = FALSE], decreasing = TRUE)[1:min(k, length(ids)-1)]]
                        #set.seed(i) oc <- sample(ids[-i], k)
                        # determine a subset of genes that show up sufficiently often
                        #fpm <- rowMeans(t(t(counts[, oc, drop = FALSE])/(ls$ls[oc])))
                        fpm <- apply(t(ca[, oc, drop = FALSE])/(ls$ls[oc]), 2, mean, trim = fpm.estimate.trim, na.rm = TRUE)
                        vi <- which(rowSums(counts[, oc] > min.count.threshold)  >=  min(ncol(oc)-1, min.nonfailed) & fpm > min.fpm)
                        df <- data.frame(count = counts[vi, ids[i]], fpm = fpm[vi])
                        plot.nb2.mixture.fit(ml[[ids[i]]], df, en = ids[i], do.par = FALSE, compressed.models = TRUE)
                    }))
                    dev.off()
                }, error = function(e) {
                    message("ERROR encountered during model fit plot outputs:")
                    message(e)
                    dev.off()
                })
        }

        return(ml)
    })


    # make a joint model matrix
    jmm <- data.frame(do.call(rbind, lapply(mll, function(tl) do.call(rbind, lapply(tl, function(m) m$model)))))
    rownames(jmm) <- unlist(lapply(mll, names))
    # reorder in the original cell order
    attr(jmm, "groups") <- rep(names(mll), unlist(lapply(mll, length)))
    return(jmm)
}


##' Normalize gene expression variance relative to transcriptome-wide expectations
##'
##' Normalizes gene expression magnitudes to ensure that the variance follows chi-squared statistics
##' with respect to its ratio to the transcriptome-wide expectation as determined by local regression
##' on expression magnitude (and optionally gene length). Corrects for batch effects.
##'
##' @param models model matrix (select a subset of rows to normalize variance within a subset of cells)
##' @param counts read count matrix
##' @param batch measurement batch (optional)
##' @param trim trim value for Winsorization (optional, can be set to 1-3 to reduce the impact of outliers, can be as large as 5 or 10 for datasets with several thousand cells)
##' @param prior expression magnitude prior
##' @param fit.genes a vector of gene names which should be used to establish the variance fit (default is NULL: use all genes). This can be used to specify, for instance, a set spike-in control transcripts such as ERCC.
##' @param plot whether to plot the results
##' @param minimize.underdispersion whether underdispersion should be minimized (can increase sensitivity in datasets with high complexity of population, however cannot be effectively used in datasets where multiple batches are present)
##' @param n.cores number of cores to use
##' @param n.randomizations number of bootstrap sampling rounds to use in estimating average expression magnitude for each gene within the given set of cells
##' @param weight.k k value to use in the final weight matrix
##' @param verbose verbosity level
##' @param weight.df.power power factor to use in determining effective number of degrees of freedom (can be increased for datasets exhibiting particularly high levels of noise at low expression magnitudes)
##' @param smooth.df degrees of freedom to be used in calculating smoothed local regression between coefficient of variation and expression magnitude (and gene length, if provided). Leave at -1 for automated guess.
##' @param max.adj.var maximum value allowed for the estimated adjusted variance (capping of adjusted variance is recommended when scoring pathway overdispersion relative to randomly sampled gene sets)
##' @param theta.range valid theta range (should be the same as was set in knn.error.models() call
##' @param gene.length optional vector of gene lengths (corresponding to the rows of counts matrix)
##'
##' @examples
##' \donttest{
##' data(pollen)
##' cd <- pollen
##' cd <- cd[,colSums(cd>0)>1.8e3]
##' cd <- cd[rowSums(cd)>10,]
##' cd <- cd[rowSums(cd>0)>5,]
##' knn <- knn.error.models(cd, k=ncol(cd)/4, n.cores=10, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
##' varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = FALSE)
##' }
##'
##' @return a list containing the following fields:
##' \itemize{
##' \item{mat} {adjusted expression magnitude values}
##' \item{matw} { weight matrix corresponding to the expression matrix}
##' \item{arv} { a vector giving adjusted variance values for each gene}
##' \item{avmodes} {a vector estimated average expression magnitudes for each gene}
##' \item{modes} {a list of batch-specific average expression magnitudes for each gene}
##' \item{prior} {estimated (or supplied) expression magnitude prior}
##' \item{edf} { estimated effective degrees of freedom}
##' \item{fit.genes} { fit.genes parameter }
##' }
##'
##' @export
pagoda.varnorm <- function(models, counts, batch = NULL, trim = 0, prior = NULL, fit.genes=NULL, plot = TRUE, minimize.underdispersion = FALSE, n.cores = detectCores(), n.randomizations = 100, weight.k = 0.9, verbose = 0, weight.df.power = 1, smooth.df = -1, max.adj.var = 10, theta.range = c(1e-2, 1e2), gene.length = NULL) {

    cd <- counts

    min.edf <- 1
    weight.k.internal <- 1
    use.mean.fpm <- FALSE
    use.expected.value <- TRUE
    cv.fit <- TRUE
    edf.damping <- 1

    # load NB extensions
    data(scde.edff, envir = environment())

    # subset cd to the cells occurring in the models
    if(verbose) { cat("checking counts ... ") }
    if(!all(rownames(models) %in% colnames(cd))) {
        stop(paste("supplied count matrix (cd) is missing data for the following cells:[", paste(rownames(models)[!rownames(models) %in% colnames(cd)], collapse = ", "), "]", sep = ""))
    }
    if(!length(rownames(models)) == length(colnames(cd)) || !all(rownames(models) == colnames(cd))) {
        cd <- cd[, match(rownames(models), colnames(cd))]
    }
    if(verbose) { cat("done\n") }

    # trim counts according to the extreme fpm values
    if(trim > 0) {
        if(verbose) { cat("Winsorizing count matrix ... ") }
        fpm <- t((t(log(cd))-models$corr.b)/models$corr.a)
        #tfpm <- log(winsorize.matrix(exp(fpm), trim = trim))
        tfpm <- winsorize.matrix(fpm, trim)
        rn <- rownames(cd)
        cn <- colnames(cd)
        cd <- round(exp(t(t(tfpm)*models$corr.a+models$corr.b)))
        cd[cd<0] <- 0
        rownames(cd) <- rn
        colnames(cd) <- cn
        rm(fpm, tfpm)
        cd <- cd[rowSums(cd) > 0, ] # omit genes without any data after Winsorization
        if(verbose) { cat("done\n") }
    }

    # check/fix batch vector
    if(verbose) { cat("checking batch ... ") }
    if(!is.null(batch)) {
        if(!is.factor(batch)) {
            batch <- as.factor(batch)
        }
        if(is.null(names(batch))) {
            if(length(batch) != nrow(models)) {
                stop("invalid batch vector supplied: length differs from nrow(models)!")
            }
            names(batch) <- rownames(models)
        } else {
            if(!all(rownames(models) %in% names(batch))) {
                stop(paste("invalid batch vector supplied: the following cell(s) are not present: [", paste(rownames(models)[!rownames(models) %in% names(batch)], collapse = ", "), "]", sep = ""))
            }
            batch <- batch[rownames(models)]
        }

        bt <- table(batch)
        min.batch.level <- 2
        if(any(bt<min.batch.level)) {
            if(verbose) { cat("omitting small batch levels [", paste(names(bt)[bt<min.batch.level], collapse = " "), "] ... ") }
            batch[batch %in% names(bt)[bt<min.batch.level]] <- names(bt)[which.max(bt)]
        }
    }
    if(verbose) { cat("ok\n") }

    # recalculate modes as needed
    if(verbose) { cat("calculating modes ... ") }
    if(is.null(prior)) {
        if(verbose) { cat("prior ") }
        prior <- scde.expression.prior(models = models, counts = cd, length.out = 400, show.plot = FALSE)
    }
    # dataset-wide mode
    if(use.mean.fpm) { # use mean fpm across cells
        avmodes <- modes <- rowMeans(exp(scde.expression.magnitude(models, cd)))
    } else { # use joint posterior mode/expected value
        jp <- scde.posteriors(models = models, cd, prior, n.cores = n.cores, return.individual.posterior.modes = TRUE, n.randomizations = n.randomizations)
        if(use.expected.value) {
            avmodes <- modes <- (jp$jp %*% as.numeric(colnames(jp$jp)))[, 1]
        } else { # use mode
            avmodes <- modes <- (as.numeric(colnames(jp$jp)))[max.col(jp$jp)]
        }
    }
    if(verbose) { cat(". ") }

    # batch-specific modes, if necessary
    if(!is.null(batch) && length(levels(batch)) > 1) {
        # calculate mode for each batch
        if(verbose) { cat("batch: [ ") }
        modes <- tapply(seq_len(nrow(models)), batch, function(ii) {
            if(verbose) { cat(as.character(batch[ii[1]]), " ") }
            if(use.mean.fpm) { # use mean fpm across cells
                modes <- rowMeans(exp(scde.expression.magnitude(models[ii, ], cd[, ii])))
            } else { # use joint posterior mode
                jp <- scde.posteriors(models = models[ii, ], cd[, ii], prior, n.cores = n.cores, return.individual.posterior.modes = TRUE, n.randomizations = n.randomizations)
                if(use.expected.value) {
                    modes <- (jp$jp %*% as.numeric(colnames(jp$jp)))[, 1]
                } else { # use mode
                    modes <- (as.numeric(colnames(jp$jp)))[max.col(jp$jp)]
                }
            }
        })
        # set dataset-wide mode
        #if(use.mean.fpm) { # use mean fpm across cells
        #  avmodes <- colMeans(do.call(rbind, modes)*as.vector(unlist(tapply(1:length(batch), batch, length))))*length(levels(batch))/length(batch)
        #jp <- scde.posteriors(models = models, cd, prior, n.cores = n.cores, return.individual.posterior.modes = TRUE, n.randomizations = n.randomizations)
        if(verbose) { cat("] ") }
    }
    if(verbose) { cat("done\n") }

    # check/calculate weights
    if(verbose) { cat("calculating weight matrix ... ") }

    # calculate default weighting scheme
    if(verbose) { cat("calculating ... ") }

    # dataset-wide version of matw (disregarding batch)
    sfp <- do.call(cbind, lapply(seq_len(ncol(cd)), function(i) ppois(cd[, i]-1, exp(models[i, "fail.r"]), lower.tail = FALSE)))
    mfp <- scde.failure.probability(models = models, magnitudes = log(avmodes))
    ofpT <- do.call(cbind, lapply(seq_len(ncol(cd)), function(i) { # for each cell
        lfpm <- log(avmodes)
        mu <- models$corr.b[i] + models$corr.a[i]*lfpm
        thetas <- get.corr.theta(models[i, ], lfpm, theta.range)
        pnbinom(1, size = thetas, mu = exp(mu), lower.tail = TRUE)
    }))
    matw <- 1-weight.k.internal*mfp*sfp # only mode failure probability
    # mode failure or NB failure
    #tmfp <- 1-(1-mfp)*(1-ofpT)
    #matw <- 1-weight.k.internal*tmfp*sfp


    # calculate batch-specific version of the weight matrix if needed
    if(!is.null(batch) && length(levels(batch)) > 1) { # with batch correction
        # save the dataset-wide one as avmatw
        # calculate mode for each batch
        if(verbose) { cat("batch: [ ") }
        bmatw <- do.call(cbind, tapply(seq_len(nrow(models)), batch, function(ii) {
            if(verbose) { cat(as.character(batch[ii[1]]), " ") }
            # set self-fail probability to p(count|background)
            # total mode failure (including overdispersion dropouts)
            #sfp <- do.call(cbind, lapply(ii, function(i) dpois(cd[, i], exp(models[i, "fail.r"]), log = FALSE)))
            sfp <- do.call(cbind, lapply(ii, function(i) ppois(cd[, i]-1, exp(models[i, "fail.r"]), lower.tail = FALSE)))

            mfp <- scde.failure.probability(models = models[ii, ], magnitudes = log(modes[[batch[ii[1]]]]))
            ofpT <- do.call(cbind, lapply(ii, function(i) { # for each cell
                lfpm <- log(modes[[batch[i]]])
                mu <- models$corr.b[i] + models$corr.a[i]*lfpm
                thetas <- get.corr.theta(models[i, ], lfpm, theta.range)
                pnbinom(1, size = thetas, mu = exp(mu), lower.tail = TRUE)
            }))

            x <- 1-weight.k.internal*mfp*sfp # only mode failure probability
            # mode failure or NB failure
            #tmfp <- 1-(1-mfp)*(1-ofpT)
            #x <- 1-weight.k.internal*tmfp*sfp
        }))
        # reorder
        bmatw <- bmatw[, rownames(models)]
        if(verbose) { cat("] ") }
    }
    if(verbose) { cat("done\n") }

    # calculate effective degrees of freedom
    # total effective degrees of freedom per gene
    if(verbose) { cat("calculating effective degrees of freedom ..") }
    ids <- 1:ncol(cd)
    names(ids) <- colnames(cd)
    # dataset-wide version
    edf.mat <- do.call(cbind, papply(ids, function(i) {
        v <- models[i, ]
        lfpm <- log(avmodes)
        mu <- exp(lfpm*v$corr.a + v$corr.b)
        # adjust very low mu levels except for those that have 0 counts (to avoid inf values)

        thetas <- get.corr.theta(v, lfpm, theta.range)
        edf <- exp(predict(scde.edff, data.frame(lt = log(thetas))))
        edf[thetas > 1e3] <- 1
        edf
    }, n.cores = n.cores))
    if(edf.damping != 1) {
        edf.mat <- ((edf.mat/ncol(edf.mat))^edf.damping) * ncol(edf.mat)
    }

    # incorporate weight into edf
    #edf.mat <- ((matw^weight.df.power)*edf.mat)
    edf.mat <- (matw*edf.mat)^weight.df.power
    #edf <- rowSums(matw*edf.mat)+1.5 # summarize eDF per gene
    edf <- rowSums(edf.mat)+1 # summarize eDF per gene
    if(verbose) { cat(".") }

    # batch-specific version if necessary
    if(!is.null(batch) && is.list(modes)) { # batch-specific mode
        bedf.mat <- do.call(cbind, papply(ids, function(i) {
            v <- models[i, ]
            lfpm <- log(modes[[batch[i]]])
            mu <- exp(lfpm*v$corr.a + v$corr.b)
            # adjust very low mu levels except for those that have 0 counts (to avoid inf values)

            thetas <- get.corr.theta(v, lfpm, theta.range)
            edf <- exp(predict(scde.edff, data.frame(lt = log(thetas))))
            edf[thetas > 1e3] <- 1
            return(edf)
        }, n.cores = n.cores))
        if(edf.damping != 1) { bedf.mat <-  ((bedf.mat/ncol(bedf.mat))^edf.damping) * ncol(edf.mat) }

        # incorporate weight into edf
        #bedf.mat <- ((bmatw^weight.df.power)*bedf.mat)
        bedf.mat <- (bmatw*bedf.mat)^weight.df.power
        bedf <- rowSums(bedf.mat)+1 # summarize eDF per gene
        if(verbose) { cat(".") }
    }

    if(verbose) { cat(" done\n") }

    if(verbose) { cat("calculating normalized expression values ... ") }
    # evaluate negative binomial deviations and effective degrees of freedom
    ids <- 1:ncol(cd)
    names(ids) <- colnames(cd)
    mat <- do.call(cbind, papply(ids, function(i) {
        v <- models[i, ]
        lfpm <- log(avmodes)
        mu <- exp(lfpm*v$corr.a + v$corr.b)
        # adjust very low mu levels except for those that have 0 counts (to avoid inf values)
        thetas <- get.corr.theta(v, lfpm, theta.range)

        #matw[, i]*edf.mat[, i]*(cd[, i]-mu)^2/(mu+mu^2/thetas)
        #x <- (cd[, i]-mu)^2/(mu+mu^2/thetas)
        #edf.mat[, i]*(cd[, i]-mu)^2/(mu+mu^2/thetas)
        # considering Poisson-nb mixture
        fail.lambda <- exp(as.numeric(v["fail.r"]))
        #edf.mat[, i]*(cd[, i]-mu)^2/(matw[, i]*(mu+mu^2/thetas) + (1-matw[, i])*((mu-fail.lambda)^2 + fail.lambda))
        edf.mat[, i]*(cd[, i]-mu)^2/(mu+mu^2/thetas +  fail.lambda)

        #edf.mat[, i]*(cd[, i]-mu)^2/(matw[, i]*mu+(mu^2)*((1-matw[, i])+matw[, i]/thetas))
    }, n.cores = n.cores))
    rownames(mat) <- rownames(cd)
    if(verbose) { cat(".") }
    # batch-specific version of mat
    if(!is.null(batch) && is.list(modes)) { # batch-specific mode
        bmat <- do.call(cbind, papply(ids, function(i) {
            v <- models[i, ]
            lfpm <- log(modes[[batch[i]]])
            mu <- exp(lfpm*v$corr.a + v$corr.b)
            # adjust very low mu levels except for those that have 0 counts (to avoid inf values)
            thetas <- get.corr.theta(v, lfpm, theta.range)

            #matw[, i]*edf.mat[, i]*(cd[, i]-mu)^2/(mu+mu^2/thetas)
            #x <- (cd[, i]-mu)^2/(mu+mu^2/thetas)
            #edf.mat[, i]*(cd[, i]-mu)^2/(mu+mu^2/thetas)
            #edf.mat[, i]*(cd[, i]-mu)^2/(matw[, i]*mu+(mu^2)*((1-matw[, i])+matw[, i]/thetas))
            fail.lambda <- exp(as.numeric(v["fail.r"]))
            #edf.mat[, i]*(cd[, i]-mu)^2/(matw[, i]*(mu+mu^2/thetas) + (1-matw[, i])*((mu-fail.lambda)^2 + fail.lambda))
            edf.mat[, i]*(cd[, i]-mu)^2/(mu+mu^2/thetas +  fail.lambda)
        }, n.cores = n.cores))
        rownames(bmat) <- rownames(cd)

        if(verbose) { cat(".") }
    }
    if(verbose) { cat(" done\n") }

    # do a model fit on the weighted standard deviation (as a function of the batch-average expression mode)
    wvar <- rowSums(mat)/rowSums(edf.mat)

    if(!is.null(batch) && is.list(modes)) { # batch-specific mode
        # estimate the ratio of the batch-specific variance to the total dataset variance
        bwvar <- rowSums(bmat)/rowSums(bedf.mat)
        bwvar.ratio <- bwvar/wvar
        wvar <- bwvar # replace wvar now that we have the ratio of
        matw <- bmatw # replace matw with the batch-specific one that will be used from here on
        # ALTERNATIVE: could adjust wvar for the bwvar.ratio here, before fitting expression dependency
      }
    fvi <- vi <- rowSums(matw) > 0 & is.finite(wvar) & wvar > 0
    if(!is.null(fit.genes)) { fvi <- fvi & rownames(mat) %in% fit.genes }
    if(!any(fvi)) { stop("unable to find a set of valid genes to establish the variance fit") }

    # s = mgcv:::s
    s = mgcv::s
    if(cv.fit) {
        #x <- gam(as.formula("cv2 ~ s(lev)"), data = df[vi, ], weights = rowSums(matw[vi, ]))
        if(is.null(gene.length)) {
            df <- data.frame(lev = log10(avmodes), cv2 = log10(wvar/avmodes^2))
            x <- mgcv::gam(cv2 ~ s(lev, k = smooth.df), data = df[fvi, ], weights = rowSums(matw[fvi, ]))
        } else {
            df <- data.frame(lev = log10(avmodes), cv2 = log10(wvar/avmodes^2), len = gene.length[rownames(cd)])
            x <- mgcv::gam(cv2 ~ s(lev, k = smooth.df) + s(len, k = smooth.df), data = df[fvi, ], weights = rowSums(matw[fvi, ]))
        }
        #x <- lm(cv2~lev, data = df[vi, ], weights = rowSums(matw[vi, ]))

        zval.m <- 10^(df$cv2[vi]-predict(x, newdata = df[vi, ]))

        if(plot) {
            par(mfrow = c(1, 2), mar = c(3.5, 3.5, 1.0, 1.0), mgp = c(2, 0.65, 0))
            #smoothScatter(df$lev[vi], log(wvar[vi]), nbin = 256, xlab = "expression magnitude (log10)", ylab = "wvar (log)") abline(h = 0, lty = 2, col = 2)
            #points(df[paste("g", diff.exp.gene.ids, sep = ""), "lev"], log(wvar[paste("g", diff.exp.gene.ids, sep = "")]), col = 2)

            smoothScatter(df$lev[vi], df$cv2[vi], nbin = 256, xlab = "expression magnitude (log10)", ylab = "cv^2 (log10)")
            lines(sort(df$lev[vi]), predict(x, newdata = df[vi, ])[order(df$lev[vi])], col = 2, pch = ".", cex = 1)
            if(!is.null(fit.genes)) { # show genes used for the fit
              points(df$lev[fvi],df$cv2[fvi],pch=".",col="green",cex=1)
            }

            #points(df[paste("g", diff.exp.gene.ids, sep = ""), "lev"], df[paste("g", diff.exp.gene.ids, sep = ""), "cv2"], col = 2)
        }

        # optional : re-weight to minimize the underdispersed points
        if(minimize.underdispersion) {
            pv <- pchisq(zval.m*(edf[vi]-1), edf[vi], log.p = FALSE, lower.tail = TRUE)
            pv[edf[vi]<= min.edf] <- 0
            pv <- p.adjust(pv)
            #x <- gam(as.formula("cv2 ~ s(lev)"), data = df[vi, ], weights = (pmin(10, -log(pv))+1)*rowSums(matw[vi, ]))
            x <- mgcv::gam(cv2 ~ s(lev, k = smooth.df), data = df[fvi, ], weights = (pmin(10, -log(pv))+1)*rowSums(matw[fvi, ]))
            zval.m <- 10^(df$cv2[vi]-predict(x,newdata=df[vi,]))
            if(plot) {
              lines(sort(df$lev[vi]), predict(x, newdata = df[vi, ])[order(df$lev[vi])], col = 4, pch = ".", cex = 1)
            }
        }
    } else {
        df <- data.frame(lev = log10(avmodes), sd = sqrt(wvar))
        #x <- gam(as.formula("sd ~ s(lev)"), data = df[vi, ], weights = rowSums(matw[vi, ]))
        x <- mgcv::gam(sd ~ s(lev, k = smooth.df), data = df[fvi, ], weights = rowSums(matw[fvi, ]))
        zval.m <- (as.numeric((df$sd[vi])/pmax(min.sd, predict(x,newdata=df[vi,]))))^2

        if(plot) {
            par(mfrow = c(1, 2), mar = c(3.5, 3.5, 1.0, 1.0), mgp = c(2, 0.65, 0))
            smoothScatter(df$lev[vi], df$sd[vi], nbin = 256, xlab = "expression magnitude", ylab = "weighted sdiv")
            lines(sort(df$lev[vi]), predict(x, newdata = df[vi, ])[order(df$lev[vi])], col = 2, pch = ".", cex = 1)
            if(!is.null(fit.genes)) { # show genes used for the fit
              points(df$lev[fvi],df$sd[fvi],pch=".",col="green",cex=1)
            }
        }

        # optional : re-weight to minimize the underdispersed points
        if(minimize.underdispersion) {
            pv <- pchisq(zval.m*(edf[vi]-1), edf[vi], log.p = FALSE, lower.tail = TRUE)
            pv[edf[vi]<= min.edf] <- 0
            pv <- p.adjust(pv)
            #x <- gam(as.formula("sd ~ s(lev)"), data = df[vi, ], weights = (pmin(20, -log(pv))+1)*rowSums(matw[vi, ]))
            x <- mgcv::gam(sd ~ s(lev, k = smooth.df), data = df[fvi, ], weights = (pmin(20, -log(pv))+1)*rowSums(matw[fvi, ]))
            zval.m <- (as.numeric((df$sd[vi])/pmax(x$fitted.values, min.sd)))^2
            if(plot) {
              lines(sort(df$lev[vi]), predict(x, newdata = df[vi, ])[order(df$lev[vi])], col = 4, pch = ".", cex = 1)
            }
        }
    }

    # adjust for inter-batch variance
    if(!is.null(batch) && is.list(modes)) { # batch-specific mode
        #zval.m <- zval.m*pmin(bwvar.ratio[vi], 1) # don't increase zval.m even if batch-specific specific variance is higher than the dataset-wide variance
        zval.m <- zval.m*pmin(bwvar.ratio[vi], 1/bwvar.ratio[vi]) # penalize for strong deviation in either direction
    }

    # calculate adjusted variance
    qv <- pchisq(zval.m*(edf[vi]-1), edf[vi], log.p = TRUE, lower.tail = FALSE)
    qv[edf[vi]<= min.edf] <- 0
    qv[abs(qv)<1e-10] <- 0
    arv <- rep(NA, length(vi))
    arv[vi] <- qchisq(qv, ncol(matw)-1, lower.tail = FALSE, log.p = TRUE)/ncol(matw)
    arv <- pmin(max.adj.var, arv)
    names(arv) <- rownames(cd)
    if(plot) {
        smoothScatter(df$lev[vi], arv[vi], xlab = "expression magnitude (log10)", ylab = "adjusted variance (log10)", nbin = 256)
        abline(h = 1, lty = 2, col = 8)
        abline(h = max.adj.var, lty = 3, col = 2)
        if(!is.null(fit.genes)) {
          points(df$lev[fvi],arv[fvi],pch=".",col="green",cex=1)
        }
        #points(df[paste("g", diff.exp.gene.ids, sep = ""), "lev"], arv[paste("g", diff.exp.gene.ids, sep = "")], col = 2)
        #points(df$lev[vi], arv[vi], col = 2, pch = ".", cex = 2)
    }

    # Wilcox score upper bound
    wsu <- function(k, n, z = qnorm(0.975)) {
        p <- k/n
        pmin(1, (2*n*p+z^2+(z*sqrt(z^2-1/n+4*n*p*(1-p)-(4*p-2)) +1))/(2*(n+z^2)))
    }

    # use milder weight matrix
    #matw <- 1-0.9*((1-matw)^2) # milder weighting for the the PCA (1-0.9*sp*mf)
    matw <- 1-weight.k*(1-matw) # milder weighting for the the PCA (1-0.9*sp*mf)
    matw <- matw/rowSums(matw)
    mat <- log10(exp(scde.expression.magnitude(models, cd))+1)

    # estimate observed variance (for scaling) before batch adjustments
    #varm <- sqrt(arv/pmax(weightedMatVar(mat, matw, batch = batch), 1e-5)) varm[varm<1e-5] <- 1e-5 mat <- mat*varm
    ov <- weightedMatVar(mat, matw)
    vr <- arv/ov
    vr[ov <=  0] <- 0

    if(!is.null(batch) && is.list(modes)) { # batch-specific mode
        # adjust proportion of zeros
        # determine lowest upper bound of non-zero measurement probability among the batch (for each gene)
        nbub <- apply(do.call(cbind, tapply(seq_len(ncol(mat)), batch, function(ii) {
            wsu(rowSums(mat[, ii] > 0), length(ii), z = qnorm(1-1e-2))
        })), 1, min)

        # decrease the batch weights for each gene to match the total
        # expectation of the non-zero measurements
        nbo <- do.call(cbind, tapply(seq_len(ncol(mat)), batch, function(ii) {
            matw[, ii]*pmin(1, ceiling(nbub*length(ii))/rowSums(mat[, ii] > 0))
        }))
        nbo <- nbo[, colnames(matw)]
        matw <- nbo

        ## # center 0 and non-0 observations between batches separately
        ## amat <- mat amat[amat == 0] <- NA
        ## amat.av <- rowMeans(amat, na.rm = TRUE) # dataset means
        ## # adjust each batch by the mean of its non-0 measurements
        ## amat <- do.call(cbind, tapply(1:ncol(amat), batch, function(ii) {
        ##   amat[, ii]-rowMeans(amat[, ii], na.rm = TRUE)
        ## }))
        ## amat <- amat[, colnames(mat)] # fix the ordering
        ## # shift up each gene by the dataset mean
        ## amat <- amat+amat.av
        ## amat[is.na(amat)] <- 0
        ## mat <- amat

        amat <- mat
        nr <- ncol(matw)/rowSums(matw)
        amat.av <- rowMeans(amat*matw)*nr # dataset means
        amat <- do.call(cbind, tapply(seq_len(ncol(amat)), batch, function(ii) {
            amat[, ii]-(rowMeans(amat[, ii]*matw[, ii]*nr, na.rm = TRUE))
        }))
        amat <- amat[, colnames(matw)]
        mat <- amat+amat.av

        # alternative: actually zero-out entries in mat
        ## nbub <- rowMin(do.call(cbind, tapply(1:ncol(amat), batch, function(ii) {
        ##   wsu(rowSums(amat[, ii] > 0), length(ii), z = qnorm(1-1e-2))
        ## })))
        ## set.seed(0)

        ## # decrease the batch weights for each gene to match the total
        ## # expectation of the non-zero measurements
        ## matm <- do.call(cbind, tapply(1:ncol(amat), batch, function(ii) {
        ##   # number of entries to zero-out per gene
        ##   nze <- rowSums(amat[, ii] > 0) - ceiling(nbub*length(ii))
        ##   # construct mat multiplier submatrix
        ##   sa <- rep(1, length(ii))
        ##   smatm <- do.call(rbind, lapply(1:length(nze), function(ri) {
        ##     if(nze[ri]<1) { return(sa) }
        ##     vi <- which(mat[ri, ii] > 0)
        ##     a <- sa a[vi[sample.int(length(vi), nze[ri])]] <- 0
        ##     a
        ##   }))
        ##   colnames(smatm) <- colnames(mat[, ii])
        ##   rownames(smatm) <- rownames(mat)
        ##   smatm
        ## }))
        ## matm <- matm[, colnames(mat)]
        ## mat <- mat*matm
        ## matw <- matw*matm
    }

    # center (no batch)
    mat <- weightedMatCenter(mat, matw)
    mat <- mat*sqrt(vr)

    if(!is.null(batch) && is.list(modes)) { # batch-specific mode
        return(list(mat = mat, matw = matw, arv = arv, modes = modes, avmodes = avmodes, prior = prior, edf = edf, batch = batch, trim = trim, bwvar.ratio = bwvar.ratio))
    } else {
        return(list(mat = mat, matw = matw, arv = arv, modes = modes, avmodes = avmodes, prior = prior, edf = edf, batch = batch, trim = trim))
    }
}


##' Control for a particular aspect of expression heterogeneity in a given population
##'
##' Similar to subtracting n-th principal component, the current procedure determines
##' (weighted) projection of the expression matrix onto a specified aspect (some pattern
##' across cells, for instance sequencing depth, or PC corresponding to an undesired process
##' such as ribosomal pathway variation) and subtracts it from the data so that it is controlled
##' for in the subsequent weighted PCA analysis.
##'
##' @param varinfo normalized variance info (from pagoda.varnorm())
##' @param aspect a vector giving a cell-to-cell variation pattern that should be controlled for (length should be corresponding to ncol(varinfo$mat))
##' @param center whether the matrix should be re-centered following pattern subtraction
##'
##' @return a modified varinfo object with adjusted expression matrix (varinfo$mat)
##'
##' @examples
##' \donttest{
##' data(pollen)
##' cd <- pollen
##' cd <- cd[,colSums(cd>0)>1.8e3]
##' cd <- cd[rowSums(cd)>10,]
##' cd <- cd[rowSums(cd>0)>5,]
##' knn <- knn.error.models(cd, k=ncol(cd)/4, n.cores=10, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
##' varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = FALSE)
##' # create go environment
##' library(org.Hs.eg.db)
##' ids <- unlist(lapply(mget(rownames(cd), org.Hs.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
##' rids <- names(ids)
##' names(rids) <- ids
##' go.env <- eapply(org.Hs.egGO2ALLEGS, function(x) as.character(na.omit(rids[x])))
##' go.env <- go.env[unlist(lapply(go.env, length))>5]
##' library(GO.db)
##' desc <- unlist(lapply(mget(names(go.env), GOTERM, ifnotfound = NA), function(x) if(is.logical(x)) { return("") } else { slot(x, "Term")}))
##' names(go.env) <- paste(names(go.env), desc)  # append description to the names
##' go.env <- list2env(go.env)  # convert to an environment
##' # subtract the pattern
##' cc.pattern <- pagoda.show.pathways(ls(go.env)[1:2], varinfo, go.env, show.cell.dendrogram = TRUE, showRowLabels = TRUE)  # Look at pattern from 2 GO annotations
##' varinfo.cc <- pagoda.subtract.aspect(varinfo, cc.pattern)
##' }
##'
##' @export
pagoda.subtract.aspect <- function(varinfo, aspect, center = TRUE) {
    if(length(aspect) != ncol(varinfo$mat)) { stop("aspect should be a numeric vector of the same length as the number of cells (i.e. ncol(varinfo$mat))") }
    v <- aspect
    v <- v-mean(v)
    v <- v/sqrt(sum(v^2))
    nr <- ((varinfo$mat * varinfo$matw) %*% v)/(varinfo$matw %*% v^2)
    mat.c <- varinfo$mat - t(v %*% t(nr))
    if(center) {
        mat.c <- weightedMatCenter(mat.c, varinfo$matw) # this commonly re-introduces some background dependency because of the matw
    }
    varinfo$mat <- mat.c
    varinfo
}


##' Run weighted PCA analysis on pre-annotated gene sets
##'
##' For each valid gene set (having appropriate number of genes) in the provided environment (setenv),
##' the method will run weighted PCA analysis, along with analogous analyses of random gene sets of the
##' same size, or shuffled expression magnitudes for the same gene set.
##'
##' @param varinfo adjusted variance info from pagoda.varinfo() (or pagoda.subtract.aspect())
##' @param setenv environment listing gene sets (contains variables with names corresponding to gene set name, and values being vectors of gene names within each gene set)
##' @param n.components number of principal components to determine for each gene set
##' @param n.cores number of cores to use
##' @param min.pathway.size minimum number of observed genes that should be contained in a valid gene set
##' @param max.pathway.size maximum number of observed genes in a valid gene set
##' @param n.randomizations number of random gene sets (of the same size) to be evaluated in parallel with each gene set (can be kept at 5 or 10, but should be increased to 50-100 if the significance of pathway overdispersion will be determined relative to random gene set models)
##' @param n.internal.shuffles number of internal (independent row shuffles) randomizations of expression data that should be evaluated for each gene set (needed only if one is interested in gene set coherence P values, disabled by default; set to 10-30 to estimate)
##' @param n.starts number of random starts for the EM method in each evaluation
##' @param center whether the expression matrix should be recentered
##' @param batch.center whether batch-specific centering should be used
##' @param proper.gene.names alternative vector of gene names (replacing rownames(varinfo$mat)) to be used in cases when the provided setenv uses different gene names
##' @param verbose verbosity level
##'
##' @return a list of weighted PCA info for each valid gene set
##'
##' @examples
##' \donttest{
##' data(pollen)
##' cd <- pollen
##' cd <- cd[,colSums(cd>0)>1.8e3]
##' cd <- cd[rowSums(cd)>10,]
##' cd <- cd[rowSums(cd>0)>5,]
##' knn <- knn.error.models(cd, k=ncol(cd)/4, n.cores=10, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
##' varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = FALSE)
##' # create go environment
##' library(org.Hs.eg.db)
##' ids <- unlist(lapply(mget(rownames(cd), org.Hs.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
##' rids <- names(ids)
##' names(rids) <- ids
##' go.env <- eapply(org.Hs.egGO2ALLEGS, function(x) as.character(na.omit(rids[x])))
##' go.env <- go.env[unlist(lapply(go.env, length))>5]
##' library(GO.db)
##' desc <- unlist(lapply(mget(names(go.env), GOTERM, ifnotfound = NA), function(x) if(is.logical(x)) { return("") } else { slot(x, "Term")}))
##' names(go.env) <- paste(names(go.env), desc)  # append description to the names
##' go.env <- list2env(go.env)  # convert to an environment
##' pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components=1, n.cores=10, n.internal.shuffles=50)
##' }
##'
##' @export
pagoda.pathway.wPCA <- function(varinfo, setenv, n.components = 2, n.cores = detectCores(), min.pathway.size = 10, max.pathway.size = 1e3, n.randomizations = 10, n.internal.shuffles = 0, n.starts = 10, center = TRUE, batch.center = TRUE, proper.gene.names = NULL, verbose = 0) {
    mat <- varinfo$mat
    matw <- varinfo$matw
    gsl <- NULL
    return.gsl <- FALSE
    smooth <- 0
    if(batch.center) { batch <- varinfo$batch } else { batch <- NULL }

    if(is.null(proper.gene.names)) { proper.gene.names <- rownames(mat) }


    if(center) {
        mat <- weightedMatCenter(mat, matw, batch = batch)
    }

    vi <- apply(mat, 1, function(x) sum(abs(diff(x))) > 0)
    vi[is.na(vi)] <- FALSE
    mat <- mat[vi, , drop = FALSE] # remove constant rows
    matw <- matw[vi, , drop = FALSE]
    proper.gene.names <- proper.gene.names[vi]

    if(is.null(gsl)) {
        gsl <- ls(envir = setenv)
        gsl.ng <- unlist(lapply(sn(gsl), function(go) sum(unique(get(go, envir = setenv)) %in% proper.gene.names)))
        gsl <- gsl[gsl.ng >= min.pathway.size & gsl.ng<= max.pathway.size]
        names(gsl) <- gsl
    }
    if(verbose) {
        message("processing ", length(gsl), " valid pathways")
    }
    if(return.gsl) return(gsl)


    # transpose mat to save a bit of calculations
    mat <- t(mat)
    matw <- t(matw)

    mcm.pc <- papply(gsl, function(x) {
        lab <- proper.gene.names %in% get(x, envir = setenv)
        if(sum(lab)<1) { return(NULL) }

        #smooth <- round(sum(lab)*smooth.fraction)
        #smooth <- max(sum(lab), smooth)

        #xp <- pca(d, nPcs = n.components, center = TRUE, scale = "none")
        #xp <- epca(mat[, lab], ncomp = n.components, center = FALSE, nstarts = n.starts)
        xp <- bwpca(mat[, lab, drop = FALSE], matw[, lab, drop = FALSE], npcs = n.components, center = FALSE, nstarts = n.starts, smooth = smooth, n.shuffles = n.internal.shuffles)

        # get standard deviations for the random samples
        ngenes <- sum(lab)
        z <- do.call(rbind, lapply(seq_len(n.randomizations), function(i) {
            si <- sample(1:ncol(mat), ngenes)
            #epca(mat[, si], ncomp = 1, center = FALSE, nstarts = n.starts)$sd
            xp <- bwpca(mat[, si, drop = FALSE], matw[, si, drop = FALSE], npcs = 1, center = FALSE, nstarts = n.starts, smooth = smooth)$sd
        }))

        # flip orientations to roughly correspond with the means
        cs <- unlist(lapply(seq_len(ncol(xp$scores)), function(i) sign(cor(xp$scores[, i], colMeans(t(mat[, lab, drop = FALSE])*abs(xp$rotation[, i]))))))

        xp$scores <- t(t(xp$scores)*cs)
        xp$rotation <- t(t(xp$rotation)*cs)

        # local normalization of each component relative to sampled PC1 sd
        avar <- pmax(0, (xp$sd^2-mean(z[, 1]^2))/sd(z[, 1]^2))
        xv <- t(xp$scores)
        xv <- xv/apply(xv, 1, sd)*sqrt(avar)
        return(list(xv = xv, xp = xp, z = z, sd = xp$sd, n = ngenes))
    }, n.cores = n.cores)
}


##' Estimate effective number of cells based on lambda1 of random gene sets
##'
##' Examines the dependency between the amount of variance explained by the first principal component
##' of a gene set and the number of genes in a gene set to determine the effective number of cells
##' for the Tracy-Widom distribution
##'
##' @param pwpca result of the pagoda.pathway.wPCA() call with n.randomizations > 1
##' @param start optional starting value for the optimization (if the NLS breaks, trying high starting values usually fixed the local gradient problem)
##'
##' @return effective number of cells
##'
##' @examples
##' \donttest{
##' data(pollen)
##' cd <- pollen
##' cd <- cd[,colSums(cd>0)>1.8e3]
##' cd <- cd[rowSums(cd)>10,]
##' cd <- cd[rowSums(cd>0)>5,]
##' knn <- knn.error.models(cd, k=ncol(cd)/4, n.cores=10, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
##' varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = FALSE)
##' pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components=1, n.cores=10, n.internal.shuffles=50)
##' pagoda.effective.cells(pwpca)
##' }
##'
##' @export
pagoda.effective.cells <- function(pwpca, start = NULL) {
    n.genes <- unlist(lapply(pwpca, function(x) rep(x$n, nrow(x$z))))
    var <- unlist(lapply(pwpca, function(x) x$z[, 1]))^2
    if(is.null(start)) { start <- nrow(pwpca[[1]]$xp$scores)*10 }

    n.cells <- nrow(pwpca[[1]]$xp$scores)
    of <- function(p, v, sp) {
        sn <- p[1]
        vfit <- (sn+sp)^2/(sn*sn+1/2) -1.2065335745820*(sn+sp)*((1/sn + 1/sp)^(1/3))/(sn*sn+1/2)
        residuals <- (v-vfit)^2
        return(sum(residuals))
    }
    x <- nlminb(objective = of, start = c(start), v = var, sp = sqrt(n.genes-1/2), lower = c(1), upper = c(n.cells))
    return((x$par)^2+1/2)
}


##' Determine de-novo gene clusters and associated overdispersion info
##'
##' Determine de-novo gene clusters, their weighted PCA lambda1 values, and random matrix expectation.
##'
##' @param varinfo varinfo adjusted variance info from pagoda.varinfo() (or pagoda.subtract.aspect())
##' @param trim additional Winsorization trim value to be used in determining clusters (to remove clusters that group outliers occurring in a given cell). Use higher values (5-15) if the resulting clusters group outlier patterns
##' @param n.clusters number of clusters to be determined (recommended range is 100-200)
##' @param cor.method correlation method ("pearson", "spearman") to be used as a distance measure for clustering
##' @param n.samples number of randomly generated matrix samples to test the background distribution of lambda1 on
##' @param n.starts number of wPCA EM algorithm starts at each iteration
##' @param n.internal.shuffles number of internal shuffles to perform (only if interested in set coherence, which is quite high for clusters by definition, disabled by default; set to 10-30 shuffles to estimate)
##' @param n.cores number of cores to use
##' @param verbose verbosity level
##' @param plot whether a plot showing distribution of random lambda1 values should be shown (along with the extreme value distribution fit)
##' @param show.random whether the empirical random gene set values should be shown in addition to the Tracy-Widom analytical approximation
##' @param n.components number of PC to calculate (can be increased if the number of clusters is small and some contain strong secondary patterns - rarely the case)
##' @param method clustering method to be used in determining gene clusters
##' @param secondary.correlation whether clustering should be performed on the correlation of the correlation matrix instead
##' @param n.cells number of cells to use for the randomly generated cluster lambda1 model
##' @param old.results optionally, pass old results just to plot the model without recalculating the stats
##'
##' @return a list containing the following fields:
##' \itemize{
##' \item{clusters} {a list of genes in each cluster values}
##' \item{xf} { extreme value distribution fit for the standardized lambda1 of a randomly generated pattern}
##' \item{tci} { index of a top cluster in each random iteration}
##' \item{cl.goc} {weighted PCA info for each real gene cluster}
##' \item{varm} {standardized lambda1 values for each randomly generated matrix cluster}
##' \item{clvlm} {a linear model describing dependency of the cluster lambda1 on a Tracy-Widom lambda1 expectation}
##' }
##'
##' @examples
##' \donttest{
##' data(pollen)
##' cd <- pollen
##' cd <- cd[,colSums(cd>0)>1.8e3]
##' cd <- cd[rowSums(cd)>10,]
##' cd <- cd[rowSums(cd>0)>5,]
##' knn <- knn.error.models(cd, k=ncol(cd)/4, n.cores=10, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
##' varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = FALSE)
##' clpca <- pagoda.gene.clusters(varinfo, trim=7.1/ncol(varinfo$mat), n.clusters=150, n.cores=10, plot=FALSE)
##' }
##'
##' @export
pagoda.gene.clusters <- function(varinfo, trim = 3.1/ncol(varinfo$mat), n.clusters = 150, n.samples = 60, cor.method = "p", n.internal.shuffles = 0, n.starts = 10, n.cores = detectCores(), verbose = 0, plot = FALSE, show.random = FALSE, n.components = 1, method = "ward.D", secondary.correlation = FALSE, n.cells = ncol(varinfo$mat), old.results = NULL) {

    smooth <- 0
    mat <- varinfo$mat
    matw <- varinfo$matw
    batch = varinfo$batch

    if(trim > 0) {
        mat <- winsorize.matrix(mat, trim = trim)
    }
    if(!is.null(batch)) {
        # center mat by batch
        mat <- weightedMatCenter(mat, matw, batch)
    }


    if(!is.null(old.results)) {
        if(verbose) { cat ("reusing old results for the observed clusters\n")}
        gcls <- old.results$clusters
        cl.goc <- old.results$cl.goc
    } else {
        if(verbose) { cat ("determining gene clusters ...")}
        # actual clusters
        vi<-which(abs(apply(mat, 1, function(x) sum(abs(diff(x))))) > 0)
        if(is.element("WGCNA", installed.packages()[, 1])) {
            gd <- as.dist(1-WGCNA::cor(t(mat)[, vi], method = cor.method, nThreads = n.cores))
        } else {
            gd <- as.dist(1-cor(t(mat)[, vi], method = cor.method))
        }

        if(secondary.correlation) {
            if(is.element("WGCNA", installed.packages()[, 1])) {
                gd <- as.dist(1-WGCNA::cor(as.matrix(gd), method = "p", nThreads = n.cores))
            } else {
                gd <- as.dist(1-cor(as.matrix(gd), method = "p"))
            }
        }

        if(is.element("fastcluster", installed.packages()[, 1])) {
            gcl <- fastcluster::hclust(gd, method = method)
        } else {
            gcl <- stats::hclust(gd, method = method)
        }
        gcll <- cutree(gcl, n.clusters)
        gcls <- tapply(rownames(mat)[vi], as.factor(gcll), I)
        names(gcls) <- paste("geneCluster", names(gcls), sep = ".")

        rm(gd, gcl)
        gc()

        # determine PC1 for the actual clusters
        if(verbose) { cat (" cluster PCA ...")}
        il <- tapply(vi, factor(gcll, levels = c(1:length(gcls))), I)
        cl.goc <- papply(il, function(ii) {
            xp <- bwpca(t(mat[ii, , drop = FALSE]), t(matw[ii, , drop = FALSE]), npcs = n.components, center = FALSE, nstarts = n.starts, smooth = smooth, n.shuffles = n.internal.shuffles)

            cs <- unlist(lapply(seq_len(ncol(xp$scores)), function(i) sign(cor(xp$scores[, i], colMeans(mat[ii, , drop = FALSE]*abs(xp$rotation[, i]))))))

            xp$scores <- t(t(xp$scores)*cs)
            xp$rotation <- t(t(xp$rotation)*cs)

            return(list(xp = xp, sd = xp$sd, n = length(ii)))
        }, n.cores = n.cores)
        names(cl.goc) <- paste("geneCluster", names(cl.goc), sep = ".")

        if(verbose) { cat ("done\n")}
    }

    # sampled variation
    if(!is.null(old.results) && !is.null(old.results$varm)) {
        if(verbose) { cat ("reusing old results for the sampled clusters\n")}
        varm <- old.results$varm } else {
            if(verbose) { cat ("generating", n.samples, "randomized samples ")}
            varm <- do.call(rbind, papply(seq_len(n.samples), function(i) { # each sampling iteration
                set.seed(i)
                # generate random normal matrix
                # TODO: use n.cells instead of ncol(matw)
                m <- matrix(rnorm(nrow(mat)*n.cells), nrow = nrow(mat), ncol = n.cells)
                #m <- weightedMatCenter(m, matw, batch = batch)

                if(show.random) {
                    full.m <- t(m) # save untrimmed version of m for random gene set controls
                }

                if(trim > 0) {
                    m <- winsorize.matrix(m, trim = trim)
                }

                vi<-which(abs(apply(m, 1, function(x) sum(diff(abs(x))))) > 0)
                if(is.element("WGCNA", installed.packages()[, 1])) {
                    gd <- as.dist(1-WGCNA::cor(t(m[vi, ]), method = cor.method, nThreads = 1))
                } else {
                    gd <- as.dist(1-cor(t(m[vi, ]), method = cor.method))
                }
                if(secondary.correlation) {
                    if(is.element("WGCNA", installed.packages()[, 1])) {
                        gd <- as.dist(1-WGCNA::cor(as.matrix(gd), method = "p", nThreads = 1))
                    } else {
                        gd <- as.dist(1-cor(as.matrix(gd), method = "p"))
                    }
                }

                if(is.element("fastcluster", installed.packages()[, 1])) {
                    gcl <- fastcluster::hclust(gd, method = method)
                } else {
                    gcl <- stats::hclust(gd, method = method)
                }
                gcll <- cutree(gcl, n.clusters)
                rm(gd, gcl)
                gc()

                # transpose to save time
                m <- t(m) # matw <- t(matw)

                sdv <- tapply(vi, gcll, function(ii) {
                    #as.numeric(bwpca(m[, ii], matw[, ii], npcs = 1, center = FALSE, nstarts = n.starts, smooth = smooth)$sd)^2
                    pcaMethods::sDev(pcaMethods::pca(m[, ii], nPcs = 1, center = FALSE))^2
                })

                pathsizes <- unlist(tapply(vi, gcll, length))
                names(pathsizes) <- pathsizes

                if(show.random) {
                    rsdv <- unlist(lapply(names(pathsizes), function(s) {
                        vi <- sample(1:ncol(full.m), as.integer(s))
                        pcaMethods::sDev(pcaMethods::pca(full.m[, vi], nPcs = 1, center = FALSE))^2
                    }))
                    if(verbose) { cat (".")}
                    return(data.frame(n = as.integer(pathsizes), var = unlist(sdv), round = i, rvar = rsdv))
                }

                if(verbose) { cat (".")}
                data.frame(n = as.integer(pathsizes), var = unlist(sdv), round = i)

            }, n.cores = n.cores))
            if(verbose) { cat ("done\n")}
        }

    # score relative to Tracey-Widom distribution
    #require(RMTstat)
    x <- RMTstat::WishartMaxPar(n.cells, varm$n)
    varm$pm <- x$centering-(1.2065335745820)*x$scaling # predicted mean of a random set
    varm$pv <- (1.607781034581)*x$scaling # predicted variance of a random set
    #clvlm <- lm(var~pm, data = varm)
    clvlm <- lm(var~0+pm+n, data = varm)
    varm$varst <- (varm$var-predict(clvlm))/sqrt(varm$pv)
    #varm$varst <- as.numeric(varm$var - (cbind(1, varm$pm) %*% coef(clvlm)))/sqrt(varm$pv)
    #varm$varst <- as.numeric(varm$var - (varm$pm* coef(clvlm)[2]))/sqrt(varm$pv)

    #varm$varst <- (varm$var-varm$pm)/sqrt(varm$pv)
    tci <- tapply(seq_len(nrow(varm)), as.factor(varm$round), function(ii) ii[which.max(varm$varst[ii])])

    #xf <- fevd(varm$varst[tci], type = "Gumbel") # fit on top clusters
    xf <- extRemes::fevd(varm$varst, type = "Gumbel") # fit on all clusters

    if(plot) {
        require(extRemes)
        par(mfrow = c(1, 2), mar = c(3.5, 3.5, 3.5, 1.0), mgp = c(2, 0.65, 0), cex = 0.9)
        smoothScatter(varm$n, varm$var, main = "simulations", xlab = "cluster size", ylab = "PC1 variance")
        if(show.random) {
            points(varm$n, varm$rvar, pch = ".", col = "red")
        }
        #pv <- predict(rsm, newdata = data.frame(n = sort(varm$n)), se.fit = TRUE)
        on <- order(varm$n, decreasing = TRUE)
        lines(varm$n[on], predict(clvlm)[on], col = 4, lty = 3)
        lines(varm$n[on], varm$pm[on], col = 2)
        lines(varm$n[on], (varm$pm+1.96*sqrt(varm$pv))[on], col = 2, lty = 2)
        lines(varm$n[on], (varm$pm-1.96*sqrt(varm$pv))[on], col = 2, lty = 2)
        legend(x = "bottomright", pch = c(1, 19, 19), col = c(1, 4, 2), legend = c("top clusters", "clusters", "random"), bty = "n")

        points(varm$n[tci], varm$var[tci], col = 1)
        extRemes::plot.fevd(xf, type = "density", main = "Gumbel fit")
        abline(v = 0, lty = 3, col = 4)
    }

    #pevd(9, loc = xf$results$par[1], scale = xf$results$par[2], lower.tail = FALSE)
    #xf$results$par

    return(list(clusters = gcls, xf = xf, tci = tci, cl.goc = cl.goc, varm = varm, clvlm = clvlm, trim = trim))
}


##' Score statistical significance of gene set and cluster overdispersion
##'
##' Evaluates statistical significance of the gene set and cluster lambda1 values, returning
##' either a text table of Z scores, etc, a structure containing normalized values of significant
##' aspects, or a set of genes underlying the significant aspects.
##'
##' @param pwpca output of pagoda.pathway.wPCA()
##' @param clpca output of pagoda.gene.clusters() (optional)
##' @param n.cells effective number of cells (if not provided, will be determined using pagoda.effective.cells())
##' @param z.score Z score to be used as a cutoff for statistically significant patterns (defaults to 0.05 P-value
##' @param return.table whether a text table showing
##' @param return.genes whether a set of genes driving significant aspects should be returned
##' @param plot whether to plot the cv/n vs. dataset size scatter showing significance models
##' @param adjust.scores whether the normalization of the aspect patterns should be based on the adjusted Z scores - qnorm(0.05/2, lower.tail = FALSE)
##' @param score.alpha significance level of the confidence interval for determining upper/lower bounds
##' @param use.oe.scale whether the variance of the returned aspect patterns should be normalized using observed/expected value instead of the default chi-squared derived variance corresponding to overdispersion Z score
##' @param effective.cells.start starting value for the pagoda.effective.cells() call
##'
##' @return if return.table = FALSE and return.genes = FALSE (default) returns a list structure containing the following items:
##' \itemize{
##' \item{xv} {a matrix of normalized aspect patterns (rows- significant aspects, columns- cells}
##' \item{xvw} { corresponding weight matrix }
##' \item{gw} { set of genes driving the significant aspects }
##' \item{df} { text table with the significance testing results }
##' }
##'
##' @examples
##' \donttest{
##' data(pollen)
##' cd <- pollen
##' cd <- cd[,colSums(cd>0)>1.8e3]
##' cd <- cd[rowSums(cd)>10,]
##' cd <- cd[rowSums(cd>0)>5,]
##' knn <- knn.error.models(cd, k=ncol(cd)/4, n.cores=10, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
##' varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = FALSE)
##' pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components=1, n.cores=10, n.internal.shuffles=50)
##' tam <- pagoda.top.aspects(pwpca, return.table = TRUE, plot=FALSE, z.score=1.96)  # top aspects based on GO only
##' }
##'
##' @export
pagoda.top.aspects <- function(pwpca, clpca = NULL, n.cells = NULL, z.score = qnorm(0.05/2, lower.tail = FALSE), return.table = FALSE, return.genes = FALSE, plot = FALSE, adjust.scores = TRUE, score.alpha = 0.05, use.oe.scale = FALSE, effective.cells.start = NULL) {
    basevar = 1

    if(is.null(n.cells)) {
        n.cells <- pagoda.effective.cells(pwpca, start = effective.cells.start)
    }


    vdf <- data.frame(do.call(rbind, lapply(seq_along(pwpca), function(i) {
        vars <- as.numeric((pwpca[[i]]$sd)^2)
        shz <- NA
        if(!is.null(pwpca[[i]]$xp$randvar)) { shz <- (vars - mean(pwpca[[i]]$xp$randvar))/sd(pwpca[[i]]$xp$randvar) }
        cbind(i = i, var = vars, n = pwpca[[i]]$n, npc = seq(1:ncol(pwpca[[i]]$xp$scores)), shz = shz)
    })))

    # fix p-to-q mistake in qWishartSpike
    qWishartSpikeFixed <- function (q, spike, ndf = NA, pdim = NA, var = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)  {
        params <- RMTstat::WishartSpikePar(spike, ndf, pdim, var, beta)
        qnorm(q, mean = params$centering, sd = params$scaling, lower.tail, log.p)
    }

    # add right tail approximation to ptw, which gives up quite early
    pWishartMaxFixed <- function (q, ndf, pdim, var = 1, beta = 1, lower.tail = TRUE) {
        params <- RMTstat::WishartMaxPar(ndf, pdim, var, beta)
        q.tw <- (q - params$centering)/(params$scaling)
        p <- RMTstat::ptw(q.tw, beta, lower.tail, log.p = TRUE)
        p[p == -Inf] <- pgamma((2/3)*q.tw[p == -Inf]^(3/2), 2/3, lower.tail = FALSE, log.p = TRUE) + lgamma(2/3) + log((2/3)^(1/3))
        p
    }


    #bi <- which.max(unlist(lapply(pwpca, function(x) x$n)))
    #vshift <- mean(pwpca[[bi]]$z[, 1]^2)/pwpca[[bi]]$n
    #ev <- ifelse(spike > 0, qWishartSpikeFixed(0.5, spike, n.cells, pwpca[[bi]]$n, var = basevar, lower.tail = FALSE), RMTstat::qWishartMax(0.5, n.cells, pwpca[[bi]]$n, var = basevar, lower.tail = FALSE))/pwpca[[bi]]$n
    #cat("vshift = ", vshift)
    vshift <- 0
    ev <- 0

    vdf$var <- vdf$var-(vshift-ev)*vdf$n

    #vdf$var[vdf$npc == 1] <- vdf$var[vdf$npc == 1]-(vshift-ev)*vdf$n[vdf$npc == 1]
    vdf$exp <- RMTstat::qWishartMax(0.5, n.cells, vdf$n, var = basevar, lower.tail = FALSE)
    #vdf$z <- qnorm(pWishartMax(vdf$var, n.cells, vdf$n, log.p = TRUE, lower.tail = FALSE, var = basevar), lower.tail = FALSE, log.p = TRUE)
    vdf$z <- qnorm(pWishartMaxFixed(vdf$var, n.cells, vdf$n, lower.tail = FALSE, var = basevar), lower.tail = FALSE, log.p = TRUE)
    vdf$cz <- qnorm(bh.adjust(pnorm(as.numeric(vdf$z), lower.tail = FALSE, log.p = TRUE), log = TRUE), lower.tail = FALSE, log.p = TRUE)
    vdf$ub <- RMTstat::qWishartMax(score.alpha/2, n.cells, vdf$n, var = basevar, lower.tail = FALSE)
    vdf$ub.stringent <- RMTstat::qWishartMax(score.alpha/nrow(vdf)/2, n.cells, vdf$n, var = basevar, lower.tail = FALSE)

    if(!is.null(clpca)) {
        clpca$xf <- extRemes::fevd(varst, data = clpca$varm, type = "Gumbel")
        #clpca$xf <- fevd(clpca$varm$varst[clpca$tci], type = "Gumbel")
        clpca$xf$results$par <- c(clpca$xf$results$par, c(shape = 0))
        #plot(xf)

        clvdf <- data.frame(do.call(rbind, lapply(seq_along(clpca$cl.goc), function(i)  {
            vars <- as.numeric((clpca$cl.goc[[i]]$sd)^2)
            shz <- NA
            if(!is.null(clpca$cl.goc[[i]]$xp$randvar)) {
                shz <- (vars - mean(clpca$cl.goc[[i]]$xp$randvar))/sd(clpca$cl.goc[[i]]$xp$randvar)
            }
            cbind(i = i, var = vars, n = clpca$cl.goc[[i]]$n, npc = seq(1:ncol(clpca$cl.goc[[i]]$xp$scores)), shz = shz)
        })))

        clvdf$var <- clvdf$var-(vshift-ev)*clvdf$n

        x <- RMTstat::WishartMaxPar(n.cells, clvdf$n)
        clvdf$pm <- x$centering-(1.2065335745820)*x$scaling # predicted mean of a random set
        clvdf$pv <- (1.607781034581)*x$scaling # predicted variance of a random set
        pvar <- predict(clpca$clvlm, newdata = clvdf)
        clvdf$varst <- (clvdf$var-pvar)/sqrt(clvdf$pv)
        clvdf$exp <- clpca$xf$results$par[1]*sqrt(clvdf$pv)+pvar
        #clvdf$varst <- (clvdf$var-clvdf$pm)/sqrt(clvdf$pv)
        #clvdf$exp <- clpca$xf$results$par[1]*sqrt(clvdf$pv)+clvdf$pm

        lp <- pgev.upper.log(clvdf$varst, clpca$xf$results$par[1], clpca$xf$results$par[2], rep(clpca$xf$results$par[3], nrow(clvdf)))
        clvdf$z <- qnorm(lp, lower.tail = FALSE, log.p = TRUE)
        clvdf$cz <- qnorm(bh.adjust(pnorm(as.numeric(clvdf$z), lower.tail = FALSE, log.p = TRUE), log = TRUE), lower.tail = FALSE, log.p = TRUE)

        # CI relative to the background
        clvdf$ub <- extRemes::qevd(score.alpha/2, loc = clpca$xf$results$par[1], scale = clpca$xf$results$par[2], shape = clpca$xf$results$par[3], lower.tail = FALSE)*sqrt(clvdf$pv) + pvar
        clvdf$ub.stringent <- extRemes::qevd(score.alpha/2/nrow(clvdf), loc = clpca$xf$results$par[1], scale = clpca$xf$results$par[2], shape = clpca$xf$results$par[3], lower.tail = FALSE)*sqrt(clvdf$pv) + pvar

    }

    if(plot) {
        par(mfrow = c(1, 1), mar = c(3.5, 3.5, 1.0, 1.0), mgp = c(2, 0.65, 0))
        un <- sort(unique(vdf$n))
        on <- order(vdf$n, decreasing = FALSE)
        pccol <- colorRampPalette(c("black", "grey70"), space = "Lab")(max(vdf$npc))
        plot(vdf$n, vdf$var/vdf$n, xlab = "gene set size", ylab = "PC1 var/n", ylim = c(0, max(vdf$var/vdf$n)), col = pccol[vdf$npc])
        lines(vdf$n[on], (vdf$exp/vdf$n)[on], col = 2, lty = 1)
        lines(vdf$n[on], (vdf$ub.stringent/vdf$n)[on], col = 2, lty = 2)

        if(!is.null(clpca)) {
            pccol <- colorRampPalette(c("darkgreen", "lightgreen"), space = "Lab")(max(clvdf$npc))
            points(clvdf$n, clvdf$var/clvdf$n, col = pccol[clvdf$npc], pch = 1)

            #clvm <- clpca$xf$results$par[1]*sqrt(pmax(1e-3, predict(vm, data.frame(n = un)))) + predict(mm, data.frame(n = un))
            on <- order(clvdf$n, decreasing = FALSE)

            lines(clvdf$n[on], (clvdf$exp/clvdf$n)[on], col = "darkgreen")
            lines(clvdf$n[on], (clvdf$ub.stringent/clvdf$n)[on], col = "darkgreen", lty = 2)
        }
        #mi<-which.max(vdf$n) sv<- (vdf$var/vdf$n)[mi] - (vdf$exp/vdf$n)[mi]
        #lines(vdf$n[on], (vdf$exp/vdf$n)[on]+sv, col = 2, lty = 3)
        #lines(vdf$n[on], (vdf$ub.stringent/vdf$n)[on]+sv, col = 2, lty = 2)
    }


    if(!is.null(clpca)) { # merge in cluster stats based on their own model

        # merge pwpca, psd and pm
        # all processing from here is common
        clvdf$i <- clvdf$i+length(pwpca) # shift cluster ids
        pwpca <- c(pwpca, clpca$cl.goc)
        vdf <- rbind(vdf, clvdf[, c("i", "var", "n", "npc", "exp", "cz", "z", "ub", "ub.stringent", "shz")])
    }

    vdf$adj.shz <- qnorm(bh.adjust(pnorm(as.numeric(vdf$shz), lower.tail = FALSE, log.p = TRUE), log = TRUE), lower.tail = FALSE, log.p = TRUE)
    #vdf$oe <- vdf$var/vdf$exp
    rs <- (vshift-ev)*vdf$n
    #rs <- ifelse(vdf$npc == 1, (vshift-ev)*vdf$n, 0)
    vdf$oe <- (vdf$var+rs)/(vdf$exp+rs)
    #vdf$oe[vdf$oe<0] <- 0
    #vdf$oec <- (vdf$var-vdf$ub.stringent+vdf$exp)/vdf$exp
    #vdf$oec <- (vdf$var-vdf$ub+vdf$exp)/vdf$exp
    #vdf$oec <- (vdf$var-vdf$ub+vdf$exp+rs)/(vdf$exp+rs)
    vdf$oec <- (vdf$var+rs)/(vdf$ub+rs)
    #vdf$oec[vdf$oec<0] <- 0
    #vdf$z[vdf$z<0] <- 0



    df <- data.frame(name = names(pwpca)[vdf$i], npc = vdf$npc, n = vdf$n, score = vdf$oe, z = vdf$z, adj.z = vdf$cz, sh.z = vdf$shz, adj.sh.z = vdf$adj.shz, stringsAsFactors = FALSE)
    if(adjust.scores) {
        vdf$valid <- vdf$cz  >=  z.score
    } else {
        vdf$valid <- vdf$z  >=  z.score
    }

    if(return.table) {
        df <- df[vdf$valid, ]
        df <- df[order(df$score, decreasing = TRUE), ]
        return(df)
    }

    # determine genes driving significant pathways
    # return genes within top 2/3rds of PC loading
    gl <- lapply(which(vdf$valid), function(i) { s <- abs(pwpca[[vdf[i, "i"]]]$xp$rotation[, vdf[i, "npc"]] )
    s[s >= max(s)/3] })
    gw <- tapply(abs(unlist(gl)), as.factor(unlist(lapply(gl, names))), max)
    if(return.genes) {
        return(gw)
    }
    # return combined data structure

    # weight
    xvw <- do.call(rbind, lapply(pwpca, function(x) {
        xm <- t(x$xp$scoreweights)
    }))
    vi <- vdf$valid
    xvw <- xvw[vi, ]/rowSums(xvw[vi, ])

    # return scaled patterns
    xmv <- do.call(rbind, lapply(pwpca, function(x) {
        xm <- t(x$xp$scores)
    }))

    if(use.oe.scale) {
        xmv <- (xmv[vi, ] -rowMeans(xmv[vi, ]))* (as.numeric(vdf$oe[vi])/sqrt(apply(xmv[vi, ], 1, var)))
    } else {
        # chi-squared
        xmv <- (xmv[vi, ]-rowMeans(xmv[vi, ])) * sqrt((qchisq(pnorm(vdf$z[vi], lower.tail = FALSE, log.p = TRUE), n.cells, lower.tail = FALSE, log.p = TRUE)/n.cells)/apply(xmv[vi, ], 1, var))
    }
    rownames(xmv) <- paste("#PC", vdf$npc[vi], "# ", names(pwpca)[vdf$i[vi]], sep = "")

    return(list(xv = xmv, xvw = xvw, gw = gw, df = df))

}


##' Collapse aspects driven by the same combinations of genes
##'
##' Examines PC loading vectors underlying the identified aspects and clusters aspects based
##' on a product of loading and score correlation (raised to corr.power). Clusters of aspects
##' driven by the same genes are determined based on the distance.threshold and collapsed.
##'
##' @param tam output of pagoda.top.aspects()
##' @param pwpca output of pagoda.pathway.wPCA()
##' @param clpca output of pagoda.gene.clusters() (optional)
##' @param plot whether to plot the resulting clustering
##' @param cluster.method one of the standard clustering methods to be used (fastcluster::hclust is used if available or stats::hclust)
##' @param distance.threshold similarity threshold for grouping interdependent aspects
##' @param corr.power power to which the product of loading and score correlation is raised
##' @param abs Boolean of whether to use absolute correlation
##' @param n.cores number of cores to use during processing
##' @param ... additional arguments are passed to the pagoda.view.aspects() method during plotting
##'
##' @return a list structure analogous to that returned by pagoda.top.aspects(), but with addition of a $cnam element containing a list of aspects summarized by each row of the new (reduced) $xv and $xvw
##'
##' @examples
##' \donttest{
##' data(pollen)
##' cd <- pollen
##' cd <- cd[,colSums(cd>0)>1.8e3]
##' cd <- cd[rowSums(cd)>10,]
##' cd <- cd[rowSums(cd>0)>5,]
##' knn <- knn.error.models(cd, k=ncol(cd)/4, n.cores=10, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
##' varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = FALSE)
##' pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components=1, n.cores=10, n.internal.shuffles=50)
##' tam <- pagoda.top.aspects(pwpca, return.table = TRUE, plot=FALSE, z.score=1.96)  # top aspects based on GO only
##' tamr <- pagoda.reduce.loading.redundancy(tam, pwpca)
##' }
##'
##' @export
pagoda.reduce.loading.redundancy <- function(tam, pwpca, clpca = NULL, plot = FALSE, cluster.method = "complete", distance.threshold = 0.01, corr.power = 4, n.cores = detectCores(), abs = TRUE, ...) {
    pclc <- pathway.pc.correlation.distance(c(pwpca, clpca$cl.goc), tam$xv, target.ndf = 100, n.cores = n.cores)
    cda <- cor(t(tam$xv))
    if(abs) {
        cda <- abs(cda)
    } else {
        cda[cda<0] <- 0
    }
    cda <- as.dist(1-cda)
    cc <- (1-sqrt((1-pclc)*(1-cda)))^corr.power

    if(is.element("fastcluster", installed.packages()[, 1])) {
        y <- fastcluster::hclust(cc, method = cluster.method)
    } else {
        y <- stats::hclust(cc, method = cluster.method)
    }
    ct <- cutree(y, h = distance.threshold)
    ctf <- factor(ct, levels = sort(unique(ct)))
    xvl <- collapse.aspect.clusters(tam$xv, tam$xvw, ct, pick.top = FALSE, scale = TRUE)

    if(plot) {
        sc <- sample(colors(), length(levels(ctf)), replace = TRUE)
        view.aspects(tam$xv, row.clustering = y, row.cols = sc[as.integer(ctf)], ...)
    }

    # collapsed names
    if(!is.null(tam$cnam)) { # already has collapsed names
        cnam <- tapply(rownames(tam$xv), ctf, function(xn) unlist(tam$cnam[xn]))
    } else {
        cnam <- tapply(rownames(tam$xv), ctf, I)
    }
    names(cnam) <- rownames(xvl$d)
    tam$xv <- xvl$d
    tam$xvw <- xvl$w
    tam$cnam <- cnam
    return(tam)
}


##' Collapse aspects driven by similar patterns (i.e. separate the same sets of cells)
##'
##' Examines PC loading vectors underlying the identified aspects and clusters aspects based on score correlation. Clusters of aspects driven by the same patterns are determined based on the distance.threshold.
##'
##' @param tamr output of pagoda.reduce.loading.redundancy()
##' @param distance.threshold similarity threshold for grouping interdependent aspects
##' @param cluster.method one of the standard clustering methods to be used (fastcluster::hclust is used if available or stats::hclust)
##' @param distance distance matrix
##' @param weighted.correlation Boolean of whether to use a weighted correlation in determining the similarity of patterns
##' @param plot Boolean of whether to show plot
##' @param top Restrict output to the top n aspects of heterogeneity
##' @param trim Winsorization trim to use prior to determining the top aspects
##' @param abs Boolean of whether to use absolute correlation
##' @param ... additional arguments are passed to the pagoda.view.aspects() method during plotting
##'
##' @return a list structure analogous to that returned by pagoda.top.aspects(), but with addition of a $cnam element containing a list of aspects summarized by each row of the new (reduced) $xv and $xvw
##'
##' @examples
##' \donttest{
##' data(pollen)
##' cd <- pollen
##' cd <- cd[,colSums(cd>0)>1.8e3]
##' cd <- cd[rowSums(cd)>10,]
##' cd <- cd[rowSums(cd>0)>5,]
##' knn <- knn.error.models(cd, k=ncol(cd)/4, n.cores=10, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
##' varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = FALSE)
##' pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components=1, n.cores=10, n.internal.shuffles=50)
##' tam <- pagoda.top.aspects(pwpca, return.table = TRUE, plot=FALSE, z.score=1.96)  # top aspects based on GO only
##' tamr <- pagoda.reduce.loading.redundancy(tam, pwpca)
##' tamr2 <- pagoda.reduce.redundancy(tamr, distance.threshold = 0.9, plot = TRUE, labRow = NA, labCol = NA, box = TRUE, margins = c(0.5, 0.5), trim = 0)
##' }
##'
##' @export
pagoda.reduce.redundancy <- function(tamr, distance.threshold = 0.2, cluster.method = "complete", distance = NULL, weighted.correlation = TRUE, plot = FALSE, top = Inf, trim = 0, abs = FALSE, ...) {
    if(is.null(distance)) {
        if(weighted.correlation) {
            distance <- .Call("matWCorr", t(tamr$xv), t(tamr$xvw), PACKAGE = "scde")
            rownames(distance) <- colnames(distance) <- rownames(tamr$xv)
            if(abs) {
                distance <- stats::as.dist(1-abs(distance), upper = TRUE)
            } else {
                distance <- stats::as.dist(1-distance, upper = TRUE)
            }
        } else {
            if(abs) {
                distance <- stats::as.dist(1-abs(cor(t(tamr$xv))))
            } else {
                distance <- stats::as.dist(1-cor(t(tamr$xv)))
            }
        }
    }
    if(is.element("fastcluster", installed.packages()[, 1])) {
        y <- fastcluster::hclust(distance, method = cluster.method)
    } else {
        y <- stats::hclust(distance, method = cluster.method)
    }

    ct <- cutree(y, h = distance.threshold)
    ctf <- factor(ct, levels = sort(unique(ct)))
    xvl <- collapse.aspect.clusters(tamr$xv, tamr$xvw, ct, pick.top = FALSE, scale = TRUE)

    if(plot) {
        sc <- sample(colors(), length(levels(ctf)), replace = TRUE)
        view.aspects(tamr$xv, row.clustering = y, row.cols = sc[as.integer(ctf)], ...)
    }

    # collapsed names
    if(!is.null(tamr$cnam)) { # already has collapsed names
        cnam <- tapply(rownames(tamr$xv), ctf, function(xn) unlist(tamr$cnam[xn]))
    } else {
        cnam <- tapply(rownames(tamr$xv), ctf, I)
    }
    names(cnam) <- rownames(xvl$d)

    if(trim > 0) { xvl$d <- winsorize.matrix(xvl$d, trim) } # trim prior to determining the top sets

    rcmvar <- apply(xvl$d, 1, var)
    vi <- order(rcmvar, decreasing = TRUE)[1:min(length(rcmvar), top)]

    tamr2 <- tamr
    tamr2$xv <- xvl$d[vi, ]
    tamr2$xvw <- xvl$w[vi, ]
    tamr2$cnam <- cnam[vi]
    return(tamr2)
}


##' Determine optimal cell clustering based on the genes driving the significant aspects
##'
##' Determines cell clustering (hclust result) based on a weighted correlation of genes
##' underlying the top aspects of transcriptional heterogeneity. Branch orientation is optimized
##' if 'cba' package is installed.
##'
##' @param tam result of pagoda.top.aspects() call
##' @param varinfo result of pagoda.varnorm() call
##' @param method clustering method ('ward.D' by default)
##' @param verbose 0 or 1 depending on level of desired verbosity
##' @param include.aspects whether the aspect patterns themselves should be included alongside with the individual genes in calculating cell distance
##' @param return.details Boolean of whether to return just the hclust result or a list containing the hclust result plus the distance matrix and gene values
##'
##' @return hclust result
##'
##' @examples
##' \donttest{
##' data(pollen)
##' cd <- pollen
##' cd <- cd[,colSums(cd>0)>1.8e3]
##' cd <- cd[rowSums(cd)>10,]
##' cd <- cd[rowSums(cd>0)>5,]
##' knn <- knn.error.models(cd, k=ncol(cd)/4, n.cores=10, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
##' varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = FALSE)
##' pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components=1, n.cores=10, n.internal.shuffles=50)
##' tam <- pagoda.top.aspects(pwpca, return.table = TRUE, plot=FALSE, z.score=1.96)  # top aspects based on GO only
##' hc <- pagoda.cluster.cells(tam, varinfo)
##' plot(hc)
##' }
##'
##' @export
pagoda.cluster.cells <- function(tam, varinfo, method = "ward.D", include.aspects = FALSE, verbose = 0, return.details = FALSE) {
    # gene clustering
    gw <- tam$gw
    gw <- gw[(rowSums(varinfo$matw)*varinfo$arv)[names(gw)] > 1]

    gw <- gw/gw
    mi <- match(names(gw), rownames(varinfo$mat))
    wgm <- varinfo$mat[mi, ]
    wgm <- wgm*as.numeric(gw)
    wgwm <- varinfo$matw[mi, ]

    if(include.aspects) {
        if(verbose) { message("clustering cells based on ", nrow(wgm), " genes and ", nrow(tam$xv), " aspect patterns")}
        wgm <- rbind(wgm, tam$xv)
        wgwm <- rbind(wgwm, tam$xvw)
    } else {
        if(verbose) { message("clustering cells based on ", nrow(wgm), " genes")}
    }

    snam <- sample(colnames(wgm))

    dm <- .Call("matWCorr", wgm, wgwm, PACKAGE = "scde")
    dm <- 1-dm
    rownames(dm) <- colnames(dm) <- colnames(wgm)
    wcord <- stats::as.dist(dm, upper = TRUE)
    hc <- hclust(wcord, method = method)

    if(is.element("cba", installed.packages()[, 1])) {
        co <- cba::order.optimal(wcord, hc$merge)
        hc$merge <- co$merge
        hc$order <- co$order
    }
    if(return.details) {
        return(list(clustering = hc, distance = wcord, genes = gw))
    } else {
        return(hc)
    }
}


##' View PAGODA output
##'
##' Create static image of PAGODA output visualizing cell hierarchy and top aspects of transcriptional heterogeneity
##'
##' @param tamr Combined pathways that show similar expression patterns. Output of \code{\link{pagoda.reduce.redundancy}}
##' @param row.clustering Dendrogram of combined pathways clustering
##' @param top Restrict output to the top n aspects of heterogeneity
##' @param ... additional arguments are passed to the \code{\link{view.aspects}} method during plotting
##'
##' @return PAGODA heatmap
##'
##' @examples
##' \donttest{
##' data(pollen)
##' cd <- pollen
##' cd <- cd[,colSums(cd>0)>1.8e3]
##' cd <- cd[rowSums(cd)>10,]
##' cd <- cd[rowSums(cd>0)>5,]
##' knn <- knn.error.models(cd, k=ncol(cd)/4, n.cores=10, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
##' varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = FALSE)
##' pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components=1, n.cores=10, n.internal.shuffles=50)
##' tam <- pagoda.top.aspects(pwpca, return.table = TRUE, plot=FALSE, z.score=1.96)  # top aspects based on GO only
##' pagoda.view.aspects(tam)
##' }
##'
##' @export
pagoda.view.aspects <- function(tamr, row.clustering = hclust(dist(tamr$xv)), top = Inf, ...) {
    if(is.finite(top)) {
        rcmvar <- apply(tamr$xv, 1, var)
        vi <- order(rcmvar, decreasing = TRUE)[1:min(length(rcmvar), top)]
        tamr$xv <- tamr$xv[vi, ]
        tamr$xvw <- tamr$xvw[vi, ]
        tamr$cnam <- tamr$cnam[vi]
    }

    view.aspects(tamr$xv, row.clustering = row.clustering, ... )
}


##' View heatmap
##'
##' Internal function to visualize aspects of transcriptional heterogeneity as a heatmap. Used by \code{\link{pagoda.view.aspects}}.
##'
##' @param mat Numeric matrix
##' @param row.clustering Row dendrogram
##' @param cell.clustering Column dendrogram
##' @param zlim Range of the normalized gene expression levels, inputted as a list: c(lower_bound, upper_bound). Values outside this range will be Winsorized. Useful for increasing the contrast of the heatmap visualizations. Default, set to the 5th and 95th percentiles.
##' @param row.cols  Matrix of row colors.
##' @param col.cols  Matrix of column colors. Useful for visualizing cell annotations such as batch labels.
##' @param cols Heatmap colors
##' @param show.row.var.colors Boolean of whether to show row variance as a color track
##' @param top Restrict output to the top n aspects of heterogeneity
##' @param ... additional arguments for heatmap plotting
##'
##' @return A heatmap
##'
view.aspects <- function(mat, row.clustering = NA, cell.clustering = NA, zlim = c(-1, 1)*quantile(mat, p = 0.95), row.cols = NULL, col.cols = NULL, cols = colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(1024), show.row.var.colors = TRUE, top = Inf, ...) {
    #row.cols, col.cols are matrices for now
    rcmvar <- apply(mat, 1, var)
    mat[mat<zlim[1]] <- zlim[1]
    mat[mat > zlim[2]] <- zlim[2]
    if(class(row.clustering) == "hclust") { row.clustering <- as.dendrogram(row.clustering) }
    if(class(cell.clustering) == "hclust") { cell.clustering <- as.dendrogram(cell.clustering) }
    if(show.row.var.colors) {
        if(is.null(row.cols)) {
            icols <- colorRampPalette(c("white", "black"), space = "Lab")(1024)[1023*(rcmvar/max(rcmvar))+1]
            row.cols <- cbind(var = icols)
        }
    }
    my.heatmap2(mat, Rowv = row.clustering, Colv = cell.clustering, zlim = zlim, RowSideColors = row.cols, ColSideColors = col.cols, col = cols, ...)
}


##' Make the PAGODA app
##'
##' Create an interactive user interface to explore output of PAGODA.
##'
##' @param tamr Combined pathways that show similar expression patterns. Output of \code{\link{pagoda.reduce.redundancy}}
##' @param tam Combined pathways that are driven by the same gene sets. Output of \code{\link{pagoda.reduce.loading.redundancy}}
##' @param varinfo Variance information. Output of \code{\link{pagoda.varnorm}}
##' @param env Gene sets as an environment variable.
##' @param pwpca Weighted PC magnitudes for each gene set provided in the \code{env}. Output of \code{\link{pagoda.pathway.wPCA}}
##' @param clpca Weighted PC magnitudes for de novo gene sets identified by clustering on expression. Output of \code{\link{pagoda.gene.clusters}}
##' @param col.cols  Matrix of column colors. Useful for visualizing cell annotations such as batch labels. Default NULL.
##' @param cell.clustering Dendrogram of cell clustering. Output of \code{\link{pagoda.cluster.cells} } . Default   NULL.
##' @param row.clustering Dendrogram of combined pathways clustering. Default NULL.
##' @param title Title text to be used in the browser label for the app. Default, set as 'pathway clustering'
##' @param zlim Range of the normalized gene expression levels, inputted as a list: c(lower_bound, upper_bound). Values outside this range will be Winsorized. Useful for increasing the contrast of the heatmap visualizations. Default, set to the 5th and 95th percentiles.
##'
##' @return PAGODA app
##'
##' @export
make.pagoda.app <- function(tamr, tam, varinfo, env, pwpca, clpca = NULL, col.cols = NULL, cell.clustering = NULL, row.clustering = NULL, title = "pathway clustering", zlim = c(-1, 1)*quantile(tamr$xv, p = 0.95)) {
    # rcm - xv
    # matvar
    if(is.null(cell.clustering)) {
        cell.clustering <- pagoda.cluster.cells(tam, varinfo)
    }
    if(is.null(row.clustering)) {
        row.clustering <- hclust(dist(tamr$xv))
        row.clustering$order <- rev(row.clustering$order)
    }

    #fct - which tam row in which tamr$xv cluster.. remap tamr$cnams
    cn <- tamr$cnam
    fct <- rep(1:length(cn), lapply(cn, length))
    names(fct) <- unlist(cn)
    fct <- fct[rownames(tam$xv)]
    rcm <- tamr$xv
    rownames(rcm) <- as.character(1:nrow(rcm))
    fres <- list(hvc = cell.clustering, tvc = row.clustering, rcm = rcm, zlim2 = zlim, matvar = apply(tam$xv, 1, sd), ct = fct, matrcmcor = rep(1, nrow(tam$xv)), cols = colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(1024), colcol = col.cols)

    # gene df
    gene.df <- data.frame(var = varinfo$arv*rowSums(varinfo$matw))
    gene.df$gene <- rownames(varinfo$mat)
    gene.df <- gene.df[order(gene.df$var, decreasing = TRUE), ]

    # prepare pathway df
    df <- tamr$df
    if(exists("myGOTERM", envir = globalenv())) {
        df$desc <- mget(df$name, get("myGOTERM", envir = globalenv()), ifnotfound = "")
    } else {
        df$desc <- ""
    }
    min.z <- -9
    df$z[df$z<min.z] <- min.z
    df$adj.z[df$adj.z<min.z] <- min.z
    df$sh.z[df$sh.z<min.z] <- min.z
    df$adj.sh.z[df$adj.sh.z<min.z] <- min.z
    df <- data.frame(id = paste("#PC", df$npc, "# ", df$name, sep = ""), npc = df$npc, n = df$n, score = df$score, Z = df$z, aZ = df$adj.z, sh.Z = df$sh.z, sh.aZ = df$adj.sh.z, name = paste(df$name, df$desc))

    df <- df[order(df$score, decreasing = TRUE), ]

    # merge go.env
    if(!is.null(clpca)) {
        set.env <- list2env(c(as.list(env), clpca$clusters))
    } else {
        set.env <- env
    }
    sa <- ViewPagodaApp$new(fres, df, gene.df, varinfo$mat, varinfo$matw, set.env, name = title, trim = 0, batch = varinfo$batch)
}

##################### Internal functions

one.sided.test.id <- function(id, nam1, nam2, ifm, dm, prior, difference.prior = 0.5, bootstrap = TRUE, n.samples = 1e3, show.plots = TRUE, return.posterior = FALSE, return.both = FALSE) {
    gr <- 10^prior$x - 1
    gr[gr<0] <- 0
    lpp <- get.rep.set.general.model.logposteriors(ifm[[nam1]], dm[rep(id, length(gr)), names(ifm[[nam1]])], data.frame(fpm = gr), grid.weight = prior$grid.weight)
    ldp <- get.rep.set.general.model.logposteriors(ifm[[nam2]], dm[rep(id, length(gr)), names(ifm[[nam2]])], data.frame(fpm = gr), grid.weight = prior$grid.weight)

    if(bootstrap) {
        pjp <- do.call(cbind, lapply(seq_along(n.samples), function(i) {
            pjp <- rowSums(lpp[, sample(1:ncol(lpp), replace = TRUE)])
            pjp <- exp(pjp-max(pjp))
            pjp <- pjp/sum(pjp)
            return(pjp)
        }))
        pjp <- rowSums(pjp)
        pjp <- log(pjp/sum(pjp))

        djp <- do.call(cbind, lapply(seq_along(n.samples), function(i) {
            djp <- rowSums(ldp[, sample(1:ncol(ldp), replace = TRUE)])
            djp <- exp(djp-max(djp))
            djp <- djp/sum(djp)
            return(djp)
        }))
        djp <- rowSums(djp)
        djp <- log(djp/sum(djp))
    } else {
        pjp <- rowSums(lpp)
        djp <- rowSums(ldp)
    }

    dpy <- exp(prior$lp+djp)
    mpgr <- sum(exp(prior$lp+pjp+log(c(0, cumsum(dpy)[-length(dpy)])))) # m1
    mpls <- sum(exp(prior$lp+pjp+log(sum(dpy)-cumsum(dpy)))) # m0
    mpls/mpgr

    pjpc <- exp(prior$lp+pjp)
    pjpc <- pjpc/sum(pjpc)
    djpc <- exp(prior$lp+djp)
    djpc <- djpc/sum(djpc)

    if(show.plots || return.posterior || return.both) {
        # calculate log-fold-change posterior
        n <- length(pjpc)
        rp <- c(unlist(lapply(n:2, function(i) sum(pjpc[1:(n-i+1)]*djpc[i:n]))), unlist(lapply(seq_along(n), function(i) sum(pjpc[i:n]*djpc[1:(n-i+1)]))))
        rv <- seq(prior$x[1]-prior$x[length(prior$x)], prior$x[length(prior$x)]-prior$x[1], length = length(prior$x)*2-1)
        fcp <- data.frame(v = rv, p = rp)
    }

    if(show.plots) {
        # show each posterior
        layout(matrix(c(1:3), 3, 1, byrow = TRUE), heights = c(2, 1, 2), widths = c(1), FALSE)
        par(mar = c(2.5, 3.5, 2.5, 3.5), mgp = c(1.5, 0.65, 0), cex = 0.9)
        jpr <- range(c(0, pjpc), na.rm = TRUE)
        pp <- exp(lpp)
        cols <- rainbow(dim(pp)[2], s = 0.8)
        plot(c(), c(), xlim = range(prior$x), ylim = range(c(0, pp)), xlab = "expression level", ylab = "individual posterior", main = nam1)
        lapply(seq_len(ncol(pp)), function(i) lines(prior$x, pp[, i], col = cols[i]))
        legend(x = ifelse(which.max(na.omit(pjpc)) > length(pjpc)/2, "topleft", "topright"), bty = "n", col = cols, legend = colnames(pp), lty = rep(1, dim(pp)[2]))
        par(new = TRUE)
        plot(prior$x, pjpc, axes = FALSE, ylab = "", xlab = "", ylim = jpr, type = 'l', col = 1, lty = 1, lwd = 2)
        axis(4, pretty(jpr, 5), col = 1)
        mtext("joint posterior", side = 4, outer = FALSE, line = 2)

        # ratio plot
        par(mar = c(2.5, 3.5, 0.5, 3.5), mgp = c(1.5, 0.65, 0), cex = 0.9)
        plot(fcp$v, fcp$p, xlab = "log10 expression ratio", ylab = "ratio posterior", type = 'l', lwd = 2, main = "")
        r.mle <- fcp$v[which.max(fcp$p)]
        r.lb <- max(which(cumsum(fcp$p)<0.025))
        r.ub <- min(which(cumsum(fcp$p) > (1-0.025)))
        polygon(c(fcp$v[r.lb], fcp$v[r.lb:r.ub], fcp$v[r.ub]), y = c(-10, fcp$p[r.lb:r.ub], -10), col = "grey90")
        abline(v = r.mle, col = 2, lty = 2)
        abline(v = c(fcp$v[r.ub], fcp$v[r.lb]), col = 2, lty = 3)
        box()
        legend(x = ifelse(r.mle > 0, "topleft", "topright"), legend = c(paste("MLE = ", round(10^(r.mle), 1), " (", round(r.mle, 2), " in log10)", sep = ""), paste("95% CI: ", round(10^(fcp$v[r.lb]), 1), " - ", round(10^(fcp$v[r.ub]), 1), sep = ""), paste(" log10: ", round(fcp$v[r.lb], 2), " - ", round(fcp$v[r.ub], 2), sep = "")), bty = "n")

        # distal plot
        dp <- exp(ldp)
        par(mar = c(2.5, 3.5, 2.5, 3.5), mgp = c(1.5, 0.65, 0), cex = 0.9)
        jpr <- range(c(0, djpc), na.rm = TRUE)
        cols <- rainbow(dim(dp)[2], s = 0.8)
        plot(c(), c(), xlim = range(prior$x), ylim = range(c(0, dp)), xlab = "expression level", ylab = "individual posterior", main = nam2)
        lapply(seq_len(ncol(dp)), function(i) lines(prior$x, dp[, i], col = cols[i]))
        legend(x = ifelse(which.max(na.omit(djpc)) > length(djpc)/2, "topleft", "topright"), bty = "n", col = cols, legend = colnames(dp), lty = rep(1, dim(dp)[2]))

        par(new = TRUE)
        plot(prior$x, djpc, axes = FALSE, ylab = "", xlab = "", ylim = jpr, type = 'l', col = 1, lty = 1, lwd = 2)
        axis(4, pretty(jpr, 5), col = 1)
        mtext("joint posterior", side = 4, outer = FALSE, line = 2)
    }

    lbf <- mpls/mpgr
    lbf <- (difference.prior*lbf)/(difference.prior*lbf+1-difference.prior)
    #return(c(equal = qnorm(ebf, lower.tail = TRUE), less = qnorm(lbf, lower.tail = TRUE)))
    if(return.both) {
        return(list(z = qnorm(lbf, lower.tail = TRUE), post = fcp))
    } else if(return.posterior) {
        return(fcp)
    } else {
        return(qnorm(lbf, lower.tail = TRUE))
    }
}

# counts - data frame with fragment counts (rows - fragments columns -experiments)
# groups - a two-level factor describing grouping of columns. Use NA for observations that should be skipped
# min.count.threshold - the number of reads used to make an initial guess for the failed component
# threshold.segmentation - use min.count.threshold to perform very quick identification of the drop-outs
# threshold.prior - prior that should be associated with threshold segmentation
calculate.crossfit.models <- function(counts, groups, min.count.threshold = 4, nrep = 1, verbose = 0, min.prior = 1e-5, n.cores = 12, save.plots = TRUE, zero.lambda = 0.1, old.cfm = NULL, threshold.segmentation = FALSE, threshold.prior = 1-1e-6, max.pairs = 1000, min.pairs.per.cell = 10) {
    names(groups) <- colnames(counts)
    # enumerate cross-fit pairs within each group
    cl <- do.call(cbind, tapply(colnames(counts), groups, function(ids) {
        cl <- combn(ids, 2)
        min.pairs.per.cell <- min(length(ids)*(length(ids)-1)/2, min.pairs.per.cell)
        if(verbose) {
            cat("number of pairs: ", ncol(cl), "\n")
        }
        if(ncol(cl) > max.pairs) {
            if(verbose) {
                cat("reducing to a random sample of ", max.pairs, " pairs\n")
            }

            # make sure there's at least min.pairs.per.cell pairs for each cell
            cl <- cl[, unique(c(sample(1:ncol(cl), max.pairs),
                                unlist(lapply(ids, function(id) sample(which(colSums(cl == id) > 0), min.pairs.per.cell)))))]
        }
        cl
    }))

    orl <- c()
    if(!is.null(old.cfm)) {
        # check which pairs have already been fitted in compared in old.cfm
        pn1 <- unlist(apply(cl, 2, function(ii) paste(ii, collapse = ".vs.")))
        pn2 <- unlist(apply(cl, 2, function(ii) paste(rev(ii), collapse = ".vs."))) ### %%% use rev() to revert element order
        vi <- (pn1 %in% names(old.cfm)) | (pn2 %in% names(old.cfm))
        cl <- cl[, !vi, drop = FALSE]
        orl <- old.cfm[names(old.cfm) %in% c(pn1, pn2)]
    }
    if(verbose) {
        cat("total number of pairs: ", ncol(cl), "\n")
    }

    if(dim(cl)[2] > 0) {
        if(verbose)  message(paste("cross-fitting", ncol(cl), "pairs:"))
        rl <- papply(seq_len(ncol(cl)), function(cii) {
            ii <- cl[, cii]
            df <- data.frame(c1 = counts[, ii[1]], c2 = counts[, ii[2]])
            vi <- which(rowSums(df) > 0, )
            if(!threshold.segmentation) {
                if(verbose) {
                    message("fitting pair [", paste(ii, collapse = " "), "]")
                }
                mo1 <- FLXMRglmCf(c1~1, family = "poisson", components = c(1), mu = log(zero.lambda))
                mo2 <- FLXMRnb2glmC(c1~1+I(log(c2+1)), components = c(2))
                mo3 <- FLXMRnb2glmC(c2~1+I(log(c1+1)), components = c(2))
                mo4 <- FLXMRglmCf(c2~1, family = "poisson", components = c(3), mu = log(zero.lambda))
                m1 <- mc.stepFlexmix(c1~1, data = df[vi, ], k = 3, model = list(mo1, mo2, mo3, mo4), control = list(verbose = verbose, minprior = min.prior), concomitant = FLXPmultinom(~I((log(c1+1)+log(c2+1))/2)+1), cluster = cbind(df$c1[vi]<= min.count.threshold, df$c1[vi] > min.count.threshold & df$c2[vi] > min.count.threshold, df$c2[vi]<= min.count.threshold), nrep = nrep)

                # reduce return size
                m1@posterior <- lapply(m1@posterior, function(m) {
                    rownames(m) <- NULL
                    return(m)
                })
                #rownames(m1@concomitant@x) <- NULL
                m1@concomitant@x <- matrix()
                m1@model <- lapply(m1@model, function(mod) {
                    mod@x <- matrix()
                    mod@y <- matrix()
                    #rownames(mod@x) <- NULL
                    #rownames(mod@y) <- NULL
                    return(mod)
                })

                #parent.env(environment(m1@components[[1]][[1]]@logLik)) <- globalenv()
                #parent.env(environment(m1@components[[1]][[2]]@logLik)) <- globalenv()
                #parent.env(environment(m1@components[[2]][[1]]@logLik)) <- globalenv()
                #parent.env(environment(m1@components[[2]][[2]]@logLik)) <- globalenv()

                names(vi) <- NULL
                pm <- posterior(m1)[, c(1, 3)]
                rownames(pm) <- NULL
                cl <- clusters(m1)
                names(cl) <- NULL
                gc()
            } else {
                # use min.count.threshold to quickly segment the points
                cl <- rep(2, length(vi))
                cl[df[vi, 1]<min.count.threshold] <- 1
                cl[df[vi, 2]<min.count.threshold] <- 3
                cl[df[vi, 1]<min.count.threshold & df[vi, 2]<min.count.threshold] <- 0
                names(cl) <- NULL
                pm <- cbind(ifelse(cl == 1, threshold.prior, 1-threshold.prior), ifelse(cl == 3, threshold.prior, 1-threshold.prior))
                rownames(pm) <- NULL
            }
            rli <- list(ii = ii, clusters = cl, posterior = pm, vi = vi)
            #message("return object size for pair [", paste(ii, collapse = " "), "] is ", round(object.size(rli)/(1024^2), 3), " MB")
            return(rli)
        }, n.cores = round(n.cores/nrep))
        #, mc.preschedule = threshold.segmentation) # mclapply function has preschedule
        names(rl) <- apply(cl, 2, paste, collapse = ".vs.")
        # clean up invalid entries
        rl <- rl[!unlist(lapply(rl, is.null))]
        rl <- rl[unlist(lapply(rl, is.list))]
        #names(rl) <- unlist(lapply(rl, function(d) paste(d$ii, collapse = ".vs.")))
    } else {
        rl <- c()
    }

    if(!is.null(old.cfm)) rl <- c(rl, orl)

    if(save.plots) {
        #require(Cairo) require(RColorBrewer)
        tapply(colnames(counts), groups, function(ids) {
            cl <- combn(ids, 2)
            group <- as.character(groups[ids[1]])
            # log-scale hist
            t.pairs.panel.hist <- function(x, i = NULL, ...) {
                usr <- par("usr")
                on.exit(par(usr))
                par(usr = c(usr[1:2], 0, 1.5) )
                vi <- x > 0
                h <- hist(x, plot = FALSE)
                breaks <- h$breaks
                nB <- length(breaks)
                y <- log10(h$counts)
                y <- y/max(y)
                rect(breaks[-nB], 0, breaks[-1], y, col = "gray60", ...)
            }
            t.pairs.smoothScatter.spearman <- function(x, y, i = NULL, j = NULL, cex = 0.8, ...) {
                vi <- x > 0 | y > 0
                smoothScatter(x[vi], y[vi], add = TRUE, useRaster = TRUE, ...)
                legend(x = "bottomright", legend = paste("sr = ", round(cor(x[vi], y[vi], method = "spearman"), 2), sep = ""), bty = "n", cex = cex)
            }
            # component assignment scatter
            t.panel.component.scatter <- function(x, y, i, j, cex = 0.8, ...) {
                if(!is.null(rl[[paste(ids[i], "vs", ids[j], sep = ".")]])) {
                    m1 <- rl[[paste(ids[i], "vs", ids[j], sep = ".")]]
                    # forward plot
                    vi <- which(x > 0 | y > 0)
                    ci <- vi[m1$clusters == 1]
                    if(length(ci) > 3) {
                        points(x[ci], y[ci], pch = ".", col = densCols(x[ci], y[ci], colramp = colorRampPalette(brewer.pal(9, "Reds")[-(1:3)])), cex = 2)
                    }

                    ci <- vi[m1$clusters == 3]
                    if(length(ci) > 3) {
                        points(x[ci], y[ci], pch = ".", col = densCols(x[ci], y[ci], colramp = colorRampPalette(brewer.pal(9, "Greens")[-(1:3)])), cex = 2)
                    }
                    ci <- vi[m1$clusters == 2]
                    if(length(ci) > 3) {
                        points(x[ci], y[ci], pch = ".", col = densCols(x[ci], y[ci], colramp = colorRampPalette(brewer.pal(9, "Blues")[-(1:3)])), cex = 2)
                    }
                    legend(x = "topleft", pch = c(19), col = "blue", legend = paste("sr = ", round(cor(x[ci], y[ci], method = "spearman"), 2), sep = ""), bty = "n", cex = cex)
                    legend(x = "bottomright", pch = c(rep(19, 3)), col = c("red", "blue", "green"), legend = paste(round(unlist(tapply(m1$clusters, factor(m1$clusters, levels = c(1, 2, 3)), length))*100/length(vi), 1), "%", sep = ""), bty = "n", cex = cex)

                } else if(!is.null(rl[[paste(ids[i], "vs", ids[j], sep = ".")]])) {
                    m1 <- rl[[paste(ids[j], "vs", ids[i], sep = ".")]]
                    # reverse plot
                    vi <- which(x > 0 | y > 0)
                    ci <- vi[m1$clusters == 3]
                    if(length(ci) > 3) {
                        points(x[ci], y[ci], pch = ".", col = densCols(x[ci], y[ci], colramp = colorRampPalette(brewer.pal(9, "Reds")[-(1:3)])), cex = 2)
                    }

                    ci <- vi[m1$clusters == 1]
                    if(length(ci) > 3) {
                        points(x[ci], y[ci], pch = ".", col = densCols(x[ci], y[ci], colramp = colorRampPalette(brewer.pal(9, "Greens")[-(1:3)])), cex = 2)
                    }
                    ci <- vi[m1$clusters == 2]
                    if(length(ci) > 3) {
                        points(x[ci], y[ci], pch = ".", col = densCols(x[ci], y[ci], colramp = colorRampPalette(brewer.pal(9, "Blues")[-(1:3)])), cex = 2)
                    }
                    legend(x = "topleft", pch = c(19), col = "blue", legend = paste("sr = ", round(cor(x[ci], y[ci], method = "spearman"), 2), sep = ""), bty = "n", cex = cex)
                    legend(x = "bottomright", pch = c(rep(19, 3)), col = c("red", "blue", "green"), legend = paste(round(unlist(tapply(m1$clusters, factor(m1$clusters, levels = c(3, 2, 1)), length))*100/length(vi), 1), "%", sep = ""), bty = "n", cex = cex)
                } else {
                    #message(paste("ERROR: unable to find model for i = ", i, "j = ", j))
                    message(paste("INFO: cross-fit plots: skipping model for i = ", i, "j = ", j, " (increase max.pairs parameter if needed"))
                }
            }
            #pdf(file = paste(group, "crossfits.pdf", sep = "."), width = 3*length(ids), height = 3*length(ids))
            CairoPNG(filename = paste(group, "crossfits.png", sep = "."), width = 250*length(ids), height = 250*length(ids))
            pairs.extended(log10(counts[, ids]+1), lower.panel = t.pairs.smoothScatter.spearman, upper.panel = t.panel.component.scatter, diag.panel = t.pairs.panel.hist, cex = 1.5)
            dev.off()
        })
    }

    return(rl)
}

# estimates library sizes based on the correlated components
# min.size.entries - minimal number of entries (genes) used to determine scaling factors for individual experiments
# counts - data frame with fragment counts (rows - fragments columns -experiments)
# groups - a two-level factor describing grouping of columns. Use NA for observations that should be skipped
# cfm - cross-fit models (return of calculate.crossfit.models())
# vil - optional binary matrix (corresponding to counts) with 0s marking likely drop-out observations
# return value - library size vector in millions of reads
estimate.library.sizes <- function(counts, cfm, groups, min.size.entries = min(nrow(counts), 2e3), verbose = 0, return.details = FALSE, vil = NULL, ...) {
    #require(edgeR)
    names(groups) <- colnames(counts)
    # determine the set fragments that were not attributed to failure in any cross-comparison
    if(is.null(vil)) {
        #x <- lapply(cfm, function(d) { ll <- list(!(1:nrow(counts)) %in% d$vi[which(d$clusters != 1)], !(1:nrow(counts)) %in% d$vi[which(d$clusters != 3)]) names(ll) <- d$ii return(ll) })
        x <- lapply(cfm, function(d) { ll <- list(!(1:nrow(counts)) %in% d$vi[which(d$clusters > 1)], !(1:nrow(counts)) %in% d$vi[which(d$clusters %% 3  != 0)])
        names(ll) <- d$ii
        return(ll)
        })
        vil <- do.call(cbind, tapply(unlist(x, recursive = FALSE), factor(unlist(lapply(x, names)), levels = colnames(counts)[!is.na(groups)]), function(l) {
            x <- rowSums(do.call(cbind, l), na.rm = FALSE) == 0
            x[is.na(x)] <- FALSE
            return(x)
        }))
    }

    # order entries by the number of non-failed experiments,
    # select entries for library size estimation
    ni <- cbind(1:nrow(counts), rowSums(vil))
    ni <- ni[order(ni[, 2], decreasing = TRUE), ]
    if(nrow(ni)<min.size.entries) {
        stop("The number of valid genes (", nrow(ni), ") is lower then the specified min.size.entries (", min.size.entries, "). Please either increase min.size.entries or lower min.nonfailed parameter to increase the number of valid genes")
    }
    if(ni[min.size.entries, 2]<ncol(vil)) {
        # if the min.size.entries -th gene has failures, take only min.size.entries genes
        gis <- ni[1:min.size.entries, 1]
    } else {
        # otherwise take all genes that have not failed in any experiment
        gis <- ni[ni[, 2] == ncol(vil), 1]
    }

    if(verbose)  message(paste("adjusting library size based on", length(gis), "entries"))
    f <- calcNormFactors(as.matrix(counts[gis, !is.na(groups)]), ...)
    f <- f/exp(mean(log(f)))
    ls <- colSums(counts[gis, !is.na(groups)])*f/1e6
    if(return.details) { return(list(ls = ls, vil = vil)) } else { return(ls) }
}

# an alternative prior estimation procedure that weights down contributions by failure probability
# and uses pre-scaled fpkm guesses for magnitude estimates
estimate.signal.prior <- function(fpkm, fail, length.out = 400, show.plot = FALSE, pseudo.count = 1, bw = 0.1, max.quantile = 0.999, max.value = NULL) {
    fpkm <- log10(exp(as.matrix(fpkm))+1)
    wts <- as.numeric(as.matrix(1-fail[, colnames(fpkm)]))
    wts <- wts/sum(wts)
    # fit density on a mirror image
    if(is.null(max.value)) {
        x <- as.numeric(fpkm)
        max.value <- as.numeric(quantile(x[x<Inf], p = max.quantile))
    }
    md <- density(c(-1*as.numeric(fpkm), as.numeric(fpkm)), bw = bw, weights = c(wts/2, wts/2), n = 2*length.out+1, from = -1*max.value, to = max.value)

    gep <- data.frame(x = md$x[-c(1:length.out)], y = md$y[-c(1:length.out)])
    gep$y[is.na(gep$y)] <- 0
    gep$y <- gep$y+pseudo.count/nrow(fpkm) # pseudo-count
    gep$y <- gep$y/sum(gep$y)
    if(show.plot) {
        par(mfrow = c(1, 1), mar = c(3.5, 3.5, 3.5, 0.5), mgp = c(2.0, 0.65, 0), cex = 0.9)
        plot(gep$x, gep$y, col = 4, panel.first = abline(h = 0, lty = 2), type = 'l', xlab = "log10( signal+1 )", ylab = "probability density", main = "signal prior")
    }
    gep$lp <- log(gep$y)

    # grid weighting (for normalization)
    gep$grid.weight <- diff(10^c(gep$x[1], gep$x+c(diff(gep$x)/2, 0))-1)

    return(gep)
    plot(x)
}

# counts - data frame with gene counts (rows - genes columns -experiments)
# groups - a two-level factor describing grouping of columns. Use NA for observations that should be skipped
# cfm - cross-fit models (return of calculate.crossfit.models())
# min.nonfailed - minimal number of non-failed observations required for a gene to be used in the final model fitting
#  A minimum of either the specified value or number of experiments -1 will be used.
calculate.individual.models <- function(counts, groups, cfm, nrep = 1, verbose = 0, n.cores = 12, min.nonfailed = 2, min.size.entries = 2e3, zero.count.threshold = 10, save.plots = TRUE, linear.fit = TRUE, return.compressed.models = FALSE,  local.theta.fit = FALSE, theta.fit.range = c(1e-2, 1e2), ...) {
    names(groups) <- colnames(counts)
    # determine library size discarding non-zero entries
    ls <- estimate.library.sizes(counts, cfm, groups, min.size.entries, verbose = verbose, return.details = TRUE)

    # fit three-component models to unique pairs within each group
    mll <- tapply(colnames(counts), groups, function(ids) {
        cl <- combn(ids, 2)
        group <- as.character(groups[ids[1]])

        # incorporate cross-fit pairs from cfm
        pn1 <- unlist(apply(cl, 2, function(ii) paste(ii, collapse = ".vs.")))
        pn2 <- unlist(apply(cl, 2, function(ii) paste(rev(ii), collapse = ".vs."))) ### %%% use rev() to revert element order
        vi <- (pn1 %in% names(cfm)) | (pn2 %in% names(cfm)) # check both reverse and forward pairing
        #if(!all(vi)) stop("unable to find cross-fit models for the following pairs : ", paste(pn1[!vi]))
        if(!all(vi)) {
            if(verbose > 0) {
                if(verbose > 1) {
                    cat(paste("WARNING: unable to find cross-fit models for the following pairs : ", paste(pn1[!vi], collapse = " ")), "\n")
                } else {
                    cat("WARNING: unable to find cross-fit models for ", sum(!vi), " out of ", length(vi), " pairs. Using a subset.\n")
                }
            }
            # use a subset
            if(sum(vi) > 3) {
                pn1 <- pn1[vi]
                pn2 <- pn2[vi]
                vi <- vi[vi]
            } else {
                stop("less than 3 valid cross-fit pairs are available! giving up.")
            }
        }

        #rl <- cfm[vi]
        vi.names<-names(cfm)[names(cfm) %in% c(pn1, pn2)] ### a similar selection was done like this in calculate.crossfit.models() function
        rl <- cfm[vi.names]  ### with this sub-selection we select only sample pairs within the current group (e.g. pairs of ES)

        # determine the set genes that were not attributed to failure in any cross-comparison
        x <- lapply(rl, function(d) {
            ll <- list(!(1:nrow(counts)) %in% d$vi[which(d$clusters > 1)], !(1:nrow(counts)) %in% d$vi[which(d$clusters %% 3  != 0)])
            names(ll) <- d$ii
            return(ll)
        })
        vil <- do.call(cbind, tapply(unlist(x, recursive = FALSE), factor(unlist(lapply(x, names)), levels = ids), function(l) {
            x <- rowSums(do.call(cbind, l), na.rm = FALSE) == 0
            x[is.na(x)] <- FALSE
            return(x)
        }))

        #x <- lapply(rl, function(d) { ll <- list((d$failures == 1), (d$failures == 2)) names(ll) <- d$ii return(ll) })
        #vil <- do.call(cbind, tapply(unlist(x, recursive = FALSE), factor(unlist(lapply(x, names)), levels = ids), function(l) { x <- rowSums(do.call(cbind, l), na.rm = FALSE) == 0 x[is.na(x)] <- FALSE return(x) }))

        t.ls <- ls$ls[ids]
        adjust <- NULL
        if(!is.null(ls$adjustments)) { ls$adjustments[[groups[ids[1]]]] }
        # fit two-NB2 mixture for each experiment
        if(verbose) { message(paste("fitting", group, "models:")) }
        gc()

        # pair cell name matrix
        nm <- do.call(rbind, lapply(rl, function(x) x$ii))

        ml <- papply(seq_along(ids), function(i) { try({
            if(verbose)  message(paste(i, ":", ids[i]))
            # determine genes with sufficient number of non-failed observations in other experiments
            vi <- which(rowSums(vil[, -i, drop = FALSE]) >= min(length(ids)-1, min.nonfailed))
            fpm <- rowMeans(t(t(counts[vi, ids[-i], drop = FALSE])/(t.ls[-i])))
            if(!is.null(adjust)) { fpm <- adjust(fpm)  } # adjust for between-group systematic differences
            df <- data.frame(count = counts[vi, ids[i]], fpm = fpm)

            # reconstruct failure prior for the cell by averaging across
            # cross-cell comparisons where the cell did participate
            cp <- exp(rowMeans(log(cbind(
                do.call(cbind, lapply(rl[which(nm[, 1] == ids[i])], function(d) {
                    ivi <- rep(NA, nrow(counts))
                    ivi[d$vi] <- 1:length(d$vi)
                    d$posterior[ivi[vi], 1]
                })),
                do.call(cbind, lapply(rl[which(nm[, 2] == ids[i])], function(d) {
                    ivi <- rep(NA, nrow(counts))
                    ivi[d$vi] <- 1:length(d$vi)
                    d$posterior[ivi[vi], 2]
                }))
            )), na.rm = TRUE))
            cp <- cbind(cp, 1-cp)

            nai <- which(is.na(cp[, 1]))
            cp[nai, 1] <- 1-(1e-10)
            cp[nai, 2] <- (1e-10)
            if(linear.fit) {
                m1 <- fit.nb2gth.mixture.model(df, prior = cp, nrep = 1, verbose = verbose, zero.count.threshold = zero.count.threshold, full.theta.range = theta.fit.range, theta.fit.range = theta.fit.range, use.constant.theta.fit = !local.theta.fit, ...)
            } else {
                m1 <- fit.nb2.mixture.model(df, prior = cp, nrep = nrep, verbose = verbose, zero.count.threshold = zero.count.threshold, ...)
            }

            if(return.compressed.models) {
                v <- get.compressed.v1.model(m1)
                cl <- clusters(m1)
                rm(m1)
                gc()
                return(list(model = v, clusters = cl))
            }

            # otherwise try to reduce the size of a full model
            # reduce return size
            #m1@posterior <- lapply(m1@posterior, function(m) { rownames(m) <- NULL return(m)})
            m1@posterior <- NULL
            #rownames(m1@concomitant@x) <- NULL
            m1@concomitant@x <- matrix()
            m1@model <- lapply(m1@model, function(mod) {
                mod@x <- matrix()
                mod@y <- matrix()
                #rownames(mod@x) <- NULL
                #rownames(mod@y) <- NULL
                return(mod)
            })

            # make a clean copy of the internal environment
            t.cleanenv <- function(comp) {
                el <- list2env(as.list(environment(comp@logLik), all.names = TRUE), parent = globalenv())
                ep <- list2env(as.list(environment(comp@predict), all.names = TRUE), parent = globalenv())
                pf <- get("predict", envir = el)
                environment(pf) <- ep
                assign("predict", pf, envir = el)
                pf <- get("predict", envir = ep)
                environment(pf) <- ep
                assign("predict", pf, envir = ep)

                pf <- get("logLik", envir = el)
                environment(pf) <- el
                assign("logLik", pf, envir = el)
                pf <- get("logLik", envir = ep)
                environment(pf) <- el
                assign("logLik", pf, envir = ep)

                environment(comp@logLik) <- el
                environment(comp@predict) <- ep
                comp
            }
            m1@components <- lapply(m1@components, function(cl) lapply(cl, t.cleanenv))

            # clean up the formula environment (was causing multithreading problems)
            rm(list = ls(env = attr(m1@concomitant@formula, ".Environment")), envir = attr(m1@concomitant@formula, ".Environment"))
            gc()
            #rm(list = ls(env = attr(m1@formula, ".Environment")), envir = attr(m1@formula, ".Environment"))
            return(m1)
        })}, n.cores = n.cores) # end cell iteration

        # check if there were errors in the multithreaded portion
        vic <- which(unlist(lapply(seq_along(ml), function(i) {
            if(class(ml[[i]]) == "try-error") {
                message("ERROR encountered in building a model for cell ", ids[i], " - skipping the cell. Error:")
                message(ml[[i]])
                #tryCatch(stop(paste("ERROR encountered in building a model for cell ", ids[i])), error = function(e) stop(e))
                return(FALSE);
            }
            return(TRUE);
        })))
        ml <- ml[vic]; names(ml) <- ids[vic];

        if(length(vic)<length(ids)) {
          message("ERROR fitting of ", (length(ids)-length(vic)), " out of ", length(ids), " cells resulted in errors reporting remaining ", length(vic), " cells")
        }

        if(save.plots && length(ml)>0) {
            # model fits
            #CairoPNG(filename = paste(group, "model.fits.png", sep = "."), width = 1024, height = 300*length(ids))
            pdf(file = paste(group, "model.fits.pdf", sep = "."), width = ifelse(linear.fit, 15, 13), height = 4)
            #l <- layout(matrix(seq(1, 4*length(ids)), nrow = length(ids), byrow = TRUE), rep(c(1, 1, 1, 0.5), length(ids)), rep(1, 4*length(ids)), FALSE)
            l <- layout(matrix(seq(1, 4), nrow = 1, byrow = TRUE), rep(c(1, 1, 1, ifelse(linear.fit, 1, 0.5)), 1), rep(1, 4), FALSE)
            par(mar = c(3.5, 3.5, 3.5, 0.5), mgp = c(2.0, 0.65, 0), cex = 0.9)
            invisible(lapply(seq_along(vic), function(j) {
                i <- vic[j];
                vi <- which(rowSums(vil[, -i, drop = FALSE]) >= min(length(ids)-1, min.nonfailed))
                df <- data.frame(count = counts[vi, ids[i]], fpm = rowMeans(t(t(counts[vi, ids[-i], drop = FALSE])/(t.ls[-i]))))
                plot.nb2.mixture.fit(ml[[j]], df, en = ids[i], do.par = FALSE, compressed.models = return.compressed.models)
            }))
            dev.off()
        }

        return(ml)

    }) # end group iteration

    if(return.compressed.models) {
        # make a joint model matrix
        jmm <- data.frame(do.call(rbind, lapply(mll, function(tl) do.call(rbind, lapply(tl, function(m) m$model)))))
        rownames(jmm) <- unlist(lapply(mll, names))
        # reorder in the original cell order
        attr(jmm, "groups") <- rep(names(mll), unlist(lapply(mll, length)))
        return(jmm)
    } else {
        return(mll)
    }
}


#######
## V1 optimized methods
#######

# gets an array summary of gam model structure (assumes a flat ifm list)
get.compressed.v1.models <- function(ifml) {
    data.frame(do.call(rbind, lapply(ifml, get.compressed.v1.model)))
}

# get a vector representation of a given model
get.compressed.v1.model <- function(m1) {
    if(class(m1@model[[2]]) == "FLXMRnb2gthC") { # linear fit model
        v <- c(m1@concomitant@coef[c(1:2), 2], get("coef", environment(m1@components[[1]][[1]]@predict)))
        names(v) <- c("conc.b", "conc.a", "fail.r")
        vth <- m1@components[[2]][[2]]@parameters$coef
        # translate mu regression from linear to log model
        v <- c(v, c("corr.b" = log(as.numeric(vth["corr.a"])), "corr.a" = 1), vth[-match("corr.a", names(vth))], "conc.a2" = m1@concomitant@coef[3, 2])
    } else { # original publication model
        v <- c(m1@concomitant@coef[, 2], get("coef", environment(m1@components[[1]][[1]]@predict)), m1@components[[2]][[2]]@parameters$coef, get("theta", environment(m1@components[[2]][[2]]@predict)))
        names(v) <- c("conc.b", "conc.a", "fail.r", "corr.b", "corr.a", "corr.theta")
    }
    v
}

# calculates posterior matrices (log scale) for a set of ifm models
calculate.posterior.matrices <- function(dat, ifm, prior, n.cores = 32, inner.cores = 4, outer.cores = round(n.cores/inner.cores)) {
    marginals <- data.frame(fpm = 10^prior$x - 1)
    marginals$fpm[marginals$fpm<0] <- 0
    lapply(ifm, function(group.ifm) {
        papply(sn(names(group.ifm)), function(nam) {
            df <- get.exp.logposterior.matrix(group.ifm[[nam]], dat[, nam], marginals, n.cores = inner.cores, grid.weight = prior$grid.weight)
            rownames(df) <- rownames(dat)
            colnames(df) <- as.character(prior$x)
            return(df)
        }, n.cores = n.cores)
    })
}

sample.posterior <- function(dat, ifm, prior, n.samples = 1, n.cores = 32) {
    marginals <- data.frame(fpm = 10^prior$x - 1)
    lapply(ifm, function(group.ifm) {
        papply(sn(names(group.ifm)), function(nam) {
            get.exp.sample(group.ifm[[nam]], dat[, nam], marginals, prior.x = prior$x, n = n.samples)
        }, n.cores = n.cores)
    })
}

# calculate joint posterior matrix for a given group of experiments
# lmatl - list of posterior matrices (log scale) for individual experiments
calculate.joint.posterior.matrix <- function(lmatl, n.samples = 100, bootstrap = TRUE, n.cores = 15) {
    if(bootstrap) {
        jpl <- papply(seq_len(n.cores), function(i) jpmatLogBoot(Matl = lmatl, Nboot = ceiling(n.samples/n.cores), Seed = i), n.cores = n.cores)
        jpl <- Reduce("+", jpl)
        jpl <- jpl/rowSums(jpl)
    } else {
        jpl <- Reduce("+", lmatl)
        jpl <- exp(jpl-log.row.sums(jpl))
    }
    rownames(jpl) <- rownames(lmatl[[1]])
    colnames(jpl) <- colnames(lmatl[[1]])
    jpl
}

# calculate joint posterior of a group defined by a composition vector
# lmatll - list of posterior matrix lists (as obtained from calculate.posterior.matrices)
# composition - a named vector, indicating the number of samples that should be drawn from each element of lmatll to compose a group
calculate.batch.joint.posterior.matrix <- function(lmatll, composition, n.samples = 100, n.cores = 15) {
    # reorder composition vector to match lmatll names
    jpl <- papply(seq_len(n.cores), function(i) jpmatLogBatchBoot(lmatll, composition[names(lmatll)], ceiling(n.samples/n.cores), i), n.cores = n.cores)
    jpl <- Reduce("+", jpl)
    jpl <- jpl/rowSums(jpl)
    #jpl <- jpmatLogBatchBoot(lmatll, composition[names(lmatll)], n.samples, n.cores)
    rownames(jpl) <- rownames(lmatll[[1]][[1]])
    colnames(jpl) <- colnames(lmatll[[1]][[1]])
    jpl
}

# calculates the likelihood of expression difference based on
# two posterior matrices (not adjusted for prior)
calculate.ratio.posterior <- function(pmat1, pmat2, prior, n.cores = 15, skip.prior.adjustment = FALSE) {
    n <- length(prior$x)
    if(!skip.prior.adjustment) {
        pmat1 <- t(t(pmat1)*prior$y)
        pmat2 <- t(t(pmat2)*prior$y)
    }

    chunk <- function(x, n) split(x, sort(rank(x) %% n))
    if(n.cores > 1) {
        x <- do.call(rbind, papply(chunk(1:nrow(pmat1), n.cores*5), function(ii) matSlideMult(pmat1[ii, , drop = FALSE], pmat2[ii, , drop = FALSE]), n.cores = n.cores))
    } else {
        x <- matSlideMult(pmat1, pmat2)
    }
    x <- x/rowSums(x)

    rv <- seq(prior$x[1]-prior$x[length(prior$x)], prior$x[length(prior$x)]-prior$x[1], length = length(prior$x)*2-1)
    colnames(x) <- as.character(rv)
    rownames(x) <- rownames(pmat1)
    return(x)
}

# quick utility function to get the difference Z score from the ratio posterior
get.ratio.posterior.Z.score <- function(rpost, min.p = 1e-15) {
    rpost <- rpost+min.p
    rpost <- rpost/rowSums(rpost)
    zi <- which.min(abs(as.numeric(colnames(rpost))))
    gs <- rowSums(rpost[, 1:(zi-1), drop = FALSE])
    zl <- pmin(0, qnorm(gs, lower.tail = FALSE))
    zg <- pmax(0, qnorm(gs+rpost[, zi, drop = FALSE], lower.tail = FALSE))
    z <- ifelse(abs(zl) > abs(zg), zl, zg)
}

# calculate a joint posterior matrix with bootstrap
jpmatLogBoot <- function(Matl, Nboot, Seed) {
    .Call("jpmatLogBoot", Matl, Nboot, Seed, PACKAGE = "scde")
}

# similar to the above, but compiles joint by sampling a pre-set
# number of different types (defined by Comp factor)
jpmatLogBatchBoot <- function(Matll, Comp, Nboot, Seed) {
    .Call("jpmatLogBatchBoot", Matll, Comp, Nboot, Seed, PACKAGE = "scde")
}

matSlideMult <- function(Mat1, Mat2) {
    .Call("matSlideMult", Mat1, Mat2, PACKAGE = "scde")
}

calculate.failure.p <- function(dat, ifm, n.cores = 32) {
    lapply(ifm, function(group.ifm) {
        lapply(sn(names(group.ifm)), function(nam) {
            get.concomitant.prob(group.ifm[[nam]], counts = dat[, nam])
        })
    })
}

# calculate failure probabilities across all cells for a given set
# of levels (lfpm - log(fpm) vector for all genes
calculate.failure.lfpm.p <- function(lfpm, ifm, n.cores = 32) {
    lapply(ifm, function(group.ifm) {
        lapply(sn(names(group.ifm)), function(nam) {
            get.concomitant.prob(group.ifm[[nam]], lfpm = lfpm)
        })
    })
}

# get expected fpm from counts
get.fpm.estimates <- function(m1, counts) {
    if(class(m1@components[[2]][[2]]) == "FLXcomponentE") {
        # gam do inverse interpolation
        b1 <- get("b1", envir = environment(m1@components[[2]][[2]]@predict))
        z <- approx(x = b1$fitted.values, y = b1$model$x, xout = counts, rule = 1:2)$y
        z[is.na(z)] <- -Inf
        z
    } else {
        # linear model
        par <- m1@components[[2]][[2]]@parameters
        if(!is.null(par[["linear"]])) {
            log((counts-par$coef[[1]])/par$coef[[2]])
        } else {
            (log(counts)-par$coef[[1]])/par$coef[[2]]
        }
    }
}


#######
## INTERNAL FUNCTIONS
#######

# clean up stale web server reference
.onAttach <- function(...) {

    if(exists("___scde.server", envir = globalenv())) {
        old.server <- get("___scde.server", envir = globalenv())
        n.apps <- length(old.server$appList)-1
        # TODO fix server rescue...
        packageStartupMessage("scde: found stale web server instance with ", n.apps, " apps. removing.")
        # remove
        rm("___scde.server", envir = globalenv())
        return(TRUE)

        if(n.apps > 0) {
            require(Rook)
            require(rjson)
            packageStartupMessage("scde: found stale web server instance with ", n.apps, " apps. restarting.")
            rm("___scde.server", envir = globalenv()) # remove old instance (apparently saved Rook servers can't just be restarted ... we'll make a new one and re-add all of the apps

            tryCatch( {
                server <- get.scde.server(ip = old.server$listenAddr, port = old.server$listenPort) # launch a new server
                if(!is.null(server)) {
                    lapply(old.server$appList[-1], function(sa) {
                        server$add(app = sa$app, name = sa$name)
                    })
                }
            }, error = function(e) message(e))

        } else {
            packageStartupMessage("scde: found stale web server instance with ", n.apps, " apps. removing.")
            # remove
            rm("___scde.server", envir = globalenv())
        }
    }
}

.onUnload <- function(libpath) {
    library.dynam.unload("scde", libpath, verbose = TRUE)
}

# rdf : count/fpm data frame
fit.nb2.mixture.model <- function(rdf, zero.count.threshold = 10, prior = cbind(rdf$count<= zero.count.threshold, rdf$count > zero.count.threshold), nrep = 3, iter = 50, verbose = 0, background.rate = 0.1, ...) {
    #mo1 <- FLXMRnb2glmC(count~1, components = c(1), theta.range = c(0.5, Inf))
    #mo1 <- FLXMRglmCf(count~1, components = c(1), family = "poisson", mu = 0.01)
    #mo1 <- FLXMRglmC(count~1, components = c(1), family = "poisson")
    mo1 <- FLXMRglmCf(count~1, family = "poisson", components = c(1), mu = log(background.rate))
    mo2 <- FLXMRnb2glmC(count~1+I(log(fpm)), components = c(2), theta.range = c(0.5, Inf))

    m1 <- mc.stepFlexmix(count~1, data = rdf, k = 2, model = list(mo1, mo2), control = list(verbose = verbose, minprior = 0, iter = iter), concomitant = FLXPmultinom(~I(log(fpm))+1), cluster = prior, nrep = nrep, ...)

    # check if the theta was underfit
    if(get("theta", envir = environment(m1@components[[2]][[2]]@logLik)) == 0.5) {
        # refit theta
        sci <- clusters(m1) == 2
        fit <- glm.nb.fit(m1@model[[2]]@x[sci, , drop = FALSE], m1@model[[2]]@y[sci], weights = rep(1, sum(sci)), offset = c(), init.theta = 0.5)
        assign("theta", value = fit$theta, envir = environment(m1@components[[2]][[2]]@logLik))
        m1@components[[2]][[2]]@parameters$coef <- fit$coefficients
        assign("coef", value = fit$coefficients, envir = environment(m1@components[[2]][[2]]@logLik))
        message("WARNING: theta was underfit, new theta = ", fit$theta)
    }

    return(m1)
}

fit.nb2gth.mixture.model <- function(rdf, zero.count.threshold = 10, prior = as.integer(rdf$count >= zero.count.threshold | rdf$fpm<median(rdf$fpm[rdf$count<zero.count.threshold]))+1, nrep = 0, verbose = 0 , full.theta.range = c(1e-2, 1e2), theta.fit.range = full.theta.range, theta.sp = 1e-2, use.constant.theta.fit = FALSE, alpha.weight.power = 1/2, iter = 50) {
    #mo1 <- FLXMRglmC(count~1, components = c(1), family = "poisson")
    #matrix(cbind(ifelse(rdf$count<= zero.count.threshold, 0.95, 0.05), ifelse(rdf$count > zero.count.threshold, 0.95, 0.05)))
    mo1 <- FLXMRglmCf(count~1, family = "poisson", components = c(1), mu = log(0.1))
    mo2 <- FLXMRnb2gthC(count~0+fpm, components = c(2), full.theta.range = full.theta.range, theta.fit.range = theta.fit.range, theta.fit.sp = theta.sp, constant.theta = use.constant.theta.fit, alpha.weight.power = alpha.weight.power)
    m1 <- mc.stepFlexmix(count~1, data = rdf, k = 2, model = list(mo1, mo2), control = list(verbose = verbose, minprior = 0, iter = iter), concomitant = FLXPmultinom(~I(log(fpm))+I(log(fpm)^2)+1), cluster = prior, nrep = nrep)
    return(m1)
}

# rdf : count/fpm data frame
# en : experiment name for plotting
# n.zero.windows - number of windows to visualize failure model fit
# m1 - fitted model
plot.nb2.mixture.fit <- function(m1, rdf, en, do.par = TRUE, n.zero.windows = 50, compressed.models = FALSE, bandwidth = 0.05) {
    #require(Cairo) require(RColorBrewer)
    if(do.par) {
        CairoPNG(filename = paste(en, "model.fit.png", sep = "."), width = 800, height = 300)
        l <- layout(matrix(c(1:4), 1, 4, byrow = TRUE), c(1, 1, 1, 0.5), rep(1, 4), FALSE)
        par(mar = c(3.5, 3.5, 3.5, 0.5), mgp = c(2.0, 0.65, 0), cex = 0.9)
    }
    smoothScatter(log10(rdf$fpm+1), log10(rdf$count+1), xlab = "expected FPM", ylab = "observed counts", main = paste(en, "scatter", sep = " : "), bandwidth = bandwidth)

    plot(c(), c(), xlim = range(log10(rdf$fpm+1)), ylim = range(log10(rdf$count+1)), xlab = "expected FPM", ylab = "observed counts", main = paste(en, "components", sep = " : "))
    if(compressed.models) {
        vpi <- m1$clusters == 1
    } else {
        vpi <- clusters(m1) == 1
    }
    if(sum(vpi) > 2){
        points(log10(rdf$fpm[vpi]+1), log10(rdf$count[vpi]+1), pch = ".", col = densCols(log10(rdf$fpm[vpi]+1), log10(rdf$count[vpi]+1), colramp = colorRampPalette(brewer.pal(9, "Reds")[-(1:3)])), cex = 2)
    }
    if(sum(!vpi) > 2){
        points(log10(rdf$fpm[!vpi]+1), log10(rdf$count[!vpi]+1), pch = ".", col = densCols(log10(rdf$fpm[!vpi]+1), log10(rdf$count[!vpi]+1), colramp = colorRampPalette(brewer.pal(9, "Blues")[-(1:3)])), cex = 2)
        # show fit
        fpmo <- order(rdf$fpm[!vpi], decreasing = FALSE)
        if(compressed.models) {
            #rf <- scde.failure.probability(data.frame(t(m1$model)), magnitudes = log(rdf$fpm))
            lines(log10(rdf$fpm[!vpi]+1)[fpmo], log10(exp(m1$model[["corr.a"]]*log(rdf$fpm[!vpi])[fpmo]+m1$model[["corr.b"]])+1))
            if("corr.ltheta.b" %in% names(m1$model)) {
                # show 95% CI for the non-constant theta fit
                xval <- range(log(rdf$fpm[!vpi]))
                xval <- seq(xval[1], xval[2], length.out = 100)
                thetas <- get.corr.theta(m1$model, xval)
                #thetas <- exp(m1$model[["corr.ltheta.i"]]+m1$model[["corr.ltheta.lfpm"]]*xval)
                #thetas <- (1+exp((m1$model["corr.ltheta.lfpm.m"] - xval)/m1$model["corr.ltheta.lfpm.s"]))/m1$model["corr.ltheta.a"]
                alpha <- 0.05
                yval <- exp(m1$model[["corr.a"]]*xval + m1$model[["corr.b"]])
                lines(log10(exp(xval)+1), log10(qnbinom(alpha/2, size = thetas, mu = yval)+1), col = 1, lty = 2)
                lines(log10(exp(xval)+1), log10(qnbinom(1-alpha/2, size = thetas, mu = yval)+1), col = 1, lty = 2)
                lines(log10(exp(xval)+1), log10(qnbinom(alpha/2, size = m1$model[["corr.theta"]], mu = yval)+1), col = 8, lty = 2)
                lines(log10(exp(xval)+1), log10(qnbinom(1-alpha/2, size = m1$model[["corr.theta"]], mu = yval)+1), col = 8, lty = 2)
            }
        } else {
            lines(log10(rdf$fpm[!vpi]+1)[fpmo], log10(m1@components[[2]][[2]]@predict(cbind(1, log(rdf$fpm[!vpi])))+1)[fpmo], col = 4)
        }
    }
    legend(x = "topleft", col = c("red", "blue"), pch = 19, legend = c("failure component", "correlated component"), bty = "n", cex = 0.9)

    # zero fit
    if(n.zero.windows > nrow(rdf)) { n.zero.windows <- nrow(rdf) }
    bw <- floor(nrow(rdf)/n.zero.windows)
    if(compressed.models) {
        rdf$cluster <- m1$clusters
    } else {
        rdf$cluster <- clusters(m1)
    }
    rdf <- rdf[order(rdf$fpm, decreasing = FALSE), ]
    fdf <- data.frame(y = rowMeans(matrix(log10(rdf$fpm[1:(n.zero.windows*bw)]+1), ncol = bw, byrow = TRUE)), zf = rowMeans(matrix(as.integer(rdf$cluster[1:(n.zero.windows*bw)] == 1), ncol = bw, byrow = TRUE)))
    plot(zf~y, fdf, ylim = c(0, 1), xlim = range(na.omit(log10(rdf$fpm+1))), xlab = "expected FPM", ylab = "fraction of failures", main = "failure model", pch = 16, cex = 0.5)
    ol <- order(rdf$fpm, decreasing = TRUE)
    if(compressed.models) {
        fp <- scde.failure.probability(data.frame(t(m1$model)), magnitudes = log(rdf$fpm))
        lines(log10(rdf$fpm[ol]+1), fp[ol], col = 2)
    } else {
        mt <- terms(m1@concomitant@formula, data = rdf)
        mf <- model.frame(delete.response(mt), data = rdf, na.action = NULL)
        cm0 <- exp(model.matrix(mt, data = mf) %*% m1@concomitant@coef)
        cm0 <- cm0/rowSums(cm0)
        lines(log10(rdf$fpm[ol]+1), cm0[ol, 1], col = 2)
    }


    # show thetas
    #tl <- c(fail = get("theta", envir = environment(m1@components[[1]][[1]]@logLik)), corr = get("theta", envir = environment(m1@components[[2]][[2]]@logLik)))
    if(compressed.models) {
        if("corr.ltheta.b" %in% names(m1$model)) {
            p <- exp(m1$model[["corr.a"]]*log(rdf$fpm[!vpi])+m1$model[["corr.b"]])
            alpha <- ((rdf$count[!vpi]/p-1)^2 - 1/p)
            trng <- log(range(c(m1$model[["corr.theta"]], thetas))) + 0.5*c(-1, 1)
            # restrict the alpha to the confines of the estimated theta values
            alpha[alpha > exp(-trng[1])] <- exp(-trng[1])
            alpha[alpha<exp(-trng[2])] <- exp(-trng[2])

            smoothScatter(log10(rdf$fpm[!vpi]+1), -log10(alpha), ylim = trng*log10(exp(1)), xlab = "FPM", ylab = "log10(theta)", main = "overdispersion", bandwidth = bandwidth)
            xval <- range(log(rdf$fpm[!vpi]))
            xval <- seq(xval[1], xval[2], length.out = 100)
            #thetas <- exp(m1$model[["corr.ltheta.i"]]+m1$model[["corr.ltheta.lfpm"]]*xval)
            thetas <- get.corr.theta(m1$model, xval)
            #plot(log10(exp(xval)+1), log(thetas), ylim = trng, type = 'l', xlab = "FPM", ylab = "log(theta)", main = "overdispersion")
            lines(log10(exp(xval)+1), log10(thetas))
            abline(h = log10(m1$model[["corr.theta"]]), col = 1, lty = 2)
        } else {
            tl <- c(fail = c(), corr = m1$model[["corr.theta"]])
            barplot(tl, beside = TRUE, las = 2, col = c("dodgerblue1", "indianred1"), ylab = "magnitude", main = "theta")
        }

    } else {
        tl <- c(fail = c(0), corr = get("theta", envir = environment(m1@components[[2]][[2]]@logLik)))
        barplot(tl, beside = TRUE, las = 2, col = c("dodgerblue1", "indianred1"), ylab = "magnitude", main = "theta")
    }
    box()
    if(do.par) {   dev.off() }
}

## from nb2.crossmodels.r
mc.stepFlexmix <- function(..., nrep = 5, n.cores = nrep, return.all = FALSE) {
    if(nrep < 2) {
        return(flexmix(...))
    } else {
        ml <- papply(seq_len(nrep), function(m) {
            x = try(flexmix(...))
        }, n.cores = n.cores)
        if(return.all) { return(ml) }
        ml <- ml[unlist(lapply(ml, function(x) !is(x, "try-error")))]
        logLiks <- unlist(lapply(ml, logLik))
        ml[[which.max(logLiks)]]
    }
}

# df: count matrix
# xr: expression level for each row in the matrix
# ml: fitted model list for a replicate
get.rep.set.posteriors <- function(xr, df, ml, rescale = TRUE) {
    pl <- do.call(cbind, lapply(seq_along(ml), function(i) {
        edf <- data.frame(y = df[, i], xr = xr)
        m1 <- ml[[i]]
        x <- FLXgetModelmatrix(m1@model[[1]], edf, m1@model[[1]]@formula)
        #cx <- FLXgetModelmatrix(m1@concomitant, edf, m1@concomitant@formula)
        cm <- (1/(1+exp(-1*(x@x %*% m1@concomitant@coef[, -1]))))
        p1 <- (1-cm)*exp(FLXdeterminePostunscaled(x, m1@components[[1]]))
        p2 <- cm*exp(FLXdeterminePostunscaled(x, m1@components[[2]]))
        tpr <- as.numeric(p1+p2)
        if(rescale) {tpr <- tpr/sum(tpr) }
        return(tpr)
    }))
    colnames(pl) <- names(ml)

    return(pl)
}

# evaluates likelihood for a list of models and a set of
# corresponding counts
# ml - model list
# counts - observed count matrix corresponding to the models
# marginals - marginal info, to which model-specific count will be appended
get.rep.set.general.model.posteriors <- function(ml, counts, marginals, grid.weight = rep(1, nrow(marginals)), rescale = TRUE, min.p = 0) {
    pl <- do.call(cbind, lapply(seq_along(ml), function(i) {
        marginals$count <- counts[, i]
        rowSums(get.component.model.lik(ml[[i]], marginals))+min.p
    }))
    if(rescale) {
        #pl <- pl*grid.weight+min.p
        pl <- t(t(pl)/colSums(pl))
    }
    colnames(pl) <- names(ml)
    return(pl)
}

get.rep.set.general.model.logposteriors <- function(ml, counts, marginals, grid.weight = rep(1, nrow(marginals)), rescale = TRUE) {
    pl <- do.call(rbind, lapply(seq_along(ml), function(i) {
        marginals$count <- counts[, i]
        log.row.sums(get.component.model.loglik(ml[[i]], marginals))
    }))
    if(rescale) {
        pl <- pl-log.row.sums(pl)
    }
    rownames(pl) <- names(ml)
    return(t(pl))
}

# evaluate likelihood on a mixed model with a binomial concomitant
# returns posterior probability for each component: rowSums(return) gives
# total likelihood. (note it's not on a log scale!)
get.component.model.lik <- function(m1, newdata) {
    # core models
    cp <- exp(do.call("+", lapply(seq_along(m1@model), function(i) {
        y <- posterior(m1@model[[i]], newdata, lapply(m1@components, "[[", i))
    })))
    # concomitant

    # no groups!
    mt <- terms(m1@concomitant@formula, data = newdata)
    mf <- model.frame(delete.response(mt), data = newdata, na.action = NULL)
    cm0 <- exp(model.matrix(mt, data = mf) %*% m1@concomitant@coef)
    cm0 <- cm0/rowSums(cm0)
    cm0[!is.finite(cm0)] <- 1
    return(cp*cm0)
}

# same as above, but keeping log resolution
get.component.model.loglik <- function(m1, newdata) {
    # core models
    cp <- do.call("+", lapply(seq_along(m1@model), function(i) {
        y <- posterior(m1@model[[i]], newdata, lapply(m1@components, "[[", i))
    }))
    cp[!is.finite(cp)] <- sign(cp[!is.finite(cp)])*.Machine$double.xmax
    # concomitant
    # no groups!
    mt <- terms(m1@concomitant@formula, data = newdata)
    mf <- model.frame(delete.response(mt), data = newdata, na.action = NULL)
    cm0 <- model.matrix(mt, data = mf) %*% m1@concomitant@coef
    cm0[is.nan(cm0)] <- 1
    cm0[!is.finite(cm0)] <- sign(cm0[!is.finite(cm0)])*.Machine$double.xmax

    cm0 <- cm0-log.row.sums(cm0)
    return(cp+cm0)
}

# returns a matrix of posterior values, with rows corresponding to genes, and
# columns to marginal values (prior fpkm grid)
# m1 - model
# counts - vector of per-gene counts for a given experiment
# marginals - fpm data frame
get.exp.posterior.matrix <- function(m1, counts, marginals, grid.weight = rep(1, nrow(marginals)), rescale = TRUE, n.cores = 32, min.p = 0) {
    uc <- unique(counts)
    #message(paste("get.exp.posterior.matrix() :", round((1-length(uc)/length(counts))*100, 3), "% savings"))
    cat(".")
    df <- do.call(rbind, papply(uc, function(x) {
        rowSums(get.component.model.lik(m1, cbind(marginals, count = rep(x, nrow(marginals)))))+min.p
    }, n.cores = n.cores))
    if(rescale) {
        #df <- t(t(df)*grid.weight)+min.p
        df <- df/rowSums(df)
    }
    df <- df[match(counts, uc), , drop = FALSE]
    rownames(df) <- names(counts)
    df
}

get.exp.logposterior.matrix <- function(m1, counts, marginals, grid.weight = rep(1, nrow(marginals)), rescale = TRUE, n.cores = 32) {
    uc <- unique(counts)
    #message(paste("get.exp.logposterior.matrix() :", round((1-length(uc)/length(counts))*100, 3), "% savings"))
    cat(".")
    df <- do.call(rbind, papply(uc, function(x) {
        log.row.sums(get.component.model.loglik(m1, cbind(marginals, count = rep(x, nrow(marginals)))))
    }, n.cores = n.cores))
    if(rescale) {
        df <- df-log.row.sums(df)
    }
    df <- df[match(counts, uc), , drop = FALSE]
    rownames(df) <- names(counts)
    df
}

# similar to get.exp.posterior.matrix(), but returns inverse ecdf list
# note that x must be supplied
get.exp.posterior.samples <- function(pmatl, prior, n.samples = 1, n.cores = 32) {
    sl <- papply(seq_along(pmatl), function(i) t(apply(pmatl[[i]], 1, function(d) approxfun(cumsum(d), prior$x, rule = 2)(runif(n.samples)))), n.cores = n.cores)
    names(sl) <- names(pmatl)
    sl
}
# similar to get.exp.posterior.matrix(), but returns inverse ecdf list
# note that x must be supplied
get.exp.sample <- function(m1, counts, marginals, prior.x, n, rescale = TRUE) {
    do.call(rbind, papply(counts, function(x) {
        tpr <- log.row.sums(get.component.model.loglik(m1, cbind(marginals, count = rep(x, nrow(marginals)))))
        if(rescale)  {
            tpr <- exp(tpr-max(tpr))
            tpr <- tpr/sum(tpr)
        }
        return(approxfun(cumsum(tpr), prior.x, rule = 2)(runif(n)))
    }, n.cores=1))
}

# gets a probability of failed detection for a given observation
# optional vector of fpm values (log) can be supplied to evaluate mixing probability
# at a point other than MLE fpm
get.concomitant.prob <- function(m1, counts = NULL, lfpm = NULL) {
    if(is.null(lfpm)) {
        lfpm <- get.fpm.estimates(m1, counts)
    }
    newdata <- data.frame(fpm = exp(lfpm))
    mt <- terms(m1@concomitant@formula, data = newdata)
    mf <- model.frame(delete.response(mt), data = newdata, na.action = NULL)
    cm0 <- exp(model.matrix(mt, data = mf) %*% m1@concomitant@coef)
    cm0[is.nan(cm0)] <- 1
    cm0 <- cm0/rowSums(cm0)
    return(as.numeric(cm0[, 1]))
}

# copied from flexmix
log.row.sums <- function(m) {
    M <- m[cbind(seq_len(nrow(m)), max.col(m, ties.method = "first"))] # "random" doesn't work!
    M + log(rowSums(exp(m - M)))
}


#######
## from nb1glm.R
#######

# nb2 glm implementation
setClass("FLXMRnb2glm", contains = "FLXMRglm", package = "flexmix")

FLXMRnb2glm <- function(formula = . ~ .,  offset = NULL, init.theta = NULL, theta.range = c(0, 1e3), ...) {
    #require(MASS)
    family <- "negative.binomial"
    glmrefit <- function(x, y, w) {
        #message("FLXRnb2glm:refit:nb2")
        fit <- c(glm.nb.fit(x, y, weights = w, offset = offset, init.theta = init.theta, theta.range = theta.range),
                 list(call = sys.call(), offset = offset, control = eval(formals(glm.fit)$control), method = "glm.fit")
        )
        fit$df.null <- sum(w) + fit$df.null - fit$df.residual - fit$rank
        fit$df.residual <- sum(w) - fit$rank
        fit$x <- x
        fit
    }

    z <- new("FLXMRnb2glm", weighted = TRUE, formula = formula,
             name = "FLXMRnb2glm", offset = offset,
             family = family, refit = glmrefit)
    z@preproc.y <- function(x) {
        if (ncol(x)  >  1)
            stop(paste("for the", family, "family y must be univariate"))
        x
    }


    z@defineComponent <- expression({
        predict <- function(x, ...) {
            dotarg = list(...)
            #message("FLXRnb2glm:predict:nb2")
            if("offset" %in% names(dotarg)) offset <- dotarg$offset
            p <- x%*%coef
            if (!is.null(offset)) p <- p + offset
            negative.binomial(theta)$linkinv(p)
        }
        logLik <- function(x, y, ...) {
            r <- dnbinom(y, size = theta, mu = predict(x, ...), log = TRUE)
            #message(paste("FLXRnb2glm:loglik:nb2", theta))
            return(r)
        }

        new("FLXcomponent",
            parameters = list(coef = coef),
            logLik = logLik, predict = predict,
            df = df)
    })

    z@fit <- function(x, y, w){
        #message("FLXRnb2glm:fit:nb2")
        w[y<= 1] <- w[y<= 1]/1e6 # focus the fit on non-failed genes
        fit <- glm.nb.fit(x, y, weights = w, offset = offset, init.theta = init.theta, theta.range = theta.range)
        # an ugly hack to restrict to non-negative slopes
        cf <- coef(fit)
        if(cf[2]<0) { cf <- c(mean(y*w)/sum(w), 0) }
        with(list(coef = cf, df = ncol(x), theta = fit$theta, offset = offset), eval(z@defineComponent))
    }

    return(z)
}

# component-specific version of the nb2glm
# nb2 glm implementation
setClass("FLXMRnb2glmC", representation(vci = "ANY"), contains = "FLXMRnb2glm", package = "flexmix")

# components is used to specify the indices of the components on which likelihood will be
# evaluated. Others will return as loglik of 0
FLXMRnb2glmC <- function(... , components = NULL) {
    #require(MASS)
    z <- new("FLXMRnb2glmC", FLXMRnb2glm(...), vci = components)
    z
}

# nb2 glm implementation
setClass("FLXMRnb2gam", contains = "FLXMRglm", package = "flexmix")

setClass("FLXcomponentE",
         representation(refitTheta = "function",
                        theta.fit = "ANY"),
         contains = "FLXcomponent", package = "flexmix")

# nb2 implementation with a simple trimmed-mean/median slope, and a gam theta fit
setClass("FLXMRnb2gth", contains = "FLXMRglm", package = "flexmix")

# get values of theta for a given set of models and expression (log-scale) magnitudes
get.corr.theta <- function(model, lfpm, theta.range = NULL) {
    if("corr.ltheta.b" %in% names(model)) {
        #th <- exp(-1*(model[["corr.ltheta.a"]]/(1+exp((model[["corr.ltheta.lfpm.m"]] - lfpm)/model[["corr.ltheta.lfpm.s"]])) + log(model[["corr.ltheta.b"]])))
        th <- exp(-1*(model[["corr.ltheta.b"]]+(model[["corr.ltheta.t"]]-model[["corr.ltheta.b"]])/(1+10^((model[["corr.ltheta.m"]]-lfpm)*model[["corr.ltheta.s"]]))^model[["corr.ltheta.r"]]))
    } else {
        if(length(lfpm) > 1) {
            th <- rep(model[["corr.theta"]], length(lfpm))
        } else {
            th <- model[["corr.theta"]]
        }
    }
    if(!is.null(theta.range)) {
        th[th<theta.range[1]] <- theta.range[1]
        th[th > theta.range[2]] <- theta.range[2]
        th[is.nan(th)] <- theta.range[1]
    }
    th
}

FLXMRnb2gth <- function(formula = . ~ .,  offset = NULL, full.theta.range = c(1e-3, 1e3), theta.fit.range = full.theta.range*c(1e-1, 1e1), theta.fit.sp = c(-1), constant.theta = FALSE, slope.mean.trim = 0.4, alpha.weight.power = 1/2, ...) {
    if(slope.mean.trim<0) { slope.mean.trim <- 0 }
    if(slope.mean.trim > 0.5) { slope.mean.trim <- 0.5 }

    family <- "negative.binomial"
    glmrefit <- function(x, y, w) {
        message("ERROR: FLXRnb2gth:glmrefit: NOT IMPLEMENTED")
        return(NULL)
    }

    z <- new("FLXMRnb2gth", weighted = TRUE, formula = formula,
             name = "FLXMRnb2gth", offset = offset,
             family = family, refit = glmrefit)
    z@preproc.y <- function(x) {
        if (ncol(x)  >  1)
            stop(paste("for the", family, "family y must be univariate"))
        x
    }

    z@defineComponent <- expression({
        predict <- function(x, ...) {
            dotarg = list(...)
            #message("FLXRnb2gth:predict:nb2")
            coef["corr.a"]*x
        }
        logLik <- function(x, y, ...) {
            dotarg = list(...)
            #message("FLXRnb2gth:logLik")
            if(constant.theta) {
                th <- coef["corr.theta"]
            } else {
                #th <- exp(coef["corr.ltheta.i"] + coef["corr.ltheta.lfpm"]*log(x))
                #th <- exp(-1*(coef["corr.ltheta.a"]/(1+exp((coef["corr.ltheta.lfpm.m"] - log(x))/coef["corr.ltheta.lfpm.s"])) + log(coef["corr.ltheta.b"])))
                th <- get.corr.theta(coef, log(x))

            }
            # restrict theta to the pre-defined range
            th[th  >  full.theta.range[2]] <- full.theta.range[2]
            th[th < full.theta.range[1]] <- full.theta.range[1]
            # evaluate NB
            r <- dnbinom(y, size = th, mu = coef["corr.a"]*x, log = TRUE)
        }

        new("FLXcomponent",
            parameters = list(coef = coef, linear = TRUE),
            logLik = logLik, predict = predict, df = df)
    })

    z@fit <- function(x, y, w){
        # message("FLXRnb2gth:fit")

        # estimate slope using weighted trimmed mean
        #w[y == 0] <- w[y == 0]/1e6
        #r <- y/x
        #ro <- order(r)
        ## cumulative weight sum along the ratio order (to figure out where to trim)
        #cs <- cumsum(w[ro])/sum(w)
        #lb <- min(which(cs > slope.mean.trim))
        #ub <- max(which(cs<(1-slope.mean.trim)))
        #ro <- ro[lb:ub]
        ## slope fit
        #a <- weighted.mean(r[ro], w[ro])

        a <- as.numeric(coef(glm(y~0+x, family = poisson(link = "identity"), start = weighted.mean(y/x, w), weights = w))[1])

        # predicted values
        p <- a*x

        #te <- p^2/((y-p)^2 - p) # theta point estimates
        #te[te<theta.fit.range[1]] <- theta.fit.range[1] te[te > theta.fit.range[2]] <- theta.fit.range[2]
        alpha <- ((y/p-1)^2 - 1/p)
        alpha[alpha<1/theta.fit.range[2]] <- 1/theta.fit.range[2]
        alpha[alpha > 1/theta.fit.range[1]] <- 1/theta.fit.range[1]

        #theta <- MASS::theta.ml(y, p, sum(w), w)
        theta <- MASS::theta.md(y, p, sum(w)-1, w)
        theta <- pmin(pmax(theta.fit.range[1], theta), theta.fit.range[2])
        if(constant.theta) {
            v <- c("corr.a" = a, "corr.theta" = theta)
        } else {
            # fit theta linear model in the log space
            #theta.l <- glm(log(te)~log(x), weights = w)
            #ac <- tryCatch( {

            mw <- w*(as.numeric(alpha)^(alpha.weight.power))
            lx <- log(x)
            lx.rng <- range(lx)
            mid.s <- (sum(lx.rng))/2
            low <- log(x) < mid.s
            lalpha <- log(alpha)
            bottom.s <- quantile(lalpha[low], 0.025, na.rm = TRUE)
            top.s <- quantile(lalpha[!low], 0.975, na.rm = TRUE)

            wsr <- function(p, x, y, w = rep(1, length(y))) {
                # borrowing nplr approach here
                bottom <- p[1]
                top <- p[2]
                xmid <- p[3]
                scal <- p[4]
                r <- p[5]
                yfit <- bottom+(top-bottom)/(1+10^((xmid-x)*scal))^r
                #exp(-1*(model[["corr.ltheta.b"]]+(model[["corr.ltheta.t"]]-model[["corr.ltheta.b"]])/(1+10^((model[["corr.ltheta.m"]]-lfpm)*model[["corr.ltheta.s"]]))^model[["corr.ltheta.r"]]))
                residuals <- (y - yfit)^2
                return(sum(w*residuals))
            }

            #po <- nlm(f = wsr, p = c(bottom.s, top.s, mid.s, s = -1, r = 0.5), x = log(x), y = lalpha, w = w*(as.numeric(alpha)^(1/2)))
            #ac <- po$estimate

            #po <- nlminb(objective = wsr, start = c(bottom.s, top.s, mid.s, s = -1, r = 0.5), x = log(x), y = lalpha, w = mw)
            po <- nlminb(objective = wsr, start = c(bottom.s, top.s, mid.s, s = -1, r = 0.5), x = log(x), y = lalpha, w = mw, lower = c(-100, -10, -100, -100, 0.1), upper = c(10, 100, 100, 0, 20))
            ac <- po$par

            #smoothScatter(log(x), log(alpha))
            #p <- ac points(log(x), p[1]+(p[2]-p[1])/(1+10^((p[3]-log(x))*p[4]))^p[5], col = 2, pch = ".")
            #browser()
            #}, error = function(e) {
            #  message("encountered error trying to fit logistic model with guessed parameters.")
            #  # fit with fewer parameters?
            #})

            v <- c(a, theta, ac)
            names(v) <- c("corr.a", "corr.theta", "corr.ltheta.b", "corr.ltheta.t", "corr.ltheta.m", "corr.ltheta.s", "corr.ltheta.r")
        }

        with(list(coef = v, full.theta.range = full.theta.range, df = ncol(x), offset = offset), eval(z@defineComponent))
    }

    return(z)
}

# component-specific version of the nb2gth
# nb2 gam implementation
setClass("FLXMRnb2gthC", representation(vci = "ANY"), contains = "FLXMRnb2gth", package = "flexmix")

# components is used to specify the indices of the components on which likelihood will be
# evaluated. Others will return as loglik of 0
FLXMRnb2gthC <- function(... , components = NULL) {
    #require(mgcv)
    z <- new("FLXMRnb2gthC", FLXMRnb2gth(...), vci = components)
    z
}

setMethod("FLXdeterminePostunscaled", signature(model = "FLXMRnb2glmC"), function(model, components, ...) {
    if(is.null(model@vci)) {
        #message("FLXMRnb2glmC:FLXdeterminePostunscaled - applying to all components")
        m <- matrix(sapply(components, function(x) x@logLik(model@x, model@y)), nrow = nrow(model@y))
    } else {
        #message(paste("FLXMRnb2glmC:FLXdeterminePostunscaled - applying to components", paste(model@vci, collapse = " ")))
        m <- matrix(do.call(cbind, lapply(seq_along(components), function(i) {
            if(i %in% model@vci) {
                components[[i]]@logLik(model@x, model@y)
            } else {
                rep(0, nrow(model@y))
            }
        })), nrow = nrow(model@y))
    }
})

setMethod("FLXmstep", signature(model = "FLXMRnb2glmC"), function(model, weights, ...) {
    # make up a dummy component return
    coef <- rep(0, ncol(model@x))
    names(coef) <- colnames(model@x)
    control <- eval(formals(glm.fit)$control)
    comp.1 <- with(list(coef = coef, df = 0, offset = NULL,
                        family = model@family), eval(model@defineComponent))

    # iterate over components
    unlist(lapply(seq_len(ncol(weights)), function(i) {
        if(i %in% model@vci) {
            #message(paste("FLXMRnb2glmC:FLXmstep - running m-step for component", i))
            FLXmstep(as(model, "FLXMRnb2glm"), weights[, i, drop = FALSE])
        } else {
            #message(paste("FLXMRnb2glmC:FLXmstep - dummy return for component", i))
            list(comp.1)
        }
    }), recursive = FALSE)
})

# same for gth
setMethod("FLXdeterminePostunscaled", signature(model = "FLXMRnb2gthC"), function(model, components, ...) {
    if(is.null(model@vci)) {
        #message("FLXMRnb2gthC:FLXdeterminePostunscaled - applying to all components")
        m <- matrix(sapply(components, function(x) x@logLik(model@x, model@y)), nrow = nrow(model@y))
    } else {
        #message(paste("FLXMRnb2gthC:FLXdeterminePostunscaled - applying to components", paste(model@vci, collapse = " ")))
        m <- matrix(do.call(cbind, lapply(seq_along(components), function(i) {
            if(i %in% model@vci) {
                components[[i]]@logLik(model@x, model@y)
            } else {
                rep(0, nrow(model@y))
            }
        })), nrow = nrow(model@y))
    }
})

setMethod("FLXmstep", signature(model = "FLXMRnb2gthC"), function(model, weights, ...) {
    # make up a dummy component return
    coef <- rep(0, ncol(model@x))
    names(coef) <- colnames(model@x)
    control <- eval(formals(glm.fit)$control)
    comp.1 <- with(list(q1 = list(coefficients = c(1)), coef = coef, df = 0, offset = NULL,
                        family = model@family), eval(model@defineComponent))

    # iterate over components
    unlist(lapply(seq_len(ncol(weights)), function(i) {
        if(i %in% model@vci) {
            #message(paste("FLXMRnb2gthC:FLXmstep - running m-step for component", i))
            FLXmstep(as(model, "FLXMRnb2gth"), weights[, i, drop = FALSE])
        } else {
            #message(paste("FLXMRnb2gthC:FLXmstep - dummy return for component", i))
            list(comp.1)
        }
    }), recursive = FALSE)
})

# component-specific version of the nb2glm
# nb2 glm implementation
setClass("FLXMRglmC", representation(vci = "ANY"), contains = "FLXMRglm", package = "flexmix")

# components is used to specify the indices of the components on which likelihood will be
# evaluated. Others will return as loglik of 0
FLXMRglmC <- function(... , components = NULL) {
    #require(MASS)
    z <- new("FLXMRglmC", FLXMRglm(...), vci = components)
    z
}

setMethod("FLXdeterminePostunscaled", signature(model = "FLXMRglmC"), function(model, components, ...) {
    if(is.null(model@vci)) {
        #message("FLXMRnb2glmC:FLXdeterminePostunscaled - applying to all components")
        m <- matrix(sapply(components, function(x) x@logLik(model@x, model@y)), nrow = nrow(model@y))
    } else {
        #message(paste("FLXMRnb2glmC:FLXdeterminePostunscaled - applying to components", paste(model@vci, collapse = " ")))
        m <- matrix(do.call(cbind, lapply(seq_along(components), function(i) {
            if(i %in% model@vci) {
                components[[i]]@logLik(model@x, model@y)
            } else {
                rep(0, nrow(model@y))
            }
        })), nrow = nrow(model@y))
    }
    #message("FLXMRnb2glmC:FLXdeterminePostunscaled : ")
    #message(m)
    #browser()
    m
})

setMethod("FLXmstep", signature(model = "FLXMRglmC"), function(model, weights, ...) {
    # make up a dummy component return
    coef <- rep(0, ncol(model@x))
    names(coef) <- colnames(model@x)
    control <- eval(formals(glm.fit)$control)
    comp.1 <- with(list(coef = coef, df = 0, offset = NULL,
                        family = model@family), eval(model@defineComponent))

    # iterate over components
    unlist(lapply(seq_len(ncol(weights)), function(i) {
        if(i %in% model@vci) {
            #message(paste("FLXMRglmC:FLXmstep - running m-step for component", i))
            FLXmstep(as(model, "FLXMRglm"), weights[, i, drop = FALSE])
        } else {
            #message(paste("FLXMRglmC:FLXmstep - dummy return for component", i))
            list(comp.1)
        }
    }), recursive = FALSE)
})

# mu-fixed version
setClass("FLXMRglmCf", representation(mu = "numeric"), contains = "FLXMRglmC", package = "flexmix")

FLXMRglmCf <- function(... , family = c("binomial", "poisson"), mu = 0) {
    #require(MASS)
    family <- match.arg(family)
    z <- new("FLXMRglmCf", FLXMRglmC(..., family = family), mu = mu)
    z
}

setMethod("FLXmstep", signature(model = "FLXMRglmCf"), function(model, weights, ...) {
    # make up a dummy component return
    coef <- c(model@mu, rep(0, ncol(model@x)-1))
    names(coef) <- colnames(model@x)
    control <- eval(formals(glm.fit)$control)
    comp.1 <- with(list(coef = coef, df = 0, offset = NULL,
                        family = model@family), eval(model@defineComponent))

    # iterate over components
    unlist(lapply(seq_len(ncol(weights)), function(i) {
        list(comp.1)
    }), recursive = FALSE)
})

# a magnitude-weighted version of the FLXPmultinom (to down-weight low-fpkm points during concomitant fit)
# alternatively: some kind of non-decreasing function could be used
setClass("FLXPmultinomW", contains = "FLXPmultinom")

FLXPmultinomW <- function(formula = ~1) {
    z <- new("FLXPmultinom", name = "FLXPmultinom", formula = formula)
    multinom.fit <- function(x, y, w, ...) {
        r <- ncol(x)
        p <- ncol(y)
        if (p < 2) stop("Multinom requires at least two components.")
        mask <- c(rep(0, r + 1), rep(c(0, rep(1, r)), p - 1))
        #if(missing(w)) w <- rep(1, nrow(y))
        #w <- round(exp(x[, 2]))
        nnet::nnet.default(x, y, w, mask = mask, size = 0,
                           skip = TRUE, softmax = TRUE, censored = FALSE,
                           rang = 0, trace = FALSE, ...)
    }
    z@fit <- function(x, y, w, ...) multinom.fit(x, y, w, ...)$fitted.values
    z@refit <- function(x, y, w, ...) {
        if (missing(w) || is.null(w)) w <- rep(1, nrow(y))
        #w <- round(exp(x[, 2]))
        fit <- nnet::multinom(y ~ 0 + x, weights = w, data = list(y = y, x = x), Hess = TRUE, trace = FALSE)
        fit$coefnames <- colnames(x)
        fit$vcoefnames <- fit$coefnames[seq_along(fit$coefnames)]
        dimnames(fit$Hessian) <- lapply(dim(fit$Hessian) / ncol(x), function(i) paste(rep(seq_len(i) + 1, each = ncol(x)), colnames(x), sep = ":"))
        fit
    }
    z
}

# variation of negative.binomial family that keeps theta value accessible
negbin.th <- function (theta = stop("'theta' must be specified"), link = "log")
{
    linktemp <- substitute(link)
    if (!is.character(linktemp))
        linktemp <- deparse(linktemp)
    if (linktemp %in% c("log", "identity", "sqrt"))
        stats <- make.link(linktemp)
    else if (is.character(link)) {
        stats <- make.link(link)
        linktemp <- link
    }
    else {
        if (inherits(link, "link-glm")) {
            stats <- link
            if (!is.null(stats$name))
                linktemp <- stats$name
        }
        else stop(linktemp, " link not available for negative binomial family available links are \"identity\", \"log\" and \"sqrt\"")
    }
    .Theta <- theta
    env <- new.env(parent = .GlobalEnv)
    assign(".Theta", theta, envir = env)
    variance <- function(mu) mu + mu^2/.Theta
    validmu <- function(mu) all(mu  >  0)
    dev.resids <- function(y, mu, wt) 2 * wt * (y * log(pmax(1,
                                                             y)/mu) - (y + .Theta) * log((y + .Theta)/(mu + .Theta)))
    aic <- function(y, n, mu, wt, dev) {
        term <- (y + .Theta) * log(mu + .Theta) - y * log(mu) +
            lgamma(y + 1) - .Theta * log(.Theta) + lgamma(.Theta) -
            lgamma(.Theta + y)
        2 * sum(term * wt)
    }
    initialize <- expression({
        if (any(y < 0)) stop("negative values not allowed for the negative binomial family")
        n <- rep(1, nobs)
        mustart <- y + (y ==  0)/6
    })
    simfun <- function(object, nsim) {
        ftd <- fitted(object)
        val <- rnegbin(nsim * length(ftd), ftd, .Theta)
    }
    environment(variance) <- environment(validmu) <- environment(dev.resids) <- environment(aic) <- environment(simfun) <- env
    famname <- paste("Negative Binomial(", format(round(theta,
                                                        4)), ")", sep = "")
    structure(list(family = famname, link = linktemp, linkfun = stats$linkfun,
                   linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids,
                   aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
                   validmu = validmu, valideta = stats$valideta, simulate = simfun, theta = theta),
              class = "family")
}


# .fit version of the glm.nb
glm.nb.fit <- function(x, y, weights = rep(1, nobs), control = list(trace = 0, maxit = 20), offset = rep(0, nobs), etastart = NULL, start = NULL, mustart = NULL, init.theta = NULL, link = "log", method = "glm.fit", intercept = TRUE, theta.range = c(0, 1e5), ...) {
    method <- "custom.glm.fit"
    #require(MASS)
    loglik <- function(n, th, mu, y, w) sum(w * (lgamma(th +
                                                            y) - lgamma(th) - lgamma(y + 1) + th * log(th) + y *
                                                     log(mu + (y ==  0)) - (th + y) * log(th + mu)))
    # link <- substitute(link)

    Call <- match.call()
    control <- do.call("glm.control", control)
    n <- length(y)
    # family for the initial guess
    fam0 <- if (missing(init.theta) | is.null(init.theta))
        do.call("poisson", list(link = link))
    else
        do.call("negative.binomial", list(theta = init.theta, link = link))

    # fit function
    if (!missing(method)) {
        #message(paste("glm.nb.fit: method = ", method))
        if (!exists(method, mode = "function"))
            stop("unimplemented method: ", sQuote(method))
        glm.fitter <- get(method)
    }
    else {
        #message("glm.nb.fit: using default glm.fit")
        method <- "glm.fit"
        glm.fitter <- stats::glm.fit
        #glm.fitter <- custom.glm.fit
    }

    if (control$trace  >  1) {
        message("Initial fit:")
    }
    fit <- glm.fitter(x = x, y = y, weights = weights, start = start, etastart = etastart, mustart = mustart, offset = offset, family = fam0, control = control, intercept = intercept)
    class(fit) <- c("glm", "lm")

    mu <- fit$fitted.values
    th <- as.vector(theta.ml(y, mu, sum(weights), weights, limit = control$maxit, trace = control$trace  >  2))
    if(!is.null(theta.range)) {
        if(th<theta.range[1]) {
            if (control$trace  >  1)
                message("adjusting theta from ", signif(th), " to ", signif(theta.range[1]), " to fit the specified range")
            th <- theta.range[1]
        } else if(th > theta.range[2]) {
            if (control$trace  >  1)
                message("adjusting theta from ", signif(th), " to ", signif(theta.range[2]), " to fit the specified range")
            th <- theta.range[2]
        }
    }
    if (control$trace  >  1)
        message("Initial value for theta:", signif(th))
    fam <- do.call("negative.binomial", list(theta = th, link = link))
    iter <- 0
    d1 <- sqrt(2 * max(1, fit$df.residual))
    d2 <- del <- 1
    g <- fam$linkfun
    Lm <- loglik(n, th, mu, y, weights)
    Lm0 <- Lm + 2 * d1
    while ((iter <- iter + 1) <=  control$maxit && (abs(Lm0 - Lm)/d1 + abs(del)/d2)  >  control$epsilon) {
        eta <- g(mu)
        fit <- glm.fitter(x = x, y = y, weights = weights, etastart = eta, offset = offset, family = fam, control = list(maxit = control$maxit*10, epsilon = control$epsilon, trace = control$trace  >  1), intercept = intercept)
        t0 <- th
        th <- theta.ml(y, mu, sum(weights), weights, limit = control$maxit, trace = control$trace  >  2)
        if(!is.null(theta.range)) {
            if(th<theta.range[1]) {
                if (control$trace  >  1)
                    message("adjusting theta from ", signif(th), " to ", signif(theta.range[1]), " to fit the specified range")
                th <- theta.range[1]
            } else if(th > theta.range[2]) {
                if (control$trace  >  1)
                    message("adjusting theta from ", signif(th), " to ", signif(theta.range[2]), " to fit the specified range")
                th <- theta.range[2]
            }
        }
        fam <- do.call("negative.binomial", list(theta = th, link = link))
        mu <- fit$fitted.values
        del <- t0 - th
        Lm0 <- Lm
        Lm <- loglik(n, th, mu, y, weights)
        if (control$trace) {
            Ls <- loglik(n, th, y, y, weights)
            Dev <- 2 * (Ls - Lm)
            message("Theta(", iter, ")  = ", signif(th), ", 2(Ls - Lm)  = ",  signif(Dev))
        }
    }
    if (!is.null(attr(th, "warn")))
        fit$th.warn <- attr(th, "warn")
    if (iter  >  control$maxit) {
        warning("alternation limit reached")
        fit$th.warn <- gettext("alternation limit reached")
    }
    if (length(offset) && intercept) {
        null.deviance <- if ("(Intercept)" %in% colnames(x))
            glm.fitter(x[, "(Intercept)", drop = FALSE], y, weights = weights, offset = offset, family = fam, control = list(maxit = control$maxit*10, epsilon = control$epsilon, trace = control$trace  >   1), intercept = TRUE)$deviance
        else
            fit$deviance
        fit$null.deviance <- null.deviance
    }
    class(fit) <- c("negbin.th", "glm", "lm")
    Call$init.theta <- signif(as.vector(th), 10)
    Call$link <- link
    fit$call <- Call
    fit$x <- x
    fit$y <- y
    fit$theta <- as.vector(th)
    fit$SE.theta <- attr(th, "SE")
    fit$twologlik <- as.vector(2 * Lm)
    fit$aic <- -fit$twologlik + 2 * fit$rank + 2
    fit$method <- method
    fit$control <- control
    fit$offset <- offset
    fit
}

custom.glm.fit <- function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL,
                            mustart = NULL, offset = rep(0, nobs), family = gaussian(),
                            control = list(), intercept = TRUE, alpha = 0)
{
    control <- do.call("glm.control", control)
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y))
        rownames(y)
    else names(y)
    conv <- FALSE
    nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars ==  0
    if (is.null(weights))
        weights <- rep.int(1, nobs)
    if (is.null(offset))
        offset <- rep.int(0, nobs)
    variance <- family$variance
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object",
             call. = FALSE)
    dev.resids <- family$dev.resids
    aic <- family$aic
    mu.eta <- family$mu.eta
    unless.null <- function(x, if.null) if (is.null(x))
        if.null
    else x
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu <- unless.null(family$validmu, function(mu) TRUE)
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!valideta(eta))
            stop("invalid linear predictor values in empty model",
                 call. = FALSE)
        mu <- linkinv(eta)
        if (!validmu(mu))
            stop("invalid fitted means in empty model", call. = FALSE)
        dev <- sum(dev.resids(y, mu, weights))
        w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
        residuals <- (y - mu)/mu.eta(eta)
        good <- rep(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric(0L)
        iter <- 0L
    }
    else {
        coefold <- NULL
        eta <- if (!is.null(etastart))
            etastart
        else if (!is.null(start))
            if (length(start)  !=  nvars)
                stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                              nvars, paste(deparse(xnames), collapse = ", ")),
                     domain = NA)
        else {
            coefold <- start
            offset + as.vector(if (NCOL(x) ==  1)
                x * start
                else x %*% start)
        }
        else family$linkfun(mustart)
        mu <- linkinv(eta)
        if (!(validmu(mu) && valideta(eta)))
            stop("cannot find valid starting values: please specify some",
                 call. = FALSE)
        devold <- sum(dev.resids(y, mu, weights))
        boundary <- conv <- FALSE
        for (iter in 1L:control$maxit) {
            good <- weights  >  0
            varmu <- variance(mu)[good]
            if (any(is.na(varmu)))
                stop("NAs in V(mu)")
            if (any(varmu ==  0))
                stop("0s in V(mu)")
            mu.eta.val <- mu.eta(eta)
            if (any(is.na(mu.eta.val[good])))
                stop("NAs in d(mu)/d(eta)")
            good <- (weights  >  0) & (mu.eta.val  !=  0)
            if (all(!good)) {
                conv <- FALSE
                warning("no observations informative at iteration ",
                        iter)
                break
            }
            z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
            #z <- (eta - offset)[good] + (family$linkfun(y[good]) - family$linkfun(mu[good]))/family$linkfun(mu.eta.val[good])

            # attempting to be robust here, trowing out fraction with highest abs(z)
            if(alpha > 0) {
                qv <- quantile(abs(z), probs = c(alpha/2, 1.0-alpha/2))
                gvi <- which(good)[which(abs(z)<qv[1] | abs(z) > qv[2])]
                good[gvi] <- FALSE
                if (all(!good)) {
                    conv <- FALSE
                    warning("no observations informative at iteration ",
                            iter)
                    break
                }
                z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
            }

            w <- mu.eta.val[good]*sqrt(weights[good]/variance(mu)[good])

            ngoodobs <- as.integer(nobs - sum(!good))
            fit <- .Fortran("dqrls", qr = x[good, ] * w, n = ngoodobs,
                            p = nvars, y = w * z, ny = 1L, tol = min(1e-07,
                                                                     control$epsilon/1000), coefficients = double(nvars),
                            residuals = double(ngoodobs), effects = double(ngoodobs),
                            rank = integer(1L), pivot = 1L:nvars, qraux = double(nvars),
                            work = double(2 * nvars)) # , PACKAGE = "base"
            #browser()
            if (any(!is.finite(fit$coefficients))) {
                conv <- FALSE
                warning(gettextf("non-finite coefficients at iteration %d",
                                 iter), domain = NA)
                break
            }
            if (nobs < fit$rank)
                stop(gettextf("X matrix has rank %d, but only %d observations",
                              fit$rank, nobs), domain = NA)
            start[fit$pivot] <- fit$coefficients
            eta <- drop(x %*% start)
            mu <- linkinv(eta <- eta + offset)
            dev <- sum(dev.resids(y, mu, weights))
            if (control$trace)
                cat("Deviance  = ", dev, "Iterations -", iter,
                    "\n")
            boundary <- FALSE
            if (!is.finite(dev)) {
                if (is.null(coefold))
                    stop("no valid set of coefficients has been found: please supply starting values",
                         call. = FALSE)
                warning("step size truncated due to divergence",
                        call. = FALSE)
                ii <- 1
                while (!is.finite(dev)) {
                    if (ii  >  control$maxit)
                        stop("inner loop 1 cannot correct step size",
                             call. = FALSE)
                    ii <- ii + 1
                    start <- (start + coefold)/2
                    eta <- drop(x %*% start)
                    mu <- linkinv(eta <- eta + offset)
                    dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE
                if (control$trace)
                    cat("Step halved: new deviance  = ", dev, "\n")
            }
            # require deviance to go down
            if ((!is.null(coefold)) & (dev - devold)/(0.1 + abs(dev))  >  3*control$epsilon) {
                warning("step size truncated due to increasing divergence", call. = FALSE)
                ii <- 1
                while ((dev - devold)/(0.1 + abs(dev))  >  3*control$epsilon) {
                    if (ii  >  control$maxit)   {
                        warning("inner loop 1 cannot correct step size", call. = FALSE)
                        break
                    }
                    ii <- ii + 1
                    start <- (start + coefold)/2
                    eta <- drop(x %*% start)
                    mu <- linkinv(eta <- eta + offset)
                    dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE
                if (control$trace)
                    cat("Step halved: new deviance  = ", dev, "\n")
            }
            if (!(valideta(eta) && validmu(mu))) {
                if (is.null(coefold))
                    stop("no valid set of coefficients has been found: please supply starting values",
                         call. = FALSE)
                warning("step size truncated: out of bounds",
                        call. = FALSE)
                ii <- 1
                while (!(valideta(eta) && validmu(mu))) {
                    if (ii  >  control$maxit)
                        stop("inner loop 2 cannot correct step size",
                             call. = FALSE)
                    ii <- ii + 1
                    start <- (start + coefold)/2
                    eta <- drop(x %*% start)
                    mu <- linkinv(eta <- eta + offset)
                }
                boundary <- TRUE
                dev <- sum(dev.resids(y, mu, weights))
                if (control$trace)
                    cat("Step halved: new deviance  = ", dev, "\n")
            }
            if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
                conv <- TRUE
                coef <- start
                break
            }

            devold <- dev
            coef <- coefold <- start
        }
        if (!conv)
            warning("glm.fit: algorithm did not converge", call. = FALSE)
        if (boundary)
            warning("glm.fit: algorithm stopped at boundary value",
                    call. = FALSE)
        eps <- 10 * .Machine$double.eps
        if (family$family ==  "binomial") {
            if (any(mu  >  1 - eps) || any(mu < eps))
                warning("glm.fit: fitted probabilities numerically 0 or 1 occurred",
                        call. = FALSE)
        }
        if (family$family ==  "poisson") {
            if (any(mu < eps))
                warning("glm.fit: fitted rates numerically 0 occurred",
                        call. = FALSE)
        }
        if (fit$rank < nvars)
            coef[fit$pivot][seq.int(fit$rank + 1, nvars)] <- NA
        xxnames <- xnames[fit$pivot]
        residuals <- (y - mu)/mu.eta(eta)
        fit$qr <- as.matrix(fit$qr)
        nr <- min(sum(good), nvars)
        if (nr < nvars) {
            Rmat <- diag(nvars)
            Rmat[1L:nr, 1L:nvars] <- fit$qr[1L:nr, 1L:nvars]
        }
        else Rmat <- fit$qr[1L:nvars, 1L:nvars]
        Rmat <- as.matrix(Rmat)
        Rmat[row(Rmat)  >  col(Rmat)] <- 0
        names(coef) <- xnames
        colnames(fit$qr) <- xxnames
        dimnames(Rmat) <- list(xxnames, xxnames)
    }
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0, nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    if (!EMPTY)
        names(fit$effects) <- c(xxnames[seq_len(fit$rank)], rep.int("",
                                                                    sum(good) - fit$rank))
    wtdmu <- if (intercept)
        sum(weights * y)/sum(weights)
    else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights ==  0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- if (EMPTY)
        0
    else fit$rank
    resdf <- n.ok - rank
    aic.model <- aic(y, length(y), mu, weights, dev) + 2 * rank
    list(coefficients = coef, residuals = residuals, fitted.values = mu,
         effects = if (!EMPTY) fit$effects, R = if (!EMPTY) Rmat,
         rank = rank, qr = if (!EMPTY) structure(fit[c("qr", "rank",
                                                       "qraux", "pivot", "tol")], class = "qr"), family = family,
         linear.predictors = eta, deviance = dev, aic = aic.model,
         null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights,
         df.residual = resdf, df.null = nulldf, y = y, converged = conv,
         boundary = boundary)
}

# copied from limma
weighted.median.scde <- function (x, w, na.rm = FALSE)
    #       Weighted median
    #       Gordon Smyth
    #       30 June 2005
{
    if (missing(w))
        w <- rep.int(1, length(x))
    else {
        if(length(w)  !=  length(x)) stop("'x' and 'w' must have the same length")
        if(any(is.na(w))) stop("NA weights not allowed")
        if(any(w<0)) stop("Negative weights not allowed")
    }
    if(is.integer(w))
        w <- as.numeric(w)
    if(na.rm) {
        w <- w[i <- !is.na(x)]
        x <- x[i]
    }
    if(all(w == 0)) {
        warning("All weights are zero")
        return(NA)
    }
    o <- order(x)
    x <- x[o]
    w <- w[o]
    p <- cumsum(w)/sum(w)
    n <- sum(p<0.5)
    if(p[n+1]  >  0.5)
        x[n+1]
    else
        (x[n+1]+x[n+2])/2
}

# FROM common.r
sn <- function(x) {
    names(x) <- x
    return(x)
}

# panel routines for pairs()
pairs.panel.hist <- function(x, i = NULL, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks
    nB <- length(breaks)
    y <- h$counts
    y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "gray70", ...)
}
pairs.panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, i = NULL, j = NULL) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, method = "pearson"))
    #r <- abs(cor(x, y, method = "spearman"))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")
    if(missing(cex.cor)) { cex <- 0.6/strwidth(txt) }
    #text(0.5, 0.5, txt, cex = cex * r)
    text(0.5, 0.5, txt, cex = cex)
}
pairs.panel.scatter <- function(x, y, i = NULL, j = NULL, ...) {
    vi <- x > 0 | y > 0
    points(x[vi], y[vi], pch = ".", col = densCols(x[vi], y[vi], colramp = colorRampPalette(brewer.pal(9, "Blues")[-(1:2)])), cex = 2)
}
pairs.panel.smoothScatter <- function(x, y, i = NULL, j = NULL, ...) {
    vi <- x > 0 | y > 0
    smoothScatter(x[vi], y[vi], add = TRUE, ...)
}

# a slight modification of pairs that passes i/j indices to the panel methods
pairs.extended <- function (x, labels, panel = points, ...,
                            lower.panel = panel, upper.panel = panel,
                            diag.panel = NULL, text.panel = textPanel,
                            label.pos = 0.5 + has.diag/3,
                            cex.labels = NULL, font.labels = 1,
                            row1attop = TRUE, gap = 1)
{
    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) {
        text(x, y, txt, cex = cex, font = font)
    }

    localAxis <- function(side, x, y, xpd, bg, col = NULL, main, oma, ...) {
        ## Explicitly ignore any color argument passed in as
        ## it was most likely meant for the data points and
        ## not for the axis.
        if(side %%2 ==  1) { Axis(x, side = side, xpd = NA, ...) }
        else { Axis(y, side = side, xpd = NA, ...) }
    }

    localPlot <- function(..., main, oma, font.main, cex.main) { plot(...) }
    localLowerPanel <- function(..., main, oma, font.main, cex.main) { lower.panel(...) }
    localUpperPanel <- function(..., main, oma, font.main, cex.main) { upper.panel(...) }
    localDiagPanel <- function(..., main, oma, font.main, cex.main) { diag.panel(...) }

    dots <- list(...)
    nmdots <- names(dots)
    if (!is.matrix(x)) {
        x <- as.data.frame(x)
        for(i in seq_along(names(x))) {
            if(is.factor(x[[i]]) || is.logical(x[[i]])) {
                x[[i]] <- as.numeric(x[[i]])
            }
            if(!is.numeric(unclass(x[[i]]))) {
                stop("non-numeric argument to 'pairs'")
            }
        }
    } else if(!is.numeric(x)) {
        stop("non-numeric argument to 'pairs'")
    }
    panel <- match.fun(panel)
    if((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) {
        lower.panel <- match.fun(lower.panel)
    }
    if((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) {
        upper.panel <- match.fun(upper.panel)
    }
    if((has.diag  <- !is.null( diag.panel)) && !missing( diag.panel)) {
        diag.panel <- match.fun( diag.panel)
    }

    if(row1attop) {
        tmp <- lower.panel
        lower.panel <- upper.panel
        upper.panel <- tmp
        tmp <- has.lower
        has.lower <- has.upper
        has.upper <- tmp
    }

    nc <- ncol(x)
    if (nc < 2) stop("only one column in the argument to 'pairs'")
    has.labs <- TRUE
    if (missing(labels)) {
        labels <- colnames(x)
        if (is.null(labels)) {
            labels <- paste("var", 1L:nc)
        }
    } else if(is.null(labels)) {
        has.labs <- FALSE
    }
    oma <- if("oma" %in% nmdots) {
        dots$oma
    } else {
        NULL
    }
    main <- if("main" %in% nmdots) {
        dots$main
    } else {
        NULL
    }
    if (is.null(oma)) {
        oma <- c(4, 4, 4, 4)
        if (!is.null(main)) {
            oma[3L] <- 6
        }
    }
    opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
    on.exit(par(opar))

    for (i in if(row1attop) 1L:nc else nc:1L)
        for (j in 1L:nc) {
            localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, type = "n", ...)
            if(i ==  j || (i < j && has.lower) || (i  >  j && has.upper) ) {
                box()
                if(i ==  1  && (!(j %% 2) || !has.upper || !has.lower ))
                    localAxis(1 + 2*row1attop, x[, j], x[, i], ...)
                if(i ==  nc && (  j %% 2  || !has.upper || !has.lower ))
                    localAxis(3 - 2*row1attop, x[, j], x[, i], ...)
                if(j ==  1  && (!(i %% 2) || !has.upper || !has.lower ))
                    localAxis(2, x[, j], x[, i], ...)
                if(j ==  nc && (  i %% 2  || !has.upper || !has.lower ))
                    localAxis(4, x[, j], x[, i], ...)
                mfg <- par("mfg")
                if(i ==  j) {
                    if (has.diag) localDiagPanel(as.vector(x[, i]), i = i, ...)
                    if (has.labs) {
                        par(usr = c(0, 1, 0, 1))
                        if(is.null(cex.labels)) {
                            l.wid <- strwidth(labels, "user")
                            cex.labels <- max(0.8, min(2, .9 / max(l.wid)))
                        }
                        text.panel(0.5, label.pos, labels[i],
                                   cex = cex.labels, font = font.labels)
                    }
                } else if(i < j)
                    localLowerPanel(as.vector(x[, j]), as.vector(x[, i]), i = i, j = j, ...)
                else
                    localUpperPanel(as.vector(x[, j]), as.vector(x[, i]), i = i, j = j, ...)
                if (any(par("mfg")  !=  mfg))
                    stop("the 'panel' function made a new plot")
            } else {
                par(new = FALSE)
            }
        }
    if (!is.null(main)) {
        font.main <- if("font.main" %in% nmdots) {
            dots$font.main
        } else {
            par("font.main")
        }
        cex.main <- if("cex.main" %in% nmdots) {
            dots$cex.main
        } else {
            par("cex.main")
        }
        mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
    }
    invisible(NULL)
}


# given a set of pdfs (columns), calculate summary statistics (mle, 95% CI, Z-score deviations from 0)
quick.distribution.summary <- function(s.bdiffp) {
    diffv <- as.numeric(colnames(s.bdiffp))
    dq <- t(apply(s.bdiffp, 1, function(p) {
        mle <- which.max(p)
        p <- cumsum(p)
        return(diffv[c(lb = max(c(1, which(p<0.025))), mle, min(c(length(p), which(p > (1-0.025)))))])
    }))/log10(2)
    colnames(dq) <- c("lb", "mle", "ub")
    cq <- rep(0, nrow(dq))
    cq[dq[, 1] > 0] <- dq[dq[, 1] > 0, 1]
    cq[dq[, 3]<0] <- dq[dq[, 3]<0, 3]
    z <- get.ratio.posterior.Z.score(s.bdiffp)
    za <- sign(z)*qnorm(p.adjust(pnorm(abs(z), lower.tail = FALSE), method = "BH"), lower.tail = FALSE)
    data.frame(dq, "ce" = as.numeric(cq), "Z" = as.numeric(z), "cZ" = as.numeric(za))
}


#######
## INTERNAL PAGODA ROUTINES
#######

# performs weighted centering of mat rows (mat - rowSums(mat*weights)/rowSums(weights))
# possibly accounting for batch effects (i.e. centering each batch separately
weightedMatCenter <- function(mat, matw, batch = NULL) {
    if(is.null(batch)) {
        return(mat-rowSums(mat*matw)/rowSums(matw))
    } else {
        cmat <- mat
        invisible(tapply(seq_len(ncol(mat)), as.factor(batch), function(ii) {
            cmat[, ii] <<- cmat[, ii, drop = FALSE] - rowSums(cmat[, ii, drop = FALSE]*matw[, ii, drop = FALSE])/rowSums(matw[, ii, drop = FALSE])
        }))
        return(cmat)
    }
}

# per-experiment/per-gene weighted variance estimate
# weight matrix should have the same dimensions as the data matrix
weightedMatVar <- function(mat, matw, batch = NULL, center = TRUE, min.weight = 0, normalize.weights = TRUE) {
    # normalize weights
    #matw <- matw/rowSums(matw)
    #matw <- matw/rowSums(matw)*ncol(matw)
    if(center) {
        mat <- weightedMatCenter(mat, matw, batch)
    }

    #weightedMatVar.Rcpp(mat, matw)
    #return(rowSums(mat*mat*matw) / (1-rowSums(matw*matw)))
    #return(rowSums(mat*mat*matw))

    #return(rowSums(mat*mat*matw) * rowSums(matw) /pmax(rowSums(matw)^2 - rowSums(matw*matw), rep(min.weight, nrow(matw))))

    v<- rowSums(mat*mat*matw)
    if(normalize.weights) { v <- v/rowSums(matw) }
    v
}

# GEV t() function
gev.t <- function(x, loc, scale, shape = rep(0, length(loc)), log = FALSE) {
    if(log) {
        pmin(0, ifelse(shape == 0, -(x-loc)/scale, (-1/shape)*log(pmax(0, 1+shape*(x-loc)/scale))))
    } else {
        pmin(1, ifelse(shape == 0, exp(-(x-loc)/scale), ((pmax(0, 1+shape*(x-loc)/scale))^(-1/shape))))
    }
}
# returns upper tail of GEV in log scale
pgev.upper.log <- function(x, loc, scale, shape = rep(0, length(loc))) {
    tv <- gev.t(x, loc, scale, shape, log = TRUE)
    tv[tv >  -5 & tv<0] <- log(-expm1(-exp(tv[tv >  -5 & tv<0])))
    tv
}

# BH P-value adjustment with a log option
bh.adjust <- function(x, log = FALSE) {
    nai <- which(!is.na(x))
    ox <- x
    x<-x[nai]
    id <- order(x, decreasing = FALSE)
    if(log) {
        q <- x[id] + log(length(x)/seq_along(x))
    } else {
        q <- x[id]*length(x)/seq_along(x)
    }
    a <- rev(cummin(rev(q)))[order(id)]
    ox[nai]<-a
    ox
}

pathway.pc.correlation.distance <- function(pcc, xv, n.cores = 10, target.ndf = NULL) {
    # all relevant gene names
    rotn <- unique(unlist(lapply(pcc[gsub("^#PC\\d+# ", "", rownames(xv))], function(d) rownames(d$xp$rotation))))
    # prepare an ordered (in terms of genes) and centered version of each component
    pl <- papply(rownames(xv), function(nam) {
        pnam <- gsub("^#PC\\d+# ", "", nam)
        pn <- as.integer(gsub("^#PC(\\d+)# .*", "\\1", nam))
        rt <- pcc[[pnam]]$xp$rotation[, pn]
        # order names/values according to increasing name match index
        mi <- match(names(rt), rotn)
        mo <- order(mi, decreasing = FALSE)
        rt <- as.numeric(rt)-mean(rt)
        return(list(i = mi[mo], v = rt[mo]))
    }, n.cores = n.cores)

    x <- .Call("plSemicompleteCor2", pl, PACKAGE = "scde")

    if(!is.null(target.ndf)) {
        r <- x$r[upper.tri(x$r)]
        n <- x$n[upper.tri(x$n)]
        suppressWarnings(tv <- r*sqrt((n-2)/(1-r^2)))
        z <- pt(tv, df = n-2, lower.tail = FALSE, log.p = TRUE)
        nr <- qt(z, df = target.ndf-2, lower.tail = FALSE, log.p = TRUE)
        nr <- nr/sqrt(target.ndf-2+nr^2)
        nr[is.nan(nr)] <- r[is.nan(nr)]

        cr <- x$r
        cr[upper.tri(cr)] <- nr
        cr[lower.tri(cr)] <- t(cr)[lower.tri(cr)]
    } else {
        cr <- x$r
    }

    rownames(cr) <- colnames(cr) <- rownames(xv)
    d <- stats::as.dist(1-abs(cr))
    d[d<0] <- 0
    d

}

collapse.aspect.clusters <- function(d, dw, ct, scale = TRUE, pick.top = FALSE) {
    xvm <- do.call(rbind, tapply(seq_len(nrow(d)), factor(ct, levels = sort(unique(ct))), function(ii) {
        if(length(ii) == 1) return(d[ii, ])
        if(pick.top) {
            return(d[ii[which.max(apply(d[ii, ], 1, var))], ])
        }
        xp <- pcaMethods::pca(t(d[ii, ]), nPcs = 1, center = TRUE, scale = "none")
        xv <- pcaMethods::scores(xp)[, 1]
        if(sum(abs(diff(xv))) > 0 && cor(xv, colMeans(d[ii, ]*abs(pcaMethods::loadings(xp)[, 1])))<0) { xv <- -1*xv }
        #set scale at top pathway?
        if(sum(abs(diff(xv))) > 0) {
            if(scale) {
                xv <- xv*sqrt(max(apply(d[ii, ], 1, var)))/sqrt(var(xv))
            }
            if(sum(abs(xv)) == 0) { xv <- abs(rnorm(length(xv), sd = 1e-6)) }
        } else {
            xv <- abs(rnorm(length(xv), sd = 1e-6))
        }
        #xv <- xv/sqrt(length(ii))
        xv
    }))
    rownames(xvm) <- unlist(tapply(seq_len(nrow(d)), factor(ct, levels = sort(unique(ct))), function(ii) {
        if(length(ii) == 1) return(rownames(d)[ii])
        return(rownames(d)[ii[which.max(apply(d[ii, ], 1, var))]])
    }))

    xvmw <- do.call(rbind, tapply(seq_len(nrow(d)), factor(ct, levels = sort(unique(ct))), function(ii) {
        w <- colSums(dw[ii, , drop = FALSE]*apply(d[ii, , drop = FALSE], 1, sd))
        w <- w/sum(w)
    }))

    return(list(d = xvm, w = xvmw))
}
# convert R color to a web hex representation
col2hex <- function(col) {
    unlist(lapply(col, function(c) {
        c <- col2rgb(c)
        sprintf("#%02X%02X%02X", c[1], c[2], c[3])
    }))
}

my.heatmap2 <- function (x, Rowv = NULL, Colv = if(symm)"Rowv" else NULL,
                         distfun = dist, hclustfun = hclust,
                         reorderfun = function(d, w) reorder(d, w),
                         add.expr, symm = FALSE, revC = identical(Colv, "Rowv"),
                         scale = c("none", "row", "column"), na.rm = TRUE,
                         margins = c(5, 5), internal.margin = 0.5, ColSideColors, RowSideColors,
                         cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc),
                         labRow = NULL, labCol = NULL,  xlab = NULL, ylab = NULL,
                         keep.dendro = FALSE,
                         grid = FALSE, grid.col = 1, grid.lwd = 1,
                         verbose = getOption("verbose"), Colv.vsize = 0.15, Rowv.hsize = 0.15, ColSideColors.unit.vsize = "0.08", RowSideColors.hsize = 0.03, lasCol = 2, lasRow = 2, respect = FALSE, box = FALSE, zlim = NULL, ...)
{
    scale <- if(symm && missing(scale)) "none" else match.arg(scale)
    if(length(di <- dim(x))  !=  2 || !is.numeric(x))
        stop("'x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if(nr < 1 || nc <=  1)
        stop("'x' must have at least one row and 2 columns")
    if(!is.numeric(margins) || length(margins)  !=  2)
        stop("'margins' must be a numeric vector of length 2")

    if(is.null(zlim)) {
        zlim <- range(x[is.finite(x)])
    } else {
        x[x<zlim[1]] <- zlim[1]
        x[x > zlim[2]] <- zlim[2]
    }

    doRdend <- !identical(Rowv, NA)
    doCdend <- !identical(Colv, NA)
    ## by default order by row/col means
    if(is.null(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
    }
    if(is.null(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
    }

    ## get the dendrograms and reordering indices

    if(doRdend) {
        if(inherits(Rowv, "dendrogram"))
            ddr <- Rowv
        else {
            hcr <- hclustfun(distfun(x))
            ddr <- as.dendrogram(hcr)
            if(!is.logical(Rowv) || Rowv)
                ddr <- reorderfun(ddr, Rowv)
        }
        if(nr  !=  length(rowInd <- order.dendrogram(ddr)))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else rowInd <- 1:nr

    if(doCdend) {
        if(inherits(Colv, "dendrogram"))
            ddc <- Colv
        else if(identical(Colv, "Rowv")) {
            if(nr  !=  nc)
                stop('Colv = "Rowv" but nrow(x)  !=  ncol(x)')
            ddc <- ddr
        }
        else {
            hcc <- hclustfun(distfun(if(symm)x else t(x)))
            ddc <- as.dendrogram(hcc)
            if(!is.logical(Colv) || Colv)
                ddc <- reorderfun(ddc, Colv)
        }
        if(nc  !=  length(colInd <- order.dendrogram(ddc)))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else colInd <- 1:nc

    ## reorder x
    x <- x[rowInd, colInd, drop = FALSE]

    labRow <-
        if(is.null(labRow))
            if(is.null(rownames(x))) (1:nr)[rowInd] else rownames(x)
    else labRow[rowInd]
    labCol <-
        if(is.null(labCol))
            if(is.null(colnames(x))) (1:nc)[colInd] else colnames(x)
    else labCol[colInd]

    if(scale ==  "row") {
        x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
        sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if(scale ==  "column") {
        x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
        sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }

    ## Calculate the plot layout
    ds <- dev.size(units = "cm")


    lmat <- rbind(c(NA, 3), 2:1)
    if(doRdend) {
        lwid <- c(if(is.character(Rowv.hsize)) Rowv.hsize else lcm(Rowv.hsize*ds[1]), 1)
    } else {
        lmat[2, 1] <- NA
        lmat[1, 2] <- 2
        lwid <- c(0, 1)
    }
    if(doCdend) {
        lhei <- c(if(is.character(Colv.vsize)) Colv.vsize else lcm(Colv.vsize*ds[2]), 1)
    } else {
        lmat[1, 2] <- NA
        lhei <- c(0, 1)
    }
    #lwid <- c(if(doRdend) lcm(Rowv.hsize*ds[1]) else "0.5 cm", 1)
    #lhei <- c((if(doCdend) lcm(Colv.vsize*ds[2]) else "0.5 cm"), 1)
    if(!missing(ColSideColors) && !is.null(ColSideColors)) { ## add middle row to layout

        if(is.matrix(ColSideColors)) {
            if(ncol(ColSideColors) != nc)
                stop("'ColSideColors' matrix must have the same number of columns as length ncol(x)")
            if(is.character(ColSideColors.unit.vsize)) {
                ww <- paste(as.numeric(gsub("(\\d+\\.?\\d*)(.*)", "\\1", ColSideColors.unit.vsize, perl = TRUE))*nrow(ColSideColors), gsub("(\\d+\\.?\\d*)(.*)", "\\2", ColSideColors.unit.vsize, perl = TRUE), sep = "")
            } else {
                ww <- lcm(ColSideColors.unit.vsize*ds[2]*nrow(ColSideColors))
            }
            lmat <- rbind(lmat[1, ]+1, c(NA, 1), lmat[2, ]+1)
            lhei <- c(lhei[1], ww, lhei[2])
        } else {
            if(!is.character(ColSideColors) || length(ColSideColors)  !=  nc)
                stop("'ColSideColors' must be a character vector of length ncol(x)")
            if(is.character(ColSideColors.unit.vsize)) {
                ww <- paste(as.numeric(gsub("(\\d+\\.?\\d*)(.*)", "\\1", ColSideColors.unit.vsize, perl = TRUE)), gsub("(\\d+\\.?\\d*)(.*)", "\\2", ColSideColors.unit.vsize, perl = TRUE), sep = "")
            } else {
                ww <- lcm(ColSideColors.unit.vsize*ds[2])
            }
            lmat <- rbind(lmat[1, ]+1, c(NA, 1), lmat[2, ]+1)
            lhei <- c(lhei[1], ww, lhei[2])
        }
    }
    if(!missing(RowSideColors) && !is.null(RowSideColors)) { ## add middle column to layout
        if(!is.character(RowSideColors) || length(RowSideColors)  !=  nr)
            stop("'RowSideColors' must be a character vector of length nrow(x)")
        lmat <- cbind(lmat[, 1]+1, c(rep(NA, nrow(lmat)-1), 1), lmat[, 2]+1)
        lwid <- c(lwid[1], if(is.character(RowSideColors.hsize)) RowSideColors.hsize else lcm(RowSideColors.hsize*ds[1]), lwid[2])
    }
    lmat[is.na(lmat)] <- 0
    if(verbose) {
        cat("layout: widths = ", lwid, ", heights = ", lhei, " lmat = \n")
        print(lmat)
    }

    ## Graphics `output' -----------------------

    op <- par(no.readonly = TRUE)
    #on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = respect)
    ## draw the side bars
    if(!missing(RowSideColors) && !is.null(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, internal.margin))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        if(box) { box() }
    }
    if(!missing(ColSideColors) && !is.null(ColSideColors)) {
        par(mar = c(internal.margin, 0, 0, margins[2]))
        if(is.matrix(ColSideColors)) {
            image(t(matrix(1:length(ColSideColors), byrow = TRUE, nrow = nrow(ColSideColors), ncol = ncol(ColSideColors))), col = as.vector(t(ColSideColors[, colInd, drop = FALSE])), axes = FALSE)
            if(box) { box() }
        } else {
            image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
            if(box) { box() }
        }
    }
    ## draw the main carpet
    par(mar = c(margins[1], 0, 0, margins[2]))
    if(!symm || scale  !=  "none")
        x <- t(x)
    if(revC) { # x columns reversed
        iy <- nr:1
        ddr <- rev(ddr)
        x <- x[, iy, drop = FALSE]
    } else iy <- 1:nr

    image(1:nc, 1:nr, x, xlim = 0.5+ c(0, nc), ylim = 0.5+ c(0, nr),
          axes = FALSE, xlab = "", ylab = "", zlim = zlim, ...)
    if(box) { box() }
    axis(1, 1:nc, labels =  labCol, las =  lasCol, line =  -0.5, tick =  0, cex.axis =  cexCol)
    if(!is.null(xlab)) mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels =  labRow, las =  lasRow, line =  -0.5, tick =  0, cex.axis =  cexRow)
    if(!is.null(ylab)) mtext(ylab, side = 4, line = margins[2] - 1.25, las = lasRow)
    if (!missing(add.expr))
        eval(substitute(add.expr))


    if(grid) {
        abline(v = c(1:nc)-0.5, col = grid.col, lwd = grid.lwd)
        abline(h = c(1:nr)-0.5, col = grid.col, lwd = grid.lwd)
        box(col = grid.col, lwd = grid.lwd)
    }

    ## the two dendrograms :
    if(doRdend) {
        par(mar = c(margins[1], internal.margin, 0, 0))
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none", xaxs = "i")
    }

    if(doCdend) {
        par(mar = c(0, 0, internal.margin, margins[2]))
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none", yaxs = "i")
    }
    invisible(list(rowInd = rowInd, colInd = colInd,
                   Rowv = if(keep.dendro && doRdend) ddr,
                   Colv = if(keep.dendro && doCdend) ddc ))
}


# rook class for browsing differential expression results

ViewDiff <- setRefClass(
    'ViewDiff',
    fields = c('gt', 'models', 'counts', 'prior', 'groups', 'batch', 'geneLookupURL'),
    methods = list(

        initialize = function(results, models, counts, prior, groups = NULL, batch = NULL, geneLookupURL = NULL) {
            if(!is.null(results$results)) {
                gt <<- results$results
            } else {
                gt <<- results
            }
            # add raw names if this wasn't a batch-corrected sample
            if("mle" %in% colnames(gt)) {
                colnames(gt) <<- paste(colnames(gt), "raw", sep = "_")
            }
            if(!is.null(results$batch.adjusted)) {
                df <- results$batch.adjusted
                colnames(df) <- paste(colnames(df), "cor", sep = "_")
                gt <<- cbind(gt, df)
            }
            if(!is.null(results$batch.effect)) {
                df <- results$batch.effect
                colnames(df) <- paste(colnames(df), "bat", sep = "_")
                gt <<- cbind(gt, df)
            }
            colnames(gt) <<- tolower(colnames(gt))

            # append expression levels to the results
            if(!is.null(results$joint.posteriors)) {
                gt$lev1 <<- log10(as.numeric(colnames(results$joint.posteriors[[1]]))[max.col(results$joint.posteriors[[1]])]+1)
                gt$lev2 <<- log10(as.numeric(colnames(results$joint.posteriors[[2]]))[max.col(results$joint.posteriors[[2]])]+1)
            }
            gt$gene <<- rownames(gt)
            gt <<- data.frame(gt)

            # guess gene lookup for common cases
            if(is.null(geneLookupURL)) {
                # human
                if( any(grepl("ENSG\\d+", gt$gene[1])) || any(c("CCLU1", "C22orf45") %in% gt$gene)) {
                    geneLookupURL <<- "http://useast.ensembl.org/Homo_sapiens/Gene/Summary?g = {0}"
                } else if( any(grepl("ENSMUSG\\d+", gt$gene[1]))) {
                    geneLookupURL <<- "http://useast.ensembl.org/Mus_musculus/Gene/Summary?g = {0}"
                } else if( any(c("Foxp2", "Sept1", "Lrrc34") %in% gt$gene)) {
                    # mouse MGI
                    geneLookupURL <<- "http://www.informatics.jax.org/searchtool/Search.do?query = {0}"
                } else if( any(grepl("FBgn\\d+", gt$gene[1])) || any(c("CG3680", "CG8290") %in% gt$gene)) {
                    # flybase
                    geneLookupURL <<- "http://flybase.org/cgi-bin/uniq.html?db = fbgn&GeneSearch = {0}&context = {1}&species = Dmel&cs = yes&caller = genejump"
                } else {
                    # default, forward to ensemble search
                    geneLookupURL <<- "http://useast.ensembl.org/Multi/Search/Results?q = {0}site = ensembl"
                }
            } else {
                geneLookupURL <<- geneLookupURL
            }



            gt <<- gt[gt$z_raw != "NA", ]
            gt <<- gt[!is.na(gt$z_raw), ]


            models <<- models
            counts <<- counts
            prior <<- prior
            if(is.null(groups)) { # recover groups from models
                groups <<- as.factor(attr(models, "groups"))
                if(is.null(groups)) stop("ERROR: groups factor is not provided, and models structure is lacking groups attribute")
                names(groups) <<- rownames(models)
            } else {
                groups <<- groups
            }
            if(length(levels(groups)) != 2) {
                stop(paste("ERROR: wrong number of levels in the grouping factor (", paste(levels(groups), collapse = " "), "), but must be two.", sep = ""))
            }

            batch <<- batch
            callSuper()
        },
        call = function(env){
            path <- env[['PATH_INFO']]
            req <- Request$new(env)
            res <- Response$new()

            switch(path,
                   # INDEX
                   '/index.html' = {
                       body <- paste('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd" >
                                     <html >
                                     <head >
                                     <meta http-equiv = "Content-Type" content = "text/html charset = iso-8859-1" >
                                     <title > SCDE: ', paste(levels(groups), collapse = " vs. "), '</title >
                                     <!-- ExtJS -- >
                                     <link rel = "stylesheet" type = "text/css" href = "http://pklab.med.harvard.edu/sde/extjs/resources/ext-theme-neptune/ext-theme-neptune-all.css" / >

                                     <!-- Shared -- >
                                     <link rel = "stylesheet" type = "text/css" href = "http://pklab.med.harvard.edu/sde/ext-4.2.1.883/examples/shared/example.css" / >

                                     <link rel = "stylesheet" type = "text/css" href = "http://pklab.med.harvard.edu/sde/additional.css" / >
                                     <!-- GC -- >

                                     <style type = "text/css" >
                                     .x-panel-framed {
                                     padding: 0
                                     }
                                     </style >
                                     <script type = "text/javascript" src = "http://pklab.med.harvard.edu/sde/ext-4.2.1.883/ext-all.js" > </script >

                                     <script type = "text/javascript" > var geneLookupURL = "', geneLookupURL, '"</script >
                                     <script type = "text/javascript" src = "http://pklab.med.harvard.edu/sde/viewembed.js" > </script >

                                     </head >
                                     <body style = "margin-top:0padding-top:10px" >
                                     <div id = "example-grid" > </div >
                                     </body >
                                     </html >
                                     ', sep = "")
                       res$header('"Content-Type": "text/html"')
                       res$write(body)
                   },
                   # GENE TABLE
                   '/genetable.json' = {
                       lgt <- gt
                       if(!is.null(req$params()$filter)) {
                           fl <- rjson::fromJSON(URLdecode(req$params()$filter))
                           for( fil in fl) {
                               lgt <- lgt[grep(fil$value, lgt[, fil$property], perl = TRUE, ignore.case = TRUE), ]
                           }
                       }
                       start <- ifelse(is.null(req$params()$start), 1, as.integer(req$params()$start)+1)
                       limit <- ifelse(is.null(req$params()$limit), 1000, as.integer(req$params()$limit))
                       dir <- ifelse(is.null(req$params()$dir), "DESC", req$params()$dir)
                       trows <- nrow(lgt)
                       if(trows > 0) {
                           if(!is.null(req$params()$sort)) {
                               if(req$params()$sort %in% colnames(lgt)) {
                                   lgt <- lgt[order(lgt[, req$params()$sort], decreasing = (dir == "DESC")), ]
                               }
                           } else { # default sort
                               if(is.null(lgt$z_cor)) { lgt <- lgt[order(abs(lgt$z_raw), decreasing = TRUE), ] } else { lgt <- lgt[order(abs(lgt$z_cor), decreasing = TRUE), ] }
                           }
                       }
                       lgt <- format(lgt[min(start, nrow(lgt)):min((start+limit), nrow(lgt)), ], nsmall = 2, digits = 2)
                       ol <- apply(lgt, 1, function(x) as.list(x))
                       names(ol) <- NULL
                       s <- rjson::toJSON(list(totalCount = trows, genes = ol))
                       res$header('"Content-Type": "application/json"')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }

                   },
                   # POSTERIOR PLOT
                   '/posterior.png' = {
                       gene <- ifelse(is.null(req$params()$gene), sample(gt$gene), req$params()$gene)
                       bootstrap <- ifelse(is.null(req$params()$bootstrap), TRUE, req$params()$bootstrap == "T")
                       show.individual.posteriors <- ifelse(is.null(req$params()$show.individual.posteriors), TRUE, req$params()$show.individual.posteriors == "true")

                       t <- tempfile()
                       #require(Cairo)
                       CairoPNG(filename = TRUE, width = 350, height = 560)
                       scde.test.gene.expression.difference(gene = gene, models = models, counts = counts, groups = groups, prior = prior, batch = batch, ratio.range = c(-10, 10), show.individual.posteriors = show.individual.posteriors, verbose = FALSE)
                       dev.off()
                       res$header('Content-type', 'image/png')
                       res$body <- t
                       names(res$body) <- 'file'
                   },
                   # GENE EXPRESSION LEVELS
                   '/elevels.html' = {
                       geneName <- ifelse(is.null(req$params()$geneName), gt$gene[[1]], req$params()$geneName)
                       gc <- counts[rownames(counts) == geneName, , drop = FALSE]
                       fpm <- exp(scde.expression.magnitude(models, counts = gc))
                       df <- rbind(FPM = gc, level = fpm)
                       df <- round(df, 2)
                       # order columns according to groups
                       df <- df[, unlist(tapply(seq_len(ncol(df)), groups, I))]
                       cell.col <- rep(c("#E9A994", "#66CCFF"), as.integer(table(groups)))

                       render.row <- function(nam, val, col) {
                           paste("<tr > ", "<th > ", nam, "</th > ", paste("<td bgcolor = ", col, " > ", val, "</td > ", sep = "", collapse = " "), "</tr > ", sep = "")
                       }

                       sh <- paste("<tr > ", paste("<th > ", c(" ", colnames(df)), "</th > ", sep = "", collapse = " "), "</tr > ")
                       #sb <- paste(render.row("cells", colnames(df), cell.col), render.row("FPKM", df[1, ], cell.col), render.row("mode", df[2, ], cell.col), collapse = "\n")
                       sb <- paste(render.row("counts", df[1, ], cell.col), render.row("FPM", df[2, ], cell.col), collapse = "\n")
                       res$header('"Content-Type": "text/html"')
                       res$write(paste("<table id = \"elevels\" > ", sh, sb, "</table > "))
                   },
{
    res$write('default')
}
                       )
            res$finish()
        }
            )
    )

t.view.pathways <- function(pathways, mat, matw, env, proper.names = rownames(mat), colcols = NULL, zlim = NULL, labRow = NA, vhc = NULL, cexCol = 1, cexRow = 1, n.pc = 1, nstarts = 50, row.order = NULL, show.Colv = TRUE, plot = TRUE, trim = 1.1/ncol(mat), bwpca = TRUE, ...) {
    # retrieve gis
    lab <- which(proper.names %in% na.omit(unlist(mget(pathways, envir = env, ifnotfound = NA))))

    if(length(lab) == 0) {
        # try genes
        lab <- which(proper.names %in% pathways)
    }
    if(length(lab) == 0)
        return(NULL)
    #t.quick.show.mat(mat[lab, ], normalize.rows = TRUE)
    #table(rownames(mat) %in% mget(pathways, envir = env))

    if(trim > 0) {
        mat <- winsorize.matrix(mat, trim = trim)
    }
    d <- mat[lab, , drop = FALSE]
    dw <- matw[lab, , drop = FALSE]
    if(length(lab) > 2) {
        xp <- bwpca(t(d), t(dw), npcs = n.pc, center = FALSE, nstarts = 3)
    } else {
        xp <- list()
    }

    d <- d-rowMeans(d)
    dd <- as.dist(1-abs(cor(t(as.matrix(d)))))
    dd[is.na(dd)] <- 1
    if(is.null(row.order)) {
        if(length(lab) > 2) {
            if(is.element("fastcluster", installed.packages()[, 1])) {
                hc <- fastcluster::hclust(dd, method = "ward.D")
            } else {
                hc <- stats::hclust(dd, method = "ward.D")
            }
            row.order <- hc$order
        } else {
            row.order <- c(seq_along(lab))
        }
    }

    if(is.null(vhc)) {
        vd <- as.dist(1-cor(as.matrix(d)))
        vd[is.na(vd)] <- 1
        if(is.element("fastcluster", installed.packages()[, 1])) {
            vhc <- fastcluster::hclust(vd, method = "ward.D")
        } else {
            vhc <- stats::hclust(vd, method = "ward.D")
        }

    }

    if(is.null(zlim)) { zlim <- quantile(d, p = c(0.01, 0.99)) }
    vmap <- d
    vmap[vmap<zlim[1]] <- zlim[1]
    vmap[vmap > zlim[2]] <- zlim[2]
    rownames(vmap) <- rownames(d)

    aval <- colSums(d*dw)/colSums(dw)
    if(!is.null(xp$scores)) {
        oc <- xp$scores[, n.pc]
        if(cor(oc, aval, method = "spearman")<0) {
            oc <- -1*oc
            xp$scores[, n.pc] <- -1*xp$scores[, n.pc]
            xp$rotation[, n.pc] <- -1*xp$rotation[, n.pc]
        }
        xp$oc <- oc
        z <- colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(100)[round(oc/max(abs(oc))*49)+50]
        ld <- xp$rotation[row.order, , drop = FALSE]
        ld <- colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(100)[round(ld/max(abs(ld))*49)+50]
    } else {
        ld <- z <- NULL
    }

    aval <- colorRampPalette(c("blue", "white", "red"), space = "Lab")(100)[round(aval/max(abs(aval))*49)+50]



    # oc2 <- xp$scores[, 2]
    # oc2 <- colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(100)[round(oc2/max(abs(oc2))*49)+50]
    # oc3 <- xp$scores[, 3]
    # oc3 <- colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(100)[round(oc3/max(abs(oc3))*49)+50]


    #z <- do.call(rbind, list(aval, oc))
    #z <- rbind(oc3, oc2, oc)

    col <- colorRampPalette(c("blue", "white", "red"), space = "Lab")(256)
    if(plot) {
        if(!is.null(colcols)) { z <- rbind(colcols, z) }
        if(show.Colv) { Colv <- as.dendrogram(vhc) } else { Colv <- NA }
        my.heatmap2(vmap[row.order, , drop = FALSE], Rowv = NA, Colv = Colv, zlim = zlim, col = col, scale = "none", RowSideColors = ld, ColSideColors = z, labRow = labRow, cexCol = cexCol, cexRow = cexRow, ...)
    }
    xp$vhc <- vhc
    xp$lab <- lab
    xp$row.order <- row.order
    xp$col <- col
    xp$oc.col <- colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(256)
    xp$vmap <- vmap
    xp$zlim <- zlim
    return(invisible(xp))
}

##' View pathway or gene weighted PCA
##'
##' Takes in a list of pathways (or a list of genes), runs weighted PCA, optionally showing the result.
##' @param pathways character vector of pathway or gene names
##' @param varinfo output of pagoda.varnorm()
##' @param goenv environment mapping pathways to genes
##' @param n.genes number of genes to show
##' @param two.sided whether the set of shown genes should be split among highest and lowest loading (T) or if genes with highest absolute loading (F) should be shown
##' @param n.pc optional integer vector giving the number of principal component to show for each listed pathway
##' @param colcols optional column color matrix
##' @param zlim optional z color limit
##' @param showRowLabels controls whether row labels are shown in the plot
##' @param cexCol column label size (cex)
##' @param cexRow row label size (cex)
##' @param nstarts number of random starts for the wPCA
##' @param cell.clustering cell clustering
##' @param show.cell.dendrogram whether cell dendrogram should be shown
##' @param plot whether the plot should be shown
##' @param box whether to draw a box around the plotted matrix
##' @param trim optional Winsorization trim that should be applied
##' @param return.details whether the function should return the matrix as well as full PCA info instead of just PC1 vector
##' @param ... additional arguments are passed to the \code{c.view.pathways}
##' @return cell scores along the first principal component of shown genes (returned as invisible)
##' @export
pagoda.show.pathways <- function(pathways, varinfo, goenv = NULL, n.genes = 20, two.sided = FALSE, n.pc = rep(1, length(pathways)), colcols = NULL, zlim = NULL, showRowLabels = FALSE, cexCol = 1, cexRow = 1, nstarts = 10, cell.clustering = NULL, show.cell.dendrogram = TRUE, plot = TRUE, box = TRUE, trim = 0, return.details = FALSE , ...) {
    labRow <- NA
    if(showRowLabels) { labRow <- NULL }
    x <- c.view.pathways(pathways, varinfo$mat, varinfo$matw, goenv, batch = varinfo$batch, n.genes = n.genes, two.sided = two.sided, n.pc = n.pc, colcols = colcols, zlim = zlim, labRow = labRow, cexCol = cexCol, cexRow = cexRow, trim = trim, show.Colv = show.cell.dendrogram, plot = plot, vhc = cell.clustering, labCol = NA, box = TRUE, ...)
    if(return.details) {
        invisible(x)
    } else {
        invisible(x$scores[, 1])
    }
}

# takes in a list of pathways with a list of corresponding PC numbers
# recalculates PCs for each individual pathway, weighting gene loading in each pathway and then by total
# pathway variance over the number of genes (rough approximation)
c.view.pathways <- function(pathways, mat, matw, goenv = NULL, batch = NULL, n.genes = 20, two.sided = TRUE, n.pc = rep(1, length(pathways)), colcols = NULL, zlim = NULL, labRow = NA, vhc = NULL, cexCol = 1, cexRow = 1, nstarts = 50, row.order = NULL, show.Colv = TRUE, plot = TRUE, trim = 1.1/ncol(mat), showPC = TRUE,  ...) {
    # are these genes or pathways being passed?
    if(!is.null(goenv)) {
        x <- pathways %in% ls(goenv)
    } else {
        x <- rep(FALSE, length(pathways))
    }
    if(sum(x) > 0) { # some pathways matched
      if(!all(x)) {
        message("WARNING: partial match to pathway names. The following entries did not match: ", paste(pathways[!x], collapse = " "))
        }
        # look up genes for each pathway
        pathways <- pathways[x]
        p.genes <- mget(pathways, goenv, ifnotfound = NA)
    } else { # try as genes
        x <- pathways %in% rownames(mat)
        if(sum(x) > 0) {
            if(!all(x)) {
                message("WARNING: partial match to gene names. The following entries did not match: ", paste(pathways[!x], collapse = " "))
            }
            p.genes <- list("genes" = pathways[x])
            pathways <- c("genes");
        } else { # neither genes nor pathways are passed
            stop("ERROR: provided names do not match either gene nor pathway names (if the pathway environment was provided)")
        }
    }
    gvi <- rownames(mat) %in% unlist(p.genes)
    if(trim > 0) {
        mat <- winsorize.matrix(mat, trim = trim)
    }
    # recalculate wPCA for each pathway
    ppca <- pagoda.pathway.wPCA(varinfo = list(mat = mat[gvi, , drop = FALSE], matw = matw[gvi, , drop = FALSE], batch = batch), setenv = list2env(p.genes), n.cores = 1, n.randomizations = 0, n.starts = 2, n.components = max(n.pc), verbose = FALSE, min.pathway.size = 0, max.pathway.size = Inf, n.internal.shuffles = 0)

    if(length(ppca) > 1) { # if more than one pathway was supplied, combine genes using appropriate loadings and use consensus PCA (1st PC) as a pattern
        # score top loading genes for each desired PC, scaling by the sd/sqrt(n) (so that ^2 is = var/n)
        scaled.gene.loadings <- unlist(lapply(seq_along(pathways), function(i) {
            gl <- ppca[[pathways[i]]]$xp$rotation[, n.pc[i], drop = TRUE]*as.numeric(ppca[[pathways[i]]]$xp$sd)[n.pc[i]]/sqrt(ppca[[pathways[i]]]$n)
            names(gl) <- rownames(ppca[[pathways[i]]]$xp$rotation)
            gl
        }))


        if(two.sided) {
            # positive
            reduced.gene.loadings <- sort(tapply(scaled.gene.loadings, as.factor(names(scaled.gene.loadings)), max), decreasing = TRUE)
            selected.genes.pos <- reduced.gene.loadings[1:min(length(reduced.gene.loadings), round(n.genes/2))]

            # negative
            reduced.gene.loadings <- sort(tapply(scaled.gene.loadings, as.factor(names(scaled.gene.loadings)), min), decreasing = FALSE)
            selected.genes.neg <- reduced.gene.loadings[1:min(length(reduced.gene.loadings), round(n.genes/2))]
            selected.genes <- c(selected.genes.pos, selected.genes.neg)
            selected.genes <- selected.genes[match(unique(names(selected.genes)), names(selected.genes))]

        } else {
            reduced.gene.loadings <- sort(tapply(abs(scaled.gene.loadings), as.factor(names(scaled.gene.loadings)), max), decreasing = TRUE)
            selected.genes <- reduced.gene.loadings[1:min(length(reduced.gene.loadings), n.genes)]
        }

        # consensus pattern
        #lab <- match(names(selected.genes), rownames(mat))
        lab <- names(selected.genes);

        if(length(lab) == 0)
            return(NULL)
        if(length(lab)<3) { return(NULL) }
        if(trim > 0) {
            rn <- rownames(mat)
            cn <- colnames(mat)
            mat <- winsorize.matrix(mat, trim = trim)
            rownames(mat) <- rn
            colnames(mat) <- cn
        }
        d <- mat[lab, , drop = FALSE]
        dw <- matw[lab, , drop = FALSE]

        #d <- d*abs(as.numeric(selected.genes))
        xp <- bwpca(t(d), t(dw), npcs = 1, center = FALSE, nstarts = 3)

        consensus.npc = 1 # use first PC as a pattern
    } else { # only one pathway was provided
        xp <- ppca[[1]]$xp
        lab <- rownames(xp$rotation)[order(abs(xp$rotation[, n.pc[1]]), decreasing = TRUE)]
        if(length(lab) > n.genes) {
            lab <- lab[1:n.genes]
            xp$rotation <- xp$rotation[lab, , drop = FALSE]
        }

        d <- mat[lab, , drop = FALSE]
        dw <- matw[lab, , drop = FALSE]
        consensus.npc = n.pc[1] # use specified PC as a pattern
    }

    d <- d-rowMeans(d)
    dd <- as.dist(1-abs(cor(t(as.matrix(d)))))
    dd[is.na(dd)] <- 1
    if(is.null(row.order)) {
        if(length(lab) > 2) {
            if(is.element("fastcluster", installed.packages()[, 1])) {
                hc <- fastcluster::hclust(dd, method = "ward.D")
            } else {
                hc <- stats::hclust(dd, method = "ward.D")
            }
            row.order <- hc$order
        } else {
            row.order <- c(seq_along(lab))
        }
    }

    if(is.null(vhc)) {
        vd <- as.dist(1-cor(as.matrix(d)))
        vd[is.na(vd)] <- 1
        if(is.element("fastcluster", installed.packages()[, 1])) {
            vhc <- fastcluster::hclust(vd, method = "ward.D")
        } else {
            vhc <- stats::hclust(vd, method = "ward.D")
        }
    }

    if(is.null(zlim)) { zlim <- quantile(d, p = c(0.01, 0.99)) }
    vmap <- d
    vmap[vmap<zlim[1]] <- zlim[1]
    vmap[vmap > zlim[2]] <- zlim[2]
    rownames(vmap) <- rownames(d)

    aval <- colSums(d*dw*as.numeric(abs(xp$rotation[, consensus.npc])))/colSums(dw)
    oc <- xp$scores[, consensus.npc]
    if(cor(oc, aval, method = "p")<0) {
        oc <- -1*oc
        xp$scores[, consensus.npc] <- -1*xp$scores[, consensus.npc]
        xp$rotation[, consensus.npc] <- -1*xp$rotation[, consensus.npc]
    }

    aval <- colorRampPalette(c("blue", "white", "red"), space = "Lab")(100)[round(aval/max(abs(aval))*49)+50]
    z <- rbind(colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(100)[round(oc/max(abs(oc))*49)+50])

    ld <- xp$rotation[lab[row.order], consensus.npc]
    ld <- colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(100)[round(ld/max(abs(ld))*49)+50]

    # oc2 <- xp$scores[, 2]
    # oc2 <- colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(100)[round(oc2/max(abs(oc2))*49)+50]
    # oc3 <- xp$scores[, 3]
    # oc3 <- colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(100)[round(oc3/max(abs(oc3))*49)+50]

    #z <- do.call(rbind, list(aval, oc))
    #z <- rbind(oc3, oc2, oc)
    if((!showPC) || length(lab)<= 1) {
        z <- NULL
    }
    col <- colorRampPalette(c("blue", "white", "red"), space = "Lab")(256)

    if(!is.null(colcols)) {
      if(is.null(z)) {
        z <- colcols;
      } else {
        z <- rbind(colcols, z)
      }
    }

    if(plot) {
        if(show.Colv) {
            my.heatmap2(vmap[row.order, , drop = FALSE], Rowv = NA, Colv = as.dendrogram(vhc), zlim = zlim, col = col, scale = "none", RowSideColors = ld, ColSideColors = z, labRow = labRow, cexCol = cexCol, cexRow = cexRow, ...)
        } else {
            my.heatmap2(vmap[row.order, vhc$order, drop = FALSE], Rowv = NA, Colv = NA, zlim = zlim, col = col, scale = "none", RowSideColors = ld, ColSideColors = z[,vhc$order], labRow = labRow, cexCol = cexCol, cexRow = cexRow, ...)
        }

    }
    xp$vhc <- vhc
    xp$lab <- lab
    xp$row.order <- row.order
    xp$oc <- oc
    xp$col <- col
    xp$oc.col <- colorRampPalette(c("darkgreen", "white", "darkorange"), space = "Lab")(256)
    xp$vmap <- vmap
    xp$zlim <- zlim
    xp$consensus.pc <- consensus.npc
    return(invisible(xp))
}

# returns enriched categories for a given gene list as compared with a given universe
# returns a list with over and under fields containing list of over and underrepresented terms
calculate.go.enrichment <- function(genelist, universe, pvalue.cutoff = 1e-3, mingenes = 3, env = go.env, subset = NULL, list.genes = FALSE, over.only = FALSE) {
    genelist <- unique(genelist)
    all.genes <- unique(ls(env))
    # determine sizes
    universe <- unique(c(universe, genelist))
    universe <- universe[universe != ""]
    genelist <- genelist[genelist != ""]
    ns <- length(intersect(genelist, all.genes))
    us <- length(intersect(universe, all.genes))
    #pv <- lapply(go.map, function(gl) { nwb <- length(intersect(universe, gl[[1]])) if(nwb<mingenes) { return(0.5)} else { p <- phyper(length(intersect(genelist, gl[[1]])), nwb, us-nwb, ns) return(ifelse(p > 0.5, 1.0-p, p)) }})

    # compile count vectors
    stab <- table(unlist(mget(as.character(genelist), env, ifnotfound = NA), recursive = TRUE))
    utab <- table(unlist(mget(as.character(universe), env, ifnotfound = NA), recursive = TRUE))
    if(!is.null(subset)) {
        stab <- stab[names(stab) %in% subset]
        utab <- utab[names(utab) %in% subset]
    }

    tabmap <- match(rownames(stab), rownames(utab))

    cv <- data.frame(cbind(utab, rep(0, length(utab))))
    names(cv) <- c("u", "s")
    cv$s[match(rownames(stab), rownames(utab))] <- as.vector(stab)
    cv <- na.omit(cv)
    cv <- cv[cv$u > mingenes, ]

    if(over.only) {
        lpr <- phyper(cv$s-1, cv$u, us-cv$u, ns, lower.tail = FALSE, log.p = TRUE)
    } else {
        pv <- phyper(cv$s, cv$u, us-cv$u, ns, lower.tail = FALSE)
        lpr <- ifelse(pv<0.5, phyper(cv$s-1, cv$u, us-cv$u, ns, lower.tail = FALSE, log.p = TRUE), phyper(cv$s+1, cv$u, us-cv$u, ns, lower.tail = TRUE, log.p = TRUE))
    }
    lpr <- phyper(cv$s-1, cv$u, us-cv$u, ns, lower.tail = FALSE, log.p = TRUE)
    lpra <- bh.adjust(lpr, log = TRUE)
    z <- qnorm(lpr, lower.tail = FALSE, log.p = TRUE)
    za <- qnorm(lpra, lower.tail = FALSE, log.p = TRUE)
    # correct for multiple hypothesis
    mg <- length(which(cv$u > mingenes))
    if(over.only) {
        if(pvalue.cutoff<1) {
            ovi <- which(lpra<= log(pvalue.cutoff))
            uvi <- c()
        } else {
            ovi <- which((lpr+mg)<= log(pvalue.cutoff))
            uvi <- c()
        }
    } else {
        if(pvalue.cutoff<1) {
            ovi <- which(pv<0.5 & lpra<= log(pvalue.cutoff))
            uvi <- which(pv > 0.5 & lpra<= log(pvalue.cutoff))
        } else {
            ovi <- which(pv<0.5 & (lpr+mg)<= log(pvalue.cutoff))
            uvi <- which(pv > 0.5 & (lpr+mg)<= log(pvalue.cutoff))
        }
    }
    ovi <- ovi[order(lpr[ovi])]
    uvi <- uvi[order(lpr[uvi])]

    #return(list(over = data.frame(t = rownames(cv)[ovi], o = cv$s[ovi], u = cv$u[ovi], p = pr[ovi]*mg), under = data.frame(t = rownames(cv)[uvi], o = cv$s[uvi], u = cv$u[uvi], p = pr[uvi]*mg)))
    if(list.genes) {
        x <- mget(as.character(genelist), env, ifnotfound = NA)
        df <- data.frame(id = rep(names(x), unlist(lapply(x, function(d) length(na.omit(d))))), go = na.omit(unlist(x)), stringsAsFactors = FALSE)
        ggl <- tapply(df$id, as.factor(df$go), I)
        ovg <- as.character(unlist(lapply(ggl[rownames(cv)[ovi]], paste, collapse = " ")))
        uvg <- as.character(unlist(lapply(ggl[rownames(cv)[uvi]], paste, collapse = " ")))
        return(list(over = data.frame(t = rownames(cv)[ovi], o = cv$s[ovi], u = cv$u[ovi], Za = za, fe = cv$s[ovi]/(ns*cv$u[ovi]/us), genes = ovg), under = data.frame(t = rownames(cv)[uvi], o = cv$s[uvi], u = cv$u[uvi], Za = za, fe = cv$s[uvi]/(ns*cv$u[uvi]/us), genes = uvg)))
    } else {
        return(list(over = data.frame(t = rownames(cv)[ovi], o = cv$s[ovi], u = cv$u[ovi], p.raw = exp(lpr[ovi]), fdr = exp(lpra)[ovi], Z = z[ovi], Za = za[ovi], fe = cv$s[ovi]/(ns*cv$u[ovi]/us), fer = cv$s[ovi]/(length(genelist)*cv$u[ovi]/length(universe))), under = data.frame(t = rownames(cv)[uvi], o = cv$s[uvi], u = cv$u[uvi], p.raw = exp(lpr[uvi]), fdr = exp(lpra)[uvi], Z = z[uvi], Za = za[uvi], fe = cv$s[uvi]/(ns*cv$u[uvi]/us))))
    }
}

##' wrapper around different mclapply mechanisms
##'
##' Abstracts out mclapply implementation, and defaults to lapply when only one core is requested (helps with debugging)
##' @param ... parameters to pass to lapply, mclapply, bplapply, etc.
##' @param n.cores number of cores. If 1 core is requested, will default to lapply
papply <- function(...,n.cores=n) {
  if(n.cores>1) {
    # bplapply implementation
    bplapply(... , BPPARAM = MulticoreParam(workers = n.cores))
  } else { # fall back on lapply
    lapply(...);
  }
}


##' A Reference Class to represent the PAGODA application
##'
##' This ROOK application class enables communication with the client-side ExtJS framework and Inchlib HTML5 canvas libraries to create the graphical user interface for PAGODA
##' Refer to the code in \code{\link{make.pagoda.app}} for usage example
##'
##' @field results Output of the pathway clustering and redundancy reduction
##' @field genes List of genes to display in the Detailed clustering panel
##' @field pathways
##' @field mat Matrix of posterior mode count estimates
##' @field matw Matrix of weights associated with each estimate in \code{mat}
##' @field goenv Gene set list as an environment
##' @field renv Global environment
##' @field name Name of the application page; for display as the page title
##' @field trim Trim quantity used for Winsorization for visualization
##' @field batch Any batch or other known confounders to be included in the visualization as a column color track
##'
ViewPagodaApp <- setRefClass(
    'ViewPagodaApp',
    fields = c('results', 'genes', 'pathways', 'mat', 'matw', 'goenv', 'renv', 'name', 'trim', 'batch'),
    methods = list(

        initialize = function(results, pathways, genes, mat, matw, goenv, batch = NULL, name = "pathway overdispersion", trim = 1.1/ncol(mat)) {
            results <<- results
            genes <<- genes
            genes$svar <<- genes$var/max(genes$var)
            genes <<- genes
            mat <<- mat
            matw <<- matw
            batch <<- batch
            goenv <<- goenv
            pathways <<- pathways
            name <<- name
            trim <<- trim
            # reverse lookup environment
            renvt <- new.env(parent = globalenv())
            xn <- ls(envir = goenv)
            xl <- mget(xn, envir = goenv)
            gel <- tapply(rep(xn, unlist(lapply(xl, length))), unlist(xl), I)
            gel <- gel[nchar(names(gel)) > 0]
            x <- lapply(names(gel), function(n) assign(n, gel[[n]], envir = renvt))
            renv <<- renvt
            rm(xn, xl, x, gel, renvt)
            gc()
            callSuper()
        },
        getgenecldata = function(genes = NULL, gcl = NULL, ltrim = 0) { # helper function to get the heatmap data for a given set of genes
            if(is.null(gcl)) {
                gcl <- t.view.pathways(genes, mat = mat, matw = matw, env = goenv, vhc = results$hvc, plot = FALSE, trim = ltrim)
            }

            matrix <- gcl$vmap[rev(gcl$row.order), results$hvc$order, drop = FALSE]
            matrix <- list(data = as.numeric(t(matrix)),
                           dim = dim(matrix),
                           rows = rownames(matrix),
                           cols = colnames(matrix),
                           colors = gcl$col,
                           domain = seq.int(gcl$zlim[1], gcl$zlim[2], length.out = length(gcl$col))
            )

            ol <- list(matrix = matrix)
            if(nrow(gcl$vmap) > 2) {
                rcmvar <- matrix(gcl$rotation[rev(gcl$row.order), , drop = FALSE], ncol = 1)
                rowcols <- list(data = as.numeric(t(rcmvar)),
                                dim = dim(rcmvar),
                                colors = gcl$oc.col,
                                domain = seq.int(-1*max(abs(rcmvar)), max(abs(rcmvar)), length.out = length(gcl$oc.col))
                )

                colcols <- matrix(gcl$oc[results$hvc$order], nrow = 1)
                colcols <- list(data = as.numeric(t(colcols)),
                                dim = dim(colcols),
                                colors = gcl$oc.col,
                                domain = seq.int(-1*max(abs(colcols)), max(abs(colcols)), length.out = length(gcl$oc.col))
                )
                ol <- c(ol, list(rowcols = rowcols, colcols = colcols))
            }
            ol
        },
        call = function(env){
            path <- env[['PATH_INFO']]
            req <- Request$new(env)
            res <- Response$new()
            switch(path,
                   # INDEX
                   '/index.html' = {
                       body <- paste('<!DOCTYPE html >
                                     <meta charset = "utf-8" >
                                     <html >
                                     <head >
                                     <title > ', name, '</title >
                                     <meta http-equiv = "Content-Type" content = "text/html charset = iso-8859-1" >
                                     <link rel = "stylesheet" type = "text/css" href = "http://pklab.med.harvard.edu/sde/extjs/resources/ext-theme-neptune/ext-theme-neptune-all.css" / >
                                     <link rel = "stylesheet" type = "text/css" href = "http://pklab.med.harvard.edu/sde/extjs/examples/shared/example.css" / >
                                     <link rel = "stylesheet" type = "text/css" href = "http://pklab.med.harvard.edu/sde/pathcl.css" / >
                                     <head profile = "http://www.w3.org/2005/10/profile" >
                                     <link rel = "icon" type = "image/png" href = "http://pklab.med.harvard.edu/sde/pagoda.png" >
                                     <script type = "text/javascript" src = "http://pklab.med.harvard.edu/sde/extjs/ext-all.js" > </script >
                                     <script type = "text/javascript" src = "http://pklab.med.harvard.edu/sde/jquery-1.11.1.min.js" > </script >
                                     <script src = "http://d3js.org/d3.v3.min.js" charset = "utf-8" > </script >
                                     <script type = "text/javascript" src = "http://pklab.med.harvard.edu/sde/pathcl.js" > </script >
                                     </head >
                                     <body > </body >
                                     </html >
                                     ', sep = "")
                       res$header('"Content-Type": "text/html"')
                       res$write(body)
                   },
                   '/pathcl.json' = { # report pathway clustering heatmap data
                       # column dendrogram
                       t <- paste(tempfile(), "svg", sep = ".")
                       svg(file = t, width = 1, height = 1) # will be rescaled later
                       par(mar = rep(0, 4), mgp = c(2, 0.65, 0), cex = 1, oma = rep(0, 4))
                       #plot(results$hvc, main = "", sub = "", xlab = "", ylab = "", axes = FALSE, labels = FALSE, xaxs = "i", yaxs = "i", hang = 0.02)
                       plot(as.dendrogram(results$hvc), axes = FALSE, yaxs = "i", xaxs = "i", xlab = "", ylab = "", sub = "", main = "", leaflab = "none")
                       dev.off()
                       x <- readLines(t)
                       treeg <- paste(x[-c(1, 2, length(x))], collapse = "")

                       matrix <- results$rcm[rev(results$tvc$order), results$hvc$order]
                       matrix <- list(data = as.numeric(t(matrix)),
                                      dim = dim(matrix),
                                      rows = rownames(matrix),
                                      cols = colnames(matrix),
                                      colors = results$cols,
                                      domain = seq.int(results$zlim2[1], results$zlim2[2], length.out = length(results$cols)),
                                      range = range(matrix)
                       )


                       icols <- colorRampPalette(c("white", "black"), space = "Lab")(256)
                       rcmvar <- matrix(apply(results$rcm[rev(results$tvc$order), , drop = FALSE], 1, var), ncol = 1)
                       rowcols <- list(data = as.numeric(t(rcmvar)),
                                       # TODO: add annotation
                                       dim = dim(rcmvar),
                                       colors = icols,
                                       domain = seq.int(0, max(rcmvar), length.out = length(icols))
                       )
                       colcols <- list(data = unlist(lapply(as.character(t(results$colcol[nrow(results$colcol):1, results$hvc$order, drop = FALSE])), col2hex)),
                                       dim = dim(results$colcol)
                       )
                       ol <- list(matrix = matrix, rowcols = rowcols, colcols = colcols, coldend = treeg, trim = trim)
                       s <- toJSON(ol)
                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }
                   },
                   '/genecl.json' = { # report heatmap data for a selected set of genes
                       selgenes <- fromJSON(req$POST()$genes)
                       ltrim <- ifelse(is.null(req$params()$trim), 1.1/ncol(mat), as.numeric(req$params()$trim))
                       ol <- getgenecldata(selgenes, ltrim = ltrim)
                       s <- toJSON(ol)
                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }
                   },
                   '/pathwaygenes.json' = { # report heatmap data for a selected set of genes
                       ngenes <- ifelse(is.null(req$params()$ngenes), 20, as.integer(req$params()$ngenes))
                       twosided <- ifelse(is.null(req$params()$twosided), FALSE, as.logical(req$params()$twosided))
                       ltrim <- ifelse(is.null(req$params()$trim), 1.1/ncol(mat), as.numeric(req$params()$trim))
                       pws <- fromJSON(req$POST()$genes)
                       n.pcs <- as.integer(gsub("^#PC(\\d+)# .*", "\\1", pws))
                       n.pcs[is.na(n.pcs)]<-1
                       x <- c.view.pathways(gsub("^#PC\\d+# ", "", pws), mat, matw, goenv = goenv, n.pc = n.pcs, n.genes = ngenes, two.sided = twosided, vhc = results$hvc, plot = FALSE, trim = ltrim, batch = batch)
                       #x <- t.view.pathways(gsub("^#PC\\d+# ", "", pws), mat, matw, env = goenv, vhc = results$hvc, plot = FALSE, trim = ltrim, n.pc = 1)
                       ##rsc <- as.vector(rowSums(matw[rownames(x$rotation), ]))*x$rotation[, 1]
                       #rsc <- x$rotation[, 1]
                       #if(twosided) {
                       #  extgenes <- unique(c(names(sort(rsc))[1:min(length(rsc), round(ngenes/2))], names(rev(sort(rsc)))[1:min(length(rsc), round(ngenes/2))]))
                       # } else {
                       #   extgenes <- names(sort(abs(rsc), decreasing = TRUE))[1:min(length(rsc), ngenes)]
                       #}
                       #ol <- getgenecldata(extgenes, ltrim = ltrim)
                       ol <- getgenecldata(genes = NULL, gcl = x, ltrim = ltrim)
                       s <- toJSON(ol)
                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }
                   },
                   '/patterngenes.json' = { # report heatmap of genes most closely matching a given pattern
                       ngenes <- ifelse(is.null(req$params()$ngenes), 20, as.integer(req$params()$ngenes))
                       twosided <- ifelse(is.null(req$params()$twosided), FALSE, as.logical(req$params()$twosided))
                       ltrim <- ifelse(is.null(req$params()$trim), 1.1/ncol(mat), as.numeric(req$params()$trim))
                       pat <- fromJSON(req$POST()$pattern)
                       # reorder the pattern back according to column clustering
                       pat[results$hvc$order] <- pat
                       patc <- .Call("matCorr", as.matrix(t(mat)), as.matrix(pat, ncol = 1) , PACKAGE = "scde")
                       if(twosided) { patc <- abs(patc) }
                       mgenes <- rownames(mat)[order(as.numeric(patc), decreasing = TRUE)[1:ngenes]]
                       ol <- getgenecldata(mgenes, ltrim = ltrim)
                       ol$pattern <- pat
                       s <- toJSON(ol)
                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }
                   },
                   '/clinfo.json' = {
                       pathcl <- ifelse(is.null(req$params()$pathcl), 1, as.integer(req$params()$pathcl))
                       ii <- which(results$ct == pathcl)
                       tpi <- order(results$matvar[ii], decreasing = TRUE)
                       #tpi <- tpi[seq(1, min(length(tpi), 15))]
                       npc <- gsub("^#PC(\\d+)#.*", "\\1", names(ii[tpi]))
                       nams <- gsub("^#PC\\d+# ", "", names(ii[tpi]))
                       if(exists("myGOTERM", envir = globalenv())) {
                           tpn <- paste(nams, mget(nams, get("myGOTERM", envir = globalenv()), ifnotfound = ""), sep = " ")
                       } else {
                           tpn <- names(ii[tpi])
                       }

                       lgt <- data.frame(do.call(rbind, lapply(seq_along(tpn), function(i) c(id = names(ii[tpi[i]]), name = tpn[i], npc = npc[i], od = as.numeric(results$matvar[ii[tpi[i]]])/max(results$matvar), sign = as.numeric(results$matrcmcor[ii[tpi[i]]]), initsel = as.integer(results$matvar[ii[tpi[i]]] >= results$matvar[ii[tpi[1]]]*0.8)))))

                       # process additional filters
                       if(!is.null(req$params()$filter)) {
                           fl <- fromJSON(URLdecode(req$params()$filter))
                           for( fil in fl) {
                               lgt <- lgt[grep(fil$value, lgt[, fil$property], perl = TRUE, ignore.case = TRUE), ]
                           }
                       }
                       start <- ifelse(is.null(req$params()$start), 1, as.integer(req$params()$start)+1)
                       limit <- ifelse(is.null(req$params()$limit), 100, as.integer(req$params()$limit))
                       dir <- ifelse(is.null(req$params()$dir), "DESC", req$params()$dir)
                       trows <- nrow(lgt)
                       if(trows > 0) {
                           if(!is.null(req$params()$sort)) {
                               if(req$params()$sort %in% colnames(lgt)) {
                                   lgt <- lgt[order(lgt[, req$params()$sort], decreasing = (dir == "DESC")), ]
                               }
                           }
                       }
                       lgt <- lgt[min(start, nrow(lgt)):min((start+limit), nrow(lgt)), ]
                       lgt$od <- format(lgt$od, nsmall = 2, digits = 2)
                       ol <- apply(lgt, 1, function(x) as.list(x))
                       names(ol) <- NULL
                       s <- toJSON(list(totalCount = trows, genes = ol))

                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }
                   },
                   '/genes.json' = {
                       lgt <- genes
                       if(!is.null(req$params()$filter)) {
                           fl <- fromJSON(URLdecode(req$params()$filter))
                           for( fil in fl) {
                               lgt <- lgt[grep(fil$value, lgt[, fil$property], perl = TRUE, ignore.case = TRUE), ]
                           }
                       }
                       start <- ifelse(is.null(req$params()$start), 1, as.integer(req$params()$start)+1)
                       limit <- ifelse(is.null(req$params()$limit), 1000, as.integer(req$params()$limit))
                       dir <- ifelse(is.null(req$params()$dir), "DESC", req$params()$dir)
                       trows <- nrow(lgt)
                       if(trows > 0) {
                           if(!is.null(req$params()$sort)) {
                               if(req$params()$sort %in% colnames(lgt)) {
                                   lgt <- lgt[order(lgt[, req$params()$sort], decreasing = (dir == "DESC")), ]
                               }
                           } else { # default sort
                               # already done
                           }
                       }
                       lgt <- format(lgt[min(start, nrow(lgt)):min((start+limit), nrow(lgt)), ], nsmall = 2, digits = 2)
                       ol <- apply(lgt, 1, function(x) as.list(x))
                       names(ol) <- NULL
                       s <- toJSON(list(totalCount = trows, genes = ol))
                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }

                   },
                   '/pathways.json' = {
                       lgt <- pathways
                       if(!is.null(req$params()$filter)) {
                           fl <- fromJSON(URLdecode(req$params()$filter))
                           for( fil in fl) {
                               lgt <- lgt[grep(fil$value, lgt[, fil$property], perl = TRUE, ignore.case = TRUE), ]
                           }
                       }
                       start <- ifelse(is.null(req$params()$start), 1, as.integer(req$params()$start)+1)
                       limit <- ifelse(is.null(req$params()$limit), 1000, as.integer(req$params()$limit))
                       dir <- ifelse(is.null(req$params()$dir), "DESC", req$params()$dir)
                       trows <- nrow(lgt)
                       if(trows > 0) {
                           if(!is.null(req$params()$sort)) {
                               if(req$params()$sort %in% colnames(lgt)) {
                                   lgt <- lgt[order(lgt[, req$params()$sort], decreasing = (dir == "DESC")), ]
                               }
                           } else { # default sort
                               # already done
                           }
                       }
                       lgt <- format(lgt[min(start, nrow(lgt)):min((start+limit), nrow(lgt)), ], nsmall = 2, digits = 2)
                       ol <- apply(lgt, 1, function(x) as.list(x))
                       names(ol) <- NULL
                       s <- toJSON(list(totalCount = trows, genes = ol))
                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }

                   },
                   '/testenr.json' = { # run an enrichment test
                       selgenes <- fromJSON(req$POST()$genes)
                       lgt <- calculate.go.enrichment(selgenes, rownames(mat), pvalue.cutoff = 0.99, env = renv, over.only = TRUE)$over
                       if(exists("myGOTERM", envir = globalenv())) {
                           lgt$nam <- paste(lgt$t, mget(as.character(lgt$t), get("myGOTERM", envir = globalenv()), ifnotfound = ""), sep = " ")
                       } else {
                           lgt$name <- lgt$t
                       }
                       lgt <- data.frame(id = paste("#PC1#", lgt$t), name = lgt$nam, o = lgt$o, u = lgt$u, Z = lgt$Z, Za = lgt$Za, fe = lgt$fe, stringsAsFactors = FALSE)

                       if(!is.null(req$params()$filter)) {
                           fl <- fromJSON(URLdecode(req$params()$filter))
                           for( fil in fl) {
                               lgt <- lgt[grep(fil$value, lgt[, fil$property], perl = TRUE, ignore.case = TRUE), ]
                           }
                       }
                       start <- ifelse(is.null(req$params()$start), 1, as.integer(req$params()$start)+1)
                       limit <- ifelse(is.null(req$params()$limit), 1000, as.integer(req$params()$limit))
                       dir <- ifelse(is.null(req$params()$dir), "DESC", req$params()$dir)
                       trows <- nrow(lgt)
                       if(trows > 0) {
                           if(!is.null(req$params()$sort)) {
                               if(req$params()$sort %in% colnames(lgt)) {
                                   lgt <- lgt[order(lgt[, req$params()$sort], decreasing = (dir == "DESC")), ]
                               }
                           }
                       }
                       lgt <- format(lgt[min(start, nrow(lgt)):min((start+limit), nrow(lgt)), ], nsmall = 2, digits = 2)
                       ol <- apply(lgt, 1, function(x) as.list(x))
                       names(ol) <- NULL
                       s <- toJSON(list(totalCount = trows, genes = ol))
                       res$header('Content-Type', 'application/javascript')
                       if(!is.null(req$params()$callback)) {
                           res$write(paste(req$params()$callback, "(", s, ")", sep = ""))
                       } else {
                           res$write(s)
                       }

                   },
                   '/celltable.txt' = {
                       matrix <- results$rcm[rev(results$tvc$order), results$hvc$order]
                       body <- paste(capture.output(write.table(round(matrix, 1), sep = "\t")), collapse = "\n")
                       res$header('Content-Type', 'text/plain')
                       #res$header('"Content-disposition": attachment')
                       res$write(body)
                   },
{
    res$header('Location', 'index.html')
    res$write('Redirecting to <a href = "index.html" > index.html</a >  for interactive browsing.')
}
                       )
            res$finish()
        }
            )
    )
