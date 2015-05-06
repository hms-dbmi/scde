#ifndef _scde_BWPCA_H
#define _scde_BWPCA_H

#include <RcppArmadillo.h>

// outer function for performing Bailey's weighted PCA
// Nshuffles - number of internal row-specific randomizations to recalculate the lambda1 (PC1 variance) on
RcppExport SEXP baileyWPCA(SEXP Mat, SEXP Matw, SEXP Npcs, SEXP Nstarts, SEXP Smooth, SEXP EMtol, SEXP EMmaxiter, SEXP Seed, SEXP Nshuffles);

void baileyWPCAround(arma::mat& m,arma::mat& mw,int nstarts,int npcs,int seed,int maxiter,double tol,int smooth,arma::colvec& smoothc,arma::mat& bestcoef, arma::mat& besteigenv);
#endif
