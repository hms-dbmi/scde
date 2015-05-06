#ifndef _scde_PAGODA_H
#define _scde_PAGODA_H

#include <RcppArmadillo.h>
RcppExport SEXP winsorizeMatrix(SEXP Mat, SEXP Trim); 
RcppExport SEXP matWCorr(SEXP Mat, SEXP Matw);
RcppExport SEXP plSemicompleteCor2(SEXP Pl);
RcppExport SEXP matCorr(SEXP X, SEXP Y);
#endif
