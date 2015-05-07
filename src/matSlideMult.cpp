#include "matSlideMult.h"

using namespace Rcpp ;

SEXP matSlideMult(SEXP Mat1, SEXP Mat2){
    arma::mat m1=Rcpp::as<arma::mat>(Mat1);
    arma::mat m2=Rcpp::as<arma::mat>(Mat2);
    int n=m1.n_cols;
    arma::mat rm(m1.n_rows,2*n-1);

    // left half
    for(int i=n;i>1;i--) {
        rm.col(n-i)=sum(m1.cols(0,n-i) % m2.cols(i-1,n-1),1);
        R_CheckUserInterrupt();
    }
    // right half
    for(int i=1;i<=n;i++) {
        rm.col(n-2+i)=sum(m1.cols(i-1,n-1) % m2.cols(0,n-i),1);
        R_CheckUserInterrupt();
    }

    return wrap(rm);
}

