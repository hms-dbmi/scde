#include "bwpca.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <numeric>


using namespace Rcpp ;


#undef DEBUG


// compare a and b positions of a data vector
bool compare_on_other(int a, int b, arma::vec& data) {
    return data[a]<data[b];
}

// independently randomize values within a column for a matrix
void set_random_matrix(arma::mat& target,arma::mat& source) {
    std::vector<int> ind(target.n_rows);
    for(int j=0;j<target.n_rows;j++) {
        ind[j]=j;
        R_CheckUserInterrupt();
    } // set up initial index (1,2,3)

    for(int i=0;i<target.n_cols;i++) {
        std::random_shuffle(ind.begin(), ind.end());
        //std::sort(ind.start(), ind.end(), std::bind(compare_on_other,  _1, _2, rv));
        for(int j=0;j<target.n_rows;j++) {
            target(j,i)=source(ind[j],i);
            //R_CheckUserInterrupt();
        }
        R_CheckUserInterrupt();
    }
}

// independently randomize values within a column for two matrices
void set_random_matrices(arma::mat& target1,arma::mat& source1,arma::mat& target2,arma::mat& source2) {
    std::vector<int> ind(target1.n_rows);
    for(int j=0;j<target1.n_rows;j++) {
        ind[j]=j;
        R_CheckUserInterrupt();
    } // set up initial index (1,2,3)

    for(int i=0;i<target1.n_cols;i++) {
        std::random_shuffle(ind.begin(), ind.end());
        //std::sort(ind.start(), ind.end(), std::bind(compare_on_other,  _1, _2, rv));
        for(int j=0;j<target1.n_rows;j++) {
            target1(j,i)=source1(ind[j],i);
            target2(j,i)=source2(ind[j],i);
            //R_CheckUserInterrupt();
        }
        R_CheckUserInterrupt();
    }
}

SEXP baileyWPCA(SEXP Mat, SEXP Matw, SEXP Npcs, SEXP Nstarts, SEXP Smooth, SEXP EMtol, SEXP EMmaxiter, SEXP Seed, SEXP Nshuffles){
    arma::mat m=Rcpp::as<arma::mat>(Mat); // can avoid copy here
    arma::mat mw=Rcpp::as<arma::mat>(Matw); // can avoid copy here


    double tol=Rcpp::as<double>(EMtol);
    int maxiter=Rcpp::as<int>(EMmaxiter);
    int nstarts=Rcpp::as<int>(Nstarts);
    int smooth=Rcpp::as<int>(Smooth);
    int npcs=Rcpp::as<int>(Npcs);
    int seed=Rcpp::as<int>(Seed);
    int nshuffles=Rcpp::as<int>(Nshuffles);

    int d=m.n_cols; // genes
    int n=m.n_rows; // cells

    if(npcs>d) { npcs=d; } // limit the number of PCs to the number of genes


    //arma::mat mwsq_colsum=sum(mw % mw);

#ifdef DEBUG
    std::cout<<"starting up"<<std::endl<<std::flush;
#endif
    // set up smoothing coefficients
    arma::colvec smoothc;
    if(smooth>0) {
        int np=smooth/2;
        int pol_degree=3; int diff_order=0;
        arma::colvec x(2*np+1); for(int i=0;i<2*np+1;i++) { x[i]=i-np; };
        arma::mat A(2*np+1,pol_degree+1);
        for(int j=0; j<pol_degree+1; j++) {
            A.col(j)=pow(x,j);
        }
        arma::mat ATA(A.t() * A);
        arma::colvec rhs(pol_degree+1,arma::fill::zeros);
        rhs[diff_order]=pow(-1,diff_order);
        smoothc=A * solve(ATA,rhs);
        //std::cout<<"smoothc:"; smoothc.print();
    }
#ifdef DEBUG
    //std::cout<<"calculated smoothing factors"<<std::endl<<std::flush;
#endif
    // multiple random starts
    arma::mat besteigenv,bestcoef;
    baileyWPCAround(m,mw,nstarts,npcs,seed,maxiter,tol,smooth,smoothc,bestcoef,besteigenv);

    // calculate amount of weighted variance explained by each estimated component
    // calculate total weighted variance
    arma::mat totvm(m % sqrt(mw));
    totvm %=totvm;
    //double totvar=accu(sum(totvm,0) / mwsq_colsum);
    double totvar=accu(totvm);
#ifdef DEBUG
    std::cout<<"total weighted variance : "<<totvar<<std::endl;
#endif

    double tvarexp=0;
    arma::vec varexp(npcs);
    arma::mat dat(m.n_rows,m.n_cols,arma::fill::zeros);
    for(int k=0; k<npcs; k++) {
        dat += (bestcoef.col(k) * besteigenv.col(k).t());
        //std::cout<<k<<"-th total reconstruction:"<<std::endl; dat.print();
#ifdef DEBUG
        //std::cout<<k<<"-th total reconstruction delta:"<<std::endl; (dat-m).print();
#endif
        arma::mat delta((dat-m) % sqrt(mw));
        delta %=delta;
        //double npres=accu(sum(delta,0) / mwsq_colsum);
        double npres=accu(delta);
#ifdef DEBUG
        std::cout<<k<<"-th component explains "<<(totvar-npres-tvarexp)<<" ("<<((totvar-npres-tvarexp)/totvar)<<") variance;  cumulative:"<<(totvar-npres)<<" ("<<((totvar-npres)/totvar)<<")"<<std::endl;
#endif
        varexp[k]=totvar-npres-tvarexp;
        tvarexp=totvar-npres;

        R_CheckUserInterrupt();
    }

    arma::mat pcw=mw * abs(besteigenv);

    if(nshuffles>0) {
        arma::vec rvars(nshuffles);
        arma::mat rm(n,d);
        arma::mat rmw(n,d);
        arma::mat reigenv,rcoef;
        for(int i=0;i<nshuffles;i++) {
            set_random_matrices(rm,m,rmw,mw);
            baileyWPCAround(rm,rmw,nstarts,npcs,seed+i,maxiter,tol,smooth,smoothc,rcoef,reigenv);
            dat.fill(0);
            dat += (rcoef.col(0) * reigenv.col(0).t());
            arma::mat delta((dat-rm) % sqrt(rmw));
            delta %=delta;
            rvars[i]=totvar-accu(delta);
        }

        return List::create(Named("rotation") = wrap(besteigenv),
                            Named("scores") = wrap(bestcoef),
                            Named("scoreweights") = wrap(pcw),
                            Named("var") = wrap(varexp),
                            Named("totvar") = wrap(totvar),
                            Named("randvar") = wrap(rvars));

    } else {
        return List::create(Named("rotation") = wrap(besteigenv),
                            Named("scores") = wrap(bestcoef),
                            Named("scoreweights") = wrap(pcw),
                            Named("var") = wrap(varexp),
                            Named("totvar") = wrap(totvar));
    }
}


// internal function performing wPCA itself with nstarts random starts
void baileyWPCAround(arma::mat& m,arma::mat& mw,int nstarts,int npcs,int seed,int maxiter,double tol,int smooth,arma::colvec& smoothc,arma::mat& bestcoef, arma::mat& besteigenv) {
    int d=m.n_cols; // genes
    int n=m.n_rows; // cells

    double bestpres(-1);

    for(int nstart=0; nstart<nstarts; nstart++) {
#ifdef DEBUG
        std::cout<<"starting iteration "<<nstart<<std::endl<<std::flush;
#endif
        // random orthonormal start for the
        arma::arma_rng::set_seed(seed+nstart);
        arma::mat X = arma::randu<arma::mat>(d,npcs);
        arma::mat start,eigenv, R;
        arma::qr_econ(start,R,X);
        eigenv=start;

        //std::cout<<"starting eigenv"<<std::endl<<std::flush;
        //eigenv.print();


        // initial solution for coefficients
        arma::mat coef(n,npcs);

        double pres(std::numeric_limits<double>::max());
        double bpres(std::numeric_limits<double>::max());
        arma::mat beigenv,bcoef; // best models of the current run

        int ii=0; // iteration counter
        while(ii<maxiter) {
            //std::cout<<"iteration"<<ii<<std::endl<<std::flush;

            // solve for coefficients
            for(int j=0; j<n; j++) {  // for each observation j
                // solving for c, m.row(j) == eigenv * c, as a weighted least squares
                // problem with weights mw.row(j)

                //std::cout<<"m:"; m.row(j).print();
                //std::cout<<"w:"; mw.row(j).print();
                arma::rowvec b=m.row(j) % mw.row(j);
                //std::cout<<"b*w:"; b.print();
                b= b * eigenv;
                //std::cout<<"b:"; b.print();
                arma::mat A(eigenv);
                A.each_col() %= mw.row(j).t();
                A=eigenv.t() * A;
                //std::cout<<"A:"; A.print();
                coef.row(j)=solve(A,b.t()).t();
            }
            //std::cout<<"updated coeff"<<std::endl<<std::flush;
            //coef.t().print();


            // solve for eigenvectors
            arma::mat dat=m;
            for(int k=0; k<npcs; k++) { // for each vector
                arma::mat cw=mw;
                cw.each_col() %=coef.col(k);
                arma::rowvec xcw=sum(dat % cw,0);
                cw.each_col() %=coef.col(k);
                xcw /= sum(cw,0);
                eigenv.col(k)=xcw.t();

                // smoothing
                if(smooth>0) {
                    // convolve with smoothing coefficients
                    int n=(smoothc.n_elem-1)/2;
                    arma::vec res=conv(eigenv.col(k),smoothc);
                    eigenv.col(k)=res.subvec(n,res.n_elem-n-1);
                    //std::cout<<"smooth res[]:"; eigenv.col(k).print();
                }

                // subtract current eigenvector vector from the data
                if(k!= npcs-1) {
                    dat -= coef.col(k) * eigenv.col(k).t()  ;
                    //std::cout<<"dat:"; dat.print();
                }
            }
            //std::cout<<"updated eigenvectors"<<std::endl<<std::flush;

            // renormalize and re-orthogonalize the eigenvectors
            eigenv.col(0) /=  sqrt(dot(eigenv.col(0),eigenv.col(0)));
            for(int k=1; k<npcs; k++) { // for each vector
                for(int kx=0; kx<k; kx++) {
                    double c = dot(eigenv.col(k) , eigenv.col(kx));
                    eigenv.col(k) -= c*eigenv.col(kx);
                }
                eigenv.col(k) /= sqrt(dot(eigenv.col(k) , eigenv.col(k)));
            }

            //std::cout<<"recalculated eigenvectors"<<std::endl<<std::flush;
            //eigenv.t().print();

            // recalculate the model fit
            arma::mat model(coef * eigenv.t());
            //std::cout<<"coeffs:"; coef.print();
            //std::cout<<"eigenv:"; eigenv.print();
            //std::cout<<"model:"; model.print();

            //std::cout<<"delta:"<<std::endl; (model-m).print();

            arma::mat delta((model-m) % sqrt(mw));
            delta %=delta;
            //double npres=accu(sum(delta,0) / mwsq_colsum);
            double npres=accu(delta);

#ifdef DEBUG
            std::cout<<"iteration "<<ii<<" pres="<<npres;
#endif
            // record model if it improved the overall precision
            if(npres<bpres) {
                bpres=npres;
                bcoef=coef;
                beigenv=eigenv;
#ifdef DEBUG
                std::cout<<", best so far";
#endif
            }

            if(tol>0 && ii>0 && (pres-npres)/npres < tol) {
                if(pres>npres) {

#ifdef DEBUG
                    std::cout<<", reached required tolerance"<<std::endl;
#endif
                    pres=npres;
                    break;
                }
                // otherwise the total variance is actually increasing, which is not a good thing
            }
#ifdef DEBUG
            std::cout<<std::endl<<std::flush;
#endif

            ii++;
            pres=npres;
        }
        if(nstart==0 || pres<bestpres) {
#ifdef DEBUG
            std::cout<<"updating the results with "<<nstart<<" calculations"<<std::endl<<std::flush;
#endif
            bestpres=bpres;
            bestcoef=bcoef;
            besteigenv=beigenv;
        }
    }

    R_CheckUserInterrupt();

}
