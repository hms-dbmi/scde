#include "jpmatLogBoot.h"
#include <cmath>

using namespace Rcpp ;

// maximum and minimum values of negative binomial theta (size) that will be considered
#define MIN_THETA 1.0e-2
#define MAX_THETA 1.0e+3


SEXP jpmatLogBoot(SEXP Matl, SEXP Nboot, SEXP Seed){
    Rcpp::List matl(Matl);
    int nrows=Rcpp::as<Rcpp::NumericMatrix>(matl[0]).nrow();
    int ncols=Rcpp::as<Rcpp::NumericMatrix>(matl[0]).ncol();
    int nmat=matl.size();
    int nboot=as<int>(Nboot);
    arma::mat jp(nrows,ncols);  // joint posterior across boostraps
    jp.zeros();
    arma::mat tjp(nrows,ncols); // for current boostrap posterior
    int seed=Rcpp::as<int>(Seed);
    srand(seed);

    for(int i=0;i<nboot;i++) {
        tjp.zeros();
        for(int j=0;j<nmat;j++) {
            int rj;
            while(nmat <= (rj= rand()/(RAND_MAX/nmat)))
                ;
            arma::mat am(Rcpp::as<Rcpp::NumericMatrix>(matl[rj]).begin(),nrows,ncols,false,true);
            tjp +=am;
            //R_CheckUserInterrupt();
        }
        arma::colvec m=max(tjp,1);
        tjp.each_col() -= m; // shift for stability
        tjp=exp(tjp);
        arma::colvec s=sum(tjp,1);
        tjp.each_col() /= s;
        jp+=tjp;
        //R_CheckUserInterrupt();
    }
    return wrap(jp);
}


// similar to above, however the joint is built not over boostrap samples
// of the same pool, but by sampling a pre-defined (by Comp vector) composition
// of several different pools (in Matll)
SEXP jpmatLogBatchBoot(SEXP Matll, SEXP Comp, SEXP Nboot, SEXP Seed){
    Rcpp::IntegerVector comp=as<Rcpp::IntegerVector>(Comp);
    Rcpp::NumericMatrix nm0=as<Rcpp::NumericMatrix>(VECTOR_ELT( VECTOR_ELT(Matll,0) ,0));
    int nrows=nm0.nrow();
    int ncols=nm0.ncol();
    int nboot=as<int>(Nboot);
    arma::mat jp(nrows,ncols);  // joint posterior across boostraps
    jp.zeros();
    arma::mat tjp(nrows,ncols); // for current boostrap posterior
    int seed=Rcpp::as<int>(Seed);
    srand(seed);

    for(int i=0;i<nboot;i++) {
        tjp.zeros();
        for(int k=0;k<comp.size();k++) { // over types
            int nsamp=comp[k];
            if(nsamp>0) {
                int nmat=LENGTH( VECTOR_ELT(Matll, k) );
                for(int j=0;j<nsamp;j++) { // over matrices
                    int rj;
                    while(nmat <= (rj= rand()/(RAND_MAX/nmat)))
                        ;
                    arma::mat am(Rcpp::as<Rcpp::NumericMatrix>(VECTOR_ELT( VECTOR_ELT(Matll,k) ,rj)).begin(),nrows,ncols,false,true);
                    tjp +=am;
                    //R_CheckUserInterrupt();
                }
            }
            //R_CheckUserInterrupt();
        }
        arma::colvec m=max(tjp,1);
        tjp.each_col() -= m; // shift for stability
        tjp=exp(tjp);
        arma::colvec s=sum(tjp,1);
        tjp.each_col() /= s;
        jp+=tjp;
        //R_CheckUserInterrupt();
    }
    return wrap(jp);
}

//calculate log joint posteriors from models and unique count lists/integer on magniude grid
// count columns must match model rows
// Models - model matrix
// Ucl - unique count list
// Uci - unique count index
// Magnitudes - marginals - (expression magnitudes)
// Nboot - number of bootstrap iterations
// Seed - random seed
// ReturnIndividualPosterior: 0 - nothing, 1 - maxima, 2 - full posterior matrices
// LocalThetaFit: 0 - constant theta; 1 - theta(fpm) fit
// SquareLogitConc: 0 - concomitant is a function of magnitude ; 1 - of magnitude and magnitude^2
// EnsembleProbability: 1 - calculate ensemble (i.e. sum) joint probability instead of a joint (i.e product)
RcppExport SEXP logBootPosterior(SEXP Models, SEXP Ucl, SEXP CountsI, SEXP Magnitudes, SEXP Nboot, SEXP Seed, SEXP ReturnIndividualPosteriors, SEXP LocalThetaFit, SEXP SquareLogitConc, SEXP EnsembleProbability) {
#define CONCB_I 0
#define CONCA_I 1
#define FAILR_I 2
#define CORRB_I 3
#define CORRA_I 4
#define CORRT_I 5
#define CORRlTB_I 6
#define CORRlTT_I 7
#define CORRlTM_I 8
#define CORRlTS_I 9
#define CORRlTR_I 10
#define CONCA2_I 11

    Rcpp::IntegerMatrix counti=Rcpp::as<Rcpp::IntegerMatrix>(CountsI);
    Rcpp::List ucl(Ucl);
    int ncells=ucl.size();
    int returnpost=Rcpp::as<int>(ReturnIndividualPosteriors);
    int localtheta=Rcpp::as<int>(LocalThetaFit);
    int squarelogitconc=Rcpp::as<int>(SquareLogitConc);
    int ensemblep=Rcpp::as<int>(EnsembleProbability);
    arma::mat models=Rcpp::as<arma::mat>(Models);
    arma::colvec magnitudes=Rcpp::as<arma::colvec>(Magnitudes);
    std::vector< arma::mat > ucposteriors;
    std::vector< std::vector < arma::uword > > ucmaxi;
    // calculate individual posteriors for each cell, for each unique count value
    //std::cout<<"individual posteriors "<<std::flush;
    double minlogprob=-1*std::numeric_limits<double>::max()/ncells/1.1;
    for(int i=0;i<ncells;i++) {
        Rcpp::IntegerVector uc(Rcpp::as<Rcpp::IntegerVector>(ucl[i]));
        int ncounts=uc.size();
        arma::mat pm(magnitudes.n_elem,ncounts);
        std::vector< arma::uword > maxi;
        arma::vec mu=magnitudes * models(i,CORRA_I);
        mu+=models(i,CORRB_I); mu=exp(mu);
        arma::vec cfp;
        if(squarelogitconc) {
            cfp=models(i,CONCA_I) + magnitudes*models(i,CONCA2_I);
            cfp%=magnitudes;
        } else {
            cfp=magnitudes * models(i,CONCA_I);
        }
        cfp+=models(i,CONCB_I);
        cfp=1/(exp(cfp)+1);
        arma::vec cfpr=1-cfp;
        cfp=log(cfp); cfpr=log(cfpr);
        double maxcfp=max(cfp);
        arma::colvec thetas;
        if(localtheta) { // non-constant theta model - prepare theta values
            thetas=-1*magnitudes + models(i,CORRlTM_I);
            thetas*=models(i,CORRlTS_I);
            thetas=exp10(thetas)+1;
            thetas=pow(thetas,models(i,CORRlTR_I));
            thetas=(models(i,CORRlTT_I) - models(i,CORRlTB_I))/thetas;
            thetas+=models(i,CORRlTB_I);
            thetas=exp(-1*thetas);

            for(unsigned int k=0;k<thetas.n_elem;k++) {
                if((!std::isfinite(thetas[k])) || (thetas[k]<MIN_THETA)) {  thetas[k]=MIN_THETA;}
                if(thetas[k]>MAX_THETA) {  thetas[k]=MAX_THETA;}
                //R_CheckUserInterrupt();
            }
        }
        //std::cout<<"cfp=["; std::copy(cfp.begin(),cfp.end(),std::ostream_iterator<double>(std::cout," ")); std::cout<<"]"<<std::endl<<std::flush;
        //std::cout<<"thetas=["; copy(thetas.begin(),thetas.end(),std::ostream_iterator<double>(std::cout," ")); std::cout<<"]"<<std::endl<<std::flush;

        for(int j=0;j<ncounts;j++) {
            // correlated prob
            arma::vec nbp(mu.n_elem);
            if(localtheta) { // linear theta model
                for(unsigned int k=0;k<mu.n_elem;k++) {
                    double muv=mu[k]; double x=uc[j];
                    // choose maximum probability when hitting the grid with the maximum
                    if((k<(mu.n_elem-1) && x>muv && x<mu[k+1]) || (k==(mu.n_elem-1) && x>muv)) { muv=x; }
                    nbp[k]=Rf_dnbinom(x,thetas[k],thetas[k]/(thetas[k]+muv),true);
                    //R_CheckUserInterrupt();
                }
            } else { // constant theta
                double theta=models(i,CORRT_I);
                for(unsigned int k=0;k<mu.n_elem;k++) {
                    double muv=mu[k]; double x=uc[j];
                    // choose maximum probability when hitting the grid with the maximum
                    if((k<(mu.n_elem-1) && x>muv && x<mu[k+1]) || (k==(mu.n_elem-1) && x>muv)) { muv=x; }
                    nbp[k]=Rf_dnbinom(x,theta,theta/(theta+muv),true);
                    //R_CheckUserInterrupt();
                }
            }
            //std::cout<<"nbp1=["; copy(nbp.begin(),nbp.end(),std::ostream_iterator<double>(std::cout," ")); std::cout<<"]"<<std::endl<<std::flush;
            nbp+=cfpr;
            // // failure probability
            double fp=Rf_dpois(uc[j],exp(models(i,FAILR_I)),true);
            double maxp=max(nbp);
            if(maxp<(maxcfp+fp)) { maxp=maxcfp+fp; }
            nbp=(exp(nbp-maxp) + exp(cfp+fp-maxp));
            nbp/=sum(nbp);
            nbp=log(nbp);

            // find max point
            if(returnpost==1 || returnpost==3) {
                arma::uword maxij;
                double maxv=nbp.max(maxij);
                maxi.push_back(maxij);
            }
            // set the lower bound to min/n.cells
            for(unsigned int k=0;k<nbp.n_elem;k++) { if(nbp[k]<minlogprob) nbp[k]=minlogprob; }
            pm.col(j)=nbp;
        }

        ucposteriors.push_back(pm);
        if(returnpost==1 || returnpost==3) { ucmaxi.push_back(maxi);}
        //std::cout<<"."<<std::flush;
    }
    //std::cout<<" done"<<std::endl;

    //std::cout<<"boostrap iterations "<<std::flush;
    // calculate joint posterior
    int ngenes=counti.nrow();
    arma::mat jp(magnitudes.n_elem,ngenes);  // joint posterior across boostraps
    jp.zeros();
    arma::mat tjp(magnitudes.n_elem,ngenes); // for current boostrap posterior
    int seed=Rcpp::as<int>(Seed);
    srand(seed);
    int nboot=as<int>(Nboot);

    if(ensemblep) {
        for(int j=0;j<ncells;j++) {
            // exponentiate and normalize ucposteriors
            arma::mat cellucpost = exp(ucposteriors[j]);
            arma::rowvec s=sum(cellucpost,0);
            cellucpost.each_row() /= s;
            for(int k=0;k<ngenes;k++) { // fill in kth gene posterior
                jp.col(k)+=cellucpost.col(counti(k,j));
                //R_CheckUserInterrupt();
            }
            //R_CheckUserInterrupt();
        }
        arma::rowvec s=sum(jp,0); // pre-adjust so that jp is normalized
        jp.each_row() /= s;
    } else {
        if(nboot==0) { // no bootstrapping
            for(int j=0;j<ncells;j++) {
                for(int k=0;k<ngenes;k++) { // fill in kth gene posterior
                    jp.col(k)+=(ucposteriors[j]).col(counti(k,j));
                }
            }
            arma::rowvec m=max(jp,0); // calculate max for each gene (column)
            jp.each_row() -= m; // shift up for stability prior to exponentiation
            jp=exp(jp);
            arma::rowvec s=sum(jp,0); // pre-adjust so that jp is normalized
            jp.each_row() /= s;
        } else {
            for(int i=0;i<nboot;i++) {
                //std::cout<<"."<<std::flush;
                tjp.zeros();
                for(int j=0;j<ncells;j++) {
                    int rj;
                    while(ncells <= (rj= rand()/(RAND_MAX/ncells)))
                        ;
                    for(int k=0;k<ngenes;k++) { // fill in kth gene posterior
                        tjp.col(k)+=(ucposteriors[rj]).col(counti(k,rj));
                        //R_CheckUserInterrupt();
                    }
                    //R_CheckUserInterrupt();
                }
                arma::rowvec m=max(tjp,0); // calculate max for each gene (column)
                tjp.each_row() -= m; // shift up for stability prior to exponentiation
                tjp=exp(tjp);
                arma::rowvec s=sum(tjp,0)*nboot; // pre-adjust for nboot so that jp is normalized
                tjp.each_row() /= s;
                jp+=tjp;
                //R_CheckUserInterrupt();
            }
        }
    }
    //std::cout<<" done"<<std::endl;
    jp=jp.t();

    if(returnpost==1) {
        // make return matrix of individual posterior maxima
        Rcpp::NumericMatrix modes(ngenes,ncells);
        for(int i=0;i<ncells;i++) {
            for(int j=0;j<ngenes;j++) {
                modes(j,i)=magnitudes[(ucmaxi[i])[counti(j,i)]];
                //R_CheckUserInterrupt();
            }
            //R_CheckUserInterrupt();
        }
        return Rcpp::List::create(Rcpp::Named("jp") = wrap(jp),
                                  Rcpp::Named("modes") = wrap(modes));
    } else if(returnpost==2) {
        // make return matrix of individual posteriors
        Rcpp::List pl(ncells);
        for(int i=0;i<ncells;i++) {
            arma::mat ipost(magnitudes.n_elem,ngenes);
            for(int j=0;j<ngenes;j++) {
                ipost.col(j)=(ucposteriors[i]).col(counti(j,i));
                //R_CheckUserInterrupt();
            }
            ipost=ipost.t();
            pl[i]=ipost;
            //R_CheckUserInterrupt();
        }
        return Rcpp::List::create(Rcpp::Named("jp") = wrap(jp),
                                  Rcpp::Named("post") = wrap(pl));
    } else if(returnpost==3) {
        // return both modes and full posteriors
        Rcpp::NumericMatrix modes(ngenes,ncells);
        for(int i=0;i<ncells;i++) {
            for(int j=0;j<ngenes;j++) {
                modes(j,i)=magnitudes[(ucmaxi[i])[counti(j,i)]];
                //R_CheckUserInterrupt();
            }
            //R_CheckUserInterrupt();
        }
        Rcpp::List pl(ncells);
        for(int i=0;i<ncells;i++) {
            arma::mat ipost(magnitudes.n_elem,ngenes);
            for(int j=0;j<ngenes;j++) {
                ipost.col(j)=(ucposteriors[i]).col(counti(j,i));
                //R_CheckUserInterrupt();
            }
            ipost=ipost.t();
            pl[i]=ipost;
            //R_CheckUserInterrupt();
        }
        return Rcpp::List::create(Rcpp::Named("jp") = wrap(jp),
                                  Rcpp::Named("modes") = wrap(modes),
                                  Rcpp::Named("post") = wrap(pl));
    }
    // returnpost==0 // return joint posteriors only
    return wrap(jp);
}


// calculate log joint posteriors from models and unique count lists/integer on magniude grid
// count columns must match model rows
// ReturnIndividualPosterior: 0 - nothing, 1 - maxima, 2 - full posterior matrices

// similar to above, however the joint is built not over boostrap samples
// of the same pool, but by sampling a pre-defined (by Composition vector) composition
// of several different pools (specified by BatchIL)
// BatchIL - a list of indecies corresponding to the models in each batch (in the order given in Models)
// Composition - a vector giving the number of cells from each batch (batches are ordered in the same way is BatchIL)
RcppExport SEXP logBootBatchPosterior(SEXP Models, SEXP Ucl, SEXP CountsI, SEXP Magnitudes, SEXP BatchIL, SEXP Composition, SEXP Nboot, SEXP Seed, SEXP ReturnIndividualPosteriors, SEXP LocalThetaFit, SEXP SquareLogitConc) {
#define CONCB_I 0
#define CONCA_I 1
#define FAILR_I 2
#define CORRB_I 3
#define CORRA_I 4
#define CORRT_I 5
#define CORRlTB_I 6
#define CORRlTT_I 7
#define CORRlTM_I 8
#define CORRlTS_I 9
#define CORRlTR_I 10
#define CONCA2_I 11

    Rcpp::IntegerMatrix counti=Rcpp::as<Rcpp::IntegerMatrix>(CountsI);
    Rcpp::List ucl(Ucl);
    int ncells=ucl.size();
    int returnpost=Rcpp::as<int>(ReturnIndividualPosteriors);
    int localtheta=Rcpp::as<int>(LocalThetaFit);
    int squarelogitconc=Rcpp::as<int>(SquareLogitConc);
    arma::mat models=Rcpp::as<arma::mat>(Models);
    arma::colvec magnitudes=Rcpp::as<arma::colvec>(Magnitudes);
    Rcpp::List batchil(BatchIL);
    Rcpp::IntegerVector comp=as<Rcpp::IntegerVector>(Composition);

    std::vector< arma::mat > ucposteriors;
    std::vector< std::vector < arma::uword > > ucmaxi;
    // calculate individual posteriors for each cell, for each unique count value
    //std::cout<<"individual posteriors "<<std::flush;
    double minlogprob=-1*std::numeric_limits<double>::max()/ncells/1.1;
    for(int i=0;i<ncells;i++) {
        Rcpp::IntegerVector uc(Rcpp::as<Rcpp::IntegerVector>(ucl[i]));
        int ncounts=uc.size();
        arma::mat pm(magnitudes.n_elem,ncounts);
        std::vector< arma::uword > maxi;
        arma::vec mu=magnitudes * models(i,CORRA_I);
        mu+=models(i,CORRB_I); mu=exp(mu);
        arma::vec cfp;
        if(squarelogitconc) {
            cfp=models(i,CONCA_I) + magnitudes*models(i,CONCA2_I);
            cfp%=magnitudes;
        } else {
            cfp=magnitudes * models(i,CONCA_I);
        }
        cfp+=models(i,CONCB_I);
        cfp=1/(exp(cfp)+1);
        arma::vec cfpr=1-cfp;
        cfp=log(cfp); cfpr=log(cfpr);
        double maxcfp=max(cfp);
        arma::colvec thetas;
        if(localtheta) { // linear theta model - prepare theta values
            thetas=-1*magnitudes + models(i,CORRlTM_I);
            thetas*=models(i,CORRlTS_I);
            thetas=exp10(thetas)+1;
            thetas=pow(thetas,models(i,CORRlTR_I));
            thetas=(models(i,CORRlTT_I) - models(i,CORRlTB_I))/thetas;
            thetas+=models(i,CORRlTB_I);
            thetas=exp(-1*thetas);

            for(unsigned int k=0;k<thetas.n_elem;k++) {
                if((!std::isfinite(thetas[k])) || (thetas[k]<MIN_THETA)) {  thetas[k]=MIN_THETA;}
                if(thetas[k]>MAX_THETA) {  thetas[k]=MAX_THETA;}
                //R_CheckUserInterrupt();
            }
        }
        //R_CheckUserInterrupt();

        for(int j=0;j<ncounts;j++) {
            // correlated prob
            arma::vec nbp(mu.n_elem);
            if(localtheta) { // linear theta model
                for(unsigned int k=0;k<mu.n_elem;k++) {
                    double muv=mu[k]; double x=uc[j];
                    // choose maximum probability when hitting the grid with the maximum
                    if((k<(mu.n_elem-1) && x>muv && x<mu[k+1]) || (k==(mu.n_elem-1) && x>muv)) { muv=x; }
                    nbp[k]=Rf_dnbinom(x,thetas[k],thetas[k]/(thetas[k]+muv),true);
                    //R_CheckUserInterrupt();
                }
            } else { // constant theta
                double theta=models(i,CORRT_I);
                for(unsigned int k=0;k<mu.n_elem;k++) {
                    double muv=mu[k]; double x=uc[j];
                    // choose maximum probability when hitting the grid with the maximum
                    if((k<(mu.n_elem-1) && x>muv && x<mu[k+1]) || (k==(mu.n_elem-1) && x>muv)) { muv=x; }
                    nbp[k]=Rf_dnbinom(x,theta,theta/(theta+muv),true);
                }
            }
            nbp+=cfpr;
            // failure probability
            double fp=Rf_dpois(uc[j],exp(models(i,FAILR_I)),true);

            // max logp to shift by
            double maxp=max(nbp);
            if(maxp<(maxcfp+fp)) { maxp=maxcfp+fp; }
            nbp=(exp(nbp-maxp) + exp(cfp+fp-maxp));
            nbp/=sum(nbp);
            nbp=log(nbp);
            // find max point
            if(returnpost==1) {
                arma::uword maxij;
		double maxv=nbp.max(maxij);
                maxi.push_back(maxij);
            }
            // set the lower bound to min/n.cells
            for(unsigned int k=0;k<nbp.n_elem;k++) {
                if(nbp[k]<minlogprob) nbp[k]=minlogprob;
                //R_CheckUserInterrupt();
            }
            pm.col(j)=nbp;
            //R_CheckUserInterrupt();
        }
        ucposteriors.push_back(pm);
        if(returnpost==1) { ucmaxi.push_back(maxi);}
        //std::cout<<"."<<std::flush;
    }
    //std::cout<<" done"<<std::endl;

    //std::cout<<"sampling iterations "<<std::flush;
    // calculate joint posterior
    int ngenes=counti.nrow();
    arma::mat jp(magnitudes.n_elem,ngenes);  // joint posterior across boostraps
    jp.zeros();
    arma::mat tjp(magnitudes.n_elem,ngenes); // for current boostrap posterior
    int seed=Rcpp::as<int>(Seed);
    srand(seed);
    int nboot=as<int>(Nboot);

    for(int i=0;i<nboot;i++) { // sampling iteration
        //std::cout<<"."<<std::flush;
        tjp.zeros();
        for(int k=0;k<comp.size();k++) { // over batches
            int nsamp=comp[k];
            if(nsamp>0) { // sample cells from k-th batch
                Rcpp::IntegerVector bi(Rcpp::as<Rcpp::IntegerVector>(batchil[k])); // indecies of cells within the current batch
                int ncells=bi.size();
                for(int j=0;j<nsamp;j++) { // over matrices
                    int rj;
                    while(ncells <= (rj= rand()/(RAND_MAX/ncells)))
                        ;
                    for(int l=0;l<ngenes;l++) { // fill in l-th gene posterior
                        tjp.col(l)+=(ucposteriors[bi[rj]]).col(counti(l,bi[rj]));
                        //R_CheckUserInterrupt();
                    }
                    //R_CheckUserInterrupt();
                }
            }
            //R_CheckUserInterrupt();
        }
        arma::rowvec m=max(tjp,0);
        tjp.each_row() -= m; // shift up for stability prior to exponentiation
        tjp=exp(tjp);
        arma::rowvec s=sum(tjp,0)*nboot;
        tjp.each_row() /= s;
        jp+=tjp;
    }
    //std::cout<<" done"<<std::endl;
    jp=jp.t();

    if(returnpost==1) {
        // make return matrix of individual posterior maxima
        Rcpp::NumericMatrix modes(ngenes,ncells);
        for(int i=0;i<ncells;i++) {
            for(int j=0;j<ngenes;j++) {
                modes(j,i)=magnitudes[(ucmaxi[i])[counti(j,i)]];
                //R_CheckUserInterrupt();
            }
            //R_CheckUserInterrupt();
        }
        return Rcpp::List::create(Rcpp::Named("jp") = wrap(jp),
                                  Rcpp::Named("modes") = wrap(modes));
    } else if(returnpost==2) {
        // make return matrix of individual posteriors
        Rcpp::List pl(ncells);
        for(int i=0;i<ncells;i++) {
            arma::mat ipost(magnitudes.n_elem,ngenes);
            for(int j=0;j<ngenes;j++) {
                ipost.col(j)=(ucposteriors[i]).col(counti(j,i));
                //R_CheckUserInterrupt();
            }
            ipost=ipost.t();
            pl[i]=ipost;
            //R_CheckUserInterrupt();
        }
        return Rcpp::List::create(Rcpp::Named("jp") = wrap(jp),
                                  Rcpp::Named("post") = wrap(pl));
    }
    // returnpost==0 // return joint posteriors only
    return wrap(jp);
}
