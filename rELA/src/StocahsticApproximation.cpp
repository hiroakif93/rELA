// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace arma;
using namespace Rcpp;

//======================================================//

//////////////////////////////////////////////////////////
//              Stchastic Approximationes               //
//////////////////////////////////////////////////////////
// -- Functions for Stochastic Approximation
arma::mat Logmodel(arma::mat m, arma::mat alpha, arma::mat beta){
  
  m *= beta;
  mat res = 1/(exp(-alpha - m)+1);
  
  return res;
}

arma::mat Logmodel_simple(arma::mat m, arma::mat alpha, arma::mat beta){
  
  m *= beta;

  mat res = 1/(exp(-alpha - m.each_row())+1);
  
  return res;
}

// One step "Heat Bath Sampling"
arma::mat OnestepHBS(arma::mat y, arma::mat lm) {
  
  // Preparation // 
  int itr1=y.n_cols;
  int itr2=y.n_rows;

  for(int r=0; r < itr2; r++){
    
    // Initial state //
    int uni=R::runif(0, itr1);
    double ne=lm(r, uni);
    
    double randomNum=randu<double>();
    
    // adjacency state//
    if(ne > randomNum){
      y(r, uni) = 1;
    }else{
      y(r, uni) = 0;
    }
    
  }
  
  return y;
}


arma::mat Logprior(arma::mat alpha, double x) {
  
  mat y = zeros(alpha.n_rows, alpha.n_cols);
  y.fill(x);
  mat lp = -tanh( (alpha/y)/2 )/y ;
  return lp;
}

//////////////////////////////////////////////////////////

//========  only Inplicit environmetanl data ===========//

// [[Rcpp::export]]
arma::mat SA_simple( arma::mat ocData, double maxInt=50000, double momentum=0.3){
  
  // ========================================= //
  // -- Initial parameter setting
  
  double nlocation = ocData.n_rows;
  double nspecies = ocData.n_cols;
  
  mat ystats=trans(ocData);
  ystats *= ocData;
  
  mat ysim = zeros(nlocation, nspecies);
  
  mat lp = ( sum(ocData, 0)+1 )/(nlocation+2);
  mat alphas = log(lp/(1-lp));

  mat beta = zeros(nspecies, nspecies);
  
  mat delalphas=zeros(1, nspecies); 
  mat delbeta = zeros(nspecies, nspecies);
  
  mat betagrad;
  mat alphasgrad;
  
  mat logmat;  mat ydif ;
  //double learningrate0=0.1;
  
  // ========================================= //
  // Main part
  double learningrate=0.1;
  double  momtm=0.3;
  for(int tt=0; tt < maxInt; tt++){
    //double learningrate=learningrate0*1000/(998+1+tt);
    //double momtm=0.9*(1-1/(0.1*tt+2));
    
    // -- Preference for co-occurence of species i and j
    logmat=Logmodel_simple(ysim, alphas, beta);
    ysim=OnestepHBS(ysim, logmat);
    
    //mat ysimstats=trans(ysim) * ysim;
    ydif = ystats - (trans(ysim) * ysim);
    
    // -- Beta gradient 
    betagrad = (ydif + Logprior(beta, 0.5)) / (nlocation * abs(eye(nspecies, nspecies)-1)) ;
    betagrad.diag().fill(0);
    
    // -- Alpha gradient
    alphasgrad = (ydif.diag().t() + Logprior(alphas, 2))/(mat(1, alphas.n_cols).fill(nlocation));
    
    // -- delta
    betagrad %= mat(betagrad.n_rows, betagrad.n_cols).fill((1-momtm)*learningrate);
    delbeta %= mat(betagrad.n_rows, betagrad.n_cols).fill(momtm);
    delbeta += betagrad;

    alphasgrad %= mat(alphasgrad.n_rows, alphasgrad.n_cols).fill((1-momtm)*learningrate) ;
    delalphas %= mat(alphasgrad.n_rows, alphasgrad.n_cols).fill(momtm);
    delalphas += alphasgrad ;
    
    beta+=delbeta; 
    alphas+=delalphas; 
  }
  
  // ========================================= //

  
  mat ystatsList=join_rows(alphas.t(),beta.t());
  return ystatsList;
}

//////////////////////////////////////////////////////////

//========  including Explicit environmetanl data ======//

// [[Rcpp::export]]
arma::mat SA( arma::mat ocData, arma::mat envData, double maxInt=50000, double momentum=0.3){
  
  // ========================================= //
  // -- Initial parameter setting
  
  double nlocation = ocData.n_rows;
  double nspecies = ocData.n_cols;
  double nenvironment = envData.n_cols;
  
  mat ystats = trans(ocData) * ocData;
  mat yenvstats = trans(envData) * ocData;
  
  mat ysim = zeros(nlocation, nspecies);
  mat alphae = zeros(nenvironment, nspecies);
  
  mat lp = ( sum(ocData, 0)+1 )/(nlocation+2);
  mat alphas = log(lp/(1-lp));
  
  mat beta = zeros(nspecies, nspecies);
  
  mat delalphas=zeros(1, nspecies); 
  mat delalphae = zeros(nenvironment, nspecies ); 
  mat delbeta = zeros(nspecies, nspecies);
  
  mat betagrad = zeros(nspecies, nspecies);
  mat alphasgrad = zeros(1, nspecies);
  mat alphaegrad = zeros(nspecies, nspecies);
  mat alpha; mat logmat; mat ydif; mat yenvdiff;
  
  //double learningrate0=0.1;
  double learningrate0=0.1;
  // ========================================= //
  // Main part
  
  for(int tt=0; tt < maxInt; tt++){
    double learningrate=learningrate0*1000/(998+1+tt);
    double momtm=0.9*(1-1/(0.1*tt+2));
    
    // -- Preference for co-occurence of species i and j
    alpha=mat(envData * alphae).each_row() + alphas;
    
    // -- 
    logmat=Logmodel(ysim, alpha, beta);
    ysim=OnestepHBS(ysim, logmat);

    ydif = ystats - (trans(ysim) * ysim);
    yenvdiff= yenvstats-(trans(envData) * ysim);
    
    // -- Beta gradient 
    betagrad = (ydif + Logprior(beta, 0.5)) / (nlocation * abs(eye(nspecies, nspecies)-1) );
    betagrad.diag().fill(0);
    
    // -- Alpha gradient
    alphasgrad = (ydif.diag().t() + Logprior(alphas, 2))/mat(1, alphas.n_cols).fill(nlocation);
    alphaegrad = (yenvdiff+Logprior(alphae,2))/mat(1, alphae.n_cols).fill(nlocation);
    
    // -- delta
    betagrad %= mat(betagrad.n_rows, betagrad.n_cols).fill((1-momtm)*learningrate);
    delbeta %= mat(betagrad.n_rows, betagrad.n_cols).fill(momtm);
    delbeta += betagrad;
      
    alphasgrad %=mat(alphasgrad.n_rows, alphasgrad.n_cols).fill((1-momtm)*learningrate) ;
    delalphas %=mat(alphasgrad.n_rows, alphasgrad.n_cols).fill(momtm);
    delalphas+=alphasgrad ;   

    alphaegrad %= mat(alphaegrad.n_rows, alphaegrad.n_cols).fill((1-momtm)*learningrate);
    delalphae %= mat(alphaegrad.n_rows, alphaegrad.n_cols).fill(momtm);
    delalphae += alphaegrad;
    
    beta+=delbeta; 
    alphas+=delalphas; 
    alphae+=delalphae;
    
  }
  
  // ========================================= //
  
  //
  
  mat ystatsList=join_rows(alphas.t(), alphae.t(), beta.t());
  
  return ystatsList;
}

//////////////////////////////////////////////////////////