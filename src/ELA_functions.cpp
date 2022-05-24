// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppArmadilloExtensions/sample.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <random>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

//======================================================//

//////////////////////////////////////////////////////////
//              Stchastic Approximationes               //
//////////////////////////////////////////////////////////
inline arma::mat Logmodel(arma::mat m,  const arma::mat& alpha, const arma::mat& beta){
  
  m *= beta;
  const arma::mat& res = 1/(exp(-alpha - m)+1);
  
  return res;
}

inline arma::mat Logmodel_simple(arma::mat m, const arma::mat& alpha, const arma::mat& beta){
  
  m *= beta;
  const arma::mat& res = 1/(exp(-alpha - m.each_row())+1);
  
  return res;
}
// One step "Heat Bath Sampling"
arma::mat OnestepHBS(arma::mat y, const arma::mat& lm) {
  
  // Preparation // 
  const int& itr1=y.n_cols;
  const int& itr2=y.n_rows;
  
  for(int r=0; r < itr2; ++r){
    
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

// One step "Heat Bath Sampling"
void OnestepHBSvoid(arma::mat& y, const arma::mat& lm) {
  
  // Preparation // 
  const int& itr1=y.n_cols;
  const int& itr2=y.n_rows;
  
  for(int r=0; r < itr2; ++r){
    
    // Initial state //
    const int& uni=R::runif(0, itr1);
    const double& ne=lm(r, uni);
    
    const double& randomNum=randu<double>();
    
    // adjacency state//
    if(ne > randomNum){
      y(r, uni) = 1;
    }else{
      y(r, uni) = 0;
    }
    
  }
  
}

inline arma::mat Logprior(const arma::mat& alpha, const double& x) {
  
  mat y = zeros(alpha.n_rows, alpha.n_cols);
  y.fill(x);
  mat lp = -tanh( (alpha/y)/2 )/y ;
  return lp;
}

//////////////////////////////////////////////////////////

//========  only Inplicit environmetanl data ===========//

// [[Rcpp::export]]
arma::mat SA_simple( const arma::mat& ocData, const double& maxInt=50000, const double& momentum=0.3){
  
  // ========================================= //
  // -- Initial parameter setting
  
  const double& nlocation = ocData.n_rows;
  const double& nspecies = ocData.n_cols;
  
  mat ystats=trans(ocData);
  ystats *= ocData;
  
  mat ysim = zeros(nlocation, nspecies);
  
  const mat& lp = ( sum(ocData, 0)+1 )/(nlocation+2);
  mat alphas = log(lp/(1-lp));
  
  mat beta = zeros(nspecies, nspecies);
  
  mat delalphas=zeros(1, nspecies); 
  mat delbeta = zeros(nspecies, nspecies);
  
  mat betagrad;
  mat alphasgrad;
  
  mat logmat;  mat ydif ;
  const double& learningrate0=0.1;
  
  const mat& betaconst=(nlocation * abs(eye(nspecies, nspecies)-1));
  const rowvec& asconst=(mat(1, nspecies).fill(nlocation));
  
  // ========================================= //
  // Main part
  //double learningrate=0.1;
  //double  momtm=0.3;
  for(int tt=0; tt < maxInt; ++tt){
    const double& learningrate=learningrate0*1000/(998+1+tt);
    const double& momtm=0.9*(1-1/(0.1*tt+2));
    
    // -- Preference for co-occurence of species i and j
    logmat=Logmodel_simple(ysim, alphas, beta);
    OnestepHBSvoid(ysim, logmat);
    
    //mat ysimstats=trans(ysim) * ysim;
    ydif = ystats - (trans(ysim) * ysim);
    
    // -- Beta gradient 
    betagrad = (ydif + Logprior(beta, 0.5)) / betaconst ;
    betagrad.diag().fill(0);

    // -- Alpha gradient
    alphasgrad = (ydif.diag().t() + Logprior(alphas, 2))/asconst;
  
    // -- delta
    betagrad %= mat(nspecies, nspecies).fill((1-momtm)*learningrate);
    delbeta %= mat(nspecies, nspecies).fill(momtm);
    delbeta += betagrad;
    
    alphasgrad %= mat(1, nspecies).fill((1-momtm)*learningrate) ;
    delalphas %= mat(1, nspecies).fill(momtm);
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
arma::mat SA( const arma::mat& ocData, const arma::mat& envData, const double& maxInt=50000, const double& momentum=0.3){
  
  // ========================================= //
  // -- Initial parameter setting
  
  const double& nlocation = ocData.n_rows;
  const double& nspecies = ocData.n_cols;
  const double& nenvironment = envData.n_cols;
  
  mat ystats = trans(ocData) * ocData;
  mat yenvstats = trans(envData) * ocData;
  
  mat ysim = zeros(nlocation, nspecies);
  mat alphae = zeros(nenvironment, nspecies);
  
  const mat& lp = ( sum(ocData, 0)+1 )/(nlocation+2);
  mat alphas = log(lp/(1-lp));
  
  mat beta = zeros(nspecies, nspecies);
  
  mat delalphas=zeros(1, nspecies); 
  mat delalphae = zeros(nenvironment, nspecies ); 
  mat delbeta = zeros(nspecies, nspecies);
  
  mat betagrad = zeros(nspecies, nspecies);
  mat alphasgrad = zeros(1, nspecies);
  mat alphaegrad = zeros(nspecies, nspecies);
  mat alpha; mat logmat; mat ydif; mat yenvdiff;
  
  const mat& betaconst=(nlocation * abs(eye(nspecies, nspecies)-1));
  const rowvec& asconst=(mat(1, nspecies).fill(nlocation));
  const rowvec& aeconst=(mat(nenvironment, nspecies).fill(nlocation));
  
  const double& learningrate0=0.1;
  // ========================================= //
  // Main part
  
  for(int tt=0; tt < maxInt; ++tt){
    const double& learningrate=learningrate0*1000/(998+1+tt);
    const double& momtm=0.9*(1-1/(0.1*tt+2));
    
    // -- Preference for co-occurence of species i and j
    alpha=mat(envData * alphae).each_row() + alphas;
    
    // -- 
    logmat=Logmodel(ysim, alpha, beta);
    OnestepHBSvoid(ysim, logmat);
    
    ydif = ystats - (trans(ysim) * ysim);
    yenvdiff= yenvstats-(trans(envData) * ysim);
    
    // -- Beta gradient 
    betagrad = (ydif + Logprior(beta, 0.5)) / betaconst;
    betagrad.diag().fill(0);
    
    // -- Alpha gradient
    alphasgrad = (ydif.diag().t() + Logprior(alphas, 2))/asconst;
    alphaegrad = (yenvdiff+Logprior(alphae,2))/aeconst;

    // -- delta
    betagrad %= mat(nspecies, nspecies).fill((1-momtm)*learningrate);
    delbeta %= mat(nspecies, nspecies).fill(momtm);
    delbeta += betagrad;
    
    alphasgrad %=mat(1, nspecies).fill((1-momtm)*learningrate) ;
    delalphas %=mat(1, nspecies).fill(momtm);
    delalphas+=alphasgrad ;   
    
    alphaegrad %= mat(nenvironment, nspecies).fill((1-momtm)*learningrate);
    delalphae %= mat(nenvironment, nspecies).fill(momtm);
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

//////////////////////////////////////////////////////////
//              Energy landscape analysis               //
//////////////////////////////////////////////////////////

// -- Functions for Energy landscape analysis 

// -- Community Energy
// [[Rcpp::export]]
inline double cEnergy(const arma::rowvec& state, const arma::rowvec& alpha, const arma::mat& beta){
  
  mat res = -state * alpha.t() - (state* (state * beta).t() ) / 2;
  return as_scalar(res);
}

//////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::mat SteepestDescent_cpp(arma::rowvec state, arma::rowvec alpha, arma::mat beta){
  
  // ================================ //
  int term=0; 
  mat energi; 
  mat y1 = conv_to<mat>::from(state);
  double y2 = cEnergy(state, alpha, beta); 
  // ================================ //
  
  do{
    // ============================= //
    mat ystat = y1;
    double yene = y2;
    // ============================= //
    
    // -- Energy
    energi = -alpha-ystat*beta;
    mat sign = (-2 * ystat +1);
    energi %=sign; energi +=yene;
    
    double minene =energi.min() ;
    if( minene < yene ){
      uword mp = energi.index_min();
      ystat(0, mp) = abs(ystat(0, mp)-1);
      
      y1.swap(ystat);
      y2 = minene;
      
    }else{
      term+=1;
    }
  } while (term==0);
  
  return join_rows(y1, mat(1,1).fill(y2));
  
}
// [[Rcpp::export]]
arma::mat SSestimate(arma::rowvec alpha, arma::mat beta, int itr=20000){
  
  mat res = zeros(itr, beta.n_cols+1);
  
  for(int i=0; i<itr; i++ ){
    
    imat intmat=randi(1, beta.n_cols, arma::distr_param(0, 1));
    rowvec state=conv_to<rowvec>::from(intmat);
    
    mat ss = SteepestDescent_cpp(state, alpha, beta);
    
    res.row(i) = ss;
  }
  return res;
}

//////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::rowvec FindingTippingpoint_cpp(arma::rowvec s1,arma:: rowvec s2, 
                                     arma::rowvec alpha, arma::mat jj,
                                     int tmax=10000){
  
  // ======================================= //
  // Distance between stable states
  uvec pd=find( abs(s1-s2)==1);
  rowvec sequ=conv_to<rowvec>::from(shuffle(pd));
  rowvec tipState(alpha.n_elem+1);
  mat samplePath;
  
  // ======================================= //
  // ||||||||||||||||||||||||||||||||||| //
  // -- Definition
  double bde=datum::inf;
  int tt=0 ;
  double minbde=datum::inf;
  uvec rpl;
  rowvec seqnew;
  double bdenew;
  int cc;
  const int& SMrow=sequ.n_elem+1;
  const int& SMcol=s1.n_elem;
  vec::fixed<2>c; 
  // ||||||||||||||||||||||||||||||||||| //
  
  do{
    tt+=1;
    // ||||||||||||||||||||||||||||||||||| //
    seqnew=sequ;
    rpl = randperm(sequ.n_elem, 2);
    
    seqnew.swap_cols(rpl(0), rpl(1));
    
    samplePath=join_rows(repmat(s1, SMrow, 1), zeros(SMrow,1));
    
    // ||||||||||||||||||||||||||||||||||| //
    for( int i=1; i<(SMrow); i++){
      
      // change column
      cc=conv_to<int>::from(seqnew.col(i-1)) ;
      // replace mat
      mat pp((SMrow-(i)), 1);
      pp.fill(conv_to<int>::from( abs(samplePath.submat(i, cc, i, cc)-1) ) );
      // reolace
      samplePath.submat(i, 0, SMrow-1, SMcol-1).col(cc)=pp;
      samplePath(i-1, SMcol)=cEnergy(conv_to<rowvec>::from(samplePath.submat(i-1, 0, i-1, SMcol-1)), alpha, jj);
      
    }
    samplePath(SMrow-1, SMcol)=cEnergy(conv_to<rowvec>::from(samplePath.submat(SMrow-1, 0, SMrow-1, SMcol-1)), alpha, jj);
    // ||||||||||||||||||||||||||||||||||| //
    
    bdenew=samplePath.col(SMcol).max();
    if(bdenew<minbde){
      minbde=std::move(bdenew);
      tipState.cols(0, SMcol-1)=samplePath.submat(samplePath.col(SMcol).index_max(), 0, samplePath.col(SMcol).index_max(), SMcol-1);
      tipState.col(tipState.n_cols-1)=minbde;
      
    }
    
    c(0)=1; c(1)= ( (exp(bde))/(exp(bdenew)) );
    if( randu<double>() < c.min()){
      sequ.swap(seqnew);
      bde=std::move(bdenew);
    }
    
    // ||||||||||||||||||||||||||||||||||| //
  } while (tt<tmax);
  
  // ||||||||||||||||||||||||||||||||||| //
  
  return tipState;
}

//////////////////////////////////////////////////////////

// To stop calculation.
int checkIdent(const arma::mat& ss, const arma::rowvec& ysim)
{
  int rows=ss.n_rows;
  vec res(rows);
  
  for(int l=0; l<rows; l++){
    res(l)=all(ysim==ss.row(l));
  }
  
  return any(res);
}

// Convert to decimal
int convDec(arma::rowvec v)
{
  int tmp;
  std::string str;
  for(int i=0; i< v.n_elem; i++){
    tmp=v(i);
    str+=std::to_string(tmp);
  }
  
  return std::stoi(str, 0, 2); 
  
}

// Convert to string
inline std::string convStr(const arma::rowvec& v)  {
  int tmp;
  std::string str;
  for(int i=0; i< v.n_elem; i++){
    tmp=v(i);
    str+=std::to_string(tmp);
  }
  
  return str; 
  
}

// Entropy numeric
inline double entropy(const arma::vec& v)
{
  vec uniq=unique(v);
  int N=uniq.n_elem;
  const double& total=v.n_elem;
  double ent=0;
  
  for(int i=0; i<N; i++){
    const double& prob=sum(v==uniq(i))/total;
    ent+=prob*log(prob);
  }
  // Rcout << probV << endl;
  
  return ent*(-1);
}

// Entropy string
inline double entropy2(const Rcpp::StringVector& v)
{
  IntegerVector tab=table(v);
  const double& total=v.length();
  //NumericVector numtab=as<NumericVector>(tab);
  
  double ent=0;
  for(int i=0; i<tab.size(); i++){
    const double& prob=tab[i]/total;
    ent+=prob*log(prob);
  }
  
  return ent*(-1);
}

// [[Rcpp::export]]
arma::mat SSentropy_cpp(arma::mat uoc, arma::mat ss,
                        arma::rowvec alpha, arma::mat beta, 
                        int seitr=1000, int convTime=10000){
  
  mat entropyres=zeros(uoc.n_rows, 3);
  for(int i=0; i < uoc.n_rows; i++){
    
    mat logmat;
    rowvec ysim=uoc.row(i);
    int tt=0;
    Rcpp::StringVector ssid(seitr);
    mat stable=zeros(seitr, 2);
    
    for(int l=0; l < seitr; l++){
      
      do{
        tt+=1;
        logmat=Logmodel_simple(ysim, alpha, beta);
        ysim=OnestepHBS(ysim, logmat);
        
        if(checkIdent(ss, ysim)==1){
          break;
        }
      } while ( tt<convTime );
      
      ssid(l)=convStr(ysim);
      stable(l, 0)=tt;
      stable(l, 1)=((tt+1)==convTime);
      
    }
    
    entropyres(i,0)=entropy2(ssid);
    entropyres(i,1)=mean(stable.col(0));
    entropyres(i,2)=sum(stable.col(1));
  }
  
  return entropyres;
}
