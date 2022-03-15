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

//////////////////////////////////////////////////////////
//              Energy landscape analysis               //
//////////////////////////////////////////////////////////

// -- Functions for Energy landscape analysis 

// -- Community Energy
// [[Rcpp::export]]
double cEnergy(arma::rowvec state, arma::rowvec alpha, arma::mat beta){
  
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
  int SMrow=sequ.n_elem+1;
  int SMcol=s1.n_elem;
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
int checkIdent(arma::mat ss, arma::rowvec ysim)
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

// Entropy
// [[Rcpp::export]]
double entropy(arma::vec v)
{
  vec uniq=unique(v);
  int N=uniq.n_elem;
  double total=v.n_elem;
  double ent;
  
  for(int i=0; i<N; i++){
    double prob=sum(v==uniq(i))/total;
    ent+=prob*log(prob);
  }
 // Rcout << probV << endl;
  
  return ent*(-1);
}

// [[Rcpp::export]]
arma::mat SSentropy(arma::mat uoc, arma::mat ss,
                    arma::rowvec alpha, arma::mat beta, 
                    int seitr=1000, int convTime=10000){
  
  mat entropyres=zeros(uoc.n_rows, 3);
  for(int i=0; i < uoc.n_rows; i++){
    
    mat logmat;
    rowvec ysim=uoc.row(i);
    int tt=0;
    mat simures=zeros(seitr, 3);
    
    for(int l=0; l < seitr; l++){
      
      do{
        tt+=1;
        logmat=Logmodel_simple(ysim, alpha, beta);
        ysim=OnestepHBS(ysim, logmat);
        
        if(checkIdent(ss=ss, ysim)==1){
          break;
        }
      } while ( tt<convTime );
      
      simures(l, 0)=convDec(ysim);
      simures(l, 1)=tt;
      simures(l, 2)=((tt+1)==convTime);
      
    }
    
    entropyres(i,0)=entropy(simures.col(0));
    entropyres(i,1)=mean(simures.col(1));
    entropyres(i,2)=sum(simures.col(2));
  }

  return entropyres;
}

