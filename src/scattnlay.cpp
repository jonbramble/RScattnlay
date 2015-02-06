/*
Taken from the scattnlay program by Pena and Pal
For original paper see
doi:10.1016/j.cpc.2009.07.010

No licence could be found in the original program code
*/

#include <Rcpp.h>
#include <iostream>
#include <complex>

extern "C" {
  #include "ucomplex.h"     //uses a custom complex number handelling code - a <complex> rewrite might be nice
  #include "nmie.h"         // c code for individual scattering calculations 
}

using namespace Rcpp;

#define MAXLAYERS 1100
#define MAXTHETA 800

// [[Rcpp::export]]
NumericVector S4_SCATTNLAY(Rcpp::S4 fullstack){
  int layer_count;
  double lambda, na, mr, mi, d;
  Rcomplex mz;
  
  Rcpp::List layers;
  
  int nmax = 0; // return value from nmie
  int nTheta;
  
  double x[MAXLAYERS],Theta[MAXLAYERS];
  complex m[MAXLAYERS];
  double Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo;
  double ti= 0.0, tf = 90.0;
  int nt = 0;  // what is nt - if you want to spec the theta angles
  complex S1[MAXLAYERS], S2[MAXLAYERS];
  
  // Not implemented in S4 yet
  // fill in the values for theta
  if(nt>1)
  {
   for (int i = 0; i < nt; i++) {
    Theta[i] = (ti + (double)i*(tf - ti)/(nt - 1))*PI/180.0;
   }
  }
    
  lambda = fullstack.slot("lambda");
  na = fullstack.slot("na");
  layers = fullstack.slot("layers");
  layer_count = layers.size();
  
  for(int i=0;i<layer_count;i++){
    S4 S4layer((SEXP)layers[i]);  // get the layer S4 object from the list
    mz = S4layer.slot("m");
    mr = mz.r;
    mi = mz.i;
    
   // some test data from one of the examples
   // for some reason, the arrays in nmie start at 1 not 0. 
  
    m[i+1].r = mr/na;             //scaled values of m
    m[i+1].i = mi/na;
     
    d = S4layer.slot("d");
    
    x[i+1] = 2*PI*na*d/lambda;    //scaled value of x
  }
  
  // call the c code here
  nmax = nMie(layer_count, x, m, nt, Theta, &Qext, &Qsca, &Qabs, &Qbk, &Qpr, &g, &Albedo, S1, S2);
  NumericVector z = NumericVector::create(Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, nmax);
  return z; 
}

  //std::cout << Qsca << std::endl;
  //convert via a std::vector - slow
  
  //std::vector<double> dataVec;
  //unsigned dataArraySize = sizeof(Theta) / sizeof(double);
  //dataVec.insert(dataVec.end(), &Theta[0], &Theta[dataArraySize]);
  
  //NumericVector y = NumericVector::create( 0.0, 1.0 ) ;
  //NumericVector NvQsca( dataVec.begin(), dataVec.end() );
  
  
