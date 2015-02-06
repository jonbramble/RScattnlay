
/*
Taken from the scattnlay program by Pena and Pal
For original paper see
doi:10.1016/j.cpc.2009.07.010

No licence could be found in the original program code
*/

#include <Rcpp.h>
#include <iostream>

extern "C" {
  #include "ucomplex.h"     //uses a custom complex number handelling code - a <complex> rewrite might be nice
  #include "nmie.h"         // c code for individual scattering calculations 
}

#define MAXLAYERS 1100
#define MAXTHETA 800

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector S4_SCATTNLAY(S4 fullstack){
  // convert the stack to component parts
  // call scatnlay with these values
  // return the set of values to the S4 class
  int layer_count;
  double lambda, na;
  Rcpp::List layers;
  
  double layer_d;
  
  lambda = fullstack.slot("lambda");
  na = fullstack.slot("na");
  layers = fullstack.slot("layers");
  layer_count = layers.size();
  
  for(int i=0;i<layer_count;i++){
    S4 S4layer((SEXP)layers[i]);// get the layer S4 object from the list
    layer_d = S4layer.slot("d");//get the data from the slots
    //m[i+1].r = 
    
  }
  
}

// [[Rcpp::export]]
NumericVector scattnlay_test(){
  
  int nmax = 0; // return value from nmie
  int L = 2; // number of layers
  int nTheta;
  
  double x[MAXLAYERS],Theta[MAXLAYERS];
  complex m[MAXLAYERS];
  double Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo;
  double ti= 0.0, tf = 90.0;
  int nt = 0;  // what is nt - if you want to spec the theta angles
  complex S1[MAXLAYERS], S2[MAXLAYERS];
  
  // fill in the values for theta
  if(nt>1)
  {
   for (int i = 0; i < nt; i++) {
    Theta[i] = (ti + (double)i*(tf - ti)/(nt - 1))*PI/180.0;
   }
  }
  
  // some test data from one of the examples
  // for some reason, the arrays in nmie start at 1 not 0. 
  
  x[1]= 0.099668871748;
  m[1].r = 1.33;
  m[1].i = 0.00;
  
  x[2]= 0.1;
  m[2].r = 1.59;
  m[2].i = 0.66;
  
  nmax = nMie(L, x, m, nt, Theta, &Qext, &Qsca, &Qabs, &Qbk, &Qpr, &g, &Albedo, S1, S2);
  
  // test data should return this:
  //0.1, +1.39348e-03, +1.12219e-05, +1.38226e-03, +1.67581e-05, +1.39346e-03, +1.64365e-03, +8.05315e-03


  NumericVector z = NumericVector::create(Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, nmax);
  return z ;
}


  //std::cout << Qsca << std::endl;
  //convert via a std::vector - slow
  
  //std::vector<double> dataVec;
  //unsigned dataArraySize = sizeof(Theta) / sizeof(double);
  //dataVec.insert(dataVec.end(), &Theta[0], &Theta[dataArraySize]);
  
  //NumericVector y = NumericVector::create( 0.0, 1.0 ) ;
  //NumericVector NvQsca( dataVec.begin(), dataVec.end() );
