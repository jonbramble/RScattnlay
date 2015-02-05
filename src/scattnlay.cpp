
/*
Taken from the scattnlay program by Pena and Pal
For original paper see
doi:10.1016/j.cpc.2009.07.010
*/

#include <Rcpp.h>
#include <iostream>

extern "C" {
  #include "ucomplex.h"     //uses a custom complex number handelling code
  #include "nmie.h" 
}

#define MAXLAYERS 1100
#define MAXTHETA 800

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector scattnlay(){
  
  int nmax = 0; // return value from nmie
  int L = 2; // number of layers
  
  double x[MAXLAYERS];
  complex m[MAXLAYERS];
  int nTheta;
  double Theta[MAXLAYERS];
  double Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo;
  double ti= 0.0, tf = 90.0;
  int nt = 0;  // what is nt?
  complex S1[MAXLAYERS], S2[MAXLAYERS];
  
  if(nt>1)
  {
   for (int i = 0; i < nt; i++) {
    Theta[i] = (ti + (double)i*(tf - ti)/(nt - 1))*PI/180.0;
   }
  }
  
  x[0]= 0.099668871748;
  m[0].r = 1.33;
  m[0].i = 0.00;
  
  x[1]= 0.1;
  m[1].r = 1.59;
  m[1].i = 0.66;
  
  nmax = nMie(L, x, m, nt, Theta, &Qext, &Qsca, &Qabs, &Qbk, &Qpr, &g, &Albedo, S1, S2);
  
  //0.1, +1.39348e-03, +1.12219e-05, +1.38226e-03, +1.67581e-05, +1.39346e-03, +1.64365e-03, +8.05315e-03

  //std::cout << Qsca << std::endl;
  //convert via a std::vector - slow
  std::vector<double> dataVec;
  unsigned dataArraySize = sizeof(Theta) / sizeof(double);
  dataVec.insert(dataVec.end(), &Theta[0], &Theta[dataArraySize]);
  
  NumericVector y = NumericVector::create( 0.0, 1.0 ) ;
  NumericVector NvQsca( dataVec.begin(), dataVec.end() );
  NumericVector z = NumericVector::create(Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, nmax);
  return z ;
}

/*List rcpp_hello_world() {

    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y ) ;

    return z ;
}*/
