/*
Taken from the scattnlay program by Pena and Pal
For original paper see
doi:10.1016/j.cpc.2009.07.010  http://www.cpc.cs.qub.ac.uk/  catalog AEEY_v1_0 
 
 O. Peña, U. Pal
 Scattnlay
 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#include <Rcpp.h>
#include <iostream>

extern "C" {
  #include "ucomplex.h"     //uses a custom complex number handelling code - a <complex> rewrite might be nice
  #include "nmie.h"         // c code for individual scattering calculations 
}

using namespace Rcpp;

#define MAXLAYERS 1100
#define MAXTHETA 800

// lots of code duplication here and scatnlay doesn't need to calc theta

// [[Rcpp::export]]
DataFrame S4_AMPL(Rcpp::S4 fullstack){
  int layer_count;
  double lambda, na, mr, mi, r;
  Rcomplex mz;
  
  Rcpp::List layers;
  
  int nmax = 0; // return value from nmie
  
  double x[MAXLAYERS];
  complex m[MAXLAYERS];
  double Theta[MAXTHETA];
  double Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo;
  double ti= 0.0, tf = 90.0;
  int nt = 0; 
  complex S1[MAXTHETA], S2[MAXTHETA];
  
  lambda = fullstack.slot("lambda");
  na = fullstack.slot("na");
  layers = fullstack.slot("layers");
  nt = fullstack.slot("nt");
  ti = fullstack.slot("ti");
  tf = fullstack.slot("tf");
  layer_count = layers.size();
  
  if(nt==1)
  {
    Theta[0] = ti*PI/180.0;
  }
  else
  {
    for (int i = 0; i < nt; i++) {
      Theta[i] = (ti + (double)i*(tf - ti)/(nt - 1))*PI/180.0;
    }
  }
  
  for(int i=0;i<layer_count;i++){
    S4 S4layer((SEXP)layers[i]);  // get the layer S4 object from the list
    mz = S4layer.slot("m");
    r = S4layer.slot("r");
    
    mr = mz.r;
    mi = mz.i;
    
    // for some reason, the arrays in nmie start at 1 not 0. 
    m[i+1].r = mr/na;             //scaled values of m
    m[i+1].i = mi/na;
    x[i+1] = 2*PI*na*r/lambda;    //scaled value of x
  }
  
  // call the c code here
  nmax = nMie(layer_count, x, m, nt, Theta, &Qext, &Qsca, &Qabs, &Qbk, &Qpr, &g, &Albedo, S1, S2);
  // return S1 and S2 real and imag parts in a matrix or dataframe
  
  NumericVector C = NumericVector::create(5.8, 9.1, 3.2);
  NumericVector D = NumericVector::create(5.0, 9.3, 3.7, 5.6, 3.7);
  
  
  double data[5] = {1.0,2.0,1.0,2.0,1.0};
  NumericVector E = NumericVector(data,data+sizeof(data)/sizeof(*data)); //constructor of numeric vector? 
  
  DataFrame G = DataFrame::create(Named("S1.r")=E,Named("S2.r")=D);
  return G;
  
}

// [[Rcpp::export]]
DoubleVector S4_SCATTNLAY(Rcpp::S4 fullstack){
  int layer_count;
  double lambda, na, mr, mi, r;
  Rcomplex mz;
  
  Rcpp::List layers;
  
  int nmax = 0; // return value from nmie
  
  double x[MAXLAYERS];
  complex m[MAXLAYERS];
  double Theta[MAXTHETA];
  double Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo;
  const int nt = 0; 
  complex S1[MAXTHETA], S2[MAXTHETA];
    
  lambda = fullstack.slot("lambda");
  na = fullstack.slot("na");
  layers = fullstack.slot("layers");
  layer_count = layers.size();

  Theta[0] = 0.0;

  for(int i=0;i<layer_count;i++){
    S4 S4layer((SEXP)layers[i]);  // get the layer S4 object from the list
    mz = S4layer.slot("m");
    r = S4layer.slot("r");
  
    mr = mz.r;
    mi = mz.i;

   // for some reason, the arrays in nmie start at 1 not 0. 
    m[i+1].r = mr/na;             //scaled values of m
    m[i+1].i = mi/na;
    x[i+1] = 2*PI*na*r/lambda;    //scaled value of x
  }
  
  // call the c code here
  nmax = nMie(layer_count, x, m, nt, Theta, &Qext, &Qsca, &Qabs, &Qbk, &Qpr, &g, &Albedo, S1, S2);
  DoubleVector z = DoubleVector::create(Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, nmax);
  return z; 
}