#include "scatterer.h"

int Scatterer::nstop(double xL) {
  int result;
  
  if (xL <= 8) {
    result = static_cast<int>(std::round(xL + 4*std::pow(xL, 1/3) + 1));
  } else if (xL <= 4200) {
    result = static_cast<int>(std::round(xL + 4.05*std::pow(xL, 1/3) + 2));
  } else {
    result = static_cast<int>(std::round(xL + 4*std::pow(xL, 1/3) + 2));
  }
  
  return result;
}

double xmr(double x, std::complex<double> m) {
  return std::abs(x*m);
}

void Scatterer::calc_nmax() {
  std::vector<double> vres;
  vres.push_back(nstop(x.back()));              // nstop on last layer - original is L+1 which is unset
  std::transform(x.begin(), x.end(), m.begin(), std::back_inserter(vres), xmr); 
  vres.push_back(0);                            // the first value must always be zero
  std::transform(x.begin(), std::prev(x.end()), std::next(m.begin()), std::back_inserter(vres), xmr);
  auto it = std::max_element(std::begin(vres), std::end(vres));
  nmax = std::round(*it) + 15;
}

std::complex<double> calc_an(int n, double xl, std::complex<double> ha, std::complex<double> ml, 
                             std::complex<double> psixl, std::complex<double> zetaxl, 
                             std::complex<double> psixlm1, std::complex<double> zetaxlm1) {
  
  std::complex<double> m = ((((ha/ml)+std::complex<double>(n/xl,0))*psixl)-psixlm1)/((((ha/ml)+std::complex<double>(n/xl,0))*zetaxl)-zetaxlm1);
  return m;
}

std::complex<double> calc_bn(int n, double xl, std::complex<double> hb, std::complex<double> ml, 
                             std::complex<double> psixl, std::complex<double> zetaxl, 
                             std::complex<double> psixlm1, std::complex<double> zetaxlm1) {
  
  std::complex<double> m = ((((hb*ml)+std::complex<double>(n/xl,0))*psixl)-psixlm1)/((((hb/ml)+std::complex<double>(n/xl,0))*zetaxl)-zetaxlm1);
  return m;
}


 // Calculates S1_n - equation (25a)
std::complex<double> calc_S1_n(int n, std::complex<double> an, std::complex<double> bn, double Pin, double Taun) {
 return ((double)(n + n + 1)/(double)(n*n + n))*(Pin*an + Taun*bn);
}
 
 // Calculates S2_n - equation (25b) (it's the same as (25a), just switches Pin and Taun)
std::complex<double> calc_S2_n(int n, std::complex<double> an, std::complex<double> bn, double Pin, double Taun) {
  return ((double)(n + n + 1)/(double)(n*n + n))*(Taun*an + Pin*bn);
}
