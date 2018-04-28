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



void Scatterer::calc_nmax() {
  int i, result;
  
  result = nstop(x[L]);
  for (i = 0; i < L; i++) {
    double xm_mag = std::abs(x[i]*m[i]);
    if (result < xm_mag) {
      result = static_cast<int>(std::round(xm_mag));
    }
    xm_mag = std::abs(x[i-1]*m[i]); // this is ackward then becuase x is not actually set!
    if (result < xm_mag) {
      result = static_cast<int>(std::round(xm_mag));
    }
  }
  nmax = result + 15;
}
