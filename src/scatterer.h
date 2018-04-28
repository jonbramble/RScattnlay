#include <complex>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

class Scatterer {
public:
  
  Scatterer(int n, std::vector<double> _x, std::vector<std::complex<double>> _m): L(n), x(_x), m(_m){
    calc_nmax();
    std::cout << "scatterer in cpp" << std::endl;
    std::cout << "n_max " << nmax << std::endl;
  }
  
private:
  const int L; //number of layers
  int nmax;
  std::vector<double> x;
  std::vector<std::complex<double>> m;
  
  int nstop(double);
  
  // calculate the number of iterations required
  void calc_nmax();
  
};