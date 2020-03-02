/// \file Ex_08.cpp
/// \brief Finds eigenvalues for the Orr-Sommerfeld operator.
///
/// Solving the Eigenvalue problem:
/// \f[
/// \left(\frac{1}{Re}\partial_{yyyy} - \left(\frac{2}{Re} + (1-y^2)\,i \right)
/// \partial_{yy} + \left(\frac{1}{Re} -2i + (1-y^2)\,i \right)\right)u \;=\;
/// \lambda
/// \;
/// (\partial_{yy} - 1) u \f]

#include <fstream>
#include <iostream>
#include <sis.hpp>

using namespace std;
typedef complex<double> Cd_t;
typedef valarray<double> Vd_t;

int main() {
  using namespace sis;
  N = 255;
  sis_setup();

  double Re = 5772.0;
  Cd_t ii(0.0, 1.0);
  Linop<double> D0(0), D1(1), D2(2), D3(3), D4(4);
  D4.coef << 1.0, 0.0, 0.0, 0.0, 0.0;
  D3.coef << 1.0, 0.0, 0.0, 0.0;
  D2.coef << 1.0, 0.0, 0.0; // for (1.0 * D2 + 0.0 * D1 + (0.0) )v
  D1.coef << 1.0, 0.0;
  D0.coef << 1.0;

  LinopMat<Cd_t> Lmat(1,1), Mmat(1,1);

  Lmat << D4 / Re  - Vd_t(2.0 / Re + (1.0 - pow(y, 2.0)))*ii*D2 +
          (D0 / Re) - (2.0 * ii*D0) + ii*(Vd_t(1.0 - pow(y, 2.0))*D0);

  Mmat << D2 - D0;

  BcMat<Cd_t> bcs(4,1);
  bcs.L << D0,//
          D0,//
          D1,//
          D1;
  bcs.eval << 1.0,//
            -1.0,//
             1.0,//
            -1.0;

  GeneralizedEigenSolver<complex<double> > eigs;
  eigs.compute(Lmat, Mmat, 6, bcs);
  std::cout << "Eigenvalues are:\n" << eigs.eigenvalues << '\n';

  return 0;
}
