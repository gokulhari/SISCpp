/// \file Ex_06.cpp
/// \brief This solves the generalized eigenvalue problem \f$ L \phi = \lambda M
/// \phi \f$.
///
/// Solves the generalized eigenvalue problem,
/// where, \f$ L = \partial_{yyyy} - 8\partial_{yy} + 16,\;\; u(\pm 1) = u''(\pm
/// 1) = 0\f$ and \f$ M = \partial_{yy} - 4,\;\; u(\pm 1) = 0 \f$.

#include <fstream> // To output data to files
#include <iostream>
#include <sis.hpp>

using namespace std;

typedef complex<double> Cd_t;
typedef valarray<double> Vd_t;

int main() {
  using namespace sis;
  N = 63; // sis::N, Defined in sis.hpp. This specifies no. of Chebyshev
          // coefficients
  sis_setup();
  ChebfunMat<complex<double> > forc(1,1);

  Linop<double> D0(0), D1(1), D2(2), D3(3), D4(4);
  D4.coef << 1.0, 0.0, 0.0, 0.0, 0.0;
  D3.coef << 1.0, 0.0, 0.0, 0.0;
  D2.coef << 1.0, 0.0, 0.0; // for (1.0 * D2 + 0.0 * D1 + (0.0) )v
  D1.coef << 1.0, 0.0;
  D0.coef << 1.0;

  LinopMat<Cd_t> Lmat(1,1), Mmat(1,1);
  Lmat << D4 - (8.0 * D2) + (16.0 * D0);
  Mmat << D2 - (4.0 * D0);
  BcMat<Cd_t> bcs(4,1);
  bcs.L << D0,//
          D2,//
          D0,//
          D2;
  bcs.eval << -1.0,//
              -1.0,//
              1.0,//
              1.0;

  GeneralizedEigenSolver<Cd_t> eigs;

  eigs.compute(Lmat, Mmat, 6, bcs);
  std::cout << "Eigenvalues: \n" << '\n';
  std::cout << eigs.eigenvalues << '\n';
  return 0;
}
