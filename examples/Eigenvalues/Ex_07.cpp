/// \file Ex_07.cpp
/// \brief This solves the generalized eigenvalue problem for a block matrix
/// operator system
///
/// The eigenvalue problem is given by:
/// \f[
/// \left[
/// \begin{array}{cc}
/// \partial_{yyyy}/16 + 2y & 1.0 \\
/// \sin (2y) & \partial_{yy}/4
/// \end{array}
/// \right]
/// \left[
/// \begin{array}{c}
/// u \\
/// v
/// \end{array}
/// \right] \;=\;
/// \lambda
/// \left[
/// \begin{array}{cc}
/// \partial_{yy}/4 & 0\\
/// 0 & 1
/// \end{array}
/// \right]
/// \left[
/// \begin{array}{c}
/// u \\
/// v
/// \end{array}
/// \right],
/// \f]
/// with boundary conditions:
/// \f[
/// u'(\pm 1) \;=\; u(\pm 1) \;=\; v(\pm 1) \;=\; 0.
/// \f]
///

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

  LinopMat<Cd_t> Lmat(2,2), Mmat(2,2);
  Lmat << D4/16.0 + Vd_t(2.0*y), 1.0,//
              sin(2.0*y), D2/4.0;

  Mmat << D2/4.0 , 0.0,//
          0.0, 1.0;
  BcMat<Cd_t> bcs(6,2);
  bcs.L << D0, 0.0,//
          D1, 0.0, //
          D0, 0.0,//
          D1, 0.0,//
          0.0, D0,//
          0.0, D0;
  bcs.eval << -1.0, 0.0,//
              -1.0, 0.0,//
              1.0, 0.0,//
              1.0, 0.0,//
              0.0, 1.0,//
              0.0, -1.0;

  GeneralizedEigenSolver<Cd_t> eigs;

  eigs.compute(Lmat, Mmat, 6, bcs);
  std::cout << "Eigenvalues: \n" << '\n';
  std::cout << eigs.eigenvalues << '\n';
  return 0;
}
