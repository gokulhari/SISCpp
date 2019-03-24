/// \file Ex_04.cpp
/// \brief This solves the problem for a block matrix operator system
///
/// The eigenvalue problem is given by:
/// \f[
/// \left[
/// \begin{array}{cc}
/// \partial_{yy}/4 + 2y & 1.0 \\
/// \sin (2y) & \partial_{yy}/4
/// \end{array}
/// \right]
/// \left[
/// \begin{array}{c}
/// u \\
/// v
/// \end{array}
/// \right] \;=\;
/// \lambda \left[
/// \begin{array}{c}
/// u \\
/// v
/// \end{array}
/// \right],
/// \f]
/// with boundary conditions:
/// \f[
/// u(\pm 1) \;=\; v(\pm 1) \;=\; 0.
/// \f]
///

#include <fstream>
#include <iostream>
#include <sis.hpp>

using namespace std;
typedef complex<double> Cd_t;
typedef valarray<double> Vd_t;


int main() {
  using namespace sis;
  N = 63;
  sis_setup();
  ChebfunMat<Cd_t> forc(1,1); // 1 by 1 ChebfunMat
  Linop<double> D2(2), D1(1), D0(0); // Linops with highest order 2, 1 and 0.
  Cd_t ii(0.0, 1.0);

  D2.coef << 1.0, 0.0, 0.0; // for (1.0 * D2 + 0.0 * D1 + (0.0) )v
  D1.coef << 1.0, 0.0;
  D0.coef << 1.0;

  LinopMat<Cd_t> Lmat(2,2);
  Lmat << D2/4.0 + Vd_t(2*y), 1.0,//
          sin(2*y), D2/4.0;

  // Specify boundary conditions using BcMat, 4 bcs on 2 variables:
  BcMat<Cd_t> bcs(4,2);
  bcs.L << D0, 0.0,//
           0.0, D0,//
           D0, 0.0,//
           0.0, D0;
  bcs.eval << -1.0, 0.0,//
               0.0, -1.0,//
               1.0, 0.0, //
               0.0, 1.0;

  GeneralizedEigenSolver<Cd_t> eigs;

  eigs.compute(Lmat, 6, bcs);
  std::cout << "Eigenvalues: \n" << '\n';
  std::cout << eigs.eigenvalues << '\n';
  return 0;
}
