/// \file Ex_10.cpp
/// \brief Solving for frequency responses of the reaction-diffusion operator,
/// \f[
/// \partial_t u(y,t) \;=\; \partial_{yy}u(y,t) -  \epsilon^2 u(y,t) + d(y) \mathrm{e},^{\! j\omega t}
///\f]
/// with homogeneous Neumann boundary conditions, \f$[\partial_y u (\cdot,t)](y = \pm 1) = 0\f$
///
/// As this system is a linear differential equation, its solution is of the form \f$u(y,t) = u(y)\,\mathrm{e},^{\! j\omega t}\f$ using this on the governing equation we arrive at
/// \f[
///  ( j\omega - \mathrm{D}^2  + \epsilon^2 \mathrm{I})\, u(y) \;=\; d(y),
/// \f]
/// where \f$\mathrm{D} = \mathrm d/\mathrm dy\f$ and \f$\mathrm{I}\f$ is the identity operator. The singular value decomposition (SVD) of the operator \f$ ( j\omega - \mathrm{D}^2  + \epsilon^2 \mathrm{I})^{ -1}\f$ finds the body force of unit \f$L^2[-1~\, 1]\f$ norm that induces the largest amplification of \f$u(y)\f$, which is the principal singular value of \f$ ( j\omega - \mathrm{D}^2  + \epsilon^2 \mathrm{I}).^{\!\! -1}\f$
/// 
/// Singular value decomposition of operators are computed using the class SingularValueDecomposition as illustrated in this example. This class supports calculation of the singular values, power spectral density (a.k.a sum of the squares of singular values of finite-trace operators), and the H-infinity norm. See our related paper for more information.
#define SIS_USE_LAPACK
#include <fstream> // To output data to files
#include <iostream>
#include <sis.hpp>

using namespace std;
typedef complex<double> Cd_t;

int main() {
  using namespace sis;
  N = 63; // sis::N, Defined in sis.hpp. This specifies no. of Chebyshev
          // coefficients
  sis_setup();

  LinopMat<Cd_t> Lmat(1, 1), Mmat(1,1);
  Linop<double> D2(2), D1(1), D0(0); // Linops with highest order 2, 1 and 0.
  Cd_t ii(0.0,1.0);

  D2.coef << 1.0, 0.0, 0.0; // for (1.0 * D2 + 0.0 * D1 + (0.0) )v
  D1.coef << 1.0, 0.0;
  D0.coef << 1.0;




  // for operator D2v - 4 v
  Lmat << D2 - 0.0001 * D0;
  Mmat << 1.0;

  double omega = 0.0;
  // Specify boundary conditions via BcMat, lbc and rbc separately:
  BcMat<Cd_t> lbcs(1, 1), rbcs(1,1);

  lbcs.L << D1;
  rbcs.L << D1;
  lbcs.eval << -1.0, // evaluate bcs.L right end
  rbcs.eval <<  1.0;          // evaluate bcs.L at 0.5, must be a value in domain [-1, 1]

  SingularValueDecomposition<Cd_t> svds;

  LinopMat<Cd_t> A(1,1), B(1,1), C(1,1);
  A = (ii*omega*Mmat) - Lmat;
  B(0,0) = 1.0;
  C(0,0) = 1.0;
  // Compute first 12 singular values
  svds.compute(A, B, C, lbcs, rbcs, 12);

  cout << "Singular Values: \n" << svds.eigenvalues << '\n';
  cout << "Power Spectral Density: \n" << svds.PowerSpectralDensity(A,B,C,lbcs,rbcs) << '\n';
  // cout << "H infinity norm: \n" << svds.HinfNorm() << "\n";
  return 0;
}
