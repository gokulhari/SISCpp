/// \file Ex_05.cpp
/// \brief Eigenvalues of the diffusion operator.
///
/// Finds the eigenvalues the system:
/// \f{align}{
/// \partial_{yy}u(y)\,-\,4\,u(y) \;=\; 0,
/// \f}
/// with Dirichlet boundary conditions.
///
/// Analytical solution: \f$-4 - n^2\pi^2/4\f$

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

  LinopMat<complex<double> > Lmat(1, 1);
  Linop<double> D2(2), D1(1), D0(0); // Linops with highest order 2, 1 and 0.

  D2.coef << 1.0, 0.0, 0.0; // for (1.0 * D2 + 0.0 * D1 + (0.0) )v
  D1.coef << 1.0, 0.0;
  D0.coef << 1.0;

  // for operator D2v - 4 v
  Lmat << D2 - 4.0 * D0;

  // Specify boundary conditions via BcMat:
  BcMat<complex<double> > bcs(
      2, 1); // Two boundary conditions in total, on one variable

  bcs.L << D0,//
          D0;
  bcs.eval << -1.0, // evaluate bcs.L right end
      1.0;          // evaluate bcs.L at 0.5, must be a value in domain [-1, 1]

  GeneralizedEigenSolver<Cd_t> eigs;

  // Compute first 6 eigenvalues
  eigs.compute(Lmat, 6, bcs);

  std::cout << "Eigenvalues: \n" << eigs.eigenvalues << '\n';

  return 0;
}
