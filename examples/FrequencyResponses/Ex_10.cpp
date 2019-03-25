/// \file Ex_10.cpp
/// \brief Solving for frequency responses of diffusion operator


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

  std::cout << "Singular Values: \n" << svds.eigenvalues << '\n';

  return 0;
}
