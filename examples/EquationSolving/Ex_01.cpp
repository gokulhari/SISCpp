/// \file Ex_01.cpp
/// \brief This simulates a simple 2nd order ODE with mixed boundary
/// conditions.
///
/// Solves:
/// \f[
/// \partial_{yy}u(y)\,-\,4\, u(y)\;=\; 1 \; - \; y^2,
/// \f]
/// with boundary conditions:
/// \f{align*}{
/// u(-1)  \,+\,  4\,u'(-1)\;&=\; 2,\\
/// u(1)  \,-\,  5\,u'(1)\;&=\; 3.
/// \f}
/// Solving by entering the forcing in physical-space, so solution
/// is returned in physical space.
///
/// \image html pics/Ex_01.svg <!--
/// --> width=500cm height=500cm

#include <fstream> // To output data to files
#include <iostream>
#include <sis.hpp>

using namespace std;

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

  // Some mixed boundary condition:
  // bcs.L << D0 + 4.0*D1,
  //         D0 - 5.0*D1;
  bcs.L << D0 + 4.0*D1,//
          D0 - 5.0*D1;
  bcs.eval << -1.0, // evaluate bcs.L right end
              0.5;          // evaluate bcs.L at 0.5, must be a value in domain [-1, 1]
  bcs.vals << 2.0,  // value at right
      3.0;          // value at left

  // Construct forcing, forc = 1 - y^2
  ChebfunMat<complex<double> > forc(1, 1); // 1 by 1 ChebfunMat
  forc << 1.0 - y * y;                     // y is defined in sis::y

  // Solve, by replacing the forcing as the solution.
  forc = linSolve(Lmat, bcs, forc);
  cout << "lbc: " << forc(1.0) << endl;
  cout << "rbc: " << forc(-1.0) << endl;

  // Write to file to plot, need only the real part:
  ofstream outf;
  outf.open("data/Ex_01.txt");
  Eigen::MatrixXd temp(N + 1, 2);
  temp << yEigen, forc(0,0).evr();
  outf << temp;
  // note that sis::yEigen is the same as y, but of type Eigen Matrix,
  outf.close();


  return 0;
}
