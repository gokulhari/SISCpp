/// \file Ex_02.cpp
/// \brief An example with a fourth order system.
///
/// Solves a fourth order system:
///\f[
/// \partial_{yyyy}u(y) \,-\, 4\,\partial_{yy}u(y) \,+\, 8\,u(y)
/// \;= \; 1 \,-\, y^2 \,+\, y^3,
/// \f]
/// with clamped boundary conditions:
/// \f[
/// \begin{array}{cll}
///  \partial_yu(1) \;&=\; u(1) \;&=\; 0 \\
///  \partial_yu(-1) \;&=\; u(-1) \;&=\; 0.
/// \end{array}
/// \f]
///
/// \image html pics/Ex_02.svg <!--
/// --> width=500cm height=500cm

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

  LinopMat<Cd_t> Lmat(1,1);
  Lmat << D4 - (4.0 * D2) + (8.0 * D0);

  BcMat<Cd_t> bcs(4,1);
  bcs.L << D0,//
          D1,//
          D0,//
          D1;
  bcs.eval << -1.0,//
              -1.0,//
              1.0,//
              1.0;
  bcs.vals << 0.0,//
              0.0,//
              0.0,//
              0.0;
  
 
  // Set the forcing:
  forc << (Cd_t(1.0,0.0) - (yc * yc) + (yc * yc * yc));
  

  // Solve and replace forcing with the boundary conditions:
  forc = linSolve(Lmat, bcs, forc);


  // Check boundaries
  cout << "lbc: " << forc(-1.0) << endl;
  cout << "rbc: " << forc(1.0) << endl;

  // Write to file to plot, need only the real part:
  ofstream outf;
  outf.open("data/Ex_02.txt");
  Eigen::MatrixXd temp(N + 1, 2);
  temp << yEigen, forc(0,0).evr();
  outf << temp;
  // note that sis::yEigen is the same as y, but of type Eigen Matrix,
  outf.close();
  return 0;
}
