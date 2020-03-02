/// \file Ex_03.cpp
/// \brief An example with nonconstant coefficients and comparissons with
/// analytical solution.
///
/// Solving the system: \f[ \partial_y u(y) + \frac{1}{y^2
/// +1}u(y) \;=\; 0, \f] with boundary condition: \f[ u(-1) = 1.0. \f] Anaytical
/// solution is given by: \f[ u(y) = \exp\left(-\mathrm{tan}^{-1}(y)\,-\,
/// \mathrm{tan}^{-1}(1) \right) \f]
///
/// \image html pics/Ex_03.svg "Comparing analytical and solutions from SIS"
/// <!--
/// --> width=500cm height=500cm
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

  D2.coef << 1.0, 0.0, 0.0; // for (1.0 * D2 + 0.0 * D1 + (0.0) )v
  D1.coef << 1.0, 0.0;
  D0.coef << 1.0;

  valarray<double> anal_sol;
  Eigen::VectorXd analSol(N+1);
  anal_sol = exp(-(atan(y) + atan(1.0)));
  cout << "Anal sol: \n";
  for (int i = 0; i < N + 1; i++) {
    cout << anal_sol[i] << "\n";
    analSol[i] = anal_sol[i];
  }

  LinopMat<Cd_t> Lmat(1,1);

  Lmat << D1 + Vd_t(1.0 / (pow(y, 2.0) + 1.0)) * D0;

  BcMat<Cd_t> bcs(1,1);

  bcs.L << D0;
  bcs.eval << -1.0;
  bcs.vals << 1.0;

  forc << Cd_t(0.0,0.0);
 

  // Replace forcing with the solution:
  forc = linSolve(Lmat,bcs,forc);

  // Check boundaries:
  std::cout << "lbc: \n" << forc(-1.0) << '\n';
  std::cout << "rbc: \n" << forc(1.0) << '\n';
  std::cout << "error from analytical solution ";
  double err = abs(real(forc(0,0).v) - anal_sol).sum();
  cout << err << endl;

  // Print to file:
  Eigen::MatrixXd outMat(N+1,3);
  outMat << yEigen,forc(0,0).evr(), analSol; 
ofstream outfile;
outfile.open("data/Ex_03.txt");
outfile << outMat;
outfile.close();
  return 0;
}
