/// \file Ex_16.cpp
/// \brief Example with discriptor approach, needs an input argument N.

#define EIGEN_USE_BLAS
#define SIS_USE_LAPACK
#define EIGEN_FAST_MATH 0
#include <fstream>
#include <iostream>
#include <sis.hpp>
using namespace std;
typedef complex<double> Cd_t;
typedef valarray<double> Vd_t;
typedef valarray<Cd_t> Vcd_t;
int main(int argc, char** argv) {
  using namespace sis;
  int bre;
  // Number of Chebyshev polynomials
  sis::N = 127;
  sscanf(argv[1], "%d", &N);
  std::cout << "N = " << N << '\n';
  sis_setup();
  valarray<double> y, U(N + 1), Uy(N + 1), Uyy(N + 1);
  // Set in cheb-point
  setChebPts(y);
  // Number of Eigenvalues to compute:
  int num_vals = 15 * (N + 1);

  /*double Re = 2310.0;
  double We = 0.001 * Re;
  double kx = 1.31;
  */
 
  double Re = 0.0;
  double We = 2;
  double kx = 1.0;
  double beta = 0.0;
  //string flowType("Poiseuille");
  string flowType("Couette");

  if (flowType.compare("Poiseuille") == 0) {
    U = 1.0 - y * y;
    Uy = -2.0 * y;
    Uyy = -2.0;
  } else if (flowType.compare("Couette") == 0) {
    U = y;
    Uy = 1.0;
    Uyy = 0.0;
  } else {
    std::cout << "Unknown flow type, in line " << __LINE__ << '\n'
              << "Exiting...\n";
    exit(1);
  }
  complex<double> ii(0.0, 1.0);
  Linop<double> Dyy, Dy, Delta, Delta2;
  LinopMat<complex<double> > Lmat(6, 6), Mmat(6, 6);

  Dyy.set(2);
  Dyy.coef << 1.0, 0.0, 0.0;

  Dy.set(1);
  Dy.coef << 1.0, 0.0;

  Delta.set(2);
  Delta.coef << 1.0, 0.0, -kx * kx;

  Delta2.n = 4;
  Delta2.set();
  Delta2.coef << 1.0, 0.0, -2 * kx * kx, 0.0, pow(kx, 4.0);

  Lmat(0, 0) = (-ii * kx * U * Re) + beta * (Delta);
  Lmat(0, 1) = -Uy * Re;
  Lmat(0, 2) = -ii * kx;
  Lmat(0, 3) = (1.0 - beta) * ii * kx;
  Lmat(0, 4) = (1.0 - beta) * Dy;
  Lmat(0, 5) = 0.0;

  Lmat(1, 0) = 0.0;
  Lmat(1, 1) = (-ii * kx * U * Re) + beta * (Delta);
  Lmat(1, 2) = -(Dy);
  Lmat(1, 3) = 0.0;
  Lmat(1, 4) = (1.0 - beta) * ii * kx;
  Lmat(1, 5) = (1.0 - beta) * Dy;

  Lmat(2, 0) = (2.0 * ii * kx / We) +
               (Vcd_t(4.0 * ii * kx * We * Uy * Uy) + 2.0 * (Vd_t(Uy) * Dy));
  Lmat(2, 1) = -4.0 * We * Uy * Uyy;
  Lmat(2, 2) = 0.0;
  Lmat(2, 3) = ((-1.0 / We) - (ii * kx * U));
  Lmat(2, 4) = 2.0 * Uy;
  Lmat(2, 5) = 0.0;

  Lmat(3, 0) = Vcd_t(ii * kx * Uy) + (Dy / We);
  Lmat(3, 1) = (ii * kx / We) + (Vcd_t(2.0 * ii * kx * We * Uy * Uy) +
                                 (Vd_t(Uy) * Dy) - Vd_t(Uyy));
  Lmat(3, 2) = 0.0;
  Lmat(3, 3) = 0.0;
  Lmat(3, 4) = ((-1.0 / We) - (ii * kx * U));
  Lmat(3, 5) = 1.0 * Uy;

  Lmat(4, 0) = 0.0;
  Lmat(4, 1) = 2.0 * ii * kx * Uy + (2.0 * Dy / We);
  Lmat(4, 2) = 0.0;
  Lmat(4, 3) = 0.0;
  Lmat(4, 4) = 0.0;
  Lmat(4, 5) = ((-1.0 / We) - (ii * kx * U));

  Lmat(5, 0) = ii * kx;
  Lmat(5, 1) = Dy;
  Lmat(5, 2) = 0.0;
  Lmat(5, 3) = 0.0;
  Lmat(5, 4) = 0.0;
  Lmat(5, 5) = 0.0;

  Mmat << Re, 0.0, 0.0, 0.0, 0.0, 0.0, //
      0.0, Re, 0.0, 0.0, 0.0, 0.0,     //
      0.0, 0.0, 0.0, 1.0, 0.0, 0.0,    //
      0.0, 0.0, 0.0, 0.0, 1.0, 0.0,    //
      0.0, 0.0, 0.0, 0.0, 0.0, 1.0,    //
      0.0, 0.0, 0.0*Dyy, 0.0, 0.0, 0.0;

  BcMat<Cd_t> bcs(8,6);
  bcs.L << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
      0.0, 1.0, 0.0, 0.0, 0.0, 0.0,      //
      0.0, Dy, 0.0, 0.0, 0.0, 0.0,       //
      0.0, 0.0, 0.0, 0.0, 0.0, 1.0,     //
      1.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
      0.0, 1.0, 0.0, 0.0, 0.0, 0.0,      //
      0.0, Dy, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 1.0;       //

  bcs.eval << -1.0, 0.0,0.0,0.0,0.0,0.0,//
              0.0, -1.0, 0.0, 0.0, 0.0, 0.0,//
              0.0, -1.0, 0.0, 0.0, 0.0, 0.0,//
              0.0, 0.0, -1.0, 0.0, 0.0,1.0,//
              1.0, 0.0, 0.0, 0.0,0.0,0.0,//
              0.0, 1.0, 0.0, 0.0, 0.0, 0.0,//
              0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, -1.0, 0.0, 0.0,-1.0;

  GeneralizedEigenSolver<complex<double> > eigs;

  eigs.computeAppend(Lmat,Mmat, 14*(N + 1), bcs);
  //eigs.removeInf();
  //eigs.sortByLargestReal();

  cout << "Eigenvalues: \n" << eigs.eigenvalues << "\n";
  ofstream outf;
  outf.open(string(string(string("data/Ex_16_Re") + int2str(int(Re)) +
            string("_We") + int2str(int(We)) + string("_beta") +
            int2str(int(beta * 100)) + string("_N") + int2str(N) + flowType +
         //   string("_") + int2str(i) +
            string(".txt"))).c_str());
  Eigen::MatrixXd temp(eigs.eigenvalues.rows(),2);
  temp << eigs.eigenvalues.real(),eigs.eigenvalues.imag();
  outf << temp;
  outf.close();

  //Eigen::MatrixXd evec(N+1,3);
  //evec << yEigen, eigs.eigenvectors(3,23).evr(), eigs.eigenvectors(3,23).evi();
  //outf.open("data/Ex_14_evec.txt");
  //outf << evec;
  //outf.close();
  return 0;
}
