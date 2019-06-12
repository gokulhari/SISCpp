/// \brief example with discriptor approach.
#define EIGEN_USE_BLAS
#define SIS_USE_LAPACK
#include <fstream>
#include <iostream>
#include <sis.hpp>
#include <string>

using namespace std;
typedef complex<double> Cd_t;
typedef valarray<double> Vd_t;

typedef valarray<complex<double> > Vcd_t;
complex<double> ii(0.0, 1.0);

int main() {
  using namespace sis;
  int bre;
  // Number of Chebyshev polynomials
  N = 127;
  sis_setup();
  valarray<double> y(N + 1), U(N + 1), Uy(N + 1), Uyy(N + 1), T11(N + 1),T12(N + 1);
  // Set in cheb-point
  setChebPts(y);
  // Number of Eigenvalues to compute:
  int num_vals = 15 * (N + 1);

  // Number of fourier modes in span-wise direction:
  int Nz = 1;
  int Nx = 1;

  valarray<double> kz(Nz);
  valarray<double> kx(Nx);

  string flowType("Poiseuille");
  // string flowType("Couette");

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

  double Re;
  double We;
  double beta;

  Linop<double> Dy(1);
  Dy.coef << 1.0, 0.0;

int k = 0;
      We = 3.96;
      Re = 3960.0;
      kx[k] = 1.15;
      beta = 0.5;
      kz[k] = 0.0;
      LinopMat<complex<double> > Lmat(10, 10), Mmat(10, 10);
      double k2 = kx[k] * kx[k] + kz[k] * kz[k];
      double k4 = k2 * k2;
      Linop<double> Delta(2), Delta2(4);

      T11 = 2.0 * We * Uy * Uy;
      T12 = Uy;

      Delta.coef << 1.0, 0.0, -k2;
      Delta2.coef << 1.0, 0.0, -2 * k2, 0.0, k4;

      Lmat(0, 0) = (-ii * kx[k] * U * Re) + beta * (Delta);
      Lmat(0, 1) = -Uy * Re;
      Lmat(0, 2) = 0.0;
      Lmat(0, 3) = -ii * kx[k];
      Lmat(0, 4) = (1.0 - beta) * ii * kx[k];
      Lmat(0, 5) = (1.0 - beta) * Dy;
      Lmat(0, 6) = (1.0 - beta) * ii * kz[k];
      Lmat(0, 7) = 0.0;
      Lmat(0, 8) = 0.0;
      Lmat(0, 9) = 0.0;

      Lmat(1, 0) = 0.0;
      Lmat(1, 1) = (-ii * kx[k] * U * Re) + beta * (Delta);
      Lmat(1, 2) = 0.0;
      Lmat(1, 3) = -(Dy);
      Lmat(1, 4) = 0.0;
      Lmat(1, 5) = (1.0 - beta) * ii * kx[k];
      Lmat(1, 6) = 0.0;
      Lmat(1, 7) = (1.0 - beta) * Dy;
      Lmat(1, 8) = (1.0 - beta) * ii * kz[k];
      Lmat(1, 9) = 0.0;

      Lmat(2, 0) = 0.0;
      Lmat(2, 1) = 0.0;
      Lmat(2, 2) = (-ii * kx[k] * U * Re) + beta * (Delta);
      Lmat(2, 3) = -ii * kz[k];
      Lmat(2, 4) = 0.0;
      Lmat(2, 5) = 0.0;
      Lmat(2, 6) = (1.0 - beta) * ii * kx[k];
      Lmat(2, 7) = 0.0;
      Lmat(2, 8) = (1.0 - beta) * Dy;
      Lmat(2, 9) = (1.0 - beta) * ii * kz[k];

      Lmat(3, 0) = ii*kx[k];
      Lmat(3, 1) = Dy;
      Lmat(3, 2) = ii* kz[k];
      Lmat(3, 3) = 0.0;
      Lmat(3, 4) = 0.0;
      Lmat(3, 5) = 0.0;
      Lmat(3, 6) = 0.0;
      Lmat(3, 7) = 0.0;
      Lmat(3, 8) = 0.0;
      Lmat(3, 9) = 0.0;

      Lmat(4, 0) = (2.0 * ii * kx[k] / We) +
                   (Vcd_t(2.0 * ii * kx[k] * T11) + 2.0 * (Vd_t(T12) * Dy));
      Lmat(4, 1) = -4.0 * We * Uy * Uyy;
      Lmat(4, 2) = 0.0;
      Lmat(4, 3) = 0.0;
      Lmat(4, 4) = ((-1.0 / We) - (ii * kx[k] * U));
      Lmat(4, 5) = 2.0 * Uy;
      Lmat(4, 6) = 0.0;
      Lmat(4, 7) = 0.0;
      Lmat(4, 8) = 0.0;
      Lmat(4, 9) = 0.0;

      Lmat(5, 0) = Vcd_t(ii * kx[k] * Uy) + (Dy / We);
      Lmat(5, 1) = (ii * kx[k] / We) + (Vcd_t(2.0 * ii * kx[k] * We * Uy * Uy) +
                                     (Vd_t(Uy) * Dy) - Vd_t(Uyy));
      Lmat(5, 2) = 0.0;
      Lmat(5, 3) = 0.0;
      Lmat(5, 4) = 0.0;
      Lmat(5, 5) = ((-1.0 / We) - (ii * kx[k] * U));
      Lmat(5, 6) = 0.0;
      Lmat(5, 7) = 1.0 * Uy;
      Lmat(5, 8) = 0.0;
      Lmat(5, 9) = 0.0;

      Lmat(6, 0) = ii * kx[k] / We;
      Lmat(6, 1) = 0.0;
      Lmat(6, 2) = (ii * kx[k] / We) + (Vcd_t(2.0 * ii * kx[k] * We * Uy * Uy) +
                                     (Vd_t(Uy) * Dy));
      Lmat(6, 3) = 0.0;
      Lmat(6, 4) = 0.0;
      Lmat(6, 5) = 0.0;
      Lmat(6, 6) = ((-1.0 / We) - (ii * kx[k] * U));
      Lmat(6, 7) = 0.0;
      Lmat(6, 8) = 1.0 * Uy;
      Lmat(6, 9) = 0.0;


      Lmat(7, 0) = 0.0;
      Lmat(7, 1) = 2.0 * ii * kx[k] * Uy + (2.0 * Dy / We);
      Lmat(7, 2) = 0.0;
      Lmat(7, 3) = 0.0;
      Lmat(7, 4) = 0.0;
      Lmat(7, 5) = 0.0;
      Lmat(7, 6) = 0.0;
      Lmat(7, 7) = ((-1.0 / We) - (ii * kx[k] * U));
      Lmat(7, 8) = 0.0;
      Lmat(7, 9) = 0.0;

      Lmat(8, 0) = 0.0;
      Lmat(8, 1) = ii * kx[k] / We;
      Lmat(8, 2) = ii * kx[k] * Uy + (Dy / We);
      Lmat(8, 3) = 0.0;
      Lmat(8, 4) = 0.0;
      Lmat(8, 5) = 0.0;
      Lmat(8, 6) = 0.0;
      Lmat(8, 7) = 0.0;
      Lmat(8, 8) = ((-1.0 / We) - (ii * kx[k] * U));
      Lmat(8, 9) = 0.0;

      Lmat(9, 0) = 0.0;
      Lmat(9, 1) = 0.0;
      Lmat(9, 2) = 2.0 * ii * kx[k] / We;
      Lmat(9, 3) = 0.0;
      Lmat(9, 4) = 0.0;
      Lmat(9, 5) = 0.0;
      Lmat(9, 6) = 0.0;
      Lmat(9, 7) = 0.0;
      Lmat(9, 8) = 0.0;
      Lmat(9, 9) = ((-1.0 / We) - (ii * kx[k] * U));


      BcMat<complex<double> > Lbc(5, 10), Rbc(5, 10),bc(10,10);
      Lbc.L << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0,//
              0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0,//
              0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0,//
              0.0, Dy, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0,//
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 , 0.0, 0.0;
      Rbc.L << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0,//
              0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0,//
              0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0,//
              0.0, Dy, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0,//
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 , 0.0, 0.0;

      Lbc.eval.setConstant(-1.0);
      Rbc.eval.setConstant(1.0);

      bc.L << Lbc.L,//
              Rbc.L;
      bc.eval << Lbc.eval,//
                Rbc.eval;
      Mmat << Re, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0,//
              0.0, Re, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0,     //
              0.0, 0.0, Re, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0,    //
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0,    //
              0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 , 0.0, 0.0,    //
              0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 , 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 , 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 , 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 1.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 1.0;

      GeneralizedEigenSolver<complex<double> > eigs;
      eigs.compute(Lmat, Mmat, (N + 1) * Lmat.r, bc);
      //eigs.removeInf();
      //eigs.sortByLargestReal();
      num_vals = eigs.eigenvalues.size();
      ofstream outf;
      outf.open(string("data/OldroydB_") + flowType + string("_beta_") +
                int2str(int(beta * 1000)) + string("_We_") + int2str(int(We)) +
                string("_Re_") + int2str(int(Re)) + string("_N_") + int2str(N) +
                string("_") + int2str(k) +
                string(".txt"));
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> out_to_file(num_vals,
                                                                        3);
      out_to_file.setConstant(0.0);
      out_to_file.col(0) = eigs.eigenvalues.real();
      out_to_file.col(1) = eigs.eigenvalues.imag();
      //out_to_file.col(2) = eigs.MPorNot.cast<double>();
      outf << out_to_file;
      outf.close();

  return 0;
}
