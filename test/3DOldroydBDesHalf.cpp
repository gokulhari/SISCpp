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
  valarray<Cd_t> U(N + 1), Uy(N + 1), Uyy(N + 1), T11(N + 1),T12(N + 1);

  // Number of Eigenvalues to compute:
  int num_vals = 15 * (N + 1);

  // Number of fourier modes in span-wise direction:
  int Nz = 1;
  int Nx = 1;

  Cd_t kz;
  Cd_t kx;

  //string flowType("Poiseuille");
   string flowType("Couette");

  if (flowType.compare("Poiseuille") == 0) {
    U = 1.0 - yc * yc;
    Uy = -2.0 * yc;
    Uyy = -2.0;
  } else if (flowType.compare("Couette") == 0) {
    U = yc;
    Uy = 1.0;
    Uyy = 0.0;
  } else {
    std::cout << "Unknown flow type, in line " << __LINE__ << '\n'
              << "Exiting...\n";
    exit(1);
  }

  Cd_t R;
  Cd_t We;
  Cd_t beta;

  Linop<double> Dy(1);
  Dy.coef << 1.0, 0.0;

int k = 0;
      We = 1.0;
      R = 0.0;
      kx = 1.0;
      beta = 0.5;
      kz = 1.0;
      LinopMat<complex<double> > Lmat(4, 4), Mmat(4, 4);
      Cd_t k2 = kx * kx + kz * kz;
      Cd_t k4 = k2 * k2;
      Linop<Cd_t> Delta(2), Delta2(4);
      Linop<Cd_t> D4(4), D3(3), D2(2), D1(1);
      Vcd_t c(N + 1), cy(N + 1), cyy(N + 1), a0(N + 1), a1(N + 1), a2(N + 1),
          a3(N + 1), a4(N + 1);

      D4.coef << 1.0, 0.0, 0.0, 0.0, 0.0;
      D3.coef << 1.0, 0.0, 0.0, 0.0;
      D2.coef << 1.0, 0.0, 0.0;
      D1.coef << 1.0, 0.0;
      // Set in cheb-point
      Cd_t omega = 0.0;
      c = (ii * omega + (1.0/We) + (ii * kx * U));
      cy = ii * kx * Uy;
      cyy = ii * kx * Uyy;
      T11 = 2.0 * We * Uy * Uy;
      T12 = Uy;

      Delta.coef << 1.0, 0.0, -k2;
      Delta2.coef << 1.0, 0.0, -2.0 * k2, 0.0, k4;

      Mmat << 1.0, 0.0, 0.0, 0.0, //
          0.0, 1.0, 0.0, 0.0,     //
          0.0, 0.0, 1.0, 0.0,     //
          0.0, 0.0, 0.0, 0.0 * Delta;

      //
      Lmat(0, 0) = Vcd_t(beta + (1.0 - beta) / (c * We)) * D2 + Vcd_t((-1.0 + beta) * (cy - Cd_t(0.0, 1.0) * kx * (Uy * (2.0 + c * We) + 2.0 * c * We * T12)) / (pow(c, 2.0) * We)) * D1 + Vcd_t(-(beta * k2) - Cd_t(0.0, 1.0) * omega * R - Cd_t(0.0, 1.0) * kx * R * U + ((-1.0 + beta) * kx * (Uy * (Cd_t(0.0, 1.0) * cy + 2.0 * kx * Uy) * We + c * (2.0 * kx + kz - Cd_t(0.0, 1.0) * Uyy * We + 2.0 * kx * We * T11)) / (pow(c, 2.0) * We)));

      Lmat(0, 1) = Vcd_t(-(((-1.0 + beta) * Uy * (2.0 + c * We)) / (pow(c, 2.0) * We))) * D2 + Vcd_t((Cd_t(0.0, -1.0) * (-1.0 + beta) * (4.0 * Uy * (Cd_t(0.0, 1.0) * cy + kx * Uy) + c * (Cd_t(0.0, -2.0) * Uyy + Uy * (Cd_t(0.0, 1.0) * cy + 4.0 * kx * Uy) * We) + pow(c, 2.0) * (kx + 2.0 * kx * pow(Uy, 2.0) * pow(We, 2.0)))) / (pow(c, 3.0) * We)) * D1 + Vcd_t(-(R * Uy) + ((-1.0 + beta) * (4.0 * kx * pow(Uy, 2.0) * (Cd_t(0.0, 1.0) * cy + kx * Uy) * We - Cd_t(0.0, 4.0) * pow(c, 2.0) * kx * Uy * Uyy * pow(We, 2.0) + c * (Cd_t(0.0, 1.0) * cy * (kx + Cd_t(0.0, 1.0) * Uyy * We + 2.0 * kx * pow(Uy, 2.0) * pow(We, 2.0)) + kx * Uy * (2.0 * kx + kz - Cd_t(0.0, 2.0) * Uyy * We + 4.0 * kx * pow(Uy, 2.0) * pow(We, 2.0))))) / (pow(c, 3.0) * We));

      Lmat(0, 2) = Vcd_t((Cd_t(0.0, -1.0) * (-1.0 + beta) * kz * Uy * (1.0 + c * We)) / (pow(c, 2.0) * We)) * D1 + Vcd_t(((-1.0 + beta) * kx * kz * (c + pow(Uy, 2.0) * We * (1.0 + 2.0 * c * We))) / (pow(c, 2.0) * We));

      Lmat(0,3) = -ii*kx;

      Lmat(1, 0) = Vcd_t((Cd_t(0.0, -1.0) * (-1.0 + beta) * kx) / (c * We)) * D1 + Vcd_t(((-1.0 + beta) * pow(kx, 2.0) * Uy) / c);

      Lmat(1, 1) = Vcd_t(beta + (2.0 - 2.0 * beta) / (c * We)) * D2 + Vcd_t(((-1.0 + beta) * (2.0 * cy - Cd_t(0.0, 1.0) * kx * Uy * (2.0 + 3.0 * c * We))) / (pow(c, 2.0) * We)) * D1 + Vcd_t(-(beta * k2) - Cd_t(0.0, 1.0) * omega * R - Cd_t(0.0, 1.0) * kx * R * U + ((-1.0 + beta) * kx * (2.0 * Uy * (Cd_t(0.0, 1.0) * cy + kx * Uy) * We + c * (kx + kz - Cd_t(0.0, 1.0) * Uyy * We + 2.0 * kx * pow(Uy, 2.0) * pow(We, 2.0)))) / (pow(c, 2.0) * We));

      Lmat(1, 2) = Vcd_t((Cd_t(0.0, -1.0) * (-1.0 + beta) * kz) / (c * We)) * D1 + Vcd_t(((-1.0 + beta) * kx * kz * Uy) / c);

      Lmat(1,3) = -D1;

      Lmat(2, 0) = ((-1.0 + beta) * pow(kx, 2.0)) / (c * We);

      Lmat(2, 1) = Vcd_t((Cd_t(0.0, -1.0) * (-1.0 + beta) * kx) / (c * We)) * D1 + Vcd_t(((-1.0 + beta) * kx * (Cd_t(0.0, 1.0) * cy + kx * Uy)) / (pow(c, 2.0) * We));

      Lmat(2, 2) = Vcd_t(beta + (1.0 - beta) / (c * We)) * D2 + Vcd_t(((-1.0 + beta) * (cy - Cd_t(0.0, 1.0) * kx * Uy * (1.0 + 2.0 * c * We))) / (pow(c, 2.0) * We)) * D1 + Vcd_t(-(beta * k2) - Cd_t(0.0, 1.0) * omega * R - Cd_t(0.0, 1.0) * kx * R * U + ((-1.0 + beta) * kx * (Uy * (Cd_t(0.0, 1.0) * cy + kx * Uy) * We + c * (kx + 2.0 * kz - Cd_t(0.0, 1.0) * Uyy * We + 2.0 * kx * pow(Uy, 2.0) * pow(We, 2.0)))) / (pow(c, 2.0) * We));

      Lmat(2,3) = -ii*kz;

      Lmat(3,0) = ii*kx;

      Lmat(3,1) = D1;

      Lmat(3,2) = ii*kz;

      Lmat(3,3) = 0.0;

      LinopMat<std::complex<double> > A, B, C;
      BcMat<std::complex<double> > Lbc(4, 4), Rbc(4, 4), bc(8, 4);

//      B.resize(4, 3);
//      C.resize(2, 4);
//      B << 1.0, 0.0, 0.0, //
//          0.0, 1.0, 0.0,  //
//          0.0, 0.0, 0.0,  //
//          0.0, 0.0, 0.0;
//      C << 1.0, 0.0, 0.0, 0.0, //
//          0.0, 1.0, 0.0, 0.0;  //
//          //0.0, 0.0, 1.0, 0.0;  //

      B.resize(4, 3);
      C.resize(3, 4);
      B << 1.0, 0.0, 0.0, //
          0.0, 1.0, 0.0,  //
          0.0, 0.0, 1.0,  //
          0.0, 0.0, 0.0;
      C << 1.0, 0.0, 0.0, 0.0, //
          0.0, 1.0, 0.0, 0.0,  //
          0.0, 0.0, 1.0, 0.0;  //
      
      Lbc.resize(4, 4);
      Rbc.resize(4, 4);
      Lbc.L << 1.0, 0.0, 0.0, 0.0, //
          0.0, 1.0, 0.0, 0.0,      //
          0.0, 0.0, 1.0, 0.0,      //
          0.0, Dy, 0.0, 0.0;
      Rbc.L << 1.0, 0.0, 0.0, 0.0, //
          0.0, 1.0, 0.0, 0.0,      //
          0.0, 0.0, 1.0, 0.0,      //
          0.0, Dy, 0.0, 0.0;

      SingularValueDecomposition<std::complex<double> > svd;
      Lbc.eval.setConstant(-1.0);
      Rbc.eval.setConstant(1.0);

      A.resize(4, 4);
      A = ((ii * omega * Mmat) - Lmat);
      svd.compute(A, B, C, Lbc, Rbc, Lbc, Rbc, 15 * (N + 1));
      std::cout << "eigenvalue: " << svd.eigenvalues[0] << "\n";
      std::cout << "eigenvalue: " << svd.eigenvalues[1] << "\n";
      std::cout << "eigenvalue: " << svd.eigenvalues[2] << "\n";

      return 0;
}
