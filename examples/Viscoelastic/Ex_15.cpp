#define EIGEN_DONT_PARALLELIZE
#define SIS_USE_LAPACK
#define SIS_DONT_SORT
#define EIGEN_FAST_MATH 0
#include <fstream>
#include <iostream>
#include <omp.h>
#include <sis.hpp>
#include <string>

using namespace std;
typedef complex<double> Cd_t;
typedef valarray<double> Vd_t;
typedef valarray<complex<double> > Vcd_t;
int main() {
  using namespace sis;
  int bre;
  valarray<double> Wes(20);
  // Number of Chebyshev polynomials
  for (int i = 0; i < 20; i++) {
    Wes[i] = double(i + 1);
  }

  N = 479;
  sis_setup();
  Vcd_t U(N + 1), Uy(N + 1), Uyy(N + 1);
  string flowType("Poiseuille");
  // string flowType("Couette");

  if (flowType.compare("Poiseuille") == 0) {
    U = Cd_t(1.0, 0.0) - yc * yc;
    Uy = Cd_t(-2.0, 0.0) * yc;
    Uyy = Cd_t(-2.0, 0.0);
  } else if (flowType.compare("Couette") == 0) {
    U = yc;
    Uy = Cd_t(1.0, 0.0);
    Uyy = Cd_t(0.0, 0.0);
  } else {
    std::cout << "Unknown flow type, in line " << __LINE__ << '\n'
              << "Exiting...\n";
    exit(1);
  }

#pragma omp parallel for
  for (int i = 0; i < 20; i++) {
    Cd_t Re = 0.0;
    Cd_t We = Wes[i];
    Cd_t kx = 1.0;
    Cd_t beta = 0.5;
    complex<double> ii(0.0, 1.0);

    Linop<complex<double> > Dyyyy(4), Dyyy(3), Dyy(2), Dy(1);
    Vcd_t c(N + 1), cy(N + 1), cyy(N + 1), a0(N + 1), a1(N + 1), a2(N + 1),
        a3(N + 1);

    Dyyyy.coef << 1.0, 0.0, 0.0, 0.0, 0.0;
    Dyyy.coef << 1.0, 0.0, 0.0, 0.0;
    Dyy.coef << 1.0, 0.0, 0.0;
    Dy.coef << 1.0, 0.0;
    // Set in cheb-point
    Cd_t omega = 0.0;
    c = (ii * omega + (1.0 / We) + (ii * kx * U));
    cy = ii * kx * Uy;
    cyy = ii * kx * Uyy;

    a3 = (2.0 * (-1.0 + beta) * (cy - ii * c * kx * Uy * We)) /
         (c * (1.0 - beta + beta * c * We));

    a2 = (-2.0 * (-1.0 + beta) *
              (pow(cy, 2.0) + Cd_t(0.0, 2.0) * cy * kx * Uy +
               2.0 * pow(kx, 2.0) * pow(Uy, 2.0)) -
          ii * pow(c, 3.0) *
              (Cd_t(0.0, -2.0) * beta * pow(kx, 2.0) + omega * Re +
               kx * Re * U) *
              We -
          (-1.0 + beta) * c *
              (-cyy - Cd_t(0.0, 2.0) * kx * Uy +
               2.0 * kx * Uy * (Cd_t(0.0, -1.0) * cy + kx * Uy) * We) +
          (-1.0 + beta) * pow(c, 2.0) * kx *
              (2.0 * kx - Cd_t(0.0, 3.0) * Uy * We +
               2.0 * kx * pow(Uy, 2.0) * pow(We, 2.0))) /
         (pow(c, 2.0) * (1.0 - beta + beta * c * We));

    a1 = (2.0 * (-1.0 + beta) * kx *
          (6.0 * cy * Uy * (Cd_t(0.0, 1.0) * cy + kx * Uy) +
           pow(c, 3.0) * kx * Uy * We * (Cd_t(0.0, 1.0) * kx + 2.0 * Uy * We) +
           pow(c, 2.0) *
               (-(cy * (kx - Cd_t(0.0, 1.0) * Uy * We)) +
                kx * Uy * We *
                    (-3.0 * Uy - Cd_t(0.0, 2.0) * kx * pow(Uy, 2.0) * We)) -
           Cd_t(0.0, 2.0) * c *
               (2.0 * Uy * (cy - Cd_t(0.0, 1.0) * kx * Uy) +
                Uy * (cyy + kx * Uy * (Cd_t(0.0, 2.0) * cy + kx * Uy) * We)))) /
         (pow(c, 3.0) * (1.0 - beta + beta * c * We));

    a0 = (kx *
          (-12.0 * (-1.0 + beta) * cy * kx * pow(Uy, 2.0) *
               (cy - ii * kx * Uy) * We +
           pow(c, 4.0) *
               (beta * pow(kx, 3.0) + ii * kx * omega * Re +
                ii * Re * (pow(kx, 2.0) * U + Uy)) *
               We -
           (-1.0 + beta) * pow(c, 3.0) * pow(kx, 2.0) *
               (kx - ii * Uy * We + 2.0 * kx * pow(Uy, 2.0) * pow(We, 2.0)) -
           (-1.0 + beta) * pow(c, 2.0) *
               (-(cyy * kx) + 2.0 * kx * pow(Uy, 2.0) * We +
                ii * Uy * (2.0 * pow(kx, 2.0) - cyy * We) +
                Cd_t(0.0, 2.0) * cy * kx * Uy * We *
                    (kx + Cd_t(0.0, 2.0) * Uy * We) +
                2.0 * kx * pow(Uy, 2.0) * We *
                    (pow(kx, 2.0) - cyy * We + Cd_t(0, 6.0) * kx * Uy * We)) +
           2.0 * (-1.0 + beta) * c *
               (2.0 * kx * pow(Uy, 2.0) * (cyy - Cd_t(0, 3.0) * kx * Uy) * We -
                pow(cy, 2.0) * (kx + ii * Uy * We +
                                2.0 * kx * pow(Uy, 2.0) * pow(We, 2.0)) +
                2.0 * cy * kx * Uy *
                    (3.0 * Uy * We +
                     ii * kx * (1.0 + 2.0 * pow(Uy, 2.0) * pow(We, 2.0)))))) /
         (pow(c, 3.0) * (1.0 - beta + beta * c * We));

    LinopMat<Cd_t> Amat(1, 1), B(1, 2), C(2, 1), Ctau(3, 1);
    Amat << (Dyyyy + (a3 * Dyyy) + (a2 * Dyy) + (a1 * Dy) + a0);

    BcMat<Cd_t> lbc(2, 1), rbc(2, 1);
    lbc.L << 1.0, //
        Dy;
    rbc.L << 1.0, //
        Dy;
    lbc.eval.setConstant(-1.0);
    rbc.eval.setConstant(1.0);
    Linop<Cd_t> tau22Tov, tau12Tov, tau11Tov;
    tau22Tov =
        (Vcd_t(2.0 / (c * We)) * Dy) + Vcd_t((Cd_t(0.0, 2.0) * kx * Uy) / c);

    tau12Tov = (Vcd_t(Cd_t(0.0, 1.0) / (c * kx * We)) * Dyy) +
               (Vcd_t((2.0 * Uy) / (pow(c, 2.0) * We)) * Dy) +
               Vcd_t((-(c * Uy * We) +
                      Cd_t(0.0, 1.0) * kx *
                          (c + 2.0 * pow(Uy, 2.0) * We * (1.0 + c * We))) /
                     (pow(c, 2.0) * We));

    tau11Tov =
        (Vcd_t((Cd_t(0.0, 2.0) * Uy * (1.0 + c * We)) /
               (pow(c, 2.0) * kx * We)) *
         Dyy) +
        (Vcd_t((-2.0 * pow(c, 2.0) -
                4.0 * pow(Uy, 2.0) * (-1.0 + pow(c, 2.0) * pow(We, 2.0))) /
               (pow(c, 3.0) * We)) *
         Dy) +
        Vcd_t((2.0 * Uy *
               (-(c * Uy * We * (1.0 + 2.0 * c * We)) +
                Cd_t(0.0, 1.0) * kx *
                    (c + 2.0 * pow(Uy, 2.0) * We * (1.0 + c * We)))) /
              (pow(c, 3.0) * We));

    Ctau << tau11Tov, //
        tau12Tov,     //
        tau22Tov;
    C << Dy / (-ii * kx), //
        1.0;
    B << Vcd_t((Cd_t(0.0, -1.0) * c * kx * We) / (1.0 - beta + beta * c * We)) *
             Dy / (-1.0),
        Vcd_t((c * pow(kx, 2.0) * We) / (1.0 - beta + beta * c * We));

    SingularValueDecomposition<complex<double> > svd;

    // For stress output
    svd.compute(Amat, B, Ctau, lbc, rbc, lbc, rbc, Amat.r * 2 * (N + 1));
    std::cout << "i = " << i << '\n';
    ofstream outf;
    string filename =
        (string("data/") + string(flowType) +
         string("_Re_") + int2str(real(Re)) + string("_We_") +
         int2str(real(We)) + string("_beta_") + int2str(real(beta) * 100.0) +
         string("_N_") + int2str(N) + string("_stress.txt"));
    outf.open(filename.c_str());
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> out_to_file(
        svd.eigenvalues.size(), 2);

    out_to_file << svd.eigenvalues.real(), svd.eigenvalues.imag();
    outf << out_to_file;
    outf.close();

    // For velocity output
    svd.compute(Amat, B, C, lbc, rbc, lbc, rbc, Amat.r * 2 * (N + 1));
    std::cout << "i = " << i << '\n';
    string filename2 =
        (string("data/") + string(flowType) +
         string("_Re_") + int2str(real(Re)) + string("_We_") +
         int2str(real(We)) + string("_beta_") + int2str(real(beta) * 100.0) +
         string("_N_") + int2str(N) + string("_vel.txt"));
    outf.open(filename2.c_str());
    out_to_file.resize(svd.eigenvalues.size(), 2);

    out_to_file << svd.eigenvalues.real(), svd.eigenvalues.imag();
    outf << out_to_file;
    outf.close();
  }
  return 0;
}
