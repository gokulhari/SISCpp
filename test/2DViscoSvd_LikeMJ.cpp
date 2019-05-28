

#define EIGEN_USE_MKL_ALL
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
  valarray<double> Wes(10);
  Eigen::VectorXd psd(10);
  // Number of Chebyshev polynomials
  for (int i = 0; i < 20; i++) {
    Wes[i] = double(i + 1);
  }

  N = 511;
  sis_setup();
  Vcd_t U(N + 1), Uy(N + 1), Uyy(N + 1);
  string flowType("Couette");
  // string flowType("Couette");
std::cout << "I'm here" << std::flush << '\n';
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

//#pragma omp parallel for
  for (int i = 0; i < 10; i++) {
    Cd_t Re = 0.0;
    Cd_t We = Wes[i];
    Cd_t kx = 1.0;
    Cd_t beta = 0.5;
    Cd_t ii(0.0, 1.0);

    Linop<Cd_t> Dyyyy(4), Dyyy(3), Dyy(2), Dy(1);
    Vcd_t c(N + 1), cy(N + 1), cyy(N + 1), a0(N + 1), a1(N + 1), a2(N + 1),
        a3(N + 1), a4(N + 1);

    Dyyyy.coef << 1.0, 0.0, 0.0, 0.0, 0.0;
    Dyyy.coef << 1.0, 0.0, 0.0, 0.0;
    Dyy.coef << 1.0, 0.0, 0.0;
    Dy.coef << 1.0, 0.0;
    // Set in cheb-point
    Cd_t omega = 0.0;
    c = (ii * omega + 1.0 + (ii * kx * We * U));
    cy = ii * kx * We * Uy;
    cyy = ii * kx * We * Uyy;

    a4 = -beta + (-1.0 + beta)/c;

    a3 = (-2.0*(-1.0 + beta)*(cy - Cd_t(0.0,1.0)*c*kx*Uy*We))/pow(c,2.0);

    a2 =  2.0*beta*pow(kx,2.0) - ((-1.0 + beta)*(-2.0*pow(cy,2.0) - 4.0*kx*Uy*We*(Cd_t(0.0,1.0)*cy + kx*Uy*We) +
         pow(c,2.0)*kx*(2.0*kx - Cd_t(0.0,3.0)*Uyy*We + 2.0*kx*pow(Uy,2.0)*pow(We,2.0)) + c*(cyy + Cd_t(0.0,2.0)*kx*We*(Uyy + Uy*(cy + Cd_t(0.0,1.0)*kx*Uy*We)))))/
     pow(c,3.0);

    a1 = (Cd_t(0.0,-2.0)*(-1.0 + beta)*kx*(6.0*cy*Uy*We*(cy - Cd_t(0.0,1.0)*kx*Uy*We) + pow(c,3.0)*kx*Uy*We*(kx - Cd_t(0.0,2.0)*Uyy*We) +
        pow(c,2.0)*(Uyy*We*(cy + Cd_t(0.0,3.0)*kx*Uy*We) + kx*(Cd_t(0.0,1.0)*cy - 2.0*kx*pow(Uy,3.0)*pow(We,3.0))) -
        2.0*c*We*(2.0*cy*(Uyy + Cd_t(0.0,1.0)*kx*pow(Uy,2.0)*We) + Uy*(cyy + kx*We*(Cd_t(0.0,-2.0)*Uyy + kx*pow(Uy,2.0)*We)))))/pow(c,4.0);

    a0 = (kx*(-(beta*pow(c,4.0)*pow(kx,3.0)) + 12.0*(-1.0 + beta)*cy*kx*pow(Uy,2.0)*pow(We,2.0)*(cy - Cd_t(0.0,1.0)*kx*Uy*We) +
        (-1.0 + beta)*pow(c,3.0)*pow(kx,2.0)*(kx - Cd_t(0.0,1.0)*Uyy*We + 2.0*kx*pow(Uy,2.0)*pow(We,2.0)) +
        (-1.0 + beta)*pow(c,2.0)*(-(cyy*kx) + Cd_t(0.0,1.0)*(-cyy + 2.0*pow(kx,2.0))*Uyy*We + 2.0*kx*pow(Uyy,2.0)*pow(We,2.0) +
           Cd_t(0.0,2.0)*cy*kx*Uy*We*(kx + Cd_t(0.0,2.0)*Uyy*We) + 2.0*kx*pow(Uy,2.0)*pow(We,2.0)*(-cyy + pow(kx,2.0) + Cd_t(0.0,6.0)*kx*Uyy*We)) +
        2.0*(-1.0 + beta)*c*(-2.0*cyy*kx*pow(Uy,2.0)*pow(We,2.0) + cy*kx*(cy - Cd_t(0.0,2.0)*kx*Uy*We)*(1.0 + 2.0*pow(Uy,2.0)*pow(We,2.0)) +
           Cd_t(0.0,1.0)*Uyy*We*(pow(cy,2.0) + Cd_t(0.0,6.0)*cy*kx*Uy*We + 6.0*pow(kx,2.0)*pow(Uy,2.0)*pow(We,2.0)))))/pow(c,4.0);

    LinopMat<Cd_t> Amat(1, 1), B(1, 2), C(2, 1), Ctau(3, 1);
    Amat << ((a4*Dyyyy) + (a3 * Dyyy) + (a2 * Dyy) + (a1 * Dy) + a0);

    BcMat<Cd_t> lbc(2, 1), rbc(2, 1);
    lbc.L << 1.0, //
        Dy;
    rbc.L << 1.0, //
        Dy;
    lbc.eval.setConstant(-1.0);
    rbc.eval.setConstant(1.0);
    Linop<Cd_t> tau22Tov, tau12Tov, tau11Tov;
    tau22Tov =
        Vcd_t((Cd_t(0.0,-2.0)*kx*Uy*We)/pow(c,2.0))*Dy + Vcd_t((2.0*pow(kx,2.0)*Uy*We)/c);
    tau12Tov = Vcd_t(1.0/c)*Dyy + Vcd_t((Cd_t(0.0,-2.0)*kx*Uy*We)/pow(c,2.0))*Dy +
      Vcd_t((kx*(Cd_t(0.0,1.0)*c*Uyy*We + kx*(c + 2.0*(1.0 + c)*pow(Uy,2.0)*pow(We,2.0))))/pow(c,2.0));

    tau11Tov = Vcd_t((2.0*(1.0 + c)*Uy*We)/pow(c,2.0))*Dyy + Vcd_t((Cd_t(0.0,2.0)*pow(c,2.0)*kx + Cd_t(0.0,4.0)*(-1.0 + pow(c,2.0))*kx*pow(Uy,2.0)*pow(We,2.0))/pow(c,3.0))*Dy +
   Vcd_t((2.0*kx*Uy*We*(Cd_t(0.0,1.0)*c*(1.0 + 2.0*c)*Uyy*We + kx*(c + 2.0*(1.0 + c)*pow(Uy,2.0)*pow(We,2.0))))/pow(c,3.0));


    Ctau << tau11Tov, //
            tau12Tov, //
            tau22Tov;
    C << Dy, //
        -ii*kx;
    B << Dy,-ii*kx;

    SingularValueDecomposition<Cd_t> svd;
    Cd_t ans = svd.PowerSpectralDensity(Amat, B, C, lbc, rbc);
    psd(i) = real(ans);


  }
  std::cout << "psd: \n" << psd << '\n';
  Eigen::MatrixXd mat(10,2);

  mat << Eigen::VectorXd::LinSpaced(10,1,10), psd;
  ofstream outfile;
  outfile.open("data/psd2.txt");
  outfile << mat;
  outfile.close();
  return 0;
}
