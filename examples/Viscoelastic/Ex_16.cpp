/// \file Visco_3D_pipe.cpp
/// \brief Finds Eigenvalues for 3D Cylindrical flow of Oldroyd-B fluid.
///
/// The velocity vector is given by \f$ \mathbf{v} \;=\; u\,\mathbf{e}_r +
/// v\,\mathbf{e}_\theta + w
/// \,\mathbf{e}_z  \f$.
/// The problem can be reduced to the following discriptor
/// form: \f[
///  \lambda \; \left[ \begin{array}{cccccccccc}
/// r^2 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0  \\
/// 0 & r^2 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
/// 0 & 0 & r^2 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
/// 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
/// 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 \\
/// 0 & 0 & 0 & 0 & 0 & r & 0 & 0 & 0 & 0 \\
/// 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
/// 0 & 0 & 0 & 0 & 0 & 0 & 0 & r & 0 & 0 \\
/// 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & r & 0 \\
/// 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 \\
/// \end{array} \right]  \left[ \begin{array}{c}
/// u\\
/// v \\
/// w \\
/// p\\
/// \tau_{rr}\\
/// \tau_{r\theta}\\
/// \tau_{rz}\\
/// \tau_{\theta \theta}\\
/// \tau_{\theta z}\\
/// \tau_{zz}
/// \end{array} \right] \;=\;
/// \left[ \begin{array}{cccccccccc}
/// \frac{r^2\,\beta}{Re}\Delta -\frac{\beta}{Re} - i\,k_z\,r^2\,
/// W(r) & -2\,i\,k_\theta\,\beta\, /Re & 0
/// & -r^2\,\partial_r & \frac{(1-\beta)\,r}{Re} (r\,\partial_r + 1)
/// & \frac{1-\beta}{Re}\,i\,k_\theta \,r &  \frac{1-\beta}{Re}\,i\,k_z \,r^2
/// & -\frac{1-\beta}{Re}\, r & 0 & 0 \\
/// 2\,i\,k_\theta\,\beta\,/Re & \frac{\beta \, r^2}{Re}\Delta -\frac{\beta}{Re}
/// -
///  i\,k_z\,r^2\, W(r) & 0 & -i\,k_\theta\, r & 0 & \frac{1-\beta}{Re}\,r\,
/// (r\,\partial_r + 2 ) & 0 & \frac{1-\beta}{Re}\,i\,k_\theta\,r
/// & \frac{1-\beta}{Re}\,i\,k_z \,r^2 & 0 \\
/// -r^2 W' & 0  &  \frac{\beta\,r^2}{Re}\Delta - i\,k_z\,r^2 W(r) &
/// -i\,k_z\,r^2 & 0 & 0 & \frac{1-\beta}{Re}\,r\,(r\,\partial_r + 1) & 0 &
/// \frac{1-\beta}{Re}
/// \,i\,k_\theta\,r & \frac{1-\beta}{Re}\,i\,k_z\,r^2\\
/// r\,\partial_r +1 & i\,k_\theta  &  r\,i\,k_z & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
/// 2\,\partial_r/We + 2\,i\,k_z\,T_{rz}
/// & 0 & 0
/// & 0 &-1/We - i\,k_z\,W  & 0 & 0 & 0 & 0 & 0\\
/// i\,k_\theta /We
/// & \frac{1}{We}\, (r\,\partial_r -1) + i\,k_z\,r\,T_{rz} & 0 & 0 & 0
/// & -r/We - i\,k_z\,W\,r & 0 & 0 & 0 & 0\\
///  T_{rz}\,\partial_r + i\,k_z/We + i\,k_z\,T_{zz} -T_{rz}'
/// & 0 & \partial_r/We + i\,k_z\,T_{rz}  & 0 & W'&0
/// & -1/We - i\,k_z\,W & 0 & 0 & 0 \\
/// 2/We & 2\,i\,k_\theta/We & 0 & 0 & 0 & 0 & 0
/// & -r/We - i\,k_z\,W\,r & 0  & 0\\
/// 0 & r\,T_{rz}\,\partial_r
///  + i\,k_z\,r/We - T_{rz} + i\,k_z\,r\,T_{zz}
/// &  i\,k_\theta/We & 0 & 0 & r\,W' & 0 & 0 & -r/We - i\,k_z\,W\,r & 0\\
/// -T_{zz}' & 0 & 2\,T_{rz}\,\partial_r + 2\,i\,k_z/We + 2\,i\,k_z\,T_{zz}
/// & 0 & 0 & 0 & 2\,W' & 0 & 0 & -1/We - i\,k_z\,W
/// \end{array} \right] \left[ \begin{array}{c}
/// u \\
/// v \\
/// w \\
/// p\\
/// \tau_{rr}\\
/// \tau_{r\theta}\\
/// \tau_{rz}\\
/// \tau_{\theta \theta}\\
/// \tau_{\theta z}\\
/// \tau_{zz}
/// \end{array} \right],
/// \f]
///
/// where, \f$\Delta = \frac{1}{r}(r\partial_{rr} + \partial_r) - k_z^2  -
/// k_\theta^2 / r^2\f$, the Laplacian for 3D Cylindrical coordinates, and \f$W
/// = 1-r^2\f$ for pipe flow.
///
/// Boundary conditions that come from no slip at the wall:
/// \f{align}
/// u(1) \;&=\; 0\\
/// v(1)  \;&=\; 0\\
/// w(1)  \;&=\; 0 \\
/// u'(1)  \;&=\; 0
/// \f}
/// The last one comes by applying continuity to the first three.
///
/// Boundary conditions in the centerline come from applying \f$\partial_\theta
/// (\mathbf v) = 0,\quad r\rightarrow 0\f$,
/// \f{align}
/// u(0)  \,+\, i\,k_\theta\,v(0)  \;&=\;  0\\
/// i\, k_{\theta} \, u(0)  \,-\, v(0) \;&=\; 0\\
/// i\,k_\theta \,w(0)  \;&=\; 0
/// \f}
///
/// Note that \f$k_\theta\f$ has to be an integer due to periodicity. There are
/// special cases boundary conditions. These boundary conditions are also valid
/// for inviscid flows.
/// For \f$k_\theta=0\f$
/// \f{align}
/// u(0) \;&=\; 0\\
/// v(0)  \;&=\; 0\\
/// u'(0)  \,+\, i\,k_z\,w(0) \;&=\; 0
/// \f}
/// Note that we need one more boundary condition. We impose the boundary
/// condition that \f$w'(0) = 0\f$, which is valid only for viscous flows
/// (due to implied symmetry). In
/// inviscid flows one will not need this boundary condition.

#define EIGEN_FAST_MATH 0
#define EIGEN_DONT_PARALLELIZE
//#define EIGEN_USE_BLAS
#define SIS_USE_LAPACK
#define SIS_DONT_SORT
//#define SIS_USE_FEAST
#include <omp.h>
#include <fstream>
#include <iostream>
#include <sis.hpp>

using namespace std;
typedef complex<double> Cd_t;
typedef valarray<double> Vd_t;
typedef valarray<complex<double> > Vcd_t;
int main(int argc, char** argv) {
  using namespace sis;
  int bre;

double Re = 50;
double kz = 1.0;
double kth = 0.0;
double beta = 0.5;
double We = 50;
N = 127;
sscanf(argv[1], "%d", &N);
Re = atof(argv[2]);
We = atof(argv[3]);
beta = atof(argv[4]);
//sscanf(argv[3], "%d", &We);

for (int i=0; i < argc; i++){
  cout << "argv["<< i << "]: " << argv[i]<< "\n";
}
std::cout << "N = " << N << '\n';
std::cout << "Re = " << Re << '\n';
std::cout << "We = " << We << '\n';
std::cout << "beta = " << beta << '\n';
sis_setup();
//exit(1);
  // Number of Chebyshev polynomials

  //N = 63;

  valarray<double> y, W(N + 1), Wr(N + 1), Wrr(N + 1), r(N + 1), Trz(N + 1),
      Tzz(N + 1), Trz_r(N + 1), Tzz_r(N + 1);
  // Set in cheb-point
  setChebPts(y);
  // Scaled y, for 0 to 1
  r = 0.5 * y + 0.5;



  W = 1.0 - r * r;
  Wr = -2.0 * r;
  Wrr = -2.0;

  Trz = Wr;
  Tzz = 2.0 * We * Wr * Wr;

  Trz_r = Wrr;
  Tzz_r = diff(Tzz);

  complex<double> ii(0.0, 1.0);

  Linop<double> Drr, Dr, Delta, I;
  LinopMat<complex<double> > Lmat(10, 10), Mmat(10, 10);
  // Note that Dyy, Dy and Delta are scaled for r = 0.5 y + 0.5, as Chebyshev
  // polynomials work in (-1,1)

  // Identity operator
  I.set(0);
  I.coef << 1.0;

  // second derivative operator.
  Drr.set(2);
  Drr.coef << 4.0, 0.0, 0.0;

  // First derivative operator.
  Dr.set(1);
  Dr.coef << 2.0, 0.0;

  // Laplacian operator.
  Delta.ncc(2);
  Delta.coefFun << 4.0 * r * r, 2.0 * r, -(kz * kz * r * r) - (kth * kth);

  // Set Lmat:
  Lmat(0, 0) =
      (beta * Delta / Re - (beta / Re)) - Vcd_t(ii * Vd_t(kz * r * r * W));
  Lmat(0, 1) = ii * (-2.0 * kth * beta / Re);
  Lmat(0, 2) = 0.0;
  Lmat(0, 3) = -(r * (r * Dr));
  Lmat(0, 4) = (1.0 - beta) * ((r * (r * Dr)) + r) / Re;
  Lmat(0, 5) = ii * Vd_t((1.0 - beta) * kth * r / Re);
  Lmat(0, 6) = ii * Vd_t((1.0 - beta) * kz * r * r / Re);
  Lmat(0, 7) = -(1.0 - beta) * r / Re;
  Lmat(0, 8) = 0.0;
  Lmat(0, 9) = 0.0;

  Lmat(1, 0) = ii * (2.0 * kth * beta / Re);
  Lmat(1, 1) =
      (beta * Delta / Re - (beta / Re)) - Vcd_t(ii * Vd_t(kz * r * r * W));
  Lmat(1, 2) = 0.0;
  Lmat(1, 3) = -ii * Vd_t(kth * r);
  Lmat(1, 4) = 0.0;
  Lmat(1, 5) = (1.0 - beta) * ((r * (r * Dr)) + Vd_t(2.0 * r)) / Re;
  Lmat(1, 6) = 0.0;
  Lmat(1, 7) = ii * Vd_t((1.0 - beta) * kth * r / Re);
  Lmat(1, 8) = ii * Vd_t((1.0 - beta) * kz * r * r / Re);
  Lmat(1, 9) = 0.0;

  Lmat(2, 0) = -r * r * Wr;
  Lmat(2, 1) = 0.0;
  Lmat(2, 2) = (beta * Delta / Re) - Vcd_t(ii * Vd_t(kz * r * r * W));
  Lmat(2, 3) = -ii * Vd_t(kz * r * r);
  Lmat(2, 4) = 0.0;
  Lmat(2, 5) = 0.0;
  Lmat(2, 6) = (1.0 - beta) * ((r * (r * Dr)) + r) / Re;
  Lmat(2, 7) = 0.0;
  Lmat(2, 8) = ii * Vd_t((1.0 - beta) * kth * r / Re);
  Lmat(2, 9) = ii * Vd_t((1.0 - beta) * kz * r * r / Re);

  Lmat(3, 0) = (r * Dr) + 1.0;
  Lmat(3, 1) = ii * kth;
  Lmat(3, 2) = ii * Vd_t(r * kz);
  Lmat(3, 3) = 0.0;
  Lmat(3, 4) = 0.0;
  Lmat(3, 5) = 0.0;
  Lmat(3, 6) = 0.0;
  Lmat(3, 7) = 0.0;
  Lmat(3, 8) = 0.0;
  Lmat(3, 9) = 0.0;

  Lmat(4, 0) = (2.0 * Dr / We) + Vcd_t(ii * Vd_t(2.0 * kz * Trz));
  Lmat(4, 1) = 0.0;
  Lmat(4, 2) = 0.0;
  Lmat(4, 3) = 0.0;
  Lmat(4, 4) = (-1.0 / We) - Vcd_t(ii * Vd_t(kz * W));
  Lmat(4, 5) = 0.0;
  Lmat(4, 6) = 0.0;
  Lmat(4, 7) = 0.0;
  Lmat(4, 8) = 0.0;
  Lmat(4, 9) = 0.0;

  Lmat(5, 0) = ii * kth / We;
  Lmat(5, 1) = (r * Dr / We - (1.0 / We)) + Vcd_t(ii * Vd_t(kz * r * Trz));
  Lmat(5, 2) = 0.0;
  Lmat(5, 3) = 0.0;
  Lmat(5, 4) = 0.0;
  Lmat(5, 5) = -(r * I / We) - Vcd_t(ii * Vd_t(kz * W * r));
  Lmat(5, 6) = 0.0;
  Lmat(5, 7) = 0.0;
  Lmat(5, 8) = 0.0;
  Lmat(5, 9) = 0.0;

  Lmat(6, 0) =
      ((Trz * Dr - Trz_r) + (ii * kz / We)) + Vcd_t(ii * Vd_t(kz * Tzz));
  Lmat(6, 1) = 0.0;
  Lmat(6, 2) = (Dr / We) + Vcd_t(ii * Vd_t(kz * Trz));
  Lmat(6, 3) = 0.0;
  Lmat(6, 4) = Wr;
  Lmat(6, 5) = 0.0;
  Lmat(6, 6) = (-1.0 / We) - Vcd_t(ii * Vd_t(kz * W));
  Lmat(6, 7) = 0.0;
  Lmat(6, 8) = 0.0;
  Lmat(6, 9) = 0.0;

  Lmat(7, 0) = 2.0 / We;
  Lmat(7, 1) = ii * (2.0 * kth / We);
  Lmat(7, 2) = 0.0;
  Lmat(7, 3) = 0.0;
  Lmat(7, 4) = 0.0;
  Lmat(7, 5) = 0.0;
  Lmat(7, 6) = 0.0;
  Lmat(7, 7) = -(r * I / We) - Vcd_t(ii * Vd_t(kz * W * r));
  Lmat(7, 8) = 0.0;
  Lmat(7, 9) = 0.0;

  Lmat(8, 0) = 0.0;
  Lmat(8, 1) = ((r * (Trz * Dr)) - Trz) +
               Vcd_t(ii * Vd_t((kz * r / We) + (kz * r * Tzz)));
  Lmat(8, 2) = ii * kth / We;
  Lmat(8, 3) = 0.0;
  Lmat(8, 4) = 0.0;
  Lmat(8, 5) = r * Wr;
  Lmat(8, 6) = 0.0;
  Lmat(8, 7) = 0.0;
  Lmat(8, 8) = -(r * I / We) - Vcd_t(ii * Vd_t(kz * W * r));
  Lmat(8, 9) = 0.0;

  Lmat(9, 0) = -Tzz_r;
  Lmat(9, 1) = 0.0;
  Lmat(9, 2) =
      (2.0 * (Trz * Dr)) + Vcd_t(ii * Vd_t((2.0 * kz / We) + (2.0 * kz * Tzz)));
  Lmat(9, 3) = 0.0;
  Lmat(9, 4) = 0.0;
  Lmat(9, 5) = 0.0;
  Lmat(9, 6) = 2.0 * Wr;
  Lmat(9, 7) = 0.0;
  Lmat(9, 8) = 0.0;
  Lmat(9, 9) = (-1.0 / We) - (ii * Vd_t(kz * W));

  Mmat << r * r, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
      0.0, r * r, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,     //
      0.0, 0.0, r * r, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,     //
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,       //
      0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,       //
      0.0, 0.0, 0.0, 0.0, 0.0, r, 0.0, 0.0, 0.0, 0.0,         //
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,       //
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, r, 0.0, 0.0,         //
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, r, 0.0,         //
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;       //

  // This is for boundary constraint from continuity:
  BcMat<complex<double> > Lbc(1, 10), Mbc(1, 10);

    Lbc.resize(10, 10);
    Lbc.L << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,      //
        Dr, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,       //
        Dr, 0.0, ii * kz, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   //
        0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,      //
        0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,      //
        0.0, 0.0, Dr, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,       //
        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,      //
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,      //
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;      //

    Lbc.eval << -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
                1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,          //
                1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,          //
                -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
                0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                      //
                0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  //
                0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  //
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, //
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;

  sis::GeneralizedEigenSolver<complex<double> > eigs;
  //feast::M0 = 2 * (N + 1);
  //feast::M0 = 100;
    eigs.compute(Lmat, Mmat, 15 * (N + 1), Lbc);
    eigs.removeInf();
    eigs.sortByLargestReal();

//  std::cout << "Eigenvalues are : \n" << eigs.eigenvalues << '\n';
  ofstream outf;
  outf.open(string(string("data/Axisymmetric3D_") + string("_Re_") + int2str(int(Re)) +
            string("_We_") + int2str(int(We)) + string("_beta_") +
            int2str(int(beta * 100)) + string("_N_") + int2str(N) +
         //   string("_") + int2str(i) +
            string(".txt")).c_str());

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> out_to_file(
      eigs.eigenvalues.size(), 2);
  out_to_file.col(0) = eigs.eigenvalues.real();
  out_to_file.col(1) = eigs.eigenvalues.imag();
  outf << out_to_file;
  outf.close();

  return 0;
}
