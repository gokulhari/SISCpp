/// \brief example with discriptor approach.
/// For pipe flow, tallying eigenvalues with those of Khorrami.
#define EIGEN_USE_BLAS
#define SIS_USE_LAPACK
#include <fstream>
#include <iostream>
#include <sis.hpp>
#include <string>

using namespace std;
typedef valarray<double> Vd_t;
typedef valarray<complex<double> > Vcd_t;
complex<double> ii(0.0, 1.0);

int main() {
  using namespace sis;
  int bre;
  // Number of Chebyshev polynomials
  N = 127;
  sis_setup();




  valarray<double> kz(1);
  valarray<double> kth(1);
  kth[0] = 0.0;
  kz[0] = 1.0;


  valarray<double> y, W(N + 1), Wr(N + 1), Wrr(N + 1), r(N + 1);
  // Set in cheb-point
  setChebPts(y);
  // Scaled y, for 0 to 1
  r = 0.5 * y + 0.5;
  // Velocity and derivative for PP flow
  W = 1.0 - r * r;
  Wr = -2.0 * r;
  Wrr = -2.0;
  double Re = 2000.0;

  Linop<double> Drr, Dr;
  // Note that Drr, Dr are scaled for r = 0.5 y + 0.5, as Chebyshev
  // polynomials work in (-1,1)
  Drr.set(2);
  Drr.coef << 4.0, 0.0, 0.0;

  Dr.set(1);
  Dr.coef << 2.0, 0.0;

  LinopMat<std::complex<double> > A, B(2, 3), C(3, 2);

  SingularValueDecomposition<std::complex<double> > svd;
  GeneralizedEigenSolver<std::complex<double> > eigs;
  svd.isSingular = SIS_SINGULAR;

  ofstream outf;

  double omega = 0.0;

  B.resize(4, 3);
  C.resize(3, 4);
  B << 1.0, 0.0, 0.0, //
      0.0, 1.0, 0.0,  //
      0.0, 0.0, 1.0,  //
      0.0, 0.0, 0.0;
  C << 1.0, 0.0, 0.0, 0.0, //
      0.0, 1.0, 0.0, 0.0,  //
      0.0, 0.0, 1.0, 0.0;  //
  // The following uses conjugate symmetry as all values have to be real in
  /// physical space. So Nz/2.

  // Below are 2D matrices stored in row-major format in a 1D array in the order
  // Ny x Nz. In the row-major format, the element (j,k) in the 2D matrix is
  // indexed by (j*Nz + k) in the 1D array.
  // These can be easily manipulated using valarray slices.

  // slice notation : slice(start,size,stride);
  // slice to access nl(:,k) in matlab notation: slice(k,N+1,Nz)
  // slice to access nl(i,:) in matlab notation: slice(i*Nz,Nz,1);
  // Replace Nz by Nz / 2 in above 3 lines for complex type due to
  // conjugate symmetry



  {
int j = 0;
    int k = 0;
    complex<double> iiomega = ii*omega;
    BcMat<std::complex<double> > Lbc(3, 4), Rbc(4, 4), LbcAd(4, 4),
        RbcAd(4, 4);
    Linop<double> Delta(2), Deltabyr(2);
    LinopMat<complex<double> > Lmat(4, 4), Mmat(4, 4);
    if (kth[k] == 0) {
      Lbc.L << 1.0, 0.0, 0.0, 0.0,               //
          0.0, 1.0, 0.0, 0.0, //
          0.0, 0.0, Dr, 0.0; //
        //  Dr, 0.0, ii * kz[j], 0.0;
      Rbc.L << 1.0, 0.0, 0.0, 0.0,                //
          0.0, 1.0, 0.0, 0.0, //
          0.0, 0.0, 1.0, 0.0, //
          Dr, 0.0, 0.0, 0.0;

      LbcAd.L << 1.0, 0.0, 0.0, 0.0,             //
          0.0, 1.0, 0.0, 0.0,//
           0.0, 0.0, Dr, 0.0, //
          Dr, 0.0, -ii * kz[j], 0.0;
      RbcAd.L << 1.0, 0.0, 0.0, 0.0,              //
          0.0, 1.0, 0.0, 0.0, //
          0.0, 0.0, 1.0, 0.0, //
          Dr, 0.0, 0.0, 0.0;
      Lbc.eval.setConstant(-1);
      Rbc.eval.setConstant(1);
      LbcAd.eval.setConstant(-1);
      RbcAd.eval.setConstant(1);

    } else if (kth[k] == 1 || kth[k] == -1){
      Lbc.L << ii * kth[k], -1.0, 0.0, 0.0,               //
          0.0, 0.0, 1.0, 0.0, //
          Dr, 0.0, 0.0, 0.0;
      Rbc.L << 1.0, 0.0, 0.0, 0.0,                //
          0.0, 1.0, 0.0, 0.0, //
          0.0, 0.0, 1.0, 0.0, //
          Dr, 0.0, 0.0, 0.0;
      LbcAd.L << -ii * kth[k], -1.0, 0.0, 0.0,               //
          0.0, 0.0, 0.0, 1.0, //
          0.0, 0.0, 1.0, 0.0, //
          Dr, 0.0, 0.0, 0.0;
      RbcAd.L << 1.0, 0.0, 0.0, 0.0,              //
          0.0, 1.0, 0.0, 0.0, //
          0.0, 0.0, 1.0, 0.0, //
          Dr, 0.0, 0.0, 0.0;
      Lbc.eval.setConstant(-1);
      Rbc.eval.setConstant(1);
      LbcAd.eval.setConstant(-1);
      RbcAd.eval.setConstant(1);
    } else {
      Lbc.L << 1.0, 0.0, 0.0, 0.0,               //
          0.0, 1.0, 0.0, 0.0, //
          0.0, 0.0, 1.0, 0.0; //
      //    Dr, 0.0, 0.0, 0.0;
      Rbc.L << 1.0, 0.0, 0.0, 0.0,                //
          0.0, 1.0, 0.0, 0.0, //
          0.0, 0.0, 1.0, 0.0, //
          Dr, 0.0, 0.0, 0.0;
      LbcAd.L << 1.0, 0.0, 0.0, 0.0,               //
          0.0, 1.0, 0.0, 0.0, //
          0.0, 0.0, 1.0, 0.0, //
          Dr, 0.0, 0.0, 0.0;
      RbcAd.L << 1.0, 0.0, 0.0, 0.0,              //
          0.0, 1.0, 0.0, 0.0, //
          0.0, 0.0, 1.0, 0.0, //
          Dr, 0.0, 0.0, 0.0;
      Lbc.eval.setConstant(-1);
      Rbc.eval.setConstant(1);
      LbcAd.eval.setConstant(-1);
      RbcAd.eval.setConstant(1);
    }
    Delta.ncc(2);
    Delta.coefFun << 4.0 * r * r, 2.0 * r,
        -(kz[j] * kz[j] * r * r) - (kth[k] * kth[k]);

    Mmat << r * r, 0.0, 0.0, 0.0, //
        0.0, r * r, 0.0, 0.0,     //
        0.0, 0.0, r * r, 0.0,     //
        0.0, 0.0, 0.0, 0.0; //

    Lmat(0, 0) =
        (Delta / Re - (1.0 / Re)) - Vcd_t(ii * Vd_t(kz[j] * r * r * W));
    Lmat(0, 1) = ii * (-2.0 * kth[k] / Re);
    Lmat(0, 2) = 0.0;
    Lmat(0, 3) = -(r * (r * Dr));

    Lmat(1, 0) = ii * (2.0 * kth[k] / Re);
    Lmat(1, 1) =
        (Delta / Re - (1.0 / Re)) - Vcd_t(ii * Vd_t(kz[j] * r * r * W));
    Lmat(1, 2) = 0.0;
    Lmat(1, 3) = -ii * Vd_t(kth[k] * r);

    Lmat(2, 0) = -r * r * Wr;
    Lmat(2, 1) = 0.0;
    Lmat(2, 2) = (Delta / Re) - Vcd_t(ii * Vd_t(kz[j] * r * r * W));
    Lmat(2, 3) = -ii * Vd_t(kz[j] * r * r);

    Lmat(3, 0) = (r * Dr) + 1.0;
    Lmat(3, 1) = ii * kth[k];
    Lmat(3, 2) = ii * Vd_t(r * kz[j]);
    Lmat(3, 3) = 0.0;

    BcMat<complex<double> > bc(7,4);
    bc.L << Lbc.L,//
            Rbc.L;
    bc.eval << Lbc.eval,//
              Rbc.eval;
    eigs.compute(Lmat, Mmat, 8*(N + 1), bc);
    //eigs.sortByLargestReal();
    cout.precision(9);

    eigs.keepConverged();
    eigs.removeInf();
    eigs.sortByLargestReal();
    std::cout << "Evals:\n " << eigs.eigenvalues << '\n';
    ofstream outf;
    outf.open("data/pipeFlowkth0Re2000.txt");
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> out_to_file(eigs.eigenvalues.rows(),2);
    out_to_file << eigs.eigenvalues.real(), eigs.eigenvalues.imag();
    outf << out_to_file;
    outf.close();
  }

  return 0;
}
