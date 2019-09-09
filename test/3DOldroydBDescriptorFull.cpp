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
  valarray<double> y(N + 1), U(N + 1), Uy(N + 1), Uyy(N + 1), T11(N + 1),
      T12(N + 1);
  // Set in cheb-point
  setChebPts(y);
  // Number of Eigenvalues to compute:
  int num_vals = 15 * (N + 1);

  // Number of fourier modes in span-wise direction:
  int Nz = 4;
  int Nx = 4;
  // Length of domain:
  valarray<double> kz(Nz);
  valarray<double> kx(Nx);
  valarray<double> z(Nz);
  valarray<double> x(Nx);
  valarray<double> omval(Nx);

  kx[0] = 1.0;
  kx[1] = 1.0;
  kx[2] = 1.0;
  kx[3] = 1.0;

  kz[0] = 1.0;
  kz[1] = 1.0;
  kz[2] = 1.0;
  kz[3] = 1.0;

  omval[0] = -1.0;
  omval[1] = -1.0;
  omval[2] = 1.0;
  omval[3] = 1.0;
  double omega = 0.385;

  omval = omega * omval;

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
  // double Re = 2000.0;
  // double We = 2.0;
  // double beta = 0.5;

  double Re = 0.0;
  double We = 1.0;
  double beta = 0.5;
  Linop<double> Dy(1);
  Dy.coef << 1.0, 0.0;

  ofstream outf;

  // Below are 3D matrices stored in row-major format in a 1D array in the order
  // Ny x Nx x Nz. In the row-major format, the element (i,j,k) in the 3D matrix
  // is indexed by (i * Nx * Nz + j * Nz + k) in a 1D array.
  // These can be easily manipulated using valarray slices.
  //
  // slice notation : slice(start, size, stride);
  //
  // slice to access nl(:,j,k) in matlab notation:
  // slice(j * Nz + k, N + 1, Nx * Nz)
  //
  // slice to access nl(i,:,:) in matlab notation: slice(i * Nx * Nz, Nx * Nz,
  // 1);
  //
  // Replace Nz by Nz / 2 in above 3 lines for complex type due to
  // conjugate symmetry as all values have to be real in
  // physical space.
  //
  // Here, I sketch the symmetry of a 8 x 8 matrix. Note that the Nyquist
  // frequency values have to be set to zero, here we denote zeros by blanks
  // "-". Numbers refer to element indices for matrices, indices start from 0 to
  // conform with C++.
  //
  // Note that we choose (x,z) - x for rows, z columns.
  // IMP: 00 must be real.
  // 00  01  02  03  -  03* 02* 01*
  // 10  11  12  13  -  73* 72* 71*
  // 20  21  22  23  -  63* 62* 61*
  // 30  31  32  33  -  53* 52* 51*
  // -   -   -   -  -   -   -   -
  // 30* 51  52  53  -  33* 32* 31*
  // 20* 61  62  63  -  23* 22* 21*
  // 10* 71  72  73  -  13* 12* 11*
  //
  // Hence we store only the elements (:,1:Nz/2) in Matlab notation. There is
  // slight repetition, that is 30*, 20* and 10* not needed to be stored.

  valarray<complex<double> > uvec(Cd_t(0.0, 0.0), (N + 1) * 4),
      vvec(Cd_t(0.0, 0.0), (N + 1) * 4), wvec(Cd_t(0.0, 0.0), (N + 1) * 4),
      pvec(Cd_t(0.0, 0.0), (N + 1) * 4), t11vec(Cd_t(0.0, 0.0), (N + 1) * 4),
      t12vec(Cd_t(0.0, 0.0), (N + 1) * 4), t13vec(Cd_t(0.0, 0.0), (N + 1) * 4),
      t22vec(Cd_t(0.0, 0.0), (N + 1) * 4), t23vec(Cd_t(0.0, 0.0), (N + 1) * 4),
      t33vec(Cd_t(0.0, 0.0), (N + 1) * 4);
  valarray<double> u3d((N + 1) * Nx * Nz), v3d((N + 1) * Nx * Nz),
      w3d((N + 1) * Nx * Nz), p3d((N + 1) * Nx * Nz), t11_3d((N + 1) * Nx * Nz),
      t12_3d((N + 1) * Nx * Nz), t13_3d((N + 1) * Nx * Nz),
      t22_3d((N + 1) * Nx * Nz), t23_3d((N + 1) * Nx * Nz),
      t33_3d((N + 1) * Nx * Nz);
  std::vector<double> omegas;
  omegas.resize(11);
  omega = 0.0;
  double domega = 0.2;

  omegas[0] = 0.0;
  for (int i = 1; i < 11; i++) {
    omegas[i] = omegas[i - 1] + domega;
  }
  for (int k = 0; k < 4; k++) {
#ifdef VALIDATE
    We = 3.96;
    Re = 3960.0;
    kx[k] = 1.15;
    kz[k] = 0.0;
#endif
    LinopMat<complex<double> > Lmat(10, 10), Mmat(10, 10);
    complex<double> iiomega = ii * omval[k];
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

    Lmat(3, 0) = ii * kx[k];
    Lmat(3, 1) = Dy;
    Lmat(3, 2) = ii * kz[k];
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

    Lmat(6, 0) = ii * kz[k] / We;
    Lmat(6, 1) = 0.0;
    Lmat(6, 2) = (ii * kx[k] / We) +
                 (Vcd_t(2.0 * ii * kx[k] * We * Uy * Uy) + (Vd_t(Uy) * Dy));
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
    Lmat(8, 1) = ii * kz[k] / We;
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
    Lmat(9, 2) = 2.0 * ii * kz[k] / We;
    Lmat(9, 3) = 0.0;
    Lmat(9, 4) = 0.0;
    Lmat(9, 5) = 0.0;
    Lmat(9, 6) = 0.0;
    Lmat(9, 7) = 0.0;
    Lmat(9, 8) = 0.0;
    Lmat(9, 9) = ((-1.0 / We) - (ii * kx[k] * U));

    BcMat<complex<double> > Lbc(5, 10), Rbc(5, 10), LbcAd(6, 10), RbcAd(6, 10),
        bc(10, 10);
    Lbc.L << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
        0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,      //
        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,      //
        0.0, Dy, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,       //
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;
    Rbc.L << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
        0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,      //
        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,      //
        0.0, Dy, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,       //
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;
    LbcAd.L << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
        0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,        //
        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,        //
        0.0, Dy, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,         //
        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,        //
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
    RbcAd.L << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
        0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,        //
        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,        //
        0.0, Dy, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,         //
        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,        //
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
    Lbc.eval.setConstant(-1.0);
    Rbc.eval.setConstant(1.0);
    LbcAd.eval.setConstant(-1.0);
    RbcAd.eval.setConstant(1.0);
    bc.L << Lbc.L, //
        Rbc.L;
    bc.eval << Lbc.eval, //
        Rbc.eval;
    Mmat << Re, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
        0.0, Re, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,     //
        0.0, 0.0, Re, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,     //
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    //
        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,    //
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

    //GeneralizedEigenSolver<complex<double> > eigs;
    sis::SingularValueDecomposition<complex<double> > svd;

    LinopMat<complex<double> > A(10, 10), B(10, 3), C(3, 10), Ctau(6, 10);
    iiomega = 0.0;
    A = ((iiomega * Mmat) - Lmat);
    B << 1.0, 0.0, 0.0, //
        0.0, 1.0, 0.0,  //
        0.0, 0.0, 1.0,  //
        0.0, 0.0, 0.0,  //
        0.0, 0.0, 0.0,  //
        0.0, 0.0, 0.0,  //
        0.0, 0.0, 0.0,  //
        0.0, 0.0, 0.0,  //
        0.0, 0.0, 0.0,  //
        0.0, 0.0, 0.0;

    C << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
        0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  //
        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;  //

    Ctau << 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,     //
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,     //
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,     //
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,     //
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
#ifdef VALIDATE
    cout << "I'm in validate ..." << endl;
    eigs.compute(Lmat, Mmat, (N + 1) * Lmat.r, bc);
    eigs.removeInf();
    eigs.sortByLargestReal();
    num_vals = eigs.eigenvalues.size();
    ofstream outf;
    outf.open(string("data/data_for_file18/OldroydB_") + flowType +
              string("_beta_") + int2str(int(beta * 1000)) + string("_We_") +
              int2str(int(We)) + string("_Re_") + int2str(int(Re)) +
              string("_N_") + int2str(N) + string("_") + int2str(k) +
              string(".txt"));
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> out_to_file(num_vals,
                                                                      3);
    out_to_file.col(0) = eigs.eigenvalues.real();
    out_to_file.col(1) = eigs.eigenvalues.imag();
    out_to_file.col(2) = eigs.MPorNot.cast<double>();
    outf << out_to_file;
    outf.close();
// exit(1);
#endif
    cout << "Began computing ... " << endl << flush;
    svd.compute(A, B, C, Lbc, Rbc, LbcAd, RbcAd, 2 * (N + 1) * A.r);
    std::cout << "eigenvalue: " << svd.eigenvalues[0] << "\n";
    num_vals = svd.eigenvalues.size();
    exit(1);
    // Store first two eigenvalues.
//    outf << kx[k] << " " << kz[k] << " " << svd.eigenvalues[0].real() << " "
//         << svd.eigenvalues[1].real() << "\n";
//
//    cout << kx[k] << " " << kz[k] << " " << svd.eigenvalues[0].real() << " "
//         << svd.eigenvalues[1].real() << "\n";
//
//    num_vals = svd.eigenvalues.size();
//    ofstream outf;
//    outf.open(string("data/data_for_file18/OldroydBkxkz1SingularVals_") +
//              flowType + string("_beta_") + int2str(int(beta * 1000)) +
//              string("_We_") + int2str(int(We)) + string("_Re_") +
//              int2str(int(Re)) + string("_N_") + int2str(N) + string("_") +
//              int2str(k) + string(".txt"));
//    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> out_to_file(num_vals,
//                                                                      3);
//    out_to_file.col(0) = svd.eigenvalues.real();
//    out_to_file.col(1) = svd.eigenvalues.imag();
//    out_to_file.col(2) = svd.MPorNot.cast<double>();
//    outf << out_to_file;
//    outf.close();
//    uvec[slice(k, N + 1, Nz)] =
//        svd.eigenvalues[0].real() * svd.eigenvectorsMat[0][0].v;
//
//    vvec[slice(k, N + 1, Nz)] =
//        svd.eigenvalues[0].real() * svd.eigenvectorsMat[0][1].v;
//
//    wvec[slice(k, N + 1, Nz)] =
//        svd.eigenvalues[0].real() * svd.eigenvectorsMat[0][2].v;
//
//    pvec[slice(k, N + 1, Nz)] =
//        svd.eigenvalues[0].real() * svd.eigenvectorsMat[0][3].v;
//
//    t11vec[slice(k, N + 1, Nz)] =
//        svd.eigenvalues[0].real() * svd.eigenvectorsMat[0][4].v;
//
//    t12vec[slice(k, N + 1, Nz)] =
//        svd.eigenvalues[0].real() * svd.eigenvectorsMat[0][5].v;
//
//    t13vec[slice(k, N + 1, Nz)] =
//        svd.eigenvalues[0].real() * svd.eigenvectorsMat[0][6].v;
//
//    t22vec[slice(k, N + 1, Nz)] =
//        svd.eigenvalues[0].real() * svd.eigenvectorsMat[0][7].v;
//    t23vec[slice(k, N + 1, Nz)] =
//        svd.eigenvalues[0].real() * svd.eigenvectorsMat[0][8].v;
//
//    t33vec[slice(k, N + 1, Nz)] =
//        svd.eigenvalues[0].real() * svd.eigenvectorsMat[0][9].v;
  }
  Eigen::VectorXd xval, zval;

  zval = Eigen::VectorXd::LinSpaced(100, -7.8, 7.8);
  xval = Eigen::VectorXd::LinSpaced(100, 0, 12.7);
  Nx = 100;
  Nz = 100;
  Vcd_t u3dc(Cd_t(0.0, 0.0), (N + 1) * Nx * Nz),
      v3dc(Cd_t(0.0, 0.0), (N + 1) * Nx * Nz),
      w3dc(Cd_t(0.0, 0.0), (N + 1) * Nx * Nz),
      p3dc(Cd_t(0.0, 0.0), (N + 1) * Nx * Nz),
      t11_3dc(Cd_t(0.0, 0.0), (N + 1) * Nx * Nz),
      t12_3dc(Cd_t(0.0, 0.0), (N + 1) * Nx * Nz),
      t13_3dc(Cd_t(0.0, 0.0), (N + 1) * Nx * Nz),
      t22_3dc(Cd_t(0.0, 0.0), (N + 1) * Nx * Nz),
      t23_3dc(Cd_t(0.0, 0.0), (N + 1) * Nx * Nz),
      t33_3dc(Cd_t(0.0, 0.0), (N + 1) * Nx * Nz);

  for (int j = 0; j < xval.size(); j++) {
    double x = xval[j];
    for (int k = 0; k < zval.size(); k++) {
      double z = zval[k];
      for (int i = 0; i < 4; i++) {
        double kx1 = kx[i];
        double kz1 = kz[i];
        u3dc[slice(j * Nz + k, N + 1, Nx * Nz)] +=
            Vcd_t(uvec[slice(i, N + 1, 4)]) *
            exp((ii * kx1 * x) + (ii * kz1 * z));
        v3dc[slice(j * Nz + k, N + 1, Nx * Nz)] +=
            Vcd_t(vvec[slice(i, N + 1, 4)]) *
            exp((ii * kx1 * x) + (ii * kz1 * z));
        w3dc[slice(j * Nz + k, N + 1, Nx * Nz)] +=
            Vcd_t(wvec[slice(i, N + 1, 4)]) *
            exp((ii * kx1 * x) + (ii * kz1 * z));
        p3dc[slice(j * Nz + k, N + 1, Nx * Nz)] +=
            Vcd_t(pvec[slice(i, N + 1, 4)]) *
            exp((ii * kx1 * x) + (ii * kz1 * z));
        t11_3dc[slice(j * Nz + k, N + 1, Nx * Nz)] +=
            Vcd_t(t11vec[slice(i, N + 1, 4)]) *
            exp((ii * kx1 * x) + (ii * kz1 * z));
        t12_3dc[slice(j * Nz + k, N + 1, Nx * Nz)] +=
            Vcd_t(t12vec[slice(i, N + 1, 4)]) *
            exp((ii * kx1 * x) + (ii * kz1 * z));
        t13_3dc[slice(j * Nz + k, N + 1, Nx * Nz)] +=
            Vcd_t(t13vec[slice(i, N + 1, 4)]) *
            exp((ii * kx1 * x) + (ii * kz1 * z));
        t22_3dc[slice(j * Nz + k, N + 1, Nx * Nz)] +=
            Vcd_t(t22vec[slice(i, N + 1, 4)]) *
            exp((ii * kx1 * x) + (ii * kz1 * z));
        t23_3dc[slice(j * Nz + k, N + 1, Nx * Nz)] +=
            Vcd_t(t23vec[slice(i, N + 1, 4)]) *
            exp((ii * kx1 * x) + (ii * kz1 * z));
        t33_3dc[slice(j * Nz + k, N + 1, Nx * Nz)] +=
            Vcd_t(t33vec[slice(i, N + 1, 4)]) *
            exp((ii * kx1 * x) + (ii * kz1 * z));
      }
    }
  }
  u3d = real(u3dc);
  v3d = real(v3dc);
  w3d = real(w3dc);
  p3d = real(p3dc);
  t11_3d = real(t11_3dc);
  t12_3d = real(t12_3dc);
  t13_3d = real(t13_3dc);
  t22_3d = real(t22_3dc);
  t23_3d = real(t23_3dc);
  t33_3d = real(t33_3dc);

  // Convert from Cheb-space to physical space:
  for (int j = 0; j < Nx; j++) {
    for (int k = 0; k < Nz; k++) {
      u3d[slice(j * Nz + k, N + 1, Nx * Nz)] =
          idct(u3d[slice(j * Nz + k, N + 1, Nx * Nz)]);
      v3d[slice(j * Nz + k, N + 1, Nx * Nz)] =
          idct(v3d[slice(j * Nz + k, N + 1, Nx * Nz)]);
      w3d[slice(j * Nz + k, N + 1, Nx * Nz)] =
          idct(w3d[slice(j * Nz + k, N + 1, Nx * Nz)]);
      p3d[slice(j * Nz + k, N + 1, Nx * Nz)] =
          idct(p3d[slice(j * Nz + k, N + 1, Nx * Nz)]);
      t11_3d[slice(j * Nz + k, N + 1, Nx * Nz)] =
          idct(t11_3d[slice(j * Nz + k, N + 1, Nx * Nz)]);
      t12_3d[slice(j * Nz + k, N + 1, Nx * Nz)] =
          idct(t12_3d[slice(j * Nz + k, N + 1, Nx * Nz)]);
      t13_3d[slice(j * Nz + k, N + 1, Nx * Nz)] =
          idct(t13_3d[slice(j * Nz + k, N + 1, Nx * Nz)]);
      t22_3d[slice(j * Nz + k, N + 1, Nx * Nz)] =
          idct(t22_3d[slice(j * Nz + k, N + 1, Nx * Nz)]);
      t23_3d[slice(j * Nz + k, N + 1, Nx * Nz)] =
          idct(t23_3d[slice(j * Nz + k, N + 1, Nx * Nz)]);
      t33_3d[slice(j * Nz + k, N + 1, Nx * Nz)] =
          idct(t33_3d[slice(j * Nz + k, N + 1, Nx * Nz)]);
    }
  }
  x.resize(xval.size());
  z.resize(zval.size());
  std::cout << "xval: \n" << xval << '\n';
  std::cout << "zval: \n" << zval << '\n';
  std::cout << "x.size(): " << xval.size() << '\n';
  std::cout << "z.size(): " << zval.size() << '\n';
  for (int i = 0; i < xval.size(); i++) {
    x[i] = xval[i];
  }
  for (int i = 0; i < zval.size(); i++) {
    z[i] = zval[i];
  }

  outf.close();
  string filename(string("data/data_for_file18/OldroydB_") + flowType +
                  string("_beta_") + int2str(int(beta * 1000)) +
                  string("_We_") + int2str(int(We)) + string("_Re_") +
                  int2str(int(Re)) + string("_N_") + int2str(N) + string("_") +
                  string("velocity"));
  string filename2(string("data/data_for_file18/OldroydB_") + flowType +
                   string("_beta_") + int2str(int(beta * 1000)) +
                   string("_We_") + int2str(int(We)) + string("_Re_") +
                   int2str(int(Re)) + string("_N_") + int2str(N) + string("_") +
                   string("stress"));
  std::cout << "z.size(): " << z.size() << '\n';

  // Export to a vtk file:
  vtkExportCartesian3D(filename, x, y, z, u3d, v3d, w3d, p3d);
  vtkExportCartesianStress3D(filename2, x, y, z, t11_3d, t12_3d, t13_3d, t22_3d,
                             t23_3d, t33_3d);

  return 0;
}
