/// \file Ex_12.cpp
/// \brief Example with evolution form LNS.

#define EIGEN_USE_BLAS
#define SIS_USE_LAPACK
#include <fstream>
#include <iostream>
#include <sis.hpp>
#include <string>

using namespace std;
typedef complex<double> Cd_t;
typedef valarray<complex<double> > Vcd_t;
complex<double> ii(0.0, 1.0);

int main() {
  using namespace sis;
  int bre;
  // Number of Chebyshev polynomials
  N = 63;
  sis_setup();
  // Number of Eigenvalues to compute:
  int num_vals = 40;

  Linop<double> D2(2), D1(1), D0(0); // Linops with highest order 2, 1 and 0.

  D2.coef << 1.0, 0.0, 0.0; // for (1.0 * D2 + 0.0 * D1 + (0.0) )v
  D1.coef << 1.0, 0.0;
  D0.coef << 1.0;

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
  kx[2] = -1.0;
  kx[3] = -1.0;

  kz[0] = 1.0;
  kz[1] = -1.0;
  kz[2] = 1.0;
  kz[3] = -1.0;

  omval[0] = -1.0;
  omval[1] = -1.0;
  omval[2] = 1.0;
  omval[3] = 1.0;
  double omega = -0.385;

  omval = omega * omval;

  valarray<double> y, Uy, U, Uyy(N + 1);
  // Set in cheb-point
  setChebPts(y);

  // Velocity and derivative for PP flow
  U = 1.0 - pow(y, 2.0);
  Uy = -2.0 * y;
  Uyy = -2.0;
  double Re = 2000.0;

  Linop<double> Dy(1);
  Dy.coef << 1.0, 0.0;


  SingularValueDecomposition<std::complex<double> > svd;

  ofstream outf;

  outf.open("data/Ex_12.txt");

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
      pvec(Cd_t(0.0, 0.0), (N + 1) * 4);
  valarray<double> u3d((N + 1) * Nx * Nz), v3d((N + 1) * Nx * Nz),
      w3d((N + 1) * Nx * Nz), p3d((N + 1) * Nx * Nz);

  for (int k = 0; k < 4; k++) {
    complex<double> iiomega = ii * omega;
    double k2 = kx[k] * kx[k] + kz[k] * kz[k];
    double k4 = k2 * k2;
    Linop<double> Delta(2), Delta2(4);
    LinopMat<complex<double> > Lmat(2, 2), Mmat(2, 2);

    Delta.coef << 1.0, 0.0, -k2;
    Delta2.coef << 1.0, 0.0, -2 * k2, 0.0, k4;

    Mmat << Delta, 0.0, //
        0.0, 1.0;

        Lmat << (-ii * kx[k] * U * Delta) + (ii * kx[k] * Uyy) + (Delta2 / Re), 0.0, //
            (-ii * kz[k] * Uy), (-ii * kx[k] * U) + (Delta / Re);

    LinopMat<Cd_t> A(2,2), B(2,3), C(3,2);
    A = ((iiomega * Mmat) - Lmat);

    B << -ii * kx[k] * Dy, -k2, -ii * kz[k] * Dy, //
        ii * kz[k], 0.0, -ii * kx[k];
    C << ii * kx[k] * Dy / k2, -ii * kz[k] / k2, //
        1.0, 0.0,                          //
        ii * kz[k] * Dy / k2, ii * kx[k] / k2;

    /*GeneralizedEigenSolver<complex<double> > eigs;

    eigs.compute(Lmat, Mmat, (N + 1) * Lmat.r, bc);
    eigs.removeInf();
    eigs.sortByLargestReal();
    ofstream outf;
    outf.open(string("data/data_for_file15/Newt") + string("_Re_") +
              int2str(int(Re)) + string("_N_") + int2str(N) + string("_") +
              int2str(k) + string(".txt"));
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> out_to_file(
        eigs.eigenvalues.size(), 3);
    out_to_file.col(0) = eigs.eigenvalues.real();
    out_to_file.col(1) = eigs.eigenvalues.imag();
    out_to_file.col(2) = eigs.MPorNot.cast<double>();
    outf << out_to_file;
    outf.close();*/
    BcMat<Cd_t> lbcs(3,2),rbcs(3,2);
    lbcs.L << D0, 0.0,//
              D1, 0.0,//
              0.0, D0;
    rbcs.L << D0, 0.0,//
              D1, 0.0,//
              0.0, D0;
    lbcs.eval.setConstant(-1.0);
    rbcs.eval.setConstant(1.0);
    svd.compute(A, B, C, lbcs, rbcs, 15 * (N + 1));
    std::cout << "eigenvalue: " << svd.eigenvalues[0] << "\n";
    num_vals = svd.eigenvalues.size();



    // Store first two eigenvalues.
    outf << kx[k] << " " << kz[k] << " " << svd.eigenvalues[0].real() << " "
         << svd.eigenvalues[1].real() << "\n";

    cout << kx[k] << " " << kz[k] << " " << svd.eigenvalues[0].real() << " "
         << svd.eigenvalues[1].real() << "\n"<< flush;

    // Store only half the eigenvectors, as the remaining correspond to the
    // optimal input.
    ChebfunMat<Cd_t> phi(2,1),output(3,1);
    phi(0,0) = svd.eigenvectors(0,0); //
    phi(1,0) = svd.eigenvectors(1,0);

    output = C(phi);

    uvec[slice(k,N+1,Nz)] =
        svd.eigenvalues[0].real() *
        output[0].v;

    vvec[slice(k,N+1,Nz)] =
        svd.eigenvalues[0].real() *
        output[1].v;

    wvec[slice(k,N+1,Nz)] =
        svd.eigenvalues[0].real() *
        output[2].v;
        
  }
  Eigen::VectorXd xval, zval;

  zval = Eigen::VectorXd::LinSpaced(100, -7.8, 7.8);
  xval = Eigen::VectorXd::LinSpaced(100, 0, 12.7);
  Nx = 100;
  Nz = 100;
  Vcd_t u3dc(Cd_t(0.0, 0.0), (N + 1) * Nx * Nz),
      v3dc(Cd_t(0.0, 0.0), (N + 1) * Nx * Nz),
      w3dc(Cd_t(0.0, 0.0), (N + 1) * Nx * Nz),
      p3dc(Cd_t(0.0, 0.0), (N + 1) * Nx * Nz);

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
      }
    }
  }
  u3d = real(u3dc);
  v3d = real(v3dc);
  w3d = real(w3dc);
  p3d = real(p3dc);

  x.resize(xval.size());
  z.resize(zval.size());
  //std::cout << "xval: \n" << xval << '\n';
  //std::cout << "zval: \n" << zval << '\n';
  //std::cout << "x.size(): " << xval.size() << '\n';
  //std::cout << "z.size(): " << zval.size() << '\n';
  for (int i = 0; i < xval.size(); i++) {
    x[i] = xval[i];
  }
  for (int i = 0; i < zval.size(); i++) {
    z[i] = zval[i];
  }

  outf.close();
  string filename("data/Ex_12");
  //std::cout << "z.size(): " << z.size() << '\n';

  // Export to a vtk file:
  vtkExportCartesian3D(filename, x, y, z, u3d, v3d, w3d, p3d);

  return 0;
}
