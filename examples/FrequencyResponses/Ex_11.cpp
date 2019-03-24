/// \brief example with discriptor approach.

#define EIGEN_USE_BLAS
#define SIS_USE_LAPACK
#include <fstream>
#include <iostream>
#include <sis.hpp>
#include <string>

using namespace std;

complex<double> ii(0.0, 1.0);

int main() {
  using namespace sis;
  int bre;
  // Number of Chebyshev polynomials
  N = 91;
  sis_setup();
  // Number of Eigenvalues to compute:
  int num_vals = 40;

  // Nyquist frequency
  double kz_Nyq = 2.0;
  double kx_Nyq = 2.0;
  // Number of fourier modes in span-wise direction:
  int Nz = 4;
  int Nx = 4;
  // Length of domain:
  double Lz = Nz * M_PI / (kz_Nyq);
  double Lx = Nx * M_PI / (kx_Nyq);
  valarray<double> kz(Nz);
  valarray<double> kx(Nx);
  valarray<double> z(Nz);
  valarray<double> x(Nx);
  for (int i = 0; i < Nz / 2; i++) {
    kz[i] = i * 1.0;
  }
  for (int i = Nz / 2; i < Nz; i++) {
    kz[i] = i - Nz;
  }
  kz = kz * (2.0 * M_PI) / Lz;

  for (int i = 0; i < Nx / 2; i++) {
    kx[i] = i * 1.0;
  }
  for (int i = Nx / 2; i < Nx; i++) {
    kx[i] = i - Nx;
  }
  kx = kx * (2.0 * M_PI) / Lx;

  double delz = Lz / Nz;
  double delx = Lz / Nx;

  for (int i = 0; i < Nz; i++) {
    z[i] = double(i) * delz;
  }
  for (int i = 0; i < Nx; i++) {
    x[i] = double(i) * delx;
  }

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

  LinopMat<std::complex<double> > A, B(2, 3), C(3, 2);
  BcMat<std::complex<double> > Lbc(4, 4), Rbc(4, 4);

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
  Lbc.eval.setConstant(-1.0);
  Rbc.eval.setConstant(1.0);

  SingularValueDecomposition<std::complex<double> > svd;

  ofstream outf;

  double omega = -0.385;

  B.resize(4, 3);
  C.resize(3, 4);
  B << 1.0, 0.0, 0.0, //
      0.0, 1.0, 0.0,  //
      0.0, 0.0, 1.0,  //
      0.0, 0.0, 0.0;
  C << 1.0, 0.0, 0.0, 0.0, //
      0.0, 1.0, 0.0, 0.0,  //
      0.0, 0.0, 1.0, 0.0;  //
  outf.open("data/file2.txt");

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

  valarray<complex<double> > u3dc((N + 1) * Nx * Nz / 2),
      v3dc((N + 1) * Nx * Nz / 2), w3dc((N + 1) * Nx * Nz / 2),
      p3dc((N + 1) * Nx * Nz / 2);
  valarray<double> u3d((N + 1) * Nx * Nz), v3d((N + 1) * Nx * Nz),
      w3d((N + 1) * Nx * Nz), p3d((N + 1) * Nx * Nz);

  for (int j = 0; j < Nx; j++) {
    for (int k = 0; k < Nz / 2; k++) {
      complex<double> iiomega = ii * omega;
      double k2 = kx[j] * kx[j] + kz[k] * kz[k];
      double k4 = k2 * k2;
      Linop<double> Delta(2), Delta2(4);
      LinopMat<complex<double> > Lmat(4, 4), Mmat(4, 4);

      Delta.coef << 1.0, 0.0, -k2;
      Delta2.coef << 1.0, 0.0, -2 * k2, 0.0, k4;

      Mmat << 1.0, 0.0, 0.0, 0.0, //
          0.0, 1.0, 0.0, 0.0,     //
          0.0, 0.0, 1.0, 0.0,     //
          0.0, 0.0, 0.0, 0.0 * Delta;

      Lmat << (-ii * kx[j] * U) + (Delta / Re), -Uy, 0.0, -ii * kx[j], //
          0.0, (-ii * kx[j] * U) + (Delta / Re), 0.0, -Dy,             //
          0.0, 0.0, (-ii * kx[j] * U) + (Delta / Re), -ii * kz[k],     //
          ii * kx[j], Dy, ii * kz[k], 0.0;

      A.resize(4, 4);
      A = ((iiomega * Mmat) - Lmat);
      svd.compute(A, B, C, Lbc, Rbc, Lbc, Rbc, 15 * (N + 1));
      std::cout << "eigenvalue: " << svd.eigenvalues[0] << "\n";
      num_vals = svd.eigenvalues.size();

      // Store first two eigenvalues.
      outf << kx[j] << " " << kz[k] << " " << svd.eigenvalues[0].real() << " "
           << svd.eigenvalues[1].real() << "\n";

      std::cout << "(" << j << "," << k << ")" << '\n';
      cout << kx[j] << " " << kz[k] << " " << svd.eigenvalues[0].real() << " "
           << svd.eigenvalues[1].real() << "\n";
      u3dc[slice(j * Nz / 2 + k, N + 1, Nx * Nz / 2)] =
          // svd.eigenvalues[0].real() *
          svd.eigenvectorsMat[0][0].v;

      v3dc[slice(j * Nz / 2 + k, N + 1, Nx * Nz / 2)] =
          // svd.eigenvalues[0].real() *
          svd.eigenvectorsMat[0][1].v;

      w3dc[slice(j * Nz / 2 + k, N + 1, Nx * Nz / 2)] =
          // svd.eigenvalues[0].real() *
          svd.eigenvectorsMat[0][2].v;

      p3dc[slice(j * Nz / 2 + k, N + 1, Nx * Nz / 2)] =
          // svd.eigenvalues[0].real() *
          svd.eigenvectorsMat[0][3].v;
    }
  }

  // take inverse fourier transforms
  for (int i = 0; i < N + 1; i++) {
    u3d[slice(i * Nx * Nz, Nx * Nz, 1)] =
        ifft2_cs(u3dc[slice(i * Nx * Nz / 2, Nx * Nz / 2, 1)], Nx, Nz);
    v3d[slice(i * Nx * Nz, Nx * Nz, 1)] =
        ifft2_cs(v3dc[slice(i * Nx * Nz / 2, Nx * Nz / 2, 1)], Nx, Nz);
    w3d[slice(i * Nx * Nz, Nx * Nz, 1)] =
        ifft2_cs(w3dc[slice(i * Nx * Nz / 2, Nx * Nz / 2, 1)], Nx, Nz);
    p3d[slice(i * Nx * Nz, Nx * Nz, 1)] =
        ifft2_cs(p3dc[slice(i * Nx * Nz / 2, Nx * Nz / 2, 1)], Nx, Nz);
  }

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
    }
  }

  outf.close();
  string filename("data/file2vec_noscale");
  std::cout << "z.size(): " << z.size() << '\n';

  // Export to a vtk file:
  vtkExportCartesian3D(filename, x, y, z, u3d, v3d, w3d, p3d);

  return 0;
}
