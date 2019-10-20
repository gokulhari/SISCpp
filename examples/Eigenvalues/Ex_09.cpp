/// \file Ex_09.cpp
/// \brief Finds eigenvalues for streamwise-constant linearized Navier-Stokes
/// equations, in two ways:
///   1. By solving for the generalized eigenvalue problem:
/// \f[
/// \left[
/// \begin{array}{cc}
/// \Delta & 0\\
/// 0 & I
/// \end{array}
/// \right]\partial_t \left[
/// \begin{array}{c}
/// v\\
/// \eta
/// \end{array}
/// \right]\;=\;\left[
/// \begin{array}{cc}
/// \frac{1}{Re}\Delta^2 & 0\\
/// -i\,k_z\,U' & \frac{1}{Re}\Delta
/// \end{array}
/// \right]\left[
/// \begin{array}{c}
/// v \\
/// \eta
/// \end{array}
/// \right]
/// \f]
/// with boundary conditions \f$v(\pm 1) \;=\; v'(\pm 1) \;=\; \eta(\pm 1) \;=\;
/// 0\f$
///
///   2. By solving the generalized eigenvalue problem:
/// \f[
/// \left[
/// \begin{array}{cccc}
/// 1 & 0 & 0 & 0 \\
/// 0 & 1 & 0 & 0 \\
/// 0 & 0 & 1 & 0 \\
/// 0 & 0 & 0 & 0 \\
/// \end{array}
/// \right]\partial_t \left[
/// \begin{array}{c}
/// u \\
/// v \\
/// w \\
/// p
/// \end{array}
/// \right] \;= \; \left[
/// \begin{array}{cccc}
/// \frac{1}{Re}\Delta & -U' & 0 & 0 \\
/// 0 & \frac{1}{Re}\Delta & 0 & -\partial_y \\
/// 0 & 0 & \frac{1}{Re}\Delta & -ik_z \\
/// 0 & \partial_y & ik_z & 0 \\
/// \end{array}
/// \right]\left[
/// \begin{array}{c}
/// u \\
/// v \\
/// w \\
/// p
/// \end{array}
/// \right]
/// \f]
/// with boundary conditions \f$u(\pm 1) \;=\; v(\pm 1) \;=\; w(\pm 1) \;=\;
/// p(-1) \;=\; 0\f$
///
/// Eigenvalues returned by both methods should have nature, (in terms of wave
/// speed and stability.), but may differ in value.
///
///
/// Note from the above form that the eigenvalues are simply the eigenvalues of
/// \f$ \Delta^2/Re\f$ and \f$\Delta/Re\f$, both of whom are real. Method 2
/// approximates the solution more accrately than method 1, as imaginary
/// parts in method 2 are zero to machine precision. Moreover, eigenfunctions
/// are directly related to velocities and pressure and not to the vorticity.

#define SIS_USE_LAPACK
#include <fstream>
#include <iostream>
#include <sis.hpp>

using namespace std;

int main() {
  using namespace sis;
  int bre;
  // Number of Chebyshev polynomials
  N = 91;
  sis_setup();

  // Number of Eigenvalues to compute:
  int num_vals = 40;

  valarray<double> y, Uy, U, Uyy(N + 1);
  // Set in cheb-point
  setChebPts(y);

  // Velocity and derivative for PP flow
  U = 1.0 - pow(y, 2.0);
  Uy = -2.0 * y;
  Uyy = -2.0;
  double Re = 2000;
  double kx = 1.0;
  double kz = 0.0;
  double k2 = kx * kx + kz * kz;
  double k4 = k2 * k2;
  complex<double> ii(0.0, 1.0);
  Linop<double> Delta, Delta2, Dy;
  LinopMat<complex<double> > Lmat(2, 2), Mmat(2, 2);

  Delta.n = 2;
  Delta.set();
  Delta.coef << 1.0, 0.0, -k2;

  Delta2.n = 4;
  Delta2.set();
  Delta2.coef << 1.0, 0.0, -2 * k2, 0.0, k4;

  Dy.n = 1;
  Dy.set();
  Dy.coef << 1.0, 0.0;

  // Linop<complex<double> > tempOp;
  // Eigen::ArrayXcd tempArray;
  // tempOp = -ii * kx * U + Delta;

  // For first type:
  Lmat << (-ii * kx * U * Delta) + (ii * kx * Uyy) + (Delta2 / Re), 0.0, //
      (-ii * kz * Uy), (-ii * kx * U) + (Delta / Re);

  Mmat << Delta, 0.0, //
      0.0, 1.0;
  BcMat<std::complex<double> > bc(6, 2);
  bc.L << 1.0, 0.0, //
      0.0, 1.0,     //
      Dy, 0.0,      //
      1.0, 0.0,     //
      0.0, 1.0,     //
      Dy, 0.0;
  bc.eval << -1.0, -1.0, //
      -1.0, -1.0,        //
      -1.0, -1.0,        //
      1.0, 1.0,          //
      1.0, 1.0,          //
      1.0, 1.0;
  
  
  GeneralizedEigenSolver<complex<double> > eigs;
  eigs.compute(Lmat, Mmat, num_vals, bc);
  eigs.keepConverged();
  eigs.sortByLargestReal();
  num_vals = eigs.eigenvalues.size();
  cout << "Eigenvalues in method 1: \n" << eigs.eigenvalues << "\n";
  ofstream outf;
  outf.open("data/Ex_15_1.txt");
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> out_to_file(num_vals,
                                                                    2);
  out_to_file.col(0) = eigs.eigenvalues.real();
  out_to_file.col(1) = eigs.eigenvalues.imag();
  outf << out_to_file;
  outf.close();
  std::cout << "Done 1 :::::::" << '\n';
  // For second type:
  Lmat.resize(4, 4);
  Mmat.resize(4, 4);
  /*
  Old Mmat:
  Mmat << 1.0, 0.0, 0.0, 0.0, //
      0.0, 1.0, 0.0, 0.0,     //
      0.0, 0.0, 1.0, 0.0,     //
      0.0, 0.0, 0.0, 0.0;
*/
  Mmat << 1.0, 0.0, 0.0, 0.0, //
      0.0, 1.0, 0.0, 0.0,     //
      0.0, 0.0, 1.0, 0.0,     //
      0.0, 0.0, 0.0, 0.0 * Delta;

  //
  Lmat << (-ii * kx * U) + (Delta / Re), -Uy, 0.0, -ii * kx, //
      0.0, (-ii * kx * U) + (Delta / Re), 0.0, -Dy,          //
      0.0, 0.0, (-ii * kx * U) + (Delta / Re), -ii * kz,     //
      ii * kx, Dy, ii * kz, 0.0;
     // ii*kx, 2.0*ii * kx * Uy + Dy, ii*kz, Delta; //Mj's suggestion
  /*
  Old bcs:
    Lmat.BcVec[0] = bc_dir;
    Lmat.BcVec[1] = bc_dir;
    Lmat.BcVec[2] = bc_dir;
    Lmat.BcVec[3] = bc_p;
  */
  BcMat<std::complex<double> > Lbc(8, 4);
  Lbc.L << 1.0, 0.0, 0.0, 0.0, //
      1.0, 0.0, 0.0, 0.0,      //
      0.0, 1.0, 0.0, 0.0,      //
      0.0, 1.0, 0.0, 0.0,      //
      0.0, Dy, 0.0, 0.0,       //
      0.0, Dy, 0.0, 0.0,       //
      0.0, 0.0, 1.0, 0.0,      //
      0.0, 0.0, 1.0, 0.0;
  Lbc.eval << 1.0, 0.0, 0.0, 0.0, //
      -1.0, 0.0, 0.0, 0.0,        //
      0.0, 1.0, 0.0, 0.0,         //
      0.0, -1.0, 0.0, 0.0,        //
      0.0, 1.0, 0.0, 0.0,         //
      0.0, -1.0, 0.0, 0.0,        //
      0.0, 0.0, 1.0, 0.0,         //
      0.0, 0.0, -1.0, 0.0;
  // std::cout << "done in " << __LINE__ << '\n';
  num_vals = 4 * (N + 1);
  eigs.compute(Lmat, Mmat, num_vals, Lbc);
  eigs.keepConverged();
  eigs.sortByLargestReal();
  cout << "Eigenvalues in method 2: \n" << eigs.eigenvalues << "\n" << flush;
  num_vals = eigs.eigenvalues.size();
  outf.open("data/Ex_15_2.txt");
  out_to_file.resize(num_vals, 2);
  out_to_file.col(0) = eigs.eigenvalues.real();
  out_to_file.col(1) = eigs.eigenvalues.imag();
  outf << out_to_file;
  outf.close();
return 0;
  // Lastly, solve with conventional boundary conditions.
  Mmat.resize(4, 4);
  Mmat << 1.0, 0.0, 0.0, 0.0, //
      0.0, 1.0, 0.0, 0.0,     //
      0.0, 0.0, 1.0, 0.0,     //
      0.0, 0.0, 0.0, 0.0;

  Lbc.resize(7, 4);
  Lbc.L << 1.0, 0.0, 0.0, 0.0, //
      1.0, 0.0, 0.0, 0.0,      //
      0.0, 1.0, 0.0, 0.0,      //
      0.0, 1.0, 0.0, 0.0,      //     //      //
      0.0, 0.0, 1.0, 0.0,      //
      0.0, 0.0, 1.0, 0.0,      //
      0.0, 0.0, 0.0, 1.0;

  Lbc.eval << 1.0, 0.0, 0.0, 0.0, //
      -1.0, 0.0, 0.0, 0.0,        //
      0.0, 1.0, 0.0, 0.0,         //
      0.0, -1.0, 0.0, 0.0,        //
      0.0, 0.0, 1.0, 0.0,         //
      0.0, 0.0, -1.0, 0.0,        //
      0.0, 0.0, 0.0, -1.0;
  num_vals = 4 * (N + 1);
  ind = 2;
  eigs.compute(Lmat, Mmat, num_vals, Lbc);
  eigs.keepConverged();
  eigs.sortByLargestReal();
  cout << "Eigenvalues in method 3: \n" << eigs.eigenvalues << "\n" << flush;
  num_vals = eigs.eigenvalues.size();
  outf.open("data/Ex_15_3.txt");
  out_to_file.resize(num_vals, 2);
  out_to_file.col(0) = eigs.eigenvalues.real();
  out_to_file.col(1) = eigs.eigenvalues.imag();
  outf << out_to_file;
  outf.close();
  //return 0;
}
