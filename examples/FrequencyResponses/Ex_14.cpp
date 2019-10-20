/// \file Ex_14.cpp
/// \brief Trying MJ's suggestions

#define SIS_USE_LAPACK
#include <fstream>
#include <iostream>
#include <sis.hpp>

using namespace std;
typedef complex<double> Cd_t;
typedef valarray<complex<double> > Vcd_t;
int main()
{
    using namespace sis;
    int bre;
    // Number of Chebyshev polynomials
    N = 127;
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
    double kz = 1.0;
    double omega = -0.385;
    double k2 = kx * kx + kz * kz;
    double k4 = k2 * k2;
    complex<double> ii(0.0, 1.0);
    Linop<double> Delta, Delta2, Dy, D0(0);
    LinopMat<complex<double> > Lmat(2, 2), Mmat(2, 2);

    Delta.n = 2;
    Delta.set();
    Delta.coef << 1.0, 0.0, -k2;
    D0.coef << 1.0;

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
    BcMat<Cd_t> lbcs(3, 2), rbcs(3, 2);
    lbcs.L << D0, 0.0, //
        Dy, 0.0,       //
        0.0, D0;
    rbcs.L << D0, 0.0, //
        Dy, 0.0,       //
        0.0, D0;
    lbcs.eval.setConstant(-1.0);
    rbcs.eval.setConstant(1.0);
    SingularValueDecomposition<std::complex<double> > svds;
    Cd_t iiomega  = ii*omega;
    LinopMat<Cd_t> A(2, 2), B(2, 3), C(3, 2);
    A = ((iiomega * Mmat) - Lmat);

    B << -ii * kx * Dy, -k2, -ii * kz * Dy, //
        ii * kz, 0.0, -ii * kx;
    C << ii * kx * Dy / k2, -ii * kz / k2, //
        1.0, 0.0,                          //
        ii * kz * Dy / k2, ii * kx / k2;
    svds.compute(A, B, C, lbcs, rbcs, 5 * (N + 1));
    std::cout << "eigenvalue: " << svds.eigenvalues << "\n";
    ofstream outf;
    
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> out_to_file(num_vals,2);
    //out_to_file.col(0) = svds.eigenvalues.real();
    //out_to_file.col(1) = svds.eigenvalues.imag();
    //outf << out_to_file;
    //outf.close();
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
       // ii*kx, 2.0 * ii * kx * Uy + Dy, ii*kz, Delta;
    /*
  Old bcs:
    Lmat.BcVec[0] = bc_dir;
    Lmat.BcVec[1] = bc_dir;
    Lmat.BcVec[2] = bc_dir;
    Lmat.BcVec[3] = bc_p;
  */

    // std::cout << "done in " << __LINE__ << '\n';
    num_vals = 4 * (N + 1);
    A.resize(4, 4);
    A = ((iiomega * Mmat) - Lmat);
    B.resize(4, 3);
    C.resize(3, 4);
    B << 1.0, 0.0, 0.0, //
        0.0, 1.0, 0.0,  //
        0.0, 0.0, 1.0,  //
        0.0, 0.0, 0.0;
    C << 1.0, 0.0, 0.0, 0.0, //
        0.0, 1.0, 0.0, 0.0,  //
        0.0, 0.0, 1.0, 0.0;  //
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
    svds.compute(A, B, C, Lbc, Rbc, Lbc, Rbc, 15 * (N + 1));
    cout << "In method 2: " << endl;
    cout << svds.eigenvalues.block(0, 0, 10, 1) << endl;
    exit(1);
    /* eigs.compute(Lmat, Mmat, num_vals, Lbc);
    eigs.keepConverged();
    eigs.sortByLargestReal();
    cout << "Eigenvalues in method 2: \n"
         << eigs.eigenvalues << "\n"
         << flush;
    num_vals = eigs.eigenvalues.size();
    outf.open("data/Ex_15_2.txt");
    out_to_file.resize(num_vals, 2);
    out_to_file.col(0) = eigs.eigenvalues.real();
    out_to_file.col(1) = eigs.eigenvalues.imag();
    outf << out_to_file;
    outf.close();
    exit(1);
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
    cout << "Eigenvalues in method 3: \n"
         << eigs.eigenvalues << "\n"
         << flush;
    num_vals = eigs.eigenvalues.size();
    outf.open("data/Ex_15_3.txt");
    out_to_file.resize(num_vals, 2);
    out_to_file.col(0) = eigs.eigenvalues.real();
    out_to_file.col(1) = eigs.eigenvalues.imag();
    outf << out_to_file;
    outf.close();*/
    return 0;
}
