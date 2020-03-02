
#define SIS_USE_LAPACK
#include <fstream>
#include <iostream>
#include <sis.hpp>

using namespace std;

int main()
{
    using namespace sis;
    int bre;
    // Number of Chebyshev polynomials
    N = 63;
    sis_setup();


    valarray<double> y, Uy, U, Uyy(N + 1), temp(N + 1);
    // Set in cheb-point
    setChebPts(y);

    // Velocity and derivative for PP flow
    U = 1.0 - pow(y, 2.0);
    Uy = -2.0 * y;
    Uyy = -2.0;
    double Re = 2000;
    double kx = 1.0;
    double kz = 1.0;
    double k2 = kx * kx + kz * kz;
    double k4 = k2 * k2;
    complex<double> ii(0.0, 1.0);
    Linop<double> Delta, Delta2, Dy;
    LinopMat<complex<double> > Lmat(4, 4), Mmat(4, 4);
    Delta.n = 2;
    Delta.set();
    Delta.coef << 1.0, 0.0, -k2;

    Delta2.n = 4;
    Delta2.set();
    Delta2.coef << 1.0, 0.0, -2 * k2, 0.0, k4;

    Dy.n = 1;
    Dy.set();
    Dy.coef << 1.0, 0.0;

    GeneralizedEigenSolver<complex<double> > eigs;
    
    Lmat.resize(4, 4);
    Mmat.resize(4, 4);
   
    Mmat << 1.0, 0.0, 0.0, 0.0, //
        0.0, 1.0, 0.0, 0.0,     //
        0.0, 0.0, 1.0, 0.0,     //
        0.0, 0.0, 0.0, 0.0 * Delta;

    Lmat << (-ii * kx * U) + (Delta / Re), -Uy, 0.0, -ii * kx, //
        0.0, (-ii * kx * U) + (Delta / Re), 0.0, -Dy,          //
        0.0, 0.0, (-ii * kx * U) + (Delta / Re), -ii * kz,     //
        ii * kx, Dy, ii * kz, 0.0;
    
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
     std::cout << "done in " << __LINE__ << '\n';
    eigs.compute(Lmat, Mmat, 1000000, Lbc);
    eigs.keepConverged();
    eigs.sortByLargestReal();
    cout << "Eigenvalues in method 2: \n"
         << eigs.eigenvalues << "\n"
         << flush;
    return 0;
}
