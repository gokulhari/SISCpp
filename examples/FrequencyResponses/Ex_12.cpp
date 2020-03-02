/// \file Ex_12.cpp
/// \brief Solving for the singular values, power spectral density and the \f$\mathcal{H}_\infty\f$ norm of the linearized Navier stokes equations. We reproduce Figure 4.10 in \cite schmid2012stability using spectral integration with the linearized Navier-Stokes equations in primitive variables.
///  
///
/// \image html pics/Ex_12.svg <!--
/// --> width=500cm height=500cm
///
/// The we plot the power spectral density:
/// \image html pics/Ex_12_2.svg <!--
/// --> width=500cm height=500cm
///
/// Lastly notice that the largest singular value is somewhere near \f$-0.5\f$; the exact value can be calculated using the fast algorithm by Bruinsma and Steinbuch \cite BRUINSMA1990287 that is implemented in our codes.


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

int main()
{
    using namespace sis;
    int bre;
    // Number of Chebyshev polynomials
    N = 31;
    sis_setup();

    // Length of domain:
    double kz = 1;
    double kx = 1;
    Eigen::VectorXd omval(20);
    Eigen::MatrixXd data(20,3); 
    Eigen::MatrixXd data2(20,2); 
    omval = Eigen::VectorXd::LinSpaced(20, -2, 0);

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
    BcMat<std::complex<double> > Lbc(4, 4), Rbc(4, 4), bc(8, 4);

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
    bc.L << Lbc.L, //
        Rbc.L;
    bc.eval << Lbc.eval, //
        Rbc.eval;

    SingularValueDecomposition<std::complex<double> > svd;

    B.resize(4, 3);
    C.resize(3, 4);
    B << 1.0, 0.0, 0.0, //
        0.0, 1.0, 0.0,  //
        0.0, 0.0, 1.0,  //
        0.0, 0.0, 0.0;
    C << 1.0, 0.0, 0.0, 0.0, //
        0.0, 1.0, 0.0, 0.0,  //
        0.0, 0.0, 1.0, 0.0;  //

    for (int k = 0; k < 20; k++)
    {
        complex<double> iiomega = ii * omval(k);
        double k2 = kx * kx + kz * kz;
        double k4 = k2 * k2;
        Linop<double> Delta(2), Delta2(4);
        LinopMat<complex<double> > Lmat(4, 4), Mmat(4, 4);

        Delta.coef << 1.0, 0.0, -k2;
        Delta2.coef << 1.0, 0.0, -2 * k2, 0.0, k4;

        Mmat << 1.0, 0.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, 0.0,     //
            0.0, 0.0, 1.0, 0.0,     //
            0.0, 0.0, 0.0, 0.0 * Delta;

        Lmat << (-ii * kx * U) + (Delta / Re), -Uy, 0.0, -ii * kx, //
            0.0, (-ii * kx * U) + (Delta / Re), 0.0, -Dy,             //
            0.0, 0.0, (-ii * kx * U) + (Delta / Re), -ii * kz,     //
            ii * kx, Dy, ii * kz, 0.0;

        A.resize(4, 4);
        A = ((iiomega * Mmat) - Lmat);

        svd.compute(A, B, C, Lbc, Rbc, Lbc, Rbc, 15 * (N + 1));

        // Store first two eigenvalues.
        cout << "First two singular values: " << svd.eigenvalues[0] << " "
             << svd.eigenvalues[1] << "\n";

        data(k,0) = omval(k); 
        data(k,1) = svd.eigenvalues[0].real();
        data(k,2) = svd.eigenvalues[1].real();

        // Compute power spectral density:
        data2(k,0) = omval(k);
        Cd_t psd;
        psd = sqrt(svd.PowerSpectralDensity(A, B, C, Lbc, Rbc, Lbc, Rbc));
        cout << "Power spectral density for omega = " << omval(k) << " is " << real(psd) << endl;
        data2(k,1) = real(psd);
    }
    // Export to a file:
    ofstream outfile("data/Ex12.txt");
    outfile << data;
    outfile.close();

    outfile.open("data/Ex_12_2.txt");
    outfile << data2;
    outfile.close();

    return 0;
}
