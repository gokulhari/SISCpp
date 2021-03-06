
#include <fstream> // To output data to files
#include <iostream>
#include <stdio.h>
#include <sis.hpp>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
using namespace std;
#ifndef CD_T
typedef complex<double> Cd_t;
#endif
//#define DEBUG

#include <lyap.h>
#include <intWts.h>

int main()
{
    int size = 4;
    Eigen::MatrixXcd AA(size, size), BB(size, size), QQ(size, size);
    AA = Eigen::MatrixXcd::Random(size, size);
    BB = Eigen::MatrixXcd::Random(size, size);
    QQ = Eigen::MatrixXcd::Random(size, size);
    Eigen::MatrixXcd X1 = dlyap(AA, BB, QQ);
    // Check norm:
    cout << "Testing dlyap: " << (AA * X1 * BB - X1 + QQ).norm() << endl;

    Eigen::MatrixXcd X2 = lyap(AA, BB, QQ);
    // Check norm:
    cout << "Testing lyap: " << (AA * X2 + X2 * BB + QQ).norm() << endl;
    //Testing OrdQz:

    OrdQz qz(AA, BB);

    
    cout << "qz info : " << qz.info << endl;
    cout << "Test qz 1: " << (qz.VSL * qz.S * qz.VSR.adjoint() - AA).norm() << endl;
    cout << "Test qz 2: " << (qz.VSL * qz.T * qz.VSR.adjoint() - BB).norm() << endl;

    // Test intWts:
    sis::N = 91;
    sis::sis_setup();
    sis::ChebfunMat<Cd_t> f(1,1);
    f << 1.0 - sis::y*sis::y;

    f.p2c();

    Eigen::MatrixXcd fm(sis::N+1,1);
    for (int i = 0 ;i < sis::N+1; i++){
        fm(i,0) = f(0,0).v[i];
    }
    
    Eigen::MatrixXd Iw = intWts(sis::N);
    cout << "Testing intWts: " << fm.adjoint() * Iw * fm << ", Actual value: " << 16.0/15.0 <<endl;

    // Testing square roots module.
    Eigen::MatrixXd sqrtIw = Iw.sqrt();
    cout << "Testing matrix square root, norm: " << (sqrtIw*sqrtIw - Iw).norm() << endl; 
     return 0;
}