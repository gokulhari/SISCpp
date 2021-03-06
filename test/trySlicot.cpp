
#include <fstream> // To output data to files
#include <iostream>
#include <stdio.h>
#include <sis.hpp>


using namespace std;

extern "C" {
    void sb04od_(char *REDUCE, char *TRANS, char *JOBD, int *M, int *N, double *A, int *LDA, double *B, int *LDB, double *C, int *LDC, double *D, int *LDD, double *E, int *LDE, double *F, int *LDF, double *SCALE, double *DIF, double *P, int *LDP, double *Q, int *LDQ, double *U, int *LDU, double *V, int *LDV, int* IWORK, double* DWORK, int* LDWORK, int* INFO);
}


typedef complex<double> Cd_t;

int main()
{
    cout<< "Hi! I am here!!" << endl;
    Eigen::MatrixXd A(4,4), B(4,4), C(4,4), D(4,4), E(4,4), F(4,4), P(4,4), Q(4,4), U(4,4), V(4,4), A1(4,4), B1(4,4), C1(4,4), D1(4,4), E1(4,4), F1(4,4);
    A = Eigen::MatrixXd::Random(4,4);
    B = Eigen::MatrixXd::Random(4, 4);
    C = Eigen::MatrixXd::Random(4, 4);
    D = Eigen::MatrixXd::Random(4, 4);
    E = Eigen::MatrixXd::Random(4, 4);
    F = Eigen::MatrixXd::Random(4, 4);

    A1 = A;
    B1 = B;
    C1 = C;
    D1 = D;
    E1 = E;
    F1 = F;


    int sizeA = A.rows();
    int lwork = -1; // set lwork to -1 to estimate workspace.
    char reduce = 'R';
    char trans = 'N';
    char jobd = 'N';
    int M = 4;
    int N = 4;
    int ldA = A.outerStride(); // ld for leading dimension
    int ldB = B.outerStride();
    int ldC = C.outerStride();
    int ldD = D.outerStride();
    int ldE = E.outerStride();
    int ldF = F.outerStride();
    int ldP = P.outerStride();
    int ldQ = Q.outerStride();
    int ldU = U.outerStride();
    int ldV = V.outerStride();
    double scale = 1.0;
    double dif;
    int* iwork;
    iwork = new int[M+N+6];
    double* dwork;
    dwork = new double[1];
    int ldwork = -1;
    int info;

    // call this to estimate workspace
    sb04od_(&reduce, &trans, &jobd, &M, &N, A.data(), &ldA, B.data(), &ldB, C.data(), &ldC, D.data(), &ldD, E.data(), &ldE, F.data(), &ldF, &scale, &dif, P.data(), &ldP, Q.data(), &ldQ, U.data(), &ldU, V.data(), &ldV, iwork, dwork, &ldwork, &info);


    ldwork = dwork[0];
    cout << "ldwork: "<< ldwork << endl;
    sb04od_(&reduce, &trans, &jobd, &M, &N, A.data(), &ldA, B.data(), &ldB, C.data(), &ldC, D.data(), &ldD, E.data(), &ldE, F.data(), &ldF, &scale, &dif, P.data(), &ldP, Q.data(), &ldQ, U.data(), &ldU, V.data(), &ldV, iwork, dwork, &ldwork, &info);

cout << "Scale = " << scale << endl;
cout << "A = " << A << endl;
cout << "B = " << B << endl;
cout << "C = " << C << endl;
cout << "E = " << E << endl;
cout << "F = " << F << endl;
cout << "P = " << P << endl;
cout << "Q = " << Q << endl;

cout << "First equation: " << A1*C - F*B1 - scale*C1 << endl;
cout << "Second equation: " << D1*C - F*E1 - scale * F1 << endl;
return 0;
}