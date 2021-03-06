
#include <fstream> // To output data to files
#include <iostream>
#include <stdio.h>
#include <sis.hpp>


using namespace std;
#ifndef CD_T
typedef complex<double> Cd_t;
#endif
extern "C" {
    void sb04od_(char *REDUCE, char *TRANS, char *JOBD, int *M, int *N, double *A, int *LDA, double *B, int *LDB, double *C, int *LDC, double *D, int *LDD, double *E, int *LDE, double *F, int *LDF, double *SCALE, double *DIF, double *P, int *LDP, double *Q, int *LDQ, double *U, int *LDU, double *V, int *LDV, int* IWORK, double* DWORK, int* LDWORK, int* INFO);
}



//#define DEBUG

#include <lyap.h>

int main()
{
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
    vector<int> iwork;
    iwork.resize(M+N+6);
    vector<double> dwork;
    dwork.resize(1);
    int ldwork = -1;
    int info;
    cout << __LINE__ << endl
         << flush;

    // call this to estimate workspace
    sb04od_(&reduce, &trans, &jobd, &M, &N, A.data(), &ldA, B.data(), &ldB, C.data(), &ldC, D.data(), &ldD, E.data(), &ldE, F.data(), &ldF, &scale, &dif, P.data(), &ldP, Q.data(), &ldQ, U.data(), &ldU, V.data(), &ldV, &iwork[0], &dwork[0], &ldwork, &info);


    ldwork = dwork[0];
    cout << "ldwork: "<< ldwork << endl;
    dwork.resize(ldwork);
    sb04od_(&reduce, &trans, &jobd, &M, &N, A.data(), &ldA, B.data(), &ldB, C.data(), &ldC, D.data(), &ldD, E.data(), &ldE, F.data(), &ldF, &scale, &dif, P.data(), &ldP, Q.data(), &ldQ, U.data(), &ldU, V.data(), &ldV, &iwork[0], &dwork[0], &ldwork, &info);
    cout << __LINE__ << endl
         << flush;

    cout << "Scale = " << scale << endl;
    cout << "A = " << A << endl;
    cout << "B = " << B << endl;
    cout << "C = " << C << endl;
    cout << "E = " << E << endl;
    cout << "F = " << F << endl;
    cout << "P = " << P << endl;
    cout << "Q = " << Q << endl;
    cout << __LINE__ << endl
         << flush;

    cout << "First equation: " << A1 * C - F * B1 - scale * C1 << endl;
    cout << "Second equation: " << D1 * C - F * E1 - scale * F1 << endl;
    cout << __LINE__ << endl
         << flush;

    char jobvl = 'N';           // Don't compute left evecs
    char jobvr = 'V';           // Compute right evecs
    std::complex<double> wkopt; // Eistimate optimum workspace
    std::complex<double> *work; // allocate optimum workspace

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> Ac(1000, 1000), Bc(1000, 1000), vl(1000, 1000), vr(1000, 1000), eigenvalues_temp(1000, 1), alpha_temp(1000, 1), beta_temp(1000, 1), alpha(1000, 1), beta(1000, 1);

    cout << __LINE__ << endl
         << flush;

    Ac = Eigen::MatrixXcd::Random(1000, 1000);
    Bc = Eigen::MatrixXcd::Random(1000, 1000);
    //cout << "Ac = " << Ac << endl;
    //cout << "Bc = " << Bc << endl;
    int ldac = Ac.outerStride();
    int ldbc = Bc.outerStride();

    // vl : left evecs, vr: right evecs.
    int ldvl = vl.outerStride();
    cout << "ldvl = " << ldvl << endl;
    int ldvr = vr.outerStride();
    double rwork[8 * 1000];

    char sort = 'S';
    int sdim = 0;
    lwork = -1;
    vector<Cd_t> workc;
    workc.resize(1);
    bool *bwork;
    bwork = new bool[1000];

    zgges_(&jobvl, &jobvr, &sort, &criteria_, &N, Ac.data(), &ldac, Bc.data(), &ldbc, &sdim, alpha.data(), beta.data(), vl.data(), &ldvl, vr.data(), &ldvr, &workc[0], &lwork, &rwork[0], bwork, &info);

    cout << workc[0];
    lwork = int(real(workc[0]));
    workc.resize(lwork);

    zgges_(&jobvl, &jobvr, &sort, &criteria_, &N, Ac.data(), &ldac, Bc.data(), &ldbc, &sdim, alpha.data(), beta.data(), vl.data(), &ldvl, vr.data(), &ldvr, &workc[0], &lwork, &rwork[0], bwork, &info);

   // cout << "Alpha: " << alpha << endl;
   // cout << "beta: " << beta << endl;
    cout << "info = " << info << endl;
    delete bwork;

    Eigen::MatrixXcd AA(4,4), BB(4,4), QQ(4,4);
    AA = Eigen::MatrixXcd::Random(4,4);
    BB = Eigen::MatrixXcd::Random(4, 4);
    QQ = Eigen::MatrixXcd::Random(4, 4);
    Eigen::MatrixXcd X1 = dlyap(AA,BB,QQ);
    // Check norm:
    cout << "Testing dlyap: " << (AA*X1*BB - X1 + QQ).norm() << endl;

    Eigen::MatrixXcd X2 = lyap(AA,BB,Q);
    // Check norm:
    cout << "Testing lyap: " << (AA*X2 + X2*BB + Q).norm() << endl;
    //Testing OrdQz:

    OrdQz qz(AA,BB);
    cout << "qz info : " << qz.info << endl;
    cout << "Test qz 1: " << (qz.VSL*qz.S*qz.VSR.adjoint() - AA).norm() << endl;
    cout << "Test qz 2: " << (qz.VSL * qz.T * qz.VSR.adjoint() - BB).norm() << endl;
    return 0;
}