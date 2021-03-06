
#ifndef CD_T
typedef complex<double> Cd_t;
#endif
extern "C"
{
    void zgees_(char *JOBVS, char *SORT, bool (*SELCTG)(Cd_t a), int *N, Cd_t *A, int *LDA, int *SDIM, Cd_t *W, Cd_t *VS, int *LDVS, Cd_t *WORK, int *LWORK, double *RWORK, bool *BWORK, int *INFO);
}

extern "C"
{
    void ztgsen_(int* IJOB,bool*  WANTQ, bool* WANTZ, bool* SELECT, int*  N, Cd_t* A,int* LDA, Cd_t* B,int* LDB, Cd_t* ALPHA, Cd_t* BETA, Cd_t* Q, int* LDQ,Cd_t* Z,int* LDZ, int* M, double* PL, double* PR, double* DIF, Cd_t* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO);
}
extern "C"
{
    bool criteria_(Cd_t alpha, Cd_t beta)
    {
        if (abs(alpha / beta) > 1e10 || (abs(beta) == 0.0) ) {
            return false;
        }
        else{
            return true;
        }
    };
}

extern "C"
{
    void zgges_(char *JOBVSL, char *JOBVSR, char *SORT, bool (*SELCTG)(Cd_t a, Cd_t b), int *N, Cd_t *A, int *LDA, Cd_t *B, int *LDB,
                int *SDIM, Cd_t *ALPHA, Cd_t *BETA, Cd_t *VSL, int *LDVSL, Cd_t *VSR, int *LDVSR, Cd_t *WORK,
                int *LWORK, double *RWORK, bool *BWORK, int *INFO);
}

class OrdQz {
    public:
    Eigen::MatrixXcd S,T,VSL,VSR,alpha,beta;
    int info;
    int M;
    OrdQz() {}
    OrdQz(Eigen::MatrixXcd A, Eigen::MatrixXcd B)
    {
        compute(A,B);
    }

    void compute(Eigen::MatrixXcd A, Eigen::MatrixXcd B){
        char jobvl = 'V';           // Compute left evecs
        char jobvr = 'V';           // Compute right evecs
        std::complex<double> wkopt; // Eistimate optimum workspace

        int r = A.rows();
        int N = r;


        S.resize(r,r); S = A;
        T.resize(r, r); T = B;
        VSL.resize(r,r);
        VSR.resize(r, r);
        alpha.resize(r,1);
        beta.resize(r, 1);

        int ldac = A.outerStride();
        int ldbc = B.outerStride();

        // vl : left evecs, vr: right evecs.
        int ldvl = VSL.outerStride();
        int ldvr = VSR.outerStride();
        double rwork[8 * r];

        char sort = 'S';
        int sdim = 0;
        int lwork = -1;
        vector<Cd_t> workc;
        workc.resize(1);
        bool bwork[N];
        
        bwork[0] = true;
        zgges_(&jobvl, &jobvr, &sort, &criteria_, &N, S.data(), &ldac, T.data(), &ldbc, &sdim, alpha.data(), beta.data(), VSL.data(), &ldvl, VSR.data(), &ldvr, &workc[0], &lwork, &rwork[0], bwork, &info);

        //cout << workc[0];
        lwork = int(real(workc[0]));
        workc.resize(lwork);

        zgges_(&jobvl, &jobvr, &sort, &criteria_, &N, S.data(), &ldac, T.data(), &ldbc, &sdim, alpha.data(), beta.data(), VSL.data(), &ldvl, VSR.data(), &ldvr, &workc[0], &lwork, &rwork[0], bwork, &info);


        int ijob = 0;
        bool wantq = true;
        bool wantz = true;
        bool select[N];
        for (int i = 0; i < alpha.size(); i++) {
            select[i] = criteria_(alpha(i, 0), beta(i, 0));
            cout << select[i] << endl;
        }
        int bre;
        cin >> bre;
        Eigen::MatrixXcd S1 = S, T1 = T;
        
        double pr;
        double pl;
        double dif;
        lwork = -1;
        int liwork = -1;
        vector<int> iwork;
        iwork.resize(1);
        ztgsen_(&ijob, &wantq, &wantz, select, &N, S.data(), &ldac, T.data(), &ldbc, alpha.data(), beta.data(), VSL.data(), &ldvl, VSR.data(), &ldvr, &M, &pl, &pr, &dif, &workc[0], &lwork, &iwork[0], &liwork, &info);
        lwork = int(real(workc[0]));
        workc.resize(lwork);
        liwork = iwork[0];
        iwork.resize(liwork);
        ztgsen_(&ijob, &wantq, &wantz, select, &N, S.data(), &ldac, T.data(), &ldbc, alpha.data(), beta.data(), VSL.data(), &ldvl, VSR.data(), &ldvr, &M, &pl, &pr, &dif, &workc[0], &lwork, &iwork[0], &liwork, &info);
        cout << "info : " << info << endl;
        cin >> bre;
        cout << "M : " << M << endl;
        cin >> bre;
        // cout << "Alpha: " << alpha << endl;
        // cout << "beta: " << beta << endl;
    //    for (int i = 0; i < N; i++)
    //    {
    //        cout << "bwork["<< i << "] = " << bwork[i] << endl;
    //    }
    //    int bre;
    //    cin >> bre;
        //delete bwork;
    }

};

bool forZgees(Cd_t a){
    return true;
}

/// \brief Solves the Lyapunov/ sylvester equation: AX + XB + Q = 0.
Eigen::MatrixXcd lyap(Eigen::MatrixXcd A, Eigen::MatrixXcd B, Eigen::MatrixXcd Q) {
    Eigen::MatrixXcd X, X1; // to be computed.
    X1.setConstant(Cd_t(0.0,0.0));
    
    // Compute the lower triangular Schur form of A. This is done by computing the Schur form of A' and then taking a conjugate transpose. Schur form is computed using LAPACK's zgees.
    Eigen::MatrixXcd T, Aad;
    Aad = A;
    Aad.adjointInPlace();
    T = Aad;
    #ifdef DEBUG
    cout << "A: " << A << endl;
    cout << "Aad: " << Aad << endl;
    #endif
    char jobvs = 'V';
    char sort = 'S';
    int N = A.rows();
//    cout << "N = " << N << endl;
    X1.resize(N,N);
    X.resize(N, N);
    int lda = Aad.outerStride();
 //   cout << "lda = " << lda << endl;
    int sdim = 1;
    int ldvs = N;
    Eigen::MatrixXcd W(N,1), Vs(N,N);
    int lwork = -1;
    vector<Cd_t> workc;
    workc.resize(1);
    bool *bwork;
    bwork = new bool[N];
    double rwork[8 * N];
    int info;
    zgees_(&jobvs, &sort, &forZgees, &N, T.data(), &lda, &sdim, W.data(), Vs.data(), &ldvs, &workc[0], &lwork, &rwork[0], bwork, &info);
    lwork = int(real(workc[0]));
    workc.resize(lwork);
    zgees_(&jobvs, &sort, &forZgees, &N, T.data(), &lda, &sdim, W.data(), Vs.data(), &ldvs, &workc[0], &lwork, &rwork[0], bwork, &info);
    // see that the norm is zero
    delete bwork;
    #ifdef DEBUG
    cout << "Error in norm: " << (Aad - Vs*T*Vs.adjoint()).norm() << endl;
    cout << "Error in norm: " << (A - Vs * T.adjoint() * Vs.adjoint()).norm() << endl;
    #endif
    // Hence T is represented by its adjoint.
    Eigen::MatrixXcd Tad = T.adjoint(), T2 = B, Us(N,N);


    // Now compute the schur decomposition of B:
    zgees_(&jobvs, &sort, &forZgees, &N, T2.data(), &lda, &sdim, W.data(), Us.data(), &ldvs, &workc[0], &lwork, &rwork[0], bwork, &info);
    #ifdef DEBUG
    cout << "Error in norm: " << (B - Us * T2 * Us.adjoint()).norm() << endl;
    #endif
    // Now loop through to get transformed X. Note that Adash is Tad and Bdash is T2 in Bartels and Stewarts algorithm
    int p = N, q = N;
    //cout << "Done" << endl
     //    << flush;
    //cout << "Tad: " << Tad << endl;
    //cout << "T2: " << T2 << endl;
    Eigen::MatrixXcd Q1 = Vs.adjoint() * Q * Us;
    //cout << "Done" << endl
     //   << flush;
     //   int bre;
    for (int k = 0; k < p; k++) {
        for (int l = 0; l < q; l++) {
            Cd_t sum1 = Cd_t(0.0, 0.0), sum2 = Cd_t(0.0, 0.0);
            for (int j = 0; j < k ; j++) {
                sum1 += (Tad(k, j) * X1(j, l));
            }
            for (int i = 0; i < l ; i++) {
                sum2 += (X1(k, i) * T2(i, l));
            }
            X1(k, l) = (-Q1(k, l) - sum1 - sum2) / (Tad(k, k) + T2(l, l));
        }
    }
    X = Vs*X1*Us.adjoint();

    #ifdef DEBUG
    cout << "Error in norm: " << (Tad*X1 + X1*T2 + Q1).norm() << endl;
    cout << "Error: " << (Tad * X1 + X1 * T2 + Q1) << endl;
    cout << "Error in norm: " << (A * X + X * B + Q).norm() << endl;
    cout << "Error: " << (A * X + X * B + Q) << endl;
#endif
    return X;
}

/// \brief Solves the Lyapunov/ sylvester equation: AXB - X + Q = 0.
Eigen::MatrixXcd dlyap(Eigen::MatrixXcd A, Eigen::MatrixXcd B, Eigen::MatrixXcd Q)
{
    Eigen::MatrixXcd X, X1; // to be computed.
    X1.setConstant(Cd_t(0.0, 0.0));

    // Compute the lower triangular Schur form of A. This is done by computing the Schur form of A' and then taking a conjugate transpose. Schur form is computed using LAPACK's zgees.
    Eigen::MatrixXcd T, Aad;
    Aad = A;
    Aad.adjointInPlace();
    T = Aad;
    //cout << "A: " << A << endl;
    //cout << "Aad: " << Aad << endl;

    char jobvs = 'V';
    char sort = 'S';
    int N = A.rows();
    //cout << "N = " << N << endl;
    X1.resize(N, N);
    X.resize(N, N);
    int lda = Aad.outerStride();
    //cout << "lda = " << lda << endl;
    int sdim = 1;
    int ldvs = N;
    Eigen::MatrixXcd W(N, 1), Vs(N, N);
    int lwork = -1;
    vector<Cd_t> workc;
    workc.resize(1);
    bool *bwork;
    bwork = new bool[N];
    double rwork[8 * N];
    int info;
    zgees_(&jobvs, &sort, &forZgees, &N, T.data(), &lda, &sdim, W.data(), Vs.data(), &ldvs, &workc[0], &lwork, &rwork[0], bwork, &info);
    lwork = int(real(workc[0]));
    workc.resize(lwork);
    zgees_(&jobvs, &sort, &forZgees, &N, T.data(), &lda, &sdim, W.data(), Vs.data(), &ldvs, &workc[0], &lwork, &rwork[0], bwork, &info);
    // see that the norm is zero
    delete bwork;
#ifdef DEBUG
    cout << "Error in norm: " << (Aad - Vs * T * Vs.adjoint()).norm() << endl;
    cout << "Error in norm: " << (A - Vs * T.adjoint() * Vs.adjoint()).norm() << endl;
#endif
    // Hence T is represented by its adjoint.
    Eigen::MatrixXcd Tad = T.adjoint(), T2 = B, Us(N, N);

    // Now compute the schur decomposition of B:
    zgees_(&jobvs, &sort, &forZgees, &N, T2.data(), &lda, &sdim, W.data(), Us.data(), &ldvs, &workc[0], &lwork, &rwork[0], bwork, &info);
#ifdef DEBUG
    cout << "Error in norm: " << (B - Us * T2 * Us.adjoint()).norm() << endl;
#endif
    // Now loop through to get transformed X. Note that Adash is Tad and Bdash is T2 in Bartels and Stewarts algorithm
    int p = N, q = N;
    //cout << "Done" << endl
    //     << flush;
    //cout << "Tad: " << Tad << endl;
    //cout << "T2: " << T2 << endl;
    Eigen::MatrixXcd Q1 = Vs.adjoint() * Q * Us;
    //cout << "Done" << endl
    //    << flush;
    //int bre;
    for (int k = 0; k < p; k++) {
        for (int l = 0; l < q; l++) {
            Cd_t sum1 = Cd_t(0.0, 0.0);
            for (int j = 0; j < k+1; j++) {
                for ( int i = 0; i < l+1; i++) {
                    sum1 += (Tad(k, j) * X1(j, i) * T2(i,l));
                }
            }
            //for (int i = 0; i < l; i++)
            //{
             //   sum2 += (X1(k, i) * T2(i, l));
            //}
            X1(k, l) = (-Q1(k, l) - sum1) / (Tad(k, k) * T2(l, l) - Cd_t(1.0,0.0));
        }
    }
    X = Vs * X1 * Us.adjoint();

#ifdef DEBUG
    cout << "Error in norm: " << (Tad * X1* T2 - X1 + Q1).norm() << endl;
    cout << "Error: " << (Tad * X1 * T2 - X1 + Q1) << endl;
    cout << "Error in norm: " << (A * X * B - X + Q).norm() << endl;
    cout << "Error: " << (A * X * B - X + Q) << endl;
#endif
    return X;
}

