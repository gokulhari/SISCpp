// function out = IntWts(N)
Eigen::MatrixXd intWts(int N){
    // IntWts: Integration weight matrix is generated for two functions in a
    // Chebyshev basis.
    // Working procedure:
    // \int_{-1}^{1} T_m(y) Tn(y) dy = \int_{0}^{\pi} cos(m \theta) cos(n
    // \theta) sin (\theta) d \theta.
    // This integral is given by:
    // z = -((1+(-1)^(m+n))/(-1 + m + n)) + (1 + (-1)^(m-n))/(1 + m - n) +
    // (1 + (-1)^(m-n))/(1 - m + n) + (1 + (-1)^(m+n))/(1 + m + n);
    // where if any term has denominator as zero, the term is treated as zero.
Eigen::MatrixXd out(N+1,N+1);
out = Eigen::MatrixXd::Zero(N + 1, N + 1);
cout << __LINE__ << " " << __FILE__ << endl << flush;
for (int m = 0; m < N+1; m++) {
    for (int n = 0; n < N+1; n++){
        double z = 0;
        if ((m + n) != 1) {
            z += -((1.0+pow((-1),(m+n)))/(-1.0 + m + n)); 
        }
        if ((m-n) != -1) {
            z += (1.0 + pow((-1),(m-n)))/(1.0 + m - n);
        }
        if ((-m+n) != -1) {
            z +=  (1.0 + pow((-1),(m-n)))/(1.0 - m + n);
        }
        if ((m+n) != -1) {
            z +=  (1.0 + pow((-1),(m+n)))/(1.0 + m + n);
        }
        out(m,n) = z/4;
    }
}
cout << __LINE__ << " " << __FILE__ << endl
     << flush;
out(0, 0) = out(0, 0) / 4.0;
out.block(1, 0, N, 1) = out.block(1, 0, N, 1)/2.0;
out.block(0, 1, 1, N) = out.block(0, 1, 1, N) / 2.0;
//out(2: end, 1) = out(2: end, 1) /2;
//out(1,2:end) = out(1,2:end)/2;
return out;
}