/// \file sis.hpp
/// \mainpage
///
/// Written by Gokul Hariharan, harih020@umn.edu.
///
///
/// \section Introduction
///
/// Spectral integral suite in C++ (SISC++) is a generic header to solve
/// two-point-boundary-value problems in the system
/// representation. SISC++ can solve for linear differential equations,
/// and compute eigenvalues, and frequency responses of linear differential systems.
/// We started this in
/// order to produce programs for direct numerical simulations of viscoelastic
/// channel flows, then we extended this as a generic solver, following in
/// the lines of Chebfun, see http://www.chebfun.org/. Chebfun is based on a
/// Matlab for which you need to purchase a license.
///
/// SISC++ aims to provide a Chebfun-like interface in C++. Just as Chebfun
/// overloads functions and operations on linear matrices to differential
/// operators,
/// SISC++ uses <a href="http://eigen.tuxfamily.org/index.php?title=Main_Page">Eigen's </a> matrix representation in C++ to
/// classes for linear differential operators.
/// For instance, one would input a normal matrix using <a href="http://eigen.tuxfamily.org/index.php?title=Main_Page">Eigen </a> in the following
/// manner,
/// \code{.cpp}
/// using namespace Eigen;
/// MatrixXd A(2,2);
/// A << 1, 2, //
///      3, 4;
/// \endcode
/// In SISC++, you can overload  the linear block-matrix operator
/// \f{align}
/// L =  \left[ \begin{array}{cc}
/// \partial_y & y \\
/// y\,\partial_{yy} & \partial_{yy}/2
/// \end{array} \right],\nonumber
/// \f}
///  as
/// \code{.cpp}
/// using namespace sis;
/// valarray<double> y;
/// setChebPts(y);
/// Linop <double> Dy(1), Dyy(2);
/// Dy.coef << 1.0, 0.0;
/// Dyy.coef << 1.0, 0.0, 0.0;
/// LinopMat<double> L(2,2);
/// L << Dy, y,//
///      y*Dyy, Dyy/2.0;
/// \endcode
/// Following that, just as you can use the EigenSolver in <a href="http://eigen.tuxfamily.org/index.php?title=Main_Page">Eigen </a> for matrix
/// eigenvalues, you can use EigenSolver in SIS for linear block-matrix
/// operators.
/// \n \n An advantage with a code in C++ is
/// you do not need a Matlab license to use it.
/// You can also optimized a code in C++ for speed by
/// using proprietary compilers, like the Intel compiler / Cray compiler.
///
/// One major difference from Chebfun is that in SISC++, we do not provide for
/// automatic collocation. Automatic collocation is an extremely useful utility
/// in Chebfun that keeps increasing the number of basis functions until the
/// solution reaches machine precision. We did not incorporate this
/// as automatic collocation adds heavily to the computational expense, which is
/// a great liability for direct numerical simulations.
///
/// For most part classes and functions are quite intuitive, but you need to know
/// at least a little bit about C++. A good place to learn C++
/// is <a href="http://www.cplusplus.com/doc/tutorial/">here </a>.
/// You'll find it amazing if you have used <a href="http://eigen.tuxfamily.org/index.php?title=Main_Page">Eigen </a>.
/// As far as the algorithm is concerned, SISC++ is based on the recent
/// spectral-integral method for non-constant by
/// Du,
///  \cite DuSI , which is a
/// well-conditioned method, compared to conventional spectral-collocation / Tau methods.
/// Note that Chebfun too has a well-conditioned scheme, the ultraspherical
/// discretization \cite UltraS. \n \n
///
/// SISC++ is well-equiped to deal with incompressible hydrodynamic
/// eigenvalue and singular value problems for Newtonian and Viscoelastic
/// fluids, both in primitive variables and the evolution form. To the best of
/// our knowledge, solving for eigenvalues of incompressible flow problems in
/// Chebfun can currently be done only by using the evolution form of the
/// governing equations, As a virtue of spectral integration, SISC++ can solve
/// for incompressible flow eigenvalue problems directly in the discriptor form
/// (in primitive variables). This can potentially save a lot of time,
/// specifically for viscoelastic fluids whose algebraic manipulations for a
/// transformation to the evolution form can become
/// cumbersome. In addition we also provide tools to handle eigenvalue boundary
/// constraints (problems where the eigenvalue appears in the boundary
/// conditions), for eigenvalue problems involving fluid-fluid interfaces.
///
/// \section Prerequisites
///
/// SIS requires two free libraries to be pre-installed,
///   1. Intel MKL (https://software.intel.com/en-us/mkl/choose-download) and
///   2. Eigen (http://eigen.tuxfamily.org/index.php?title=Main_Page)
///
/// We strogly recommend downloading these from the respective websites and
/// installing them. For <a href="http://eigen.tuxfamily.org/index.php?title=Main_Page">Eigen's </a> in specific, we ask that Mac users refrain from
/// using brew to install either of these. Of course, if the user is familiar
/// about linking these libraries correctly while using C++, then the method of
/// installing these libraries should not matter.
///
/// You must also set the path for Intel MKL. For linux, you can do this by
/// saying
/// \code{.sh}
/// source  /opt/intel/compilers_and_libraries_2019/linux/mkl/bin/mklvars.sh
/// intel64 \endcode
/// **in every new terminal**. Else, you can add this to you ~/.bashrc file.
/// For other platforms, take a look at setting environment variables
/// <a
/// href="https://software.intel.com/en-us/articles/intel-math-kernel-library-intel-mkl-2019-getting-started">here</a>
///
/// No part of SIS uses anything that is OS specific,
/// nonetheless, However, I have tried this on Mac and Linux.
///
/// Most routines in SIS use default <a href="http://eigen.tuxfamily.org/index.php?title=Main_Page">Eigen's </a>'s routines for solutions and
/// eigenvalue solver. At the time of writing this program, <a href="http://eigen.tuxfamily.org/index.php?title=Main_Page">Eigen's </a> does not have
/// eigenvalue solver for complex generalized systems, of the form
///  \f$ L\,\phi= \lambda \, M \, \phi\f$. As SIS uses a
/// well-conditioned discretization, in most cases both \f$L\f$ and \f$M\f$ are
/// well conditioned, so either can be inverted while the other is
/// singular to compute the eigenvalues.
///
/// We provide an option to use the macro SIS_USE_LAPACK. This will use LAPACK's
/// complex generalized eigenvalue solver and also replace all other places of
/// the codes with LAPACK's counter-part. Intel MKL also has LAPACK in it,
/// so if you have linked Intel MKL correctly, LAPACK must be available in the
/// same path as Intel MKL (you will not have to do anything extra). However,
/// as we use fortran code in C++, gfortran needs to be linked during compilation (implying that
/// gfortran must be installed, in Mac say "brew install gfortran", and in
/// linux "sudo apt-get install gfortran" in the terminal). In summary is you
/// need to use SIS_USE_LAPACK, make sure that gfortran is installed.
///
///
///
/// \section Installation
///
/// SIS does not need any installation and can be directly used as a header. You
/// can copy sis.hpp file into /usr/local/includes/ in either Mac or Linux,
/// as this is a default search path for most C++ compilers.
///
/// \section secgets Getting Started
///
/// If <a href="http://eigen.tuxfamily.org/index.php?title=Main_Page">Eigen's </a> and Intel MKL are installed in their default locations, then you
/// can begin by creating two folders, named bin and data. Then open a terminal
/// in the current folder and say \code{.sh} make all \endcode This will compile
/// all the examples in the directory examples, and the executables will be
/// placed in the directory bin. Then run the executables, by saying, say
/// \code{.sh}
/// ./bin/Ex1
/// \endcode
/// in the terminal.
/// Data generated will be placed in the folder data.
/// To solve your own problems, go through the examples and modify them.
/// \section citeus Cite Us
/// If you use this work, please cite us:
///
#define lapack_complex_double std::complex<double>
#define lapack_complex_float std::complex<float>
#ifndef SIS_HPP
#define SIS_HPP
#endif
#ifndef _STL_VECTOR_H
#include <vector>
#endif
#ifndef _GLIBCXX_VALARRAY
#include <valarray>
#endif
#ifndef _GLIBCXX_IOSTREAM
#include <iostream>
#endif
#if !defined _GLIBCXX_CMATH || !defined _GLIBCXX_MATH_H
#include <cmath>
#endif
#ifndef _GLIBCXX_COMPLEX
#include <complex>
#endif
#ifndef _GLIBCXX_NUMERIC_LIMITS
#include <limits>
#endif
#ifndef _GLIBCXX_ALGORITHM
#include <algorithm>
#endif
#ifndef _MKL_DFTI_H_
#include <mkl_dfti.h>
#endif
#ifndef EIGEN_CORE_H
#include <eigen3/Eigen/Eigen>
#endif
#ifndef _TIME_H
#include <time.h>
#endif
#ifdef SIS_USE_FEAST
extern "C" {
void feastinit_(int *feastparam);
}
extern "C" {
void zfeast_gegv_(int *N, std::complex<double> *A, int *LDA,
                  std::complex<double> *B, int *LDB, int *feastparam,
                  double *epsout, int *loop, double *Emid, double *r, int *M0,
                  std::complex<double> *lambda, std::complex<double> *q,
                  int *mode, double *res, int *info);
}
#ifndef _STDIO_H_
#include <stdio.h>
#endif
#ifndef _STDLIB_H
#include <stdlib.h>
#endif
#endif
#ifdef SIS_USE_LAPACK
extern "C" {
void dgeev_(char *jobvl, char *jobvr, int *n, double *a, int *lda, double *wr,
            double *wi, double *vl, int *ldvl, double *vr, int *ldvr,
            double *work, int *lwork, double *rwork, int *info);
}
extern "C" {
void dggev_(char *jobvl, char *jobvr, int *n, double *a, int *lda, double *b,
            int *ldb, double *alphar, double *alphai, double *beta, double *vl,
            int *ldvl, double *vr, int *ldvr, double *wkopt, int *lwork,
            double *rwork, int *info);
}
extern "C" {
void zgeev_(char *jobvl, char *jobvr, int *n, std::complex<double> *a, int *lda,
            std::complex<double> *w, std::complex<double> *vl, int *ldvl,
            std::complex<double> *vr, int *ldvr, std::complex<double> *work,
            int *lwork, double *rwork, int *info);
}
extern "C" {
void zggev_(char *jobvl, char *jobvr, int *n, std::complex<double> *a, int *lda,
            std::complex<double> *b, int *ldb, std::complex<double> *alpha,
            std::complex<double> *beta, std::complex<double> *vl, int *ldvl,
            std::complex<double> *vr, int *ldvr, std::complex<double> *wkopt,
            int *lwork, double *rwork, int *info);
}
#ifndef _STDIO_H_
#include <stdio.h>
#endif
#ifndef _STDLIB_H
#include <stdlib.h>
#endif
#endif

int ind = 0;
// I define a few useful operator overloads to std::valarray that aren't
// inherently defined, like multiplying a complex number with a real valarray.
// I place these functions within the same namespace: std.
namespace std {

/// \brief Use this to make a complex valarray out of two real valarrays
template <class T>
std::valarray<std::complex<T> > dou2com(const std::valarray<T> &a,
                                        const std::valarray<T> &b) {
  std::valarray<std::complex<T> > temp;

  temp.resize(a.size());
  for (int i = 0; i < a.size(); i++) {
    temp[i] = std::complex<T>(a[i], b[i]);
  }
  return temp;
};

/// \brief real part of a complex valarray
template <class T> valarray<T> real(const valarray<complex<T> > &in) {
  valarray<T> temp;

  temp.resize(in.size());
  for (int i = 0; i < temp.size(); i++) {
    temp[i] = real(in[i]);
  }
  return temp;
};

/// \brief imaginary part of a complex valarray
template <class T> valarray<T> imag(const valarray<complex<T> > &in) {
  valarray<T> temp;
  temp.resize(in.size());
  for (int i = 0; i < temp.size(); i++) {
    temp[i] = imag(in[i]);
  }
  return temp;
};

/// \brief Multiplying a complex number with a real valarray
template <class T>
std::valarray<std::complex<T> > operator*(std::complex<T> left,
                                          std::valarray<T> right) {
  std::valarray<T> re(right.size()), im(right.size());
  re = real(left) * right;
  im = imag(left) * right;
  return std::dou2com(re, im);
}

/// \brief Multiplying a real valarray with a complex number.
template <class T>
std::valarray<std::complex<T> > operator*(const std::valarray<T> &left,
                                          const std::complex<T> &right) {
  std::valarray<T> re(left.size()), im(left.size());
  re = left * real(right);
  im = left * imag(right);
  return std::dou2com(re, im);
}

/// \brief Adding a complex number to a real valarray
template <class T>
std::valarray<std::complex<T> > operator+(std::complex<T> left,
                                          std::valarray<T> right) {
  std::valarray<T> re(right.size()), im(right.size());
  re = real(left) + right;
  im = imag(left);

  return std::dou2com(re, im);
}

/// \brief Multiplying a real valarray with a complex number.
template <class T>
std::valarray<std::complex<T> > operator+(const std::valarray<T> &left,
                                          std::complex<T> right) {
  std::valarray<T> re(left.size()), im(left.size());
  re.resize(left.size());
  im.resize(left.size());
  re = left + real(right);
  im = imag(right);
  return std::dou2com(re, im);
}
/// \brief Subtracting a real valarray from a complex number.
template <class T>
std::valarray<std::complex<T> > operator-(std::complex<T> left,
                                          std::valarray<T> right) {
  std::valarray<T> re(right.size()), im(right.size());
  re.resize(right.size());
  im.resize(right.size());
  re = real(left) - right;
  im = imag(left);
  return std::dou2com(re, im);
}

/// \brief Subtracting a complex number from a real valarray
template <class T>
std::valarray<std::complex<T> > operator-(const std::valarray<T> &left,
                                          const std::complex<T> &right) {
  std::valarray<T> re(left.size()), im(left.size());
  re = left - real(right);
  im = imag(right);
  return std::dou2com(re, im);
}

// \brief Multiplying a complex valarray with a real number
template <class T>
std::valarray<std::complex<T> >
operator*(T left, std::valarray<std::complex<T> > right) {
  std::valarray<T> re(right.size()), im(right.size());
  im = left * imag(right);
  re = left * real(right);
  return std::dou2com(re, im);
}

/// \brief Multiplying a real number with a complex valarray.
template <class T>
std::valarray<std::complex<T> >
operator*(const std::valarray<std::complex<T> > &left, const T &right) {
  std::valarray<T> re(left.size()), im(left.size());
  re = real(left) * right;
  im = imag(left) * right;
  return std::dou2com(re, im);
}

/// \brief Adding a complex valarray to a real number
template <class T>
std::valarray<std::complex<T> >
operator+(T left, std::valarray<std::complex<T> > right) {
  std::valarray<T> re(right.size()), im(right.size());
  re = left + real(right);
  im = imag(right);
  return std::dou2com(re, im);
}

/// \brief Adding a real number with a complex valarray.
template <class T>
std::valarray<std::complex<T> >
operator+(const std::valarray<std::complex<T> > &left, const T &right) {
  std::valarray<T> re(left.size()), im(left.size());
  re = real(left) + right;
  im = imag(left);
  return std::dou2com(re, im);
}
/// \brief Subtracting a real number from a complex valarray.
template <class T>
std::valarray<std::complex<T> >
operator-(T left, std::valarray<std::complex<T> > right) {
  std::valarray<T> re(right.size()), im(right.size());
  re = left - real(right);
  im = -imag(right);
  return std::dou2com(re, im);
}

/// \brief Subtracting a complex valarray from a real number
template <class T>
std::valarray<std::complex<T> >
operator-(const std::valarray<std::complex<T> > &left, const T &right) {
  std::valarray<T> re(left.size()), im(left.size());
  re = real(left) - right;
  im = imag(left);
  return std::dou2com(re, im);
}

/// \brief Multiplying a complex valarray to a real valarray:
template <class T>
std::valarray<std::complex<T> >
operator*(const std::valarray<std::complex<T> > &left,
          const std::valarray<T> &right) {
  std::valarray<T> re(right.size()), im(right.size());
  re = real(left) * right;
  im = imag(left) * right;
  return std::dou2com(re, im);
}

/// \brief Multiplying a real valarray to a complex valarray:
template <class T>
std::valarray<std::complex<T> >
operator*(const std::valarray<T> &left,
          const std::valarray<std::complex<T> > &right) {
  std::valarray<T> re(left.size()), im(left.size());
  re = left * real(right);
  im = left * real(right);
  return std::dou2com(re, im);
}

string int2str(int i) {
  string out;
  stringstream temp;
  temp << i;
  return temp.str();
};

template <class T>
valarray<complex<T> > pow(valarray<complex<T> > base, T power) {
  complex<T> power_ = complex<T>(power, 0.0);
  valarray<complex<T> > out = pow(base, power_);
  return out;
}

}; // namespace std

namespace sis {

#ifndef PI
#define PI 3.141592653589793
#endif
#define SIS_SINGULAR 1
#define SIS_SVD 0
#define SIS_SVD_LEFT 1
#define SIS_SVD_RIGHT 2
#define SIS_PHYS_SPACE 1
#define SIS_CHEB_SPACE 0
int N = 31; ///< Specifies number of Chebyshev polynomials, default N = 31.

#ifdef SIS_USE_FEAST
/// \brief All variables for feast is stored in this namespace
namespace feast {
int M0 = 10; /// number of eigenvalues to search
int info = 0;
std::complex<double> center(0.0, 0.0);
double radius = 10.0;
int feastparam[64];
void feast_init() { feastinit_(feastparam); }
void display() {
  if (info == 202) {
    std::cout << "Error: Problem with size of system N" << '\n';
  } else if (info == 201) {
    std::cout << "Error: Problem with subspace M0" << '\n';
  } else if (info == 200) {
    std::cout << "Error: Problem with Emin, Emax, Emid, r" << '\n';
  } else if (info == 6) {
    std::cout << "Warning: FEAST converges but subspace is not bi-orthogonal"
              << '\n';
  } else if (info == 5) {
    std::cout << "Warning: Only stochastic estimation of #eigenvalues returned "
                 "fpm(14)=2"
              << '\n';
  } else if (info == 4) {
    std::cout << "Warning: Only the subspace has been returned using fpm(14)=1"
              << '\n';
  } else if (info == 3) {
    std::cout << "Warning: Size of the subspace M0 is too small (M0<=M)"
              << '\n';
  } else if (info == 2) {
    std::cout << "Warning: No Convergence (#iteration loops>fpm(4))" << '\n';
  } else if (info == 1) {
    std::cout << "Warning: No Eigenvalue found in the search interval" << '\n';
  } else if (info == 0) {
    std::cout << "Error: Successful exit" << '\n';
  } else if (info == -1) {
    std::cout << "Error: Internal error for allocation memory" << '\n';
  } else if (info == -2) {
    std::cout << "Error: Internal error of the inner system solver in FEAST "
                 "predefined interfaces"
              << '\n';
  } else if (info == -3) {
    std::cout << "Error: Internal error of the reduced eigenvalue solver. "
              << "Possible cause for Hermitian problem: matrix B may not be "
                 "positive definite"
              << '\n';
  } else {
    std::cout << "Error: Problem with " << info - 100
              << "th parameter of feast parameter" << '\n';
  }
}
} // namespace feast
#endif

template <class T> class LinopMat;
template <class T> class BcMat;
template <class T> class MatGen;
template <class T> class SingularValueDecomposition;
template <class T> class Discretize;

#ifndef SIS_TYPE
#define SIS_TYPE double
#endif

std::valarray<std::complex<SIS_TYPE> > half_shift(N + 1), rev_half_shift(N + 1),
    yc(N + 1);
std::valarray<SIS_TYPE> y(N + 1);
Eigen::Matrix<std::complex<SIS_TYPE>, Eigen::Dynamic, Eigen::Dynamic> ycEigen;
Eigen::Matrix<SIS_TYPE, Eigen::Dynamic, Eigen::Dynamic> yEigen;

template <class T> std::valarray<std::complex<T> > fft(std::valarray<T> in1) {
  int bre;
  int n = in1.size();
  std::valarray<std::complex<T> > out(n);
  std::valarray<T> zvec(n);
  zvec = 0.0;
  out[0] = std::complex<T>(in1[0], 0.0);
  DFTI_DESCRIPTOR_HANDLE descriptor;
  MKL_LONG status;

  status = DftiCreateDescriptor(&descriptor, DFTI_DOUBLE, DFTI_REAL, 1,
                                n); // Specify size and precision
  // cout<<status<<endl;
  status = DftiSetValue(descriptor, DFTI_PLACEMENT,
                        DFTI_NOT_INPLACE); // Out of place FFT
  status = DftiSetValue(descriptor, DFTI_CONJUGATE_EVEN_STORAGE,
                        DFTI_COMPLEX_COMPLEX);

  // cout<<status<<endl;
  status = DftiCommitDescriptor(descriptor); // Finalize the descriptor
  // cout<<status<<endl;
  status = DftiComputeForward(descriptor, &in1[0],
                              &out[0]); // Compute the Forward FFT
  // cout<<status<<endl;
  status = DftiFreeDescriptor(&descriptor); // Free the descriptor
  // cout<<status<<endl;
  for (int j = 0; j < n / 2; j++)
    out[n / 2 + j] = std::conj(out[n / 2 - j]);

  return out;
};

template <class T>
std::valarray<std::complex<T> > fft(std::slice_array<T> in1) {
  return fft(std::valarray<T>(in1));
}

template <class T>
std::valarray<T> ifft_cs(std::valarray<std::complex<T> > in)
// For conjugate symmetric
{
  int bre;
  int n = in.size();
  std::valarray<T> out(2 * n);
  out[0] = real(in[0]);
  std::valarray<std::complex<T> > in1(std::complex<T>(0.0, 0.0), n + 1);
  in1[std::slice(0, n, 1)] = in;
  DFTI_DESCRIPTOR_HANDLE descriptor;
  MKL_LONG status;

  status = DftiCreateDescriptor(&descriptor, DFTI_DOUBLE, DFTI_REAL, 1,
                                2 * n); // Specify size and precision
  // std::cout<<status<<endl;
  status = DftiSetValue(descriptor, DFTI_PLACEMENT,
                        DFTI_NOT_INPLACE); // Out of place FFT
  status = DftiSetValue(descriptor, DFTI_CONJUGATE_EVEN_STORAGE,
                        DFTI_COMPLEX_COMPLEX);
  // cout<<status<<endl;
  status = DftiCommitDescriptor(descriptor); // Finalize the descriptor
  // cout<<status<<endl;
  status = DftiComputeBackward(descriptor, &in1[0],
                               &out[0]); // Compute the Backward FFT
  // cout<<status<<endl;
  status = DftiFreeDescriptor(&descriptor); // Free the descriptor
                                            // std::cout<<status<<endl;
  out = out / T(2 * n);                     // need to manually scale
  return out;
};

template <class T>
std::valarray<T> ifft_cs(std::slice_array<std::complex<T> > in) {
  return ifft_cs(std::valarray<std::complex<T> >(in));
}
template <class T> std::valarray<T> dct(const std::valarray<T> &x) {
  int n = x.size();
  int bre;
  std::valarray<T> y(2 * n), v(n);
  std::valarray<std::complex<T> > V(n);
  for (int i = 0; i < n; i++)
    y[i] = x[i];

  for (int i = n; i < 2 * n; i++)
    y[i] = x[2 * n - i - 1];

  for (int i = 0; i < n; i++)
    v[i] = y[2 * i];
  V = fft(v);
  return 2.0 * real(std::valarray<std::complex<T> >(half_shift * V)) / T(n);
};

template <class T> std::valarray<T> dct(const std::slice_array<T> &x) {
  return dct(std::valarray<T>(x));
};

template <class T> std::valarray<T> idct(const std::valarray<T> &u) {
  int n = u.size();
  int bre;
  std::valarray<T> temp_u(n + 1);
  std::valarray<T> x(n), v(n);
  temp_u = 0.0;
  temp_u[std::slice(0, n, 1)] = u;
  std::valarray<std::complex<T> > V(n);
  for (int i = 0; i < n; i++) {
    V[i] =
        std::complex<T>(temp_u[i], 0.0) - std::complex<T>(0.0, temp_u[n - i]);
  }
  V = std::complex<T>(0.5, 0.0) * rev_half_shift * V;

  v = ifft_cs(std::valarray<std::complex<T> >(V[std::slice(0, n / 2, 1)]));

  for (int i = 0; i < (n + 1) / 2; i++)
    x[2 * i] = v[i];
  for (int i = (n + 1) / 2; i < n; i++)
    x[2 * n - 2 * i - 1] = v[i];
  return x * T(n);
};

template <class T> std::valarray<T> idct(const std::slice_array<T> &u) {
  return idct(std::valarray<T>(u));
};

template <class T>
std::valarray<std::complex<T> > dct(const std::valarray<std::complex<T> > &in) {
  int n = in.size() - 1;
  std::valarray<T> inr(n + 1), ini(n + 1);
  inr = real(in);
  ini = imag(in);
  inr = dct(inr);
  ini = dct(ini);
  return dou2com(inr, ini);
};

template <class T>
std::valarray<std::complex<T> >
dct(const std::slice_array<std::complex<T> > &in) {
  return dct(std::valarray<std::complex<T> >(in));
};

template <class T>
std::valarray<std::complex<T> >
idct(const std::valarray<std::complex<T> > &in) {
  int n = in.size() - 1;

  std::valarray<T> inr(in.size()), ini(in.size());
  inr = std::real(in);
  ini = std::imag(in);
  inr = idct(inr);
  ini = idct(ini);
  return dou2com(inr, ini);
};

template <class T>
std::valarray<std::complex<T> >
idct(const std::slice_array<std::complex<T> > &in) {
  return idct(std::valarray<std::complex<T> >(in));
};

/// \brief fft2 for a 2D matrix stored in the row-major format, Nx and Nz denote
/// dimensions in x and z. The return value is of size Nx * Nz / 2 of
/// complex type, which implicitly assumes conjugate symmetry and also that
/// the values at the Nyquist frequency are zero.
template <class T>
std::valarray<std::complex<T> > fft2(std::valarray<T> in1, int Nx, int Nz) {
  int bre;
  std::valarray<std::complex<T> > result(Nx * (Nz / 2 + 1)), out(Nx * Nz / 2);
  result[0] = std::complex<T>(in1[0], 0.0);
  DFTI_DESCRIPTOR_HANDLE descriptor;

  MKL_LONG status = 0;
  DFTI_DESCRIPTOR_HANDLE hand = 0;
  MKL_LONG n[2];
  n[0] = Nx;
  n[1] = Nz;
  status = DftiCreateDescriptor(&hand, DFTI_DOUBLE, DFTI_REAL, 2, n);
  status = DftiSetValue(hand, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
  status =
      DftiSetValue(hand, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
  MKL_LONG rs[3];
  rs[0] = 0;
  rs[1] = Nz;
  rs[2] = 1;
  status = DftiSetValue(hand, DFTI_INPUT_STRIDES, rs);
  MKL_LONG cs[3];
  cs[0] = 0;
  cs[1] = Nz / 2 + 1;
  cs[2] = 1;
  status = DftiSetValue(hand, DFTI_OUTPUT_STRIDES, cs);
  status = DftiCommitDescriptor(hand);
  status = DftiComputeForward(hand, &in1[0], &result[0]);
  DftiFreeDescriptor(&hand);
  for (int j = 0; j < Nx; j++) {
    for (int k = 0; k < Nz / 2; k++) {
      out[j * (Nz / 2) + k] = result[j * (Nz / 2 + 1) + k];
    }
  }
  // Set values at Nyquist to zero:
  {
    int j = Nx / 2;
    for (int k = 0; k < Nz / 2; k++) {
      out[j * (Nz / 2) + k] = 0.0;
    }
  }
  //  for (int j = 0; j < Nx; j ++){
  //    for (int k = 0; k < Nz /2; k ++){
  //      std::cout << out[j*(Nz /2) + k] << " ";
  //    }
  //    std::cout << '\n';
  //  }
  return out;
};

/// \brief ifft2 for a 2D matrix stored in the row-major format, Nx and Nz
/// denote
/// dimensions in x and z. The input is of size (Nx x Nz / 2) of type complex
/// double, which implicitly assumes conjugate symmetry and also that
/// the values at the Nyquist frequency are zero.
/// The return value is of size Nx * Nz of
/// double type.
template <class T>
std::valarray<T> ifft2_cs(std::valarray<std::complex<T> > in1, int Nx, int Nz) {
  int bre;
  std::valarray<std::complex<T> > in(Nx * (Nz / 2 + 1));
  std::valarray<T> out(Nx * Nz);
  in = std::complex<T>(0.0, 0.0);
  for (int j = 0; j < Nx; j++) {
    for (int k = 0; k < Nz / 2; k++) {
      in[j * (Nz / 2 + 1) + k] = in1[j * (Nz / 2) + k];
    }
  }

  // std::cout << "showing the before ifft" << '\n';
  //  for (int j = 0; j < Nx; j ++){
  //      for (int k = 0; k < Nz /2 + 1; k ++){
  //      std::cout << in[j * (Nz/2 + 1) + k] << " " ;
  //    }
  //    std::cout << "\n" << '\n';
  //  }
  // std::cin >> bre;
  DFTI_DESCRIPTOR_HANDLE descriptor;

  MKL_LONG status = 0;
  DFTI_DESCRIPTOR_HANDLE hand = 0;
  MKL_LONG n[2];
  n[0] = Nx;
  n[1] = Nz;
  status = DftiCreateDescriptor(&hand, DFTI_DOUBLE, DFTI_REAL, 2, n);
  status = DftiSetValue(hand, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
  status =
      DftiSetValue(hand, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
  MKL_LONG cs[3];
  cs[0] = 0;
  cs[1] = Nz / 2 + 1;
  cs[2] = 1;
  status = DftiSetValue(hand, DFTI_INPUT_STRIDES, cs);
  MKL_LONG rs[3];
  rs[0] = 0;
  rs[1] = Nz;
  rs[2] = 1;
  status = DftiSetValue(hand, DFTI_OUTPUT_STRIDES, rs);
  status = DftiCommitDescriptor(hand);
  status = DftiComputeBackward(hand, &in[0], &out[0]);
  DftiFreeDescriptor(&hand);
  out = out / double(Nx * Nz);
  //  for (int j = 0; j < Nx; j ++){
  //    for (int k = 0; k < Nz; k ++){
  //      std::cout << out[j*(Nz) + k] << " ";
  //    }
  //    std::cout << '\n';
  //  }
  return out;
};

template <class T>
std::valarray<std::complex<T> > fft2(std::slice_array<T> in1, int Nx, int Nz) {
  return fft2(std::valarray<T>(in1), Nx, Nz);
}

template <class T>
std::valarray<T> ifft2_cs(std::slice_array<std::complex<T> > in, int Nx,
                          int Nz) {
  return ifft2_cs(std::valarray<std::complex<T> >(in), Nx, Nz);
}

/// \brief This function is useful to see size of Eigen matrices. Returns a
/// complex number, where the real part indicates the number of rows and
/// imaginary part the number of columns.

template <class T>
std::complex<T>
size(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &in) {
  return std::complex<T>(in.rows(), in.cols());
}

template <class T>
std::complex<T>
size(const Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> &in) {
  return std::complex<T>(in.rows(), in.cols());
}

/// \brief Prints a valarray to terminal.
template <class T> void disp(std::valarray<T> in) {
  for (int i = 0; i < in.size(); i++) {
    std::cout << in[i] << "\n";
  }
};

/// \brief This function sets points to evaluate a function so that a DCT will
/// give represent the same function in a Chebyshev basis
template <class T> void setChebPts(std::valarray<T> &in) {
  in.resize(N + 1);
  for (int i = 0; i < N + 1; i++)
    in[i] = cos(M_PI * (i + 0.5) / (N + 1.0));
};

/// \brief This function sets points to evaluate a function so that a DCT will
/// give represent the same function in a Chebyshev basis, overloaded to complex
/// type.
template <class T> void setChebPts(std::valarray<std::complex<T> > &in) {
  in.resize(N + 1);
  for (int i = 0; i < N + 1; i++)
    in[i] = std::complex<T>(cos(PI * (i + 0.5) / (N + 1.0)), 0.0);
};

/// \brief This function sets points to evaluate a function so that a DCT will
/// give represent the same function in a Chebyshev basis, overloaded to Eigen
/// array class.
template <class T> void setChebPts(Eigen::Array<T, Eigen::Dynamic, 1> &in) {
  in.resize(N + 1);
  for (int i = 0; i < N + 1; i++)
    in[i] = cos(PI * (i + 0.5) / (N + 1.0));
};

/// \brief This function sets points to evaluate a function so that a DCT will
/// give represent the same function in a Chebyshev basis, overloaded to complex
/// Eigen array class.
template <class T>
void setChebPts(Eigen::Array<std::complex<T>, Eigen::Dynamic, 1> &in) {
  in.resize(N + 1);
  for (int i = 0; i < N + 1; i++)
    in[i] = std::complex<T>(cos(PI * (i + 0.5) / (N + 1.0)), 0.0);
};

void sis_setup() {
  std::complex<SIS_TYPE> ii;
  ii = std::complex<SIS_TYPE>(0.0, 1.0);
  half_shift.resize(N + 1);
  rev_half_shift.resize(N + 1);
  yEigen.resize(N + 1, 1);
  ycEigen.resize(N + 1, 1);
  y.resize(N + 1);
  yc.resize(N + 1);
  setChebPts(y);
  setChebPts(yc);
  for (int i = 0; i < N + 1; i++) {
    half_shift[i] =
        exp(-ii * std::complex<SIS_TYPE>(M_PI * i / (2.0 * (N + 1)), 0.0));
    rev_half_shift[i] =
        exp(ii * std::complex<SIS_TYPE>(M_PI * i / (2.0 * (N + 1)), 0.0));
  }

  for (int i = 0; i < N + 1; i++) {
    yEigen(i, 0) = y[i];
    ycEigen(i, 0) = yc[i];
  }

  // IMPORTANT: The values assigned to half_shift and rev_half_shift, must
  // never be changed in any header/ the main program. These are global
  // variables ONLY for reading, not for writing into.
  // Rewriting on these can cause data races with unknown behavior.
}

/// \brief Chebyshev differentiation operator for a vector of Chebyshev
/// coefficients
template <class T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
diff(const Eigen::Matrix<T, Eigen::Dynamic, 1> &u) {
  int n = u.size() - 1;
  Eigen::Matrix<T, Eigen::Dynamic, 1> du(n + 1);
  du[n - 1] = 2.0 * n * u[n];
  // cout<<n-1<<"\t"<<du[n-1]<<endl;
// std::cout << "here" << '\n' << std::flush;
  for (int i = 1; i < (n + 1) / 2; i++) {
    du[n - (2 * i + 1)] =
        du[n - (2 * i - 1)] + 2.0 * (n - 2 * i) * u[n - 2 * i];
  }
  du[n] = 0.0;
  for (int i = 1; i < (n + 1) / 2; i++) {
    du[n - 2 * i] =
        du[n - 2 * (i - 1)] + 2.0 * (n - (2 * i - 1)) * u[n + 1 - 2 * i];
  }
  return du;
};

/// \brief Chebyshev differentiation operator for a vector of Chebyshev
/// coefficients
template <class T> std::valarray<T> diff(const std::valarray<T> &u) {
  int n = u.size() - 1;
  std::valarray<T> du(n + 1);
  du[n - 1] = 2.0 * n * u[n];
  // cout<<n-1<<"\t"<<du[n-1]<<endl;

  for (int i = 1; i < (n + 1) / 2; i++) {
    du[n - (2 * i + 1)] =
        du[n - (2 * i - 1)] + 2.0 * (n - 2 * i) * u[n - 2 * i];
  }
  du[n] = 0.0;
  for (int i = 1; i < (n + 1) / 2; i++) {
    du[n - 2 * i] =
        du[n - 2 * (i - 1)] + 2.0 * (n - (2 * i - 1)) * u[n + 1 - 2 * i];
  }
  return du;
};

/// \brief Chebyshev integration operator for a vector of Chebyshev
/// coefficients.
///
/// Gives the indefinite integral.
template <class T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
integ(const Eigen::Matrix<T, Eigen::Dynamic, 1> &u) {
  int n = u.size() - 1;
  int i, neve, nodd;
  Eigen::Matrix<T, Eigen::Dynamic, 1> I1v(u.size());
  I1v.setConstant(0.0);
  neve = (n + 1) / 2;
  nodd = neve;

  // split into eve and odd parts
  Eigen::Matrix<T, Eigen::Dynamic, 1> uodd(nodd), ueve(neve);
  uodd.setConstant(0.0);
  ueve.setConstant(0.0);
  for (i = 0; i < neve; i++)
    ueve[i] = u[2 * i];

  for (i = 0; i < nodd; i++)
    uodd[i] = u[2 * i + 1];

  I1v[0] = uodd[0] / 4.0;

  for (i = 1; i < neve; i++)
    I1v[2 * i] = (uodd[i - 1] - uodd[i]) / (4.0 * i);

  for (i = 0; i < nodd - 1; i++)
    I1v[2 * i + 1] = (ueve[i] - ueve[i + 1]) / (4 * i + 2.0);

  I1v[2 * i + 1] = ueve[i] / (4.0 * i + 2.0);

  return I1v;
};

/// \brief Chebyshev integration operator for a vector of Chebyshev
/// coefficients.
///
/// Gives the indefinite integral.
template <class T> std::valarray<T> integ(const std::valarray<T> &u) {
  int n = u.size() - 1;
  int i, neve, nodd;
  std::valarray<T> I1v(u.size());
  I1v = 0.0;
  neve = (n + 1) / 2;
  nodd = neve;

  // split into eve and odd parts
  std::valarray<T> uodd(nodd), ueve(neve);
  uodd = 0.0;
  ueve = 0.0;
  for (i = 0; i < neve; i++)
    ueve[i] = u[2 * i];

  for (i = 0; i < nodd; i++)
    uodd[i] = u[2 * i + 1];

  I1v[0] = uodd[0] / 4.0;

  for (i = 1; i < neve; i++)
    I1v[2 * i] = (uodd[i - 1] - uodd[i]) / (4.0 * i);

  for (i = 0; i < nodd - 1; i++)
    I1v[2 * i + 1] = (ueve[i] - ueve[i + 1]) / (4 * i + 2.0);

  I1v[2 * i + 1] = ueve[i] / (4.0 * i + 2.0);

  return I1v;
};

template <class T>
Eigen::Array<T, Eigen::Dynamic, 1>
integ(const Eigen::Array<T, Eigen::Dynamic, 1> &u) {
  int n = u.size() - 1;
  int i, neve, nodd;
  Eigen::Matrix<T, Eigen::Dynamic, 1> I1v(u.size());
  I1v.setConstant(0.0);
  neve = (n + 1) / 2;
  nodd = neve;

  // split into eve and odd parts
  Eigen::Matrix<T, Eigen::Dynamic, 1> uodd(nodd), ueve(neve);
  uodd.setConstant(0.0);
  ueve.setConstant(0.0);
  for (i = 0; i < neve; i++)
    ueve[i] = u[2 * i];

  for (i = 0; i < nodd; i++)
    uodd[i] = u[2 * i + 1];

  I1v[0] = uodd[0] / 4.0;

  for (i = 1; i < neve; i++)
    I1v[2 * i] = (uodd[i - 1] - uodd[i]) / (4.0 * i);

  for (i = 0; i < nodd - 1; i++)
    I1v[2 * i + 1] = (ueve[i] - ueve[i + 1]) / (4 * i + 2.0);

  I1v[2 * i + 1] = ueve[i] / (4.0 * i + 2.0);

  return I1v;
};

/// \brief Use this to make a complex array out of two Eigen real valarrays
template <class T>
Eigen::Array<std::complex<T>, Eigen::Dynamic, 1>
dou2com(const Eigen::Array<T, Eigen::Dynamic, 1> &a,
        const Eigen::Array<T, Eigen::Dynamic, 1> &b) {
  Eigen::Array<std::complex<T>, Eigen::Dynamic, 1> temp;
  temp.resize(a.size());
  for (int i = 0; i < a.size(); i++) {
    temp[i] = std::complex<T>(a[i], b[i]);
  }
  return temp;
};

/// \brief This is a chebfun analogue.
/// Chebfun will represent both values in physical space or an
/// array of Chebyshev-coefficients. Values can be in either space for
/// calculations. To convert between physical values or Chebyshev
/// coefficients, one can call Chebfun.c2p() and Chebfun.p2c().

template <class T> class Chebfun {
private:
public:
  int dct_flag;
  /// \brief Stores a flag to denote if function is in Cheb-space or
  /// physical space.
  std::valarray<T> v;
  /// \brief Null constructor
  ///
  /// By default the function is f(y) =y, with values in physical space.
  Chebfun() : v(y), dct_flag(SIS_PHYS_SPACE) {
    dct_flag = SIS_PHYS_SPACE;
    v.resize(N + 1);
    for (int i = 0; i < N + 1; i++)
      v[i] = cos(PI * (i + 0.5) / (N + 1.0));
  };
  /// \brief constructor by an input valarray
  Chebfun(const std::valarray<T> &in) {
    if (in.size() != N + 1) {
      std::cout << "Input size for Chebfun has to be sis::N +1. If you set "
                   "sis::N = 11, then value input valarray size must be N+1.\n "
                   "Exiting ...\n ";
      exit(1);
    }
    v = in;
    dct_flag = SIS_PHYS_SPACE;
  };
  /// \brief constructor by an input Eigen::Array<T,Eigen::Dynamic,1>
  Chebfun(const Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &in) {
    int bre;

    v.resize(in.size());
    if (in.size() != N + 1) {
      std::cout << "Input size for Chebfun has to be sis::N +1. If you set "
                   "sis::N = 11, then value input valarray size must be N+1.\n "
                   "Exiting ...\n ";
      exit(1);
    }
    for (int i = 0; i < in.size(); i++) {
      v[i] = in(i, 0);
    }
    dct_flag = SIS_PHYS_SPACE;
  };

  /*  /// \brief constructor by an input Eigen::Matrix<T,Eigen::Dynamic,1>
    Chebfun(const Eigen::Matrix<T, Eigen::Dynamic, 1> &in) {
      int bre;
      if (in.size() != N + 1) {
        std::cout << "Input size for Chebfun has to be sis::N +1. If you set "
                     "sis::N = 11, then value input valarray size must be N+1.\n
    " "Exiting ...\n "; exit(1);
      }
      for (int i = 0; i < in.size(); i++) {
        v[i] = in[i];
      }
      dct_flag = SIS_PHYS_SPACE;
      pf = fftw_plan_r2r_1d(N + 1, &v[0], &v[0], FFTW_REDFT10, FFTW_ESTIMATE);
      pb = fftw_plan_r2r_1d(N + 1, &v[0], &v[0], FFTW_REDFT01, FFTW_ESTIMATE);
    };*/
  /// \brief Copy conststructor.
  Chebfun(const Chebfun &in) {
    v = in.v;
    dct_flag = in.dct_flag;
  }

  /// \brief Converts a Chebfun from values in physical-space to coefficients of
  /// Chebyshev polynomials
  void p2c() {

    if (dct_flag == SIS_PHYS_SPACE) {
      v = dct(v);
    } else {
      std::cout << "In Cheb-space. Can't move to Cheb-space\n";
      exit(1);
    }
    dct_flag = SIS_CHEB_SPACE;
  };

  /// \brief Converts a Chebfun from Chebyshev coefficients to values in
  /// physical space
  void c2p() {
    if (dct_flag == SIS_CHEB_SPACE) {
      v = idct(v);
    } else {
      std::cout << "In Physical-space. Can't move to Physical-space\n";
      exit(1);
    }
    dct_flag = SIS_PHYS_SPACE;
  };

  /// \brief Returns property v in terms of Eigen::Array.
  Eigen::Array<T, Eigen::Dynamic, 1> ev() {
    Eigen::Array<T, Eigen::Dynamic, 1> temp;
    temp.resize(v.size());
    for (int i = 0; i < v.size(); i++) {
      temp[i] = v[i];
    }
    return temp;
  }

  /// \brief This function returns the Toeplitz + almost Hankel matrix needed to
  /// express multiplication of two Chebyshev series. This is used in solving
  /// for nonconstant coefficients. \brief To assign a Chebfun to a constant
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MultMat() {
    int bre;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> toep(N + 1, N + 1),
        hank(N + 1, N + 1);
    toep.setConstant(0.0);
    hank.setConstant(0.0);
    if (dct_flag == SIS_PHYS_SPACE) {
      p2c();
    }
    Eigen::Array<T, Eigen::Dynamic, 1> v_arr;
    v_arr.resize(v.size());
    for (int i = 0; i < v.size(); i++) {
      v_arr[i] = v[i];
    }
    Eigen::Matrix<T, Eigen::Dynamic, 1> vec(N + 1);
    vec.setConstant(0.0);
    vec.block(0, 0, N + 1, 1) = v_arr.matrix();
    for (int i = 0; i < N + 1; i++) {
      toep.block(i, i, 1, N + 1 - i) =
          vec.block(0, 0, N + 1 - i, 1).transpose();
    }

    toep = toep.eval() + toep.transpose().eval();
    toep.diagonal() = (toep.diagonal() / 2.0).eval();

    hank.setConstant(0.0);
    for (int i = 1; i < N + 1; i++) {
      hank.block(i, 0, 1, N - i) = vec.block(i, 0, N - i, 1).transpose();
    }

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> out;
    out = 0.5 * (toep + hank);
    return out.transpose();
  }

  /// \brief Equating the Chebfun to a constant
  void operator=(T right) {
    if (dct_flag == SIS_CHEB_SPACE) {
      v = 0.0;
      v[0] = 2.0 * right;
    } else {
      v = right;
    }
  }

  /// \brief To assign a Chebfun through another Chebfun
  void operator=(const Chebfun &right) {
    v = right.v;
    dct_flag = right.dct_flag;
  }

  /// \brief To assign a Chebfun through a vector
  void operator=(const Eigen::Matrix<T, Eigen::Dynamic, 1> &right) {
    for (int i = 0; i < N + 1; i++) {
      v[i] = right[i];
    }
  }

  /// \brief To assign a Chebfun through a valarray
  void operator=(const std::valarray<T> &right) { v = right; }

  //////////////////////
  /// \brief Adding constant to itself
  void operator+=(const T &right) {
    if (this->dct_flag == SIS_CHEB_SPACE) {
      v[0] += 2.0 * right;
    } else {
      v += right;
    }
  }

  /// \brief Add a Chebfun to itself. If both are in same space, return will be
  /// in same space, else the lhs will retain space.
  void operator+=(Chebfun right) {
    if (dct_flag == right.dct_flag) {
      v += right.v;
    } else {
      if (dct_flag == SIS_PHYS_SPACE) {
        right.c2p();
        v += right.v;
        right.p2c();
      } else {
        right.p2c();
        v += right.v;
        right.c2p();
      }
    }
  }
  /// \brief To assign a Chebfun through a Array
  //  void operator=(const Eigen::Array<T, Eigen::Dynamic, 1> &right) {
  //  this->v = right;
  //}
  /// \brief Evaluate the Chebfun at a point in the domain, i.e., -1 to 1.
  T operator()(const T &a) {
    if (a > 1 || a < -1) {
      std::cout << "Error: Cannot evaluate a Chebfun outside the domain"
                << '\n';
      exit(1);
    }
    std::valarray<T> k(N + 1);
    for (int i = 0; i < N + 1; i++) {
      k[i] = double(i);
    }
    k = cos(k * acos(a));
    k[0] = 0.5 * k[0];
    if (dct_flag == SIS_PHYS_SPACE) {
      p2c();
      k = k * v;
      c2p();
      return (k.sum());
    } else {
      k = k * v;
      return k.sum();
    }
  }

  /// \brief Returns a Chebfun that represents the indefinite integral of the
  /// Chebfun
  Chebfun<T> cumsum() {
    Chebfun<T> temp;

    if (dct_flag == SIS_PHYS_SPACE) {
      p2c();
      temp.dct_flag = SIS_CHEB_SPACE;
      temp.v = integ(v);
      c2p();
      temp.c2p();
      return temp;
    } else {
      v = integ(v);
      temp.dct_flag = SIS_CHEB_SPACE;
      return temp;
    }
  }
  /// \brief Returns the L2norm of the function.
  T L2norm() {
    Chebfun<T> temp;
    temp.dct_flag = SIS_PHYS_SPACE;
    if (dct_flag == SIS_PHYS_SPACE) {
      temp.v = v * v;
    } else {
      c2p();
      temp.v = v * v;
      p2c();
    }
    temp.p2c();
    temp.v = integ(temp.v);

    T ans = (temp.operator()(1) - temp.operator()(-1)) / 2.0;
    return ans;
  };

  bool isMachinePrecision() {
    if (abs(v[N - 1] + v[N]) / 2 < abs(std::numeric_limits<T>::epsilon())) {
      return true;
    } else {
      return false;
    }
  }
  bool isMachinePrecisionHalf() {
    if (abs(v[N / 2] + v[N / 2 + 1]) / 2 <
        abs(std::numeric_limits<T>::epsilon())) {
      return true;
    } else {
      return false;
    }
  }
  /// \brief Provides the truncation defined in terms of the abs(last +
  /// last_but_one)/2.0
  T trunc() { return abs(v[N] + v[N - 1]) / 2.0; }

  //~Chebfun(){
  //    v.~valarray();
  //}
};

/// \brief Complex conjugate of a Chebfun
template <class T> Chebfun<T> conj(Chebfun<T> in) { return in; };

/// \brief Chebfun overload to complex type.
template <class T> class Chebfun<std::complex<T> > {
private:
public:
  int dct_flag;
  /// \brief Stores a flag to denote if function is in Cheb-space or
  /// physical space.
  std::valarray<std::complex<T> > v;
  /// \brief Null constructor
  ///
  /// By default the function is f(y) =y, with values in physical space.
  Chebfun() : v(yc), dct_flag(SIS_PHYS_SPACE) {
    dct_flag = SIS_PHYS_SPACE;
    v.resize(N + 1);
    for (int i = 0; i < N + 1; i++) {
      v[i] = cos(PI * (i + 0.5) / (N + 1.0));
    }
  };
  /// \brief constructor by an input valarray
  Chebfun(const std::valarray<std::complex<T> > &in) {
    if (in.size() != N + 1) {
      std::cout << "Input size for Chebfun has to be sis::N +1. If you set "
                   "sis::N = 11, then value input valarray size must be N+1.\n "
                   "Exiting ...\n ";
      exit(1);
    }
    v = in;
    dct_flag = SIS_PHYS_SPACE;
  };
  /// \brief constructor by an input Eigen::Array<T,Eigen::Eigen::Dynamic,1>
  Chebfun(const Eigen::Array<std::complex<T>, Eigen::Dynamic, 1> &in) {
    v.resize(in.size());
    if (in.size() != N + 1) {
      std::cout << "Input size for Chebfun has to be sis::N +1. If you set "
                   "sis::N = 11, then value input valarray size must be N+1.\n "
                   "Exiting ...\n ";
      exit(1);
    }
    for (int i = 0; i < in.size(); i++) {
      v[i] = in[i];
    }
    dct_flag = SIS_PHYS_SPACE;
  };

  /// \brief constructor by an input Eigen::Matrix<T,Eigen::Eigen::Dynamic,1>
  Chebfun(const Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> &in) {
    v.resize(N + 1);
    if (in.size() != N + 1) {
      std::cout << "Input size for Chebfun has to be sis::N +1. If you set "
                   "sis::N = 11, then value input valarray size must be N+1.\n "
                   "Exiting ...\n ";
      exit(1);
    }
    for (int i = 0; i < in.size(); i++) {
      v[i] = in[i];
    }
    dct_flag = SIS_PHYS_SPACE;
  };

  /// \brief Copy constructor.
  Chebfun(const Chebfun<std::complex<T> > &in) {
    v = in.v;
    dct_flag = in.dct_flag;
  }

  /// \brief Converts a Chebfun from values in physical-space to coefficients of
  /// Chebyshev polynomials
  void p2c() {
    if (dct_flag == SIS_PHYS_SPACE) {
      v = dct(v);
    } else {
      std::cout << "In Cheb-space. Can't move to Cheb-space\n";
      exit(1);
    }
    dct_flag = SIS_CHEB_SPACE;
  };

  /// \brief Converts a Chebfun from Chebyshev coefficients to values in
  /// physical space
  void c2p() {
    if (dct_flag == SIS_CHEB_SPACE) {
      v = idct(v);
    } else {
      std::cout << "In Physical-space. Can't move to Physical-space\n";
      exit(1);
    }
    dct_flag = SIS_PHYS_SPACE;
  };

  /// \brief Returns property vr in terms of Eigen::Array.
  Eigen::Array<T, Eigen::Dynamic, 1> evr() {
    Eigen::Array<T, Eigen::Dynamic, 1> temp;
    temp.resize(v.size());
    for (int i = 0; i < v.size(); i++) {
      temp[i] = std::real(v[i]);
    }
    return temp;
  }

  /// \brief Returns property v in terms of Eigen::Array.
  Eigen::Array<T, Eigen::Dynamic, 1> evi() {
    Eigen::Array<T, Eigen::Dynamic, 1> temp;
    temp.resize(v.size());
    for (int i = 0; i < v.size(); i++) {
      temp[i] = std::imag(v[i]);
    }
    return temp;
  }

  /// \brief Returns property complex array vr + i vi in terms of Eigen::Array.
  Eigen::Array<std::complex<T>, Eigen::Dynamic, 1> evc() {
    Eigen::Array<std::complex<T>, Eigen::Dynamic, 1> temp;
    temp.resize(v.size());
    for (int i = 0; i < v.size(); i++) {
      temp[i] = v[i];
    }
    return temp;
  }
  /// \brief This function returns the Toeplitz + almost Hankel matrix needed to
  /// express multiplication of two Chebyshev series. This is used in solving
  /// for nonconstant coefficients. \brief To assign a Chebfun to a constant
  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> MultMat() {
    int bre;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> toep(N + 1,
                                                                        N + 1),
        hank(N + 1, N + 1);
    toep.setConstant(0.0);
    hank.setConstant(0.0);
    if (dct_flag == SIS_PHYS_SPACE) {
      p2c();
    }
    Eigen::Array<T, Eigen::Dynamic, 1> vr_arr, vi_arr;
    vr_arr.resize(v.size());
    vi_arr.resize(v.size());
    for (int i = 0; i < v.size(); i++) {
      vr_arr[i] = std::real(v[i]);
      vi_arr[i] = std::imag(v[i]);
    }
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> vec(N + 1);
    vec.setConstant(0.0);
    vec.block(0, 0, N + 1, 1) = dou2com(vr_arr, vi_arr).matrix();
    for (int i = 0; i < N + 1; i++) {
      toep.block(i, i, 1, N + 1 - i) =
          vec.block(0, 0, N + 1 - i, 1).transpose();
    }

    toep = toep.eval() + toep.transpose().eval();
    toep.diagonal() = (toep.diagonal() / 2.0).eval();

    hank.setConstant(0.0);
    for (int i = 1; i < N + 1; i++) {
      hank.block(i, 0, 1, N - i) = vec.block(i, 0, N - i, 1).transpose();
    }

    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> out;
    out = std::complex<double>(0.5, 0.0) * (toep + hank);
    return out.transpose();
  }

  /// \brief assignment operator of Chebfun
  void operator=(const Chebfun<std::complex<T> > &in) {
    v = in.v;
    dct_flag = in.dct_flag;
  }
  /// \brief assignment operator of Chebfun
  void operator=(const Chebfun<T> &in) {
    v.resize(in.v.size());
    for (int i = 0; i < v.size(); i++)
      v[i] = in.v[i];

    dct_flag = in.dct_flag;
  }

  /// \brief To assign a Chebfun to a constant
  void operator=(const std::complex<T> &right) {
    if (dct_flag == SIS_CHEB_SPACE) {
      v = std::complex<T>(0.0, 0.0);
      v[0] = 2.0 * right;
    } else {
      v = right;
    }
  }

  /// \brief To assign a Chebfun through a vector
  void
  operator=(const Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> &right) {
    v.resize(right.size());
    for (int i = 0; i < right.size(); i++) {
      v[i] = right[i];
    }
  }

  /// \brief To assign a Chebfun through a valarray
  void operator=(const std::valarray<std::complex<T> > &right) { v = right; }
  void operator=(const std::valarray<T> &right) {
    v = std::complex<T>(right, 0.0);
  }

  ////////////////////
  /// \brief Add complex chebfun to self.
  void operator+=(Chebfun<std::complex<T> > b) {
    if (dct_flag == b.dct_flag) {
      v += b.v;
    } else {
      if (dct_flag == SIS_PHYS_SPACE) {
        b.c2p();
        v += b.v;
        b.p2c();
      } else {
        b.p2c();
        v += b.v;
        b.c2p();
      }
    }
  }
  /// \brief Add real Chebfun to self
  void operator+=(Chebfun<T> b) {
    if (dct_flag == b.dct_flag) {
      for (int i = 0; i < v.size(); i++)
        v[i] += b.v[i];
    } else {
      if (dct_flag == SIS_PHYS_SPACE) {
        b.c2p();
        for (int i = 0; i < v.size(); i++)
          v[i] += b.v[i];
        b.p2c();
      } else {
        b.p2c();
        for (int i = 0; i < v.size(); i++)
          v[i] += b.v[i];
        b.c2p();
      }
    }
  }

  /// \brief add a constant to self
  void operator+=(std::complex<T> right) {
    if (dct_flag == SIS_CHEB_SPACE) {
      v[0] += 2.0 * right;
    } else {
      v += right;
    }
  }

  /// \brief To assign a Chebfun through a array
  // void
  // operator=(const Eigen::Array<std::complex<T>, Eigen::Dynamic, 1> &right) {
  //  this->vr = right.real();
  //  this->vi = right.imag();
  //  }
  /// \brief Evaluate the Chebfun at a point in the domain, i.e., -1 to 1.
  std::complex<T> operator()(const double &a) {
    if (a > 1.0 || a < -1.0) {
      std::cout << "Error: Cannot evaluate a Chebfun outside the domain"
                << '\n';
      exit(1);
    }
    std::valarray<std::complex<T> > k(N + 1);
    for (int i = 0; i < N + 1; i++) {
      k[i] = std::complex<T>(i, 0.0);
    }
    k = cos(k * acos(std::complex<T>(a, 0.0)));
    k[0] = std::complex<T>(0.5, 0.0) * k[0];
    if (dct_flag == SIS_PHYS_SPACE) {
      p2c();
      k = k * v;
      this->c2p();
      return (k.sum());
    } else {
      k = k * v;
      return k.sum();
    }
  }
  Chebfun<T> real() {
    Chebfun<T> fun;
    fun.v = real(v);
    fun.dct_flag = dct_flag;
    return fun;
  }
  Chebfun<T> imag() {
    Chebfun<T> fun;
    fun.v = imag(v);
    fun.dct_flag = dct_flag;
    return fun;
  };

  Chebfun<std::complex<T> > cumsum() {
    Chebfun<std::complex<T> > temp;
    if (dct_flag == SIS_PHYS_SPACE) {
      p2c();
      temp.v = integ(v);
      temp.dct_flag = SIS_CHEB_SPACE;
      c2p();
      temp.c2p();
      return temp;
    } else {
      temp.v = integ(v);
      temp.dct_flag = SIS_CHEB_SPACE;
      return temp;
    }
  }

  T L2norm() {
    Chebfun<T> temp;
    temp.dct_flag = SIS_PHYS_SPACE;
    if (dct_flag == SIS_PHYS_SPACE) {
      temp.v =
          std::real(std::valarray<std::complex<T> >(v * v.apply(std::conj)));
    } else {
      c2p();
      temp.v =
          std::real(std::valarray<std::complex<T> >(v * v.apply(std::conj)));
      p2c();
    }

    temp.p2c();
    temp.v = integ(temp.v);

    T ans = (temp.operator()(1) - temp.operator()(-1)) / 2.0;
    return ans;
  };
  T trunc() { return std::abs(v[N] + v[N - 1]) / 2.0; }
  //~Chebfun(){
  //    v.~valarray();
  //}
  /// \relates
};
/// \brief Complex conjugate of a Chebfun
template <class T>
Chebfun<std::complex<T> > conj(Chebfun<std::complex<T> > in) {
  in.v = in.v.apply(std::conj);
  return in;
};

/// Chebyshev differentiation operator, to differentiate n times.
template <class T> Chebfun<T> diff(Chebfun<T> in, int n) {
  if (n == 0) {
    return in;
  } else {
    if (in.dct_flag == SIS_PHYS_SPACE) {
      in.p2c();
      in.v = diff(in.v);
      in.c2p();
      return diff(in, n - 1);
    } else {
      in.v = diff(in.v);
      return diff(in, n - 1);
    }
  }
};

/// \brief This function overloads the cout<< operator to display the chebfun
template <class T>
std::ostream &operator<<(std::ostream &stream, Chebfun<T> a) {
  std::valarray<T> y(N + 1);
  for (int i = 0; i < N + 1; i++) {
    y[i] = std::cos(PI * (i + 0.5) / (N + 1.0));
  }
  stream << "A Chebfun in ";
  if (a.dct_flag == SIS_CHEB_SPACE) {
    stream << "Chebyshev-space. Values of coefficients are: \n";
    for (int i = 0; i < a.v.size(); i++) {
      stream << a.v[i] << "\n";
    }
  } else {
    stream << "physical-space. Values at set points are: \n";
    stream << "y\t\t f(y)\n";
    for (int i = 0; i < N + 1; i++) {
      stream << y[i] << "\t\t" << a.v[i] << "\n";
    }
  }
  return stream;
};

// template <class T> class ChebfunMat<std::complex<T> >;
/// \brief This class holds a matrix of Chebfuns.
template <class T> class ChebfunMat {
  friend class ChebfunMat<std::complex<T> >;
  friend class LinopMat<T>;

private:
  int count;

public:
  /// \brief Number of rows.
  int r;
  /// \brief Number of columns.
  int c;
  /// \brief This 1D vector of Chebfuns holds all Chebfuns of the ChebfunMat in
  /// the row-major format.
  std::valarray<Chebfun<T> > ChebfunVec;
  /// \brief Null constructor
  ChebfunMat() : ChebfunVec(), r(0), c(0), count(0) {}

  /// \brief Initializes a Chebfun Matrix of r_ rows and c_ columns.
  ChebfunMat(int r_, int c_) {
    r = r_;
    c = c_;
    count = 0;
    ChebfunVec.resize(r * c);
  };

  /// \brief Copy constructor
  ChebfunMat(const ChebfunMat<T> &in) {
    r = in.r;
    c = in.c;
    resize(r, c);
    count = in.count;
    for (int i = 0; i < r * c; i++) {
      ChebfunVec[i] = in.ChebfunVec[i];
    }
  }
  /// \brief Assignment operator
  void operator=(const ChebfunMat &in) {
    r = in.r;
    c = in.c;
    resize(r, c);
    count = in.count;
    for (int i = 0; i < r * c; i++) {
      ChebfunVec[i].v = in.ChebfunVec[i].v;
      ChebfunVec[i].dct_flag = in.ChebfunVec[i].dct_flag;
    }
  }

  /// \brief Assigning a ChebfunMat using a Eigen-matrix. Eigen-Matrix should be
  /// of size (r*(N+1), c), else an error is thrown. This implies that
  /// one must have called the resize(r,c) function before using this.
  void operator=(const Eigen::Matrix<T, Eigen::Dynamic,
                                     Eigen::Dynamic> &in) {

    if (in.rows() != r * (N + 1) || in.cols() != c) {
      std::cout
          << "Error in assigning a ChebfunMat through an Eigen Matrix. In line "
          << __LINE__ << ". Exiting ..." << '\n';
      exit(1);
    }
    for (int i = 0; i < r; i++) {
      for (int j = 0; j < c; j++) {
        operator()(i, j) = in.block(i * (N + 1), j, N + 1, 1);
      }
    }
  }

  /// \brief Use this to input multiple Chebfuns to the ChebfunMat using comma
  /// separators. Input type of constant.
  ChebfunMat<T> &operator<<(T b) {
    int bre;
    ChebfunVec[count].v = b;
    count++;
    return *this;
  };

  /// \brief Use this to to refer to Chebfuns. For 1D Matrices, we use the Row
  /// major format when refering to [i]. See wiki for row-major format.
  Chebfun<T> &operator[](int i) {
    Chebfun<T> &temp = ChebfunVec[i];
    return temp;
  };

  /// \brief Use this to to refer to Chebfuns. For 2D Matrices, refer to matrix
  /// element Chebfun by (i,j) with the row and column index. Indices start from
  /// 0.
  Chebfun<T> &operator()(int i, int j) { return ChebfunVec[j + c * i]; };

  /// \brief Converts every element to Cheb-space, if not already in Cheb-space.
  void p2c(){
    for (int i = 0; i < r; i++){
      for (int j = 0; j < c; j++){
        if (operator()(i,j).dct_flag == SIS_PHYS_SPACE){
          operator()(i,j).p2c();
        }
      }
    }
  }

  /// \brief Converts every element to phys-space, if not already in phys-space.
  void c2p(){
    for (int i = 0; i < r; i++){
      for (int j = 0; j < c; j++){
        if (operator()(i,j).dct_flag == SIS_CHEB_SPACE){
          operator()(i,j).c2p();
        }
      }
    }
  }


  /// \brief Use this to evaluate the whole ChebfunMat at a given point. Point
  /// must lie in the domain, [-1,1].
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic >
  operator()(T p) {
    if (p <=1 && p >= -1 ){
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > temp(r,c);
      for (int i = 0 ; i < r; i ++){
        for (int j = 0; j < c; j++){
          temp(i,j) = operator()(i,j)(p);
        }
      }
      return temp;
    } else {
      std::cout << "Point of evaluation not in [-1,1]. In " << __LINE__
                << '\n';
      exit(1);
    }
  };

  /// \brief Use this to input multiple Chebfuns to the ChebfunMat using comma
  /// separators. Input type Eigen array.
  ChebfunMat<T> &operator<<(Eigen::Array<T, Eigen::Dynamic, 1> b) {
    int bre;
    for (int i = 0; i < N + 1; i++) {
      ChebfunVec[count].v[i] = b[i];
    }
    count++;
    return *this;
    /// separators. Input type Chebfun.
  };

  /// \brief Use this to input multiple Chebfuns to the ChebfunMat using comma
  /// separators. Input type valarray.
  ChebfunMat<T> &operator<<(std::valarray<T> b) {
    int bre;
    ChebfunVec[count].v = b;
    count++;
    return *this;
    /// separators. Input type Chebfun.
  };

  /// \brief Use this to input multiple Chebfuns to the ChebfunMat using comma
  ChebfunMat<T> &operator<<(Chebfun<T> b) {
    int bre;
    ChebfunVec[count] = b;
    count++;
    return *this;
  };

  /// \brief This clears all contents in the ChebfunMat, and then creates a
  /// fresh ChebfunMat of size r_ x c_.
  void resize(int r_, int c_) {
    r = r_;
    c = c_;
    ChebfunVec.resize(r * c);
    count = 0;
  }

  /// \brief Overloads comma separator to input Chebfuns into a ChebfunMat.
  /// Input type constant.
  ChebfunMat<T> &operator,(T b) {
    ChebfunVec[count].v = b;
    count++;
    return *this;
  }
  /// \brief Overloads comma separator to input Chebfuns into a ChebfunMat.
  /// Input type Eigen array.
  ChebfunMat<T> &operator,(Eigen::Array<T, Eigen::Dynamic, 1> b) {
    for (int i = 0; i < N + 1; i++) {
      ChebfunVec[count].v[i] = b[i];
    }
    count++;
    return *this;
  }

  /// \brief Overloads comma separator to input Chebfuns into a ChebfunMat.
  /// Input type valarray.
  ChebfunMat<T> &operator,(std::valarray<T> b) {
    ChebfunVec[count].v = b;
    count++;
    return *this;
  }
  /// \brief Overloads comma separator to input Chebfuns into a ChebfunMat.
  /// Input type Chebfun.
  ChebfunMat<T> &operator,(Chebfun<T> b) {
    ChebfunVec[count] = b;
    count++;
    return *this;
  }

  /// \brief sets the all the Chebfuns in the ChebfunMat to constant
  void setConstant(double in) {
    for (int i = 0; i < r; i++) {
      for (int j = 0; j < c; j++) {
        operator()(i, j) = in;
      }
    }
  }

  ChebfunMat<T> cTranspose() {
    sis::ChebfunMat<T> out;
    out.resize(c, r);
    for (int i = 0; i < r; i++) {
      for (int j = 0; j < c; j++) {
        out(j, i) = conj(operator()(i, j));
      }
    }
    return out;
  };

  /// \brief Returns an Eigen matrix representing the ChebfunMat.
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
  ChebfunMat2EigenMat() {
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> out(
        r * (N + 1), c);
    int bre;
    // std::cout << "r,c = " << std::complex<int>(r, c) << '\n';
    for (int i = 0; i < r; i++) {
      for (int j = 0; j < c; j++) {
        // std::cout << "in " << __LINE__ << '\n';
        if (operator()(i, j).dct_flag == SIS_CHEB_SPACE) {
          Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> temp;
          out.block(i * (N + 1), j, N + 1, 1) =
          // temp =
          operator()(i, j).ev();
          // std::cout << "in " << __LINE__ << '\n';
        } else {
          operator()(i, j).p2c();
          out.block(i * (N + 1), j, N + 1, 1) = operator()(i, j).ev();
          operator()(i, j).c2p();
        }
      }
    }
  }
  //  ~ChebfunMat(){
  //    ChebfunVec.~valarray();
  //  }
};

/// \brief ChebfunMat overloaded to complex type.
template <class T> class ChebfunMat<std::complex<T> > {
private:
  int count;

public:
  /// \brief Number of rows.
  int r;
  /// \brief Number of columns.
  int c;
  /// \brief This 1D vector of Chebfuns holds all Chebfuns of the ChebfunMat in
  /// the row-major format.
  std::valarray<Chebfun<std::complex<T> > > ChebfunVec;
  /// \brief Null constructor
  ChebfunMat() : ChebfunVec(), r(0), c(0), count(0) {}

  /// \brief Initializes a Chebfun Matrix of r_ rows and c_ columns.
  ChebfunMat(int r_, int c_) {
    r = r_;
    c = c_;
    ChebfunVec.resize(r * c);
  };

  /// \brief Copy constructor
  ChebfunMat(const ChebfunMat<std::complex<T> > &in) {
    r = in.r;
    c = in.c;
    resize(r, c);
    count = in.count;
    for (int i = 0; i < r * c; i++) {
      ChebfunVec[i] = in.ChebfunVec[i];
    }
  }
  /// \brief Assignment operator
  void operator=(const ChebfunMat<std::complex<T> > &in) {
    resize(in.r, in.c);
    count = in.count;
    for (int i = 0; i < r * c; i++) {
      ChebfunVec[i].v = in.ChebfunVec[i].v;
      ChebfunVec[i].dct_flag = in.ChebfunVec[i].dct_flag;
    }
  }

  /// \brief Assignment operator
  void operator=(const ChebfunMat<T> &in) {
    r = in.r;
    c = in.c;
    resize(r, c);
    count = in.count;
    for (int i = 0; i < r * c; i++) {
      ChebfunVec[i] = in.ChebfunVec[i];
    }
  }

  /// \brief Use this to input multiple Chebfuns to the ChebfunMat using comma
  /// separators. Input type of complex constant.
  ChebfunMat<std::complex<T> > &operator<<(std::complex<T> b) {
    int bre;
    ChebfunVec[count].v = b;
    count++;
    return *this;
  };

  /// \brief Use this to to refer to Chebfuns. For 1D Matrices, we use the Row
  /// major format when refering to [i]. See wiki for row-major format.
  Chebfun<std::complex<T> > &operator[](int i) {
    Chebfun<std::complex<T> > &temp = ChebfunVec[i];
    return temp;
  };

  /// \brief Converts every element to Cheb-space, if not already in Cheb-space.
  void p2c(){
    for (int i = 0; i < r; i++){
      for (int j = 0; j < c; j++){
        if (operator()(i,j).dct_flag == SIS_PHYS_SPACE){
          operator()(i,j).p2c();
        }
      }
    }
  }

  /// \brief Converts every element to phys-space, if not already in phys-space.
  void c2p(){
    for (int i = 0; i < r; i++){
      for (int j = 0; j < c; j++){
        if (operator()(i,j).dct_flag == SIS_CHEB_SPACE){
          operator()(i,j).c2p();
        }
      }
    }
  }
  /// \brief Use this to to refer to Chebfuns. For 2D Matrices, refer to matrix
  /// element Chebfun by (i,j) with the row and column index. Indices start from
  /// 0.
  Chebfun<std::complex<T> > &operator()(int i, int j) {

    return ChebfunVec[j + c * i];
  };

  /// \brief Use this to evaluate the whole ChebfunMat at a given point. Point
  /// must lie in the domain, [-1,1].
  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic >
  operator()(T p) {
    if (p <=1 && p >= -1 ){
      Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic > temp(r,c);
      for (int i = 0 ; i < r; i ++){
        for (int j = 0; j < c; j++){
          temp(i,j) = operator()(i,j)(p);
        }
      }
      return temp;
    } else {
      std::cout << "Point of evaluation not in [-1,1]. In " << __LINE__
                << '\n';
      exit(1);
    }
  };

  /// \brief Use this to input multiple Chebfuns to the ChebfunMat using comma
  /// separators. Input type complex Eigen array.
  ChebfunMat<std::complex<T> > &
  operator<<(Eigen::Array<std::complex<T>, Eigen::Dynamic, 1> b) {
    for (int i = 0; i < N + 1; i++) {
      ChebfunVec[count].v[i] = b[i];
    }
    count++;
    return *this;
  };

  // \brief Use this to input multiple Chebfuns to the ChebfunMat using comma
  /// separators. Input type complex valarray.
  ChebfunMat<std::complex<T> > &operator<<(std::valarray<std::complex<T> > b) {
    ChebfunVec[count].v = b;
    count++;
    return *this;
  };

  // \brief Use this to input multiple Chebfuns to the ChebfunMat using comma
  /// separators. Input type complex valarray.
  ChebfunMat<std::complex<T> > &operator<<(std::valarray<T> b) {
    int bre;
    for (int i = 0; i < N + 1; i++)
      ChebfunVec[count].v[i] = b[i];
    count++;
    return *this;
  };

  /// \brief Use this to input multiple Chebfuns to the ChebfunMat using comma
  /// separators. Input type Chebfun.
  ChebfunMat<std::complex<T> > &operator<<(Chebfun<std::complex<T> > b) {
    int bre;
    ChebfunVec[count] = b;
    count++;
    return *this;
  };

  /// \brief This clears all contents in the ChebfunMat, and then creates a
  /// fresh ChebfunMat of size r_ x c_.
  void resize(int r_, int c_) {
    r = r_;
    c = c_;
    ChebfunVec.resize(r * c);
    count = 0;
  }

  /// \brief Overloads comma separator to input Chebfuns into a ChebfunMat.
  /// Input type complex constant.
  ChebfunMat<std::complex<T> > &operator,(std::complex<T> b) {
    ChebfunVec[count].v = b;
    count++;
    return *this;
  }

  /// \brief Overloads comma separator to input Chebfuns into a ChebfunMat.
  /// Input type complex Eigen array.
  ChebfunMat<std::complex<T> > &operator,(
      Eigen::Array<std::complex<T>, Eigen::Dynamic, 1> b) {
    for (int i = 0; i < N + 1; i++) {
      ChebfunVec[count].v[i] = b[i];
    }
    count++;
    return *this;
  }
  /// \brief Overloads comma separator to input Chebfuns into a ChebfunMat.
  /// Input type complex valarray.
  ChebfunMat<std::complex<T> > &operator,(std::valarray<std::complex<T> > b) {
    ChebfunVec[count].v = b;
    count++;
    return *this;
  }

  /// \brief Returns an Eigen matrix representing the ChebfunMat.
  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
  ChebfunMat2EigenMat() {
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> out(
        r * (N + 1), c);
    int bre;
    // std::cout << "r,c = " << std::complex<int>(r, c) << '\n';
    for (int i = 0; i < r; i++) {
      for (int j = 0; j < c; j++) {
        // std::cout << "in " << __LINE__ << '\n';
        if (operator()(i, j).dct_flag == SIS_CHEB_SPACE) {
          Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> temp;
          out.block(i * (N + 1), j, N + 1, 1) =
          // temp =
          operator()(i, j).evc();
          // std::cout << "in " << __LINE__ << '\n';
        } else {
          operator()(i, j).p2c();
          out.block(i * (N + 1), j, N + 1, 1) = operator()(i, j).evc();
          operator()(i, j).c2p();
        }
      }
    }
    // std::cout << "out is " << out  << '\n';
    // std::cin >> bre;
    return out;
  }

  /// \brief Assigning a ChebfunMat using a Eigen-matrix. Eigen-Matrix should be
  /// of size (r*(N+1), c), else an error is thrown. This implies that
  /// one must have called the resize(r,c) function before using this.
  void operator=(const Eigen::Matrix<std::complex<T>, Eigen::Dynamic,
                                     Eigen::Dynamic> &in) {

    if (in.rows() != r * (N + 1) || in.cols() != c) {
      std::cout
          << "Error in assigning a ChebfunMat through an Eigen Matrix. In line "
          << __LINE__ << ". Exiting ..." << '\n';
      exit(1);
    }
    for (int i = 0; i < r; i++) {
      for (int j = 0; j < c; j++) {
        operator()(i, j) = in.block(i * (N + 1), j, N + 1, 1);
      }
    }
  }

  /// \brief sets the all the Chebfuns in the ChebfunMat to constant
  void setConstant(T in) {
    for (int i = 0; i < r; i++) {
      for (int j = 0; j < c; j++) {
        operator()(i, j) = in;
      }
    }
  }

  /// \brief sets the all the Chebfuns in the ChebfunMat to constant
  void setConstant(std::complex<T> in) {
    for (int i = 0; i < r; i++) {
      for (int j = 0; j < c; j++) {
        operator()(i, j) = in;
      }
    }
  }

  ChebfunMat<std::complex<T> > cTranspose() {
    sis::ChebfunMat<std::complex<T> > out;
    out.resize(c, r);
    for (int i = 0; i < r; i++) {
      for (int j = 0; j < c; j++) {
        out(j, i) = conj(operator()(i, j));
      }
    }
    return out;
  };
  //  ~ChebfunMat(){
  //    ChebfunVec.~valarray();
  //  }
};

/// \brief Complex Conjugate of a ChebfunMat. Not the conjugate transpose.
/// See cTranspose() in ChebfunMat for complex conjugate transpose.
template <class T> sis::ChebfunMat<T> conj(sis::ChebfunMat<T> in) {
  sis::ChebfunMat<T> out;
  out.resize(in.r, in.c);
  for (int i = 0; i < in.r; i++) {
    for (int j = 0; j < in.c; j++) {
      out(i, j) = conj(in(i, j));
    }
  }
  return out;
};

/// \brief Chebyshev differentiation operator, for a ChebfunMat, differeniates
/// every Chebfun in the ChebfunMat n times.
template <class T> ChebfunMat<T> diff(ChebfunMat<T> in, int n) {
  ChebfunMat<T> out;
  int bre;
  out.resize(in.r, in.c);
  for (int i = 0; i < in.r; i++) {
    for (int j = 0; j < in.c; j++) {
      out(i, j) = diff(in(i, j), n);
    }
  }
  return out;
};

/** \brief Linop
 * This class creates a Linear operator to solve TPBVPs.
 */
template <class T> class Linop {
private:
  int i, j;

public:
  int NCC; ///< Flag to set (value 1) for enabling nonconstant coefficients
  int n;   ///< The order of the Linear differential operator

  Eigen::Matrix<T, Eigen::Dynamic, 1>
      coef; ///< \brief Stores the coefficients in the differential equation.
  /// \brief \anchor coefFun Use this to input all non-constant coefficients.
  ///
  /// For instance, for the differential
  /// operator \f$2y^2\,\partial_y + 3y\f$, coefFun\f$ [0]= 2y^2\f$
  /// and coef\f$[1] = 3y\f$, where y is a array whose cheb-points have been
  /// set. Note that this can only be called after calling Linop.ncc(), which
  /// initiates non-constant coefficients. One can also use the Eigen format to
  /// input coefficients. For example, to set coefficients of above differential
  /// equation, we need to do the following (assuming double precision):
  /// \code{.cpp}
  /// Eigen::ArrayXd y;
  /// setChebPts(y);
  /// Linop <double> L(2); // Second order differential operator
  /// L.ncc();
  /// L.coefFun[0] = 2.0*pow(y,2.0);
  /// L.coefFun[1] = 3.0*y;
  /// \endcode
  /// Alternatively, in the eigen-input format, one can do the following:
  ///\code{.cpp}
  /// Eigen::ArrayXd y;
  /// setChebPts(y);
  /// Linop <double> L(2); // Second order differential operator
  /// L.ncc();
  /// L.coefFun << 2.0*pow(y,2.0), 3.0*y;
  ///\endcode
  ///
  /// Note that when ncc() is called coef is no more used in any function in
  /// Linop. In other words, DO NOT to mix coef and coefFun and expect it to
  /// work. If there is even one non-constant coefficient in the differential
  /// equation, then you must set coefficients in coefFun, and not in coef.
  /// This also means that constant coefficients problem can be solved as a
  /// non-constant coefficient problem, as in consider the differential operator
  /// \f$2\,\partial_y + 3\f$. This can input as a constant coefficient problem
  /// in the following manner:
  /// \code{.cpp}
  /// Linop <double> L(1);
  /// L.set();
  /// L.coef << 2.0, 3.0;
  /// \endcode
  /// or as a non-constant coefficient problem as in:
  /// \code{.cpp}
  /// Linop <double> L(1);
  /// L.ncc();
  /// L.coefFun << 2.0, 3.0;
  /// \endcode
  /// Ofcourse, constant input to a coefFun automatically assigns values at all
  /// cheb-points to that constant.
  ///
  /// However, the constant coefficient solver is far more efficient in terms of
  /// speed and memory than a non-constant solver, and hence we provide separate
  /// solvers for each case.
  ChebfunMat<T> coefFun;
  /// The solution here involves storing the Chebyshev coeefficients of the
  /// highest derivative in the differential equation, and also the constants of
  /// integration. This can be used to construct all other lower derivatives
  /// using integration.
  ///
  /// First N+1 elements of sol hold the Chebyshev coefficients and the
  /// remaining elements store the integration constants
  Eigen::Matrix<T, Eigen::Dynamic, 1> solution; ///< Stores the solution.

  /// \brief A null-constructor for Linop.
  Linop() { NCC = 0; };

  /// \brief Makes a Linop by specifying the order n_ of the differential
  /// equation.
  Linop(const int &n_) {
    n = n_;
    coef.resize(n + 1);
    NCC = 0;
  };

  /// \brief Copy constructor.
  Linop(const Linop<T> &in) {
    n = in.n;
    NCC = in.NCC;
    coef = in.coef;
    if (NCC == 1) {
      coefFun = in.coefFun;
    }
    solution = in.solution;
  }

  /// \brief Assignment operator
  void operator=(const Linop<T> &in) {
    n = in.n;
    NCC = in.NCC;
    coef = in.coef;
    coefFun = in.coefFun;
    solution = in.solution;
  }

  /// \brief Apply the operator on a Chebfun.
  Chebfun<T> operator()(Chebfun<T> in){
    Chebfun<T> out;
    out.v = 0.0;
    // if non constant coefficients:
    if (NCC==1){
      for (int i = 0; i < n + 1; i ++){
        out = out + (coefFun[n-i]*diff(in,i));
      }
    } else {
      for (int i = 0; i < n + 1; i ++){
        out = out + (coef[n-i]*diff(in,i));
      }
    }
    return out;
  }

  /// \brief Add to self:
  void operator+=(const Linop<T> &in) {
    Linop<T> temp;
    temp = *this + in;
    operator=(temp);
  }

  /// \brief Assignment operator for constant input
  void operator=(const T &in) {
    n = 0;
    NCC = 0;
    set();
    coef[0] = in;
  }

  /// \brief Assignment operator for function input, through Eigen array
  void operator=(const Eigen::Array<T, Eigen::Dynamic, 1> &in) {
    n = 0;
    ncc();
    coefFun << in;
  }

  /// \brief Assignment operator for function input, through valarray
  void operator=(const std::valarray<T> &in) {
    n = 0;
    ncc();
    coefFun << in;
  }
  /// \brief Assignment through a Chebfun:
  void operator=(Chebfun<T> in) {
    if (in.dct_flag == SIS_PHYS_SPACE) {
      operator=(in.v);
    } else {
      in.c2p();
      operator=(in.v);
      in.p2c();
    }
  };
  /// \brief Dividing Linop by scalar.
  Linop<T> operator/(T in) {
    Linop<T> out = *this;
    if (out.NCC == 0) {
      out.coef = out.coef / in;
    } else {
      for (int i = 0; i < out.n + 1; i++) {
        out.coefFun[i].v = out.coefFun[i].v / in;
      }
    }
    return out;
  }

  /// \brief Sets the the size of coefficient array based on the order of
  /// differential equation. Suppose order is 3, then the size of array needed
  /// is 4.
  void set() {
    coef.resize(n + 1);
    NCC = 0;
  };

  /// \brief Sets the the size of coefficient array based on the order of
  /// differential equation. Suppose order is 3, then the size of array needed
  /// is 4.
  void set(int n_) {
    n = n_;
    coef.resize(n_ + 1);
    NCC = 0;
  };

  /// \brief Sets the the size of non-constant coefficient array based on the
  /// order of differential equation. Suppose order is 3, then the size of array
  /// needed is 4. This must be called before assigning non-constant
  /// coefficients
  void ncc() {
    NCC = 1;
    coefFun.resize(1, n + 1);
  }

  /// \brief Sets the the size of non-constant coefficient array based on the
  /// order of differential equation. Suppose order is 3, then the size of array
  /// needed is 4. This must be called before assigning non-constant
  /// coefficients
  void ncc(int n_) {
    NCC = 1;
    n = n_;
    coefFun.resize(1, n + 1);
  }
};

/// \brief Overloads the Linop class to complex type.
template <class T> class Linop<std::complex<T> > {
private:
  int i, j;

public:
  int NCC; ///< Flag to set (value 1) for enabling nonconstant coefficients
  int n;   ///< the order of the equation
  /// For instance, for the differential
  /// operator \f$2i\,\partial_y + 3\f$, coef\f$ [0]= 2i\f$
  /// and coef\f$[1] = 3\f$.
  /// One can also use the Eigen format to input coefficients:
  /// \code{.cpp}
  /// coef << std::complex<double>(0.0,2.0), std::complex<double>(3,0.0);
  /// \endcode
  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1>
      coef; ///< \brief Stores the coefficients in the differential equation

  /// The solution here involves storing the Chebyshev coeefficients of the
  /// highest derivative in the differential equation, and also the constants of
  /// integration. This can be used to construct all other lower derivatives
  /// using integration.
  ///
  /// First N+1 elements of sol hold the Chebyshev coefficients and the
  /// remaining elements store the integration constants

  /// \brief Use this to input all non-constant coefficients. See \ref coefFun
  /// for details.
  ///
  /// Additional details are as follows. Although C++ is type specific, coefFun
  /// allows inputs to be real, by assuming that the imaginary part is zero. in
  /// addition to what is written in \ref coefFun, consider the following
  /// operator: \f$2i\partial_{yy}u(y) \,+\, 4\partial_yu(y)\,+\,y^2u(y)\f$,
  /// then the following code can be used to input the nonconstant coefficients:
  /// \code{.cpp}
  /// int main(){
  /// typedef std::complex<double> dc;
  /// Linop<dc> L(2); // Linop of 2nd order.
  /// Eigen::ArrayXd y;
  /// setChebPts(y);
  /// L.ncc();
  /// dc ii = dc(0.0,1);
  /// L.coefFun << 2.0*ii, 4.0, pow(y,2.0);
  /// return 0;
  /// }
  /// \endcode
  /// Note that in the above input to coefFun has mixed types, like double and
  /// complex<double>. This feature is allowed only coefFun. For coef, one has
  /// to be type specific, that is all inputs have to be complex<double> type,
  /// or double type, depending on what template is used to intiate the Linop.

  ChebfunMat<std::complex<T> > coefFun;
  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1>
      solution; ///< Stores the solution.

  /// \brief A null-constructor for Linop.
  Linop() { NCC = 0; };

  /// Makes a Linop by specifying the order n_ of the differential equation.
  Linop(const int &n_) {
    n = n_;
    coef.resize(n + 1);
    NCC = 0;
  };

  /// \brief Copy constructor.
  Linop(const Linop<std::complex<T> > &in) {
    n = in.n;
    NCC = in.NCC;
    coef = in.coef;
    if (NCC == 1) {
      coefFun = in.coefFun;
    }
    solution = in.solution;
  }
  void operator=(const Linop<std::complex<T> > &in) {
    n = in.n;
    NCC = in.NCC;
    coef = in.coef;
    if (NCC == 1) {
      coefFun = in.coefFun;
    }
    solution = in.solution;
  }

  void operator=(const Linop<T> &in) {
    n = in.n;
    NCC = in.NCC;
    coef = in.coef;
    coefFun = in.coefFun;
    solution = in.solution;
  }

  void operator+=(const Linop<T> &in) {
    Linop<std::complex<T> > temp;
    temp = *this + in;
    operator=(temp);
  }

  /// \brief Apply the operator on a Chebfun.
  Chebfun<std::complex<T> > operator()(Chebfun<std::complex<T> > in){
    Chebfun<std::complex<T> > out;
    out.v = std::complex<T>(0.0, 0.0);
    // if non constant coefficients:
    if (NCC==1){
      for (int i = 0; i < n + 1; i ++){
        out = out + (coefFun[n-i]*diff(in,i));
      }
    } else {
      for (int i = 0; i < n + 1; i ++){
        out = out + (coef[n-i]*diff(in,i));
      }
    }
    return out;
  }


  void operator+=(const Linop<std::complex<T> > &in) {
    Linop<std::complex<T> > temp;
    temp = *this + in;
    operator=(temp);
  }

  /// \brief Assignment operator for constant input
  void operator=(const std::complex<T> &in) {
    n = 0;
    NCC = 0;
    set();
    coef[0] = in;
  }

  /// \brief Assignment operator for constant input
  void operator=(T in) {
    n = 0;
    NCC = 0;
    set();
    coef[0] = in;
  }

  /// \brief Assignment operator for function input, through Eigen array
  void operator=(const Eigen::Array<std::complex<T>, Eigen::Dynamic, 1> &in) {
    n = 0;
    ncc();
    coefFun << in;
  }

  /// \brief Assignment operator for function input, through valarray
  void operator=(const std::valarray<T> &in) {
    n = 0;
    ncc();
    coefFun << in;
  }

  void operator=(const std::valarray<std::complex<T> > &in) {
    n = 0;
    ncc();
    coefFun << in;
  }

  void operator=(Chebfun<T> in) {
    if (in.dct_flag == SIS_PHYS_SPACE) {
      operator=(in.v);
    } else {
      in.c2p();
      operator=(in.v);
      in.p2c();
    }
  };

  void operator=(Chebfun<std::complex<T> > in) {
    if (in.dct_flag == SIS_PHYS_SPACE) {
      operator=(in.v);
    } else {
      in.c2p();
      operator=(in.v);
      in.p2c();
    }
  };

  /// \brief Dividing Linop by scalar.
  Linop<std::complex<T> > operator/(T in) {
    Linop<std::complex<T> > out = *this;
    if (out.NCC == 0) {
      out.coef = out.coef / in;
    } else {
      for (int i = 0; i < out.n + 1; i++) {
        out.coefFun[i].v = out.coefFun[i].v / std::complex<T>(in, 0.0);
      }
    }
    return out;
  }

  /// \brief Dividing Linop by scalar.
  Linop<std::complex<T> > operator/(std::complex<T> in) {
    Linop<std::complex<T> > out = *this;
    if (out.NCC == 0) {
      out.coef = out.coef / in;
    } else {
      for (int i = 0; i < out.n + 1; i++) {
        out.coefFun[i].v = out.coefFun[i].v / in;
      }
    }
    return out;
  }

  /// \brief Multiplying Linop by scalar.
  Linop<std::complex<T> > operator*(T in) {
    Linop<std::complex<T> > out = *this;
    if (out.NCC == 0) {
      out.coef = out.coef * in;
    } else {
      for (int i = 0; i < out.n + 1; i++) {
        out.coefFun[i].v = out.coefFun[i].v * std::complex<T>(in, 0.0);
      }
    }
    return out;
  }
  /*  /// \brief Multiplying Linop by ArrayXd.
    Linop<std::complex<T> > operator*(Eigen::Array<T, Eigen::Dynamic, 1> in) {
      Linop<std::complex<T> > out = *this;
      if (out.NCC == 0) {
        ncc();
        for (int i = 0; i < n + 1; i++) {
          out.coefFun[i].v = out.coef[i] * in;
        }
      } else {
        for (int i = 0; i < out.n + 1; i++) {
          out.coefFun[i].vr = out.coefFun[i].vr * in;
          out.coefFun[i].vi = out.coefFun[i].vi * in;
        }
      }
      return out;
    }

    /// \brief Multiplying Linop by ArrayXd.
    Linop<std::complex<T> >
    operator*(Eigen::Array<std::complex<T>, Eigen::Dynamic, 1> in) {
      Linop<std::complex<T> > out = *this;
      if (out.NCC == 0) {
        ncc();
        for (int i = 0; i < n + 1; i++) {
          out.coefFun[i].vr = (out.coef[i] * in).real();
          out.coefFun[i].vi = (out.coef[i] * in).imag();
        }
      } else {
        for (int i = 0; i < out.n + 1; i++) {
          out.coefFun[i].vr =
              out.coefFun[i].vr * in.real() - out.coefFun[i].vi * in.imag();
          out.coefFun[i].vi =
              out.coefFun[i].vi * in.real() + out.coefFun[i].vr * in.imag();
        }
      }
      return out;
    }
    */
  /// \brief Multiplying Linop by scalar.
  Linop<std::complex<T> > operator*(std::complex<T> in) {
    Linop<std::complex<T> > out = *this;
    if (out.NCC == 0) {
      out.coef = out.coef * in;
    } else {
      for (int i = 0; i < out.n + 1; i++) {
        out.coefFun[i].v = out.coefFun[i].v * in;
      }
    }
    return out;
  }
  /// \brief Sets the the size of coefficient array based on the order of
  /// differential equation. Suppose order is 3, then the size of array needed
  /// is 4.
  void set() {
    coef.resize(n + 1);
    NCC = 0;
  };

  /// \brief Sets the the size of coefficient array based on the order of
  /// differential equation. Suppose order is 3, then the size of array needed
  /// is 4.
  void set(int n_) {
    n = n_;
    coef.resize(n_ + 1);
    NCC = 0;
  };
  /// \brief Sets the the size of non-constant coefficient array based on the
  /// order of differential equation. Suppose order is 3, then the size of array
  /// needed is 4. This must be called before assigning non-constant
  /// coefficients
  void ncc() {
    NCC = 1;
    coefFun.resize(1, n + 1);
  }

  /// \brief Sets the the size of non-constant coefficient array based on the
  /// order of differential equation. Suppose order is 3, then the size of array
  /// needed is 4. This must be called before assigning non-constant
  /// coefficients
  void ncc(int n_) {
    NCC = 1;
    n = n_;
    coefFun.resize(1, n + 1);
  }
  sis::Linop<T> real() {
    sis::Linop<T> L;
    L.n = this->n;
    if (NCC == 0) {
      L.set();
      L.coef = coef.real();
    } else {
      L.ncc();
      for (int i = 0; i < n + 1; i++)
        L.coefFun[i].v = std::real(coefFun[i].v);
    }

    return L;
  };

  sis::Linop<T> imag() {
    sis::Linop<T> L;
    L.n = this->n;
    if (NCC == 0) {
      L.set();
      L.coef = coef.imag();
    } else {
      L.ncc();
      for (int i = 0; i < n + 1; i++)
        L.coefFun[i].v = std::imag(coefFun[i].v);
    }

    return L;
  };
  /// \brief This function solves the differential equation. See \ref solve
  /// for details.
  Chebfun<std::complex<T> > solve(Chebfun<std::complex<T> > in){};
};

/// Differentiation operator for Linop, to differentiate n times.
template <class T> Linop<T> diff(Linop<T> in, int n) {
  int bre;
  if (n == 0) {
    return in;
  } else {
    if (in.NCC == 0) { // if constant coefficients
      Linop<T> out(in.n + 1);
      out.coef << in.coef, 0.0;
      return diff(out, n - 1);
    } else {
      Linop<T> out;
      out.ncc(in.n + 1);
      out.coefFun[0] = in.coefFun[0];
      out.coefFun[out.n] = diff(in.coefFun[in.n], 1);
      for (int i = 1; i < out.n; i++) {
        out.coefFun[i] = diff(in.coefFun[i - 1], 1) + in.coefFun[i];
      }
      return diff(out, n - 1);
    }
  }
};

/// \brief This class stores functions and values needed to sort Eigenvalues.
template <class T> class EigenSorter {
private:
public:
  T epsilon; ///< Machine zero, needs to be set based on Matrix computation
             /// errors
  bool Mp; ///< \brief Stores true if machine precision is reached, else false.
  int numMp; ///< \brief Stores number of Chebyshev coefficients in which MP is
             /// reached, if not reached, it is set to N+1.
  T minVal;  ///< \brief Stores the smallest value in the Chebyshev series
  /// in case Machine precision is not reached. If Machine
  /// precision is reached, it is set to zero.
  bool isInf;

  /// \brief Null constructor
  EigenSorter() {
    numMp = 0;
    minVal = 0;
    epsilon = 1e-11;
    epsilon = abs(epsilon);
  }

  /// \brief Computes Mp and numMp for Linops
  void compute(Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> in) {
    std::vector<T> temp;
    temp.reserve(in.rows());
    for (int i = 0; i < in.rows(); i++) {
      temp.push_back(std::abs(in[i]));
    }
    if (floor(N / 2) < N / 2) {
      numMp = 0;
      while (numMp < N && ((temp[numMp] + temp[numMp + 1]) / 2.0) > epsilon) {
        numMp = numMp + 2;
      }
    } else {
      numMp = 1;
      while (numMp < N && ((temp[numMp] + temp[numMp + 1]) / 2.0) > epsilon) {
        numMp = numMp + 2;
      }
    }
    Chebfun<std::complex<T> > tempfun(in);
    tempfun.dct_flag = SIS_CHEB_SPACE;
    // Eigenfunction has to be nonzero. If zero, set the numMp to N+1
    if (tempfun.L2norm() < 1e-11) {
      numMp = N + 1;
    }
    if (std::isnan(tempfun.L2norm())) {
      numMp = N + 1;
    }
    if (numMp >= N) {
      minVal = (temp[N] + temp[N - 1]) / 2.0;
      Mp = false;
    } else {
      Mp = true;
    }
  }

  /// \brief Computes Mp and numMp for LinopMats
  void compute(Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> in, int c) {
    int bre;
    std::vector<T> temp;
    temp.reserve(N + 1);
    ChebfunMat<std::complex<T> > tempfuns;
    tempfuns.resize(c, 1);
    std::vector<std::vector<double> > L2norms;
    double ave_L2norm = 0.0;
    L2norms.resize(c);
    for (int i = 0; i < c; i++) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> tempMatr, tempMati;
      tempMatr = in.block(i * (N + 1), 0, N + 1, 1).real();
      tempMati = in.block(i * (N + 1), 0, N + 1, 1).imag();
      for (int j = 0; j < N + 1; j++) {
        tempfuns(i, 0).v[j] = std::complex<T>(tempMatr[j], tempMati[j]);
      }
      tempfuns(i, 0).dct_flag = SIS_CHEB_SPACE;
      L2norms[i].resize(2);
      L2norms[i][0] = (tempfuns(i, 0).L2norm());
      L2norms[i][1] = i;
      ave_L2norm += L2norms[i][0];
    }

    // Use the Chebfun with the largest L2norm as the one that ascertains
    // sorting criteria
    std::sort(L2norms.begin(), L2norms.end());
    ave_L2norm = ave_L2norm / double(c);

    if (ave_L2norm > 1e-11) {
      for (int i = int(L2norms[c - 1][1]) * (N + 1);
           i < int(L2norms[c - 1][1] + 1) * (N + 1); i++) {
        temp.push_back(std::abs(in[i]));
      }
      // This checks if N is odd or even:
      if (floor(N / 2) < N / 2) {
        numMp = 0;
        while (numMp < N && ((temp[numMp] + temp[numMp + 1]) / 2.0) > epsilon) {
          numMp = numMp + 2;
        }
      } else {
        numMp = 1;
        while (numMp < N && ((temp[numMp] + temp[numMp + 1]) / 2.0) > epsilon) {
          numMp = numMp + 2;
        }
      }
    } // Eigenfunction-vector has to be nonzero. If zero, set the numMp to N+1
    else {
      numMp = N + 1;
    }

    if (std::isnan(ave_L2norm)) {
      numMp = N + 1;
    }

    if (numMp >= N) {
      minVal = (temp[N] + temp[N - 1]) / 2.0;
      Mp = false;
    } else {
      Mp = true;
    }
  }
};

/// \brief This class will solve the generalized eigenvalue problem for two
/// linear operators. One of them can be singular.
///
/// Solves for \f[ \boldmath{L}\phi = \lambda \boldmath{M} \phi \f].
///
///

template <class T> class GeneralizedEigenSolver {
private:
public:
  std::vector<ChebfunMat<std::complex<T> > > eigenvectorsMat;
  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> eigenvalues;
  /// \brief Null constructor
  GeneralizedEigenSolver(){};
/// \brief Call this with an input Linear operator to solve for eigenvalues
/// and vectors. The number of Eigen values/vectors is num_vals, num_vals
/// has to be less than N.
#ifdef SIS_USE_LAPACK
  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> alpha, beta;
  int info; // gives info on lapack solver.
#endif
  void compute(Linop<T> L, Linop<T> M, int num_vals){};

  /// \brief This will use do the same work as compute(), but will overide all
  /// boundary conditions specified through BcVec and use the input BcMat
  /// instead. Read about class BcMat to see how this works, also see examples
  /// example/Ex_16.cpp and test/Visco_3D_pipe.cpp of how this is applied.
  /// This function is useful when boundary conditions are mixed between
  /// variables. Also, if you have boundary conditions that have an associated
  /// eigenvalue, see compute_with_constraints.
  void compute(const LinopMat<T> &Lmat_, const LinopMat<T> &Mmat_, int num_vals,
               const BcMat<T> &Lbc_) {
    LinopMat<T> Lmat = Lmat_;
    LinopMat<T> Mmat = Mmat_;
    BcMat<T> Lbc = Lbc_;
    int bre;
    int total_of_all_orders = 0;
    int total_boundary_conditions = 0;
    if (Lmat.r != Lmat.c) {
      std::cout << "Solution only possible with square LinopMats. Exiting ..."
                << '\n';
      exit(1);
    }
    if (Mmat.r != Mmat.c) {
      std::cout << "Solution only possible with square LinopMats. Exiting ..."
                << '\n';
      exit(1);
    }
    if (Mmat.c != Lmat.c) {
      std::cout << "Both matrices have to be of same size. Exiting ..." << '\n';
      exit(1);
    }
    int r = Lmat.r, c = Lmat.c;
    // Find the highest derivative in each column. To do this create a vector
    // highest_each_column, and a temp_int_vec that will hold all values of a
    // given column to and the maximum will be stored in highest_each_column
    std::vector<int> highest_each_columnL, highest_each_columnM,
        highest_each_column;
    highest_each_columnL.resize(c);
    highest_each_columnM.resize(c);
    highest_each_column.resize(c);

    std::vector<int> temp_vec_int;
    temp_vec_int.resize(r);
    for (int j = 0; j < c; j++) {
      for (int i = 0; i < r; i++) {
        temp_vec_int[i] = Lmat(i, j).n;
      }
      highest_each_columnL[j] =
          *std::max_element(temp_vec_int.begin(), temp_vec_int.end());
    }
    for (int j = 0; j < c; j++) {
      for (int i = 0; i < r; i++) {
        temp_vec_int[i] = Mmat(i, j).n;
      }
      highest_each_columnM[j] =
          *std::max_element(temp_vec_int.begin(), temp_vec_int.end());
    }
    for (int i = 0; i < c; i++) {
      total_of_all_orders += (highest_each_columnL[i] > highest_each_columnM[i])
                                 ? highest_each_columnL[i]
                                 : highest_each_columnM[i];
      highest_each_column[i] =
          (highest_each_columnL[i] > highest_each_columnM[i])
              ? highest_each_columnL[i]
              : highest_each_columnM[i];
    }
    total_boundary_conditions = Lbc.m;

    // total_of_all_orders has to be equal to total number of boundary
    // conditions, else the problem is ill-posed, if ill-posed, cout the same
    // and exit.
    if (total_of_all_orders != total_boundary_conditions) {
      std::cout << "The problem is ill-posed, the total of the highest "
                   "orders of all "
                   "dependent variables has to be equal to the total number of "
                   "boundary conditions specified."
                   "\n "
                   "Total no. of boundary conditions: "
                << total_boundary_conditions
                << "\n"
                   "Total of all orders: "
                << total_of_all_orders << "\n Exiting ...\n";
      exit(1);
    }
    // Declare the master matrix L:
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> masterL(
        r * (N + 1), c * (N + 1) + total_boundary_conditions);
    // Declare the master matrix M:
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> masterM(
        r * (N + 1), c * (N + 1) + total_boundary_conditions);

    masterL.setConstant(0.0);
    masterM.setConstant(0.0);
    MatGen<T> Mat;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> subs_mat;
    subs_mat.resize(total_boundary_conditions, (N + 1) * c);
    subs_mat.setConstant(0.0);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat_temp(
        total_boundary_conditions, total_boundary_conditions),
        constraints(total_boundary_conditions,
                    c * (N + 1) + total_of_all_orders);

    mat_temp.setConstant(0.0);
    int row_counter = 0, col_counter = 0;

    for (int i = 0; i < Lbc.m; i++) {
      for (int j = 0; j < Lbc.n; j++) {
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> temp =
            Lbc(i, j, highest_each_column[j]);
        constraints.block(i, col_counter, 1, temp.cols()) = temp;
        col_counter += temp.cols();
      }
      col_counter = 0;
    }

    Eigen::ColPivHouseholderQR<
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >
        qr;
    qr.compute(constraints);
    if (qr.rank() != constraints.rows()) {
      std::cout << "The boundary conditions supplied are not "
                << "  linearly independent." << '\n';
      std::cout << "qr.rank = " << qr.rank()
                << ", no. bcs: " << total_boundary_conditions << ". Exiting ..."
                << '\n';
      exit(1);
    }

    // Permutation matrix
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P;
    P = qr.colsPermutation();
    // Permute constraints
    constraints = constraints * P;
    mat_temp = constraints.block(0, 0, constraints.rows(), constraints.rows());

    Eigen::ColPivHouseholderQR<
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >
        consSolver(mat_temp);
    subs_mat = -consSolver.inverse() *
               constraints.block(0, constraints.rows(), constraints.rows(),
                                 constraints.cols() - constraints.rows());

    row_counter = 0;
    col_counter = 0;
    int master_row_counter = 0;
    int master_col_counter = 0;
    for (int j = 0; j < c; j++) {
      int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                  ? highest_each_columnL[j]
                  : highest_each_columnM[j];
      Mat.compute(n);
      for (int i = 0; i < r; i++) {
        int diffn = n - Lmat(i, j).n;
        if (Lmat(i, j).NCC == 0) {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            masterL.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Lmat(i, j).coef[k] *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, N + 1 + n));
          }
        } else {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            masterL.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Lmat(i, j).coefFun[k].MultMat().block(0, 0, N + 1, N + 1) *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, N + 1 + n));
          }
        }
        diffn = n - Mmat(i, j).n;
        if (Mmat(i, j).NCC == 0) {
          for (int k = 0; k < Mmat(i, j).n + 1; k++) {
            masterM.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Mmat(i, j).coef[k] *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, N + 1 + n));
          }
        } else {
          for (int k = 0; k < Mmat(i, j).n + 1; k++) {
            masterM.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Mmat(i, j).coefFun[k].MultMat().block(0, 0, N + 1, N + 1) *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, N + 1 + n));
          }
        }
        master_row_counter += N + 1;
      }
      master_row_counter = 0;
      master_col_counter += N + 1 + n;
    }
    // Permute columns of M and L:
    masterL = masterL * P;
    masterM = masterM * P;

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> masterL2, masterM2;

    masterL2 =
        masterL.block(0, constraints.rows(), c * (N + 1), c * (N + 1)) +
        (masterL.block(0, 0, c * (N + 1), constraints.rows()) * subs_mat);
    masterM2 =
        masterM.block(0, constraints.rows(), c * (N + 1), c * (N + 1)) +
        (masterM.block(0, 0, c * (N + 1), constraints.rows()) * subs_mat);
#ifndef SIS_USE_LAPACK
    Eigen::ColPivHouseholderQR<
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >
        solver(masterM2);
    Eigen::EigenSolver<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > eigs;
    bool Is_M_Invertible;
    if (solver.isInvertible()) {
      //  std::cout << "M is invertible" << '\n';
      //  std::cin >> bre;
      Is_M_Invertible = true;
      eigs.compute(solver.inverse() * masterL2);
      // std::cout << "Eigenvalues :\n" << eigs.eigenvalues() << '\n';
    } else {
      // std::cout << "M is not invertible." << '\n';
      // std::cin >> bre;
      Is_M_Invertible = false;
      solver.compute(masterL2);
      eigs.compute(solver.inverse() * masterM2);
      // std::cout << "Eigenvalues :\n" << 1 / eigs.eigenvalues().array() <<
      // '\n';
    }
    // std::cin >> bre;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
        evecMat_master(c * (N + 1), c * (N + 1));
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> temp_vec;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> fir_n_mat =
        subs_mat * eigs.eigenvectors();
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> evecMat(
        c * (N + 1) + total_boundary_conditions, c * (N + 1));
    evecMat << fir_n_mat, eigs.eigenvectors();

    // Unpermute eveMat:
    evecMat = P * evecMat;
    row_counter = 0;
    col_counter = 0;
    for (int i = 0; i < c; i++) {
      int n = (highest_each_columnL[i] > highest_each_columnM[i])
                  ? highest_each_columnL[i]
                  : highest_each_columnM[i];
      Mat.compute(n);
      evecMat_master.block(i * (N + 1), 0, N + 1, c * (N + 1)) =
          Mat.mats2[n] * evecMat.block(row_counter, 0, N + 1 + n, c * (N + 1));
      row_counter += N + 1 + n;
    }
    std::vector<std::vector<T> > eigenval_trunc_sorter;

    eigenval_trunc_sorter.resize(c * (N + 1));
    for (int i = 0; i < c * (N + 1); i++) {
      eigenval_trunc_sorter[i].resize(2);
    }

    Eigen::Matrix<int, Eigen::Dynamic, 1> MPcount(c * (N + 1));
    MPcount.setConstant(0);
    std::vector<EigenSorter<T> > sorter;
    sorter.resize(c * (N + 1));
    for (int i = 0; i < c * (N + 1); i++) {
      sorter[i].compute(evecMat_master.block(0, i, c * (N + 1), 1), c);
      if (!Is_M_Invertible) {
        // If M is not invertible, there will be zero eigenvalues, so that 1/0 =
        // inf eigenvalue, remove them.
        if (abs(eigs.eigenvalues()[i]) < 1e-11) {
          sorter[i].numMp = N + 1;
          sorter[i].Mp = false;
          sorter[i].minVal = 1e11;
        }
      }
      // sorter[i].compute(evecMat_master.block(0, i, (N + 1), 1));
      if (sorter[i].Mp == true) {
        MPcount[i] = 1;
      }
    }
    int MPcount_sum = MPcount.sum();
    eigenvectorsMat.resize(num_vals);
    for (int i = 0; i < num_vals; i++) {
      eigenvectorsMat[i].resize(c, 1);
    }
    eigenvalues.resize(num_vals);
    if (MPcount_sum >= num_vals) {
      for (int i = 0; i < c * (N + 1); i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
        // std::cout << eigs.eigenvalues()[i] << " " << sorter[i].numMp <<
        // "\n"; std::cin > > bre;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < num_vals; i++) {
        if (Is_M_Invertible) {
          eigenvalues[i] = eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
        } else {
          eigenvalues[i] =
              1.0 / eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
        }
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];

        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < num_vals; i++) {
          eigenvectorsMat[i][j] = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
        }
      }
    } else {
      std::cout << "Last " << num_vals - MPcount_sum
                << " eigenvectors are not resolved to machine precision."
                << '\n';
      for (int i = 0; i < c * (N + 1); i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < MPcount_sum; i++) {
        if (Is_M_Invertible) {
          eigenvalues[i] = eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
        } else {
          eigenvalues[i] =
              1.0 / eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
        }
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];
        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < MPcount_sum; i++) {
          eigenvectorsMat[i][j] = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
        }
      }

      std::vector<std::vector<T> > eigenval_trunc_sorter2;
      eigenval_trunc_sorter2.resize(c * (N + 1) - MPcount_sum);
      for (int i = MPcount_sum; i < c * (N + 1); i++) {
        eigenval_trunc_sorter2[i - MPcount_sum].resize(2);
        eigenval_trunc_sorter2[i - MPcount_sum][0] =
            sorter[int(eigenval_trunc_sorter[i][1])].minVal;
        eigenval_trunc_sorter2[i - MPcount_sum][1] =
            eigenval_trunc_sorter[i][1];
      }
      std::sort(eigenval_trunc_sorter2.begin(), eigenval_trunc_sorter2.end());
      for (int i = MPcount_sum; i < num_vals; i++) {
        if (Is_M_Invertible) {
          eigenvalues[i] = eigs.eigenvalues()[int(
              eigenval_trunc_sorter2[i - MPcount_sum][1])];
        } else {
          eigenvalues[i] =
              1.0 / eigs.eigenvalues()[int(
                        eigenval_trunc_sorter2[i - MPcount_sum][1])];
        }
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];
        Mat.compute(n);
        temp_vec.resize(N + 1 + n, 1);
        for (int i = MPcount_sum; i < num_vals; i++) {
          eigenvectorsMat[i][j] = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter2[i - MPcount_sum][1]),
              N + 1, 1);
        }
      }
    }
#else
    std::cout << "Using lapack routine..." << '\n';
    char jobvl = 'N';                 // Don't compute left evecs
    char jobvr = 'V';                 // Compute right evecs
    std::complex<double> wkopt;       // Eistimate optimum workspace
    std::complex<double> *work;       // allocate optimum workspace
    alpha.resize(masterL2.rows(), 1); // alpha for gen. eig. prob.
    beta.resize(masterL2.rows(), 1);  // beta for gen. eig. prob.

    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> masterL2_,
        masterM2_, vl(masterL2.rows(), masterL2.rows()),
        vr(masterL2.rows(), masterL2.rows()),
        eigenvalues_temp(masterL2.rows(), 1), alpha_temp(masterL2.rows(), 1),
        beta_temp(masterL2.rows(), 1);

    masterL2_ = masterL2;
    masterM2_ = masterM2;
    // vl : left evecs, vr: right evecs.
    int ldL = masterL2.outerStride(); // ld for leading dimension
    int ldM = masterM2.outerStride();
    int ldvl = vl.outerStride();
    int ldvr = vr.outerStride();
    int sizeL = masterL2.rows();
    int lwork = -1; // set lwork to -1 to estimate workspace.
    double rwork[8 * sizeL];

    // call this to estimate workspace
    zggev_(&jobvl, &jobvr, &sizeL, masterL2_.data(), &ldL, masterM2_.data(),
           &ldM, alpha_temp.data(), beta_temp.data(), vl.data(), &ldvl,
           vr.data(), &ldvr, &wkopt, &lwork, rwork, &info);

    // Now allocate workspace:
    lwork = (int)real(wkopt);
    work = (std::complex<double> *)malloc(lwork * sizeof(std::complex<double>));

    // Solve eigenvalue problem:
    zggev_(&jobvl, &jobvr, &sizeL, masterL2_.data(), &ldL, masterM2_.data(),
           &ldM, alpha_temp.data(), beta_temp.data(), vl.data(), &ldvl,
           vr.data(), &ldvr, work, &lwork, rwork, &info);

    // Free workspace.
    free((void *)work);

    eigenvalues_temp = alpha_temp.array() / beta_temp.array();

    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
        evecMat_master(c * (N + 1), c * (N + 1));
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> temp_vec;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> fir_n_mat =
        subs_mat * vr;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> evecMat(
        c * (N + 1) + total_boundary_conditions, c * (N + 1));
    evecMat << fir_n_mat, vr;
    // Unpermute eveMat:
    evecMat = P * evecMat;
    row_counter = 0;
    col_counter = 0;
    for (int i = 0; i < c; i++) {
      int n = (highest_each_columnL[i] > highest_each_columnM[i])
                  ? highest_each_columnL[i]
                  : highest_each_columnM[i];
      Mat.compute(n);
      evecMat_master.block(i * (N + 1), 0, N + 1, c * (N + 1)) =
          Mat.mats2[n] * evecMat.block(row_counter, 0, N + 1 + n, c * (N + 1));
      row_counter += N + 1 + n;
    }
    std::vector<std::vector<T> > eigenval_trunc_sorter;

    eigenval_trunc_sorter.resize(c * (N + 1));
    for (int i = 0; i < c * (N + 1); i++) {
      eigenval_trunc_sorter[i].resize(2);
    }

    Eigen::Matrix<int, Eigen::Dynamic, 1> MPcount(c * (N + 1));
    MPcount.setConstant(0);
    std::vector<EigenSorter<T> > sorter;
    sorter.resize(c * (N + 1));
    for (int i = 0; i < c * (N + 1); i++) {
      sorter[i].compute(evecMat_master.block(0, i, c * (N + 1), 1), c);
      // sorter[i].compute(evecMat_master.block(0, i, (N + 1), 1));
      if (sorter[i].Mp == true) {
        MPcount[i] = 1;
      }
    }
    if (num_vals > masterL2.rows()) {
      std::cout << "Only " << masterL.rows()
                << " eigenvalues can be calculated."
                << "Storing only that many." << '\n';
      num_vals = masterL2.rows();
    }

    int MPcount_sum = MPcount.sum();
    eigenvectorsMat.resize(num_vals);
    for (int i = 0; i < num_vals; i++) {
      eigenvectorsMat[i].resize(c, 1);
    }
    eigenvalues.resize(num_vals);
    beta.resize(num_vals);
    alpha.resize(num_vals);
    if (MPcount_sum >= num_vals) {
      for (int i = 0; i < c * (N + 1); i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
        // std::cout << eigs.eigenvalues()[i] << " " << sorter[i].numMp <<
        // "\n"; std::cin > > bre;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < num_vals; i++) {
        eigenvalues[i] = eigenvalues_temp(int(eigenval_trunc_sorter[i][1]), 0);
        alpha[i] = alpha_temp(int(eigenval_trunc_sorter[i][1]), 0);
        beta[i] = beta_temp(int(eigenval_trunc_sorter[i][1]), 0);
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];

        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < num_vals; i++) {
          eigenvectorsMat[i][j] = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }
    } else {
      std::cout << "Last " << num_vals - MPcount_sum
                << " eigenvectors are not resolved to machine precision."
                << '\n';
      for (int i = 0; i < c * (N + 1); i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < MPcount_sum; i++) {
        eigenvalues[i] = eigenvalues_temp(int(eigenval_trunc_sorter[i][1]), 0);
        alpha[i] = alpha_temp(int(eigenval_trunc_sorter[i][1]), 0);
        beta[i] = beta_temp(int(eigenval_trunc_sorter[i][1]), 0);
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];
        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < MPcount_sum; i++) {
          eigenvectorsMat[i][j] = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }

      std::vector<std::vector<T> > eigenval_trunc_sorter2;
      eigenval_trunc_sorter2.resize(c * (N + 1) - MPcount_sum);
      for (int i = MPcount_sum; i < c * (N + 1); i++) {
        eigenval_trunc_sorter2[i - MPcount_sum].resize(2);
        eigenval_trunc_sorter2[i - MPcount_sum][0] =
            sorter[int(eigenval_trunc_sorter[i][1])].minVal;
        eigenval_trunc_sorter2[i - MPcount_sum][1] =
            eigenval_trunc_sorter[i][1];
      }
      std::sort(eigenval_trunc_sorter2.begin(), eigenval_trunc_sorter2.end());
      for (int i = MPcount_sum; i < num_vals; i++) {
        eigenvalues[i] = eigenvalues_temp(
            int(eigenval_trunc_sorter2[i - MPcount_sum][1]), 0);
        alpha[i] =
            alpha_temp(int(eigenval_trunc_sorter2[i - MPcount_sum][1]), 0);
        beta[i] = beta_temp(int(eigenval_trunc_sorter2[i - MPcount_sum][1]), 0);
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];
        Mat.compute(n);
        temp_vec.resize(N + 1 + n, 1);
        for (int i = MPcount_sum; i < num_vals; i++) {
          eigenvectorsMat[i][j] = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter2[i - MPcount_sum][1]),
              N + 1, 1);
          eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }
    }

#endif
  }


  /// Use this to solve eigenvalue problems where eigenvalues appears in the
  /// boundary condition, example, in eigenvalue problems with fluid-fluid
  /// interfaces. Again
  /// solvability must hold. That total of all the highest orders of every
  /// independent variable has to be equal to the total number of boundary
  /// conditions specified + the number of constraints.
  /// Further, please ensure that the boundary conditions and constraints are
  /// linearly independent, else this will throw an error.
  void compute_with_constraints(const LinopMat<T> &Lmat_,
                                const LinopMat<T> &Mmat_, int num_vals,
                                const BcMat<T> &Lbc_, const BcMat<T> &Mbc_) {

    LinopMat<T> Lmat = Lmat_;
    LinopMat<T> Mmat = Mmat_;
    BcMat<T> Lbc = Lbc_;
    BcMat<T> Mbc = Mbc_;
    int bre;
    int total_of_all_orders = 0;
    int total_boundary_conditions = 0;
    if (Lmat.r != Lmat.c) {
      std::cout << "Solution only possible with square LinopMats. Exiting ..."
                << '\n';
      exit(1);
    }
    if (Mmat.r != Mmat.c) {
      std::cout << "Solution only possible with square LinopMats. Exiting ..."
                << '\n';
      exit(1);
    }
    if (Mmat.c != Lmat.c) {
      std::cout << "Both matrices have to be of same size. Exiting ..." << '\n';
      exit(1);
    }
    int r = Lmat.r, c = Lmat.c;
    // Find the highest derivative in each column. To do this create a vector
    // highest_each_column, and a temp_int_vec that will hold all values of a
    // given column to and the maximum will be stored in highest_each_column
    std::vector<int> highest_each_columnL, highest_each_columnM,
        highest_each_column, num_bc_each_var;
    highest_each_columnL.resize(c);
    highest_each_columnM.resize(c);
    highest_each_column.resize(c);
    num_bc_each_var.resize(c);
    std::vector<int> temp_vec_int;
    temp_vec_int.resize(r);
    for (int j = 0; j < c; j++) {
      for (int i = 0; i < r; i++) {
        temp_vec_int[i] = Lmat(i, j).n;
      }
      highest_each_columnL[j] =
          *std::max_element(temp_vec_int.begin(), temp_vec_int.end());
    }
    for (int j = 0; j < c; j++) {
      for (int i = 0; i < r; i++) {
        temp_vec_int[i] = Mmat(i, j).n;
      }
      highest_each_columnM[j] =
          *std::max_element(temp_vec_int.begin(), temp_vec_int.end());
    }
    for (int i = 0; i < c; i++) {
      total_of_all_orders += (highest_each_columnL[i] > highest_each_columnM[i])
                                 ? highest_each_columnL[i]
                                 : highest_each_columnM[i];
      highest_each_column[i] =
          (highest_each_columnL[i] > highest_each_columnM[i])
              ? highest_each_columnL[i]
              : highest_each_columnM[i];
    }
    for (int i = 0; i < c; i++) {
      total_boundary_conditions += Lmat.BcVec[i].nbc();
      num_bc_each_var[i] = Lmat.BcVec[i].nbc();
    }

    if (Lbc.m != Mbc.m) {
      std::cout << "The Lbc and Mbc have to be of same dimensions" << '\n'
                << "Exiting in " << __LINE__ << "...\n";

      exit(1);
    }
    // total_of_all_orders has to be equal to total number of boundary
    // conditions, else the problem is ill-posed, if ill-posed, cout the same
    // and exit.
    if (total_of_all_orders != (total_boundary_conditions + Lbc.m)) {
      std::cout << "The problem is ill-posed, the total of the highest "
                   "orders of all "
                   "dependent variables has to be equal to the total number of "
                   "boundary conditions plus number of constraints specified."
                   "\n "
                   "Total no. of boundary conditions: "
                << total_boundary_conditions
                << "\n"
                   "Total of all orders: "
                << total_of_all_orders << "\n Number of Constraints : " << Lbc.m
                << "\n Exiting in " << __LINE__ << " ...\n ";
      exit(1);
    }
    // Declare the master matrix L:
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> masterL(
        r * (N + 1) + Lbc.m, c * (N + 1) + Lbc.m);
    // Declare the master matrix M:
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> masterM(
        r * (N + 1) + Lbc.m, c * (N + 1) + Lbc.m);

    masterL.setConstant(0.0);
    masterM.setConstant(0.0);
    MatGen<T> Mat;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> subs_mat;

    subs_mat.resize(total_boundary_conditions, (N + 1) * c + Lbc.m);

    subs_mat.setConstant(0.0);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat_temp(
        total_boundary_conditions, total_boundary_conditions);
    mat_temp.setConstant(0.0);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> subs_mat_col_mul(
        c * (N + 1), total_boundary_conditions);
    subs_mat_col_mul.setConstant(0.0);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> subs_mat_col_mul2(
        c * (N + 1), total_boundary_conditions);
    subs_mat_col_mul2.setConstant(0.0);
    int subs_mat_row_counter = 0;
    int subs_mat_col_counter = 0;
    int row_counter = 0, col_counter = 0;
    for (int j = 0; j < c; j++) {
      // overall order is the one which is greater in each column:
      int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                  ? highest_each_columnL[j]
                  : highest_each_columnM[j];
      // Store the LinopMat with the larger degree to setup boundary conditions
      // LinopMat<T > Tmat =
      //    (highest_each_columnL[j] > highest_each_columnM[j]) ? Lmat : Mmat;

      Mat.compute(n);
      std::vector<Eigen::Matrix<T, 1, Eigen::Dynamic> > lbc_vecs(
          Lmat.BcVec[j].nl, Eigen::Matrix<T, 1, Eigen::Dynamic>(N + 1 + 2 * n));
      std::vector<Eigen::Matrix<T, 1, Eigen::Dynamic> > rbc_vecs(
          Lmat.BcVec[j].nr, Eigen::Matrix<T, 1, Eigen::Dynamic>(N + 1 + 2 * n));
      std::vector<Eigen::Matrix<T, 1, Eigen::Dynamic> > lbc_con_vecs(
          Lmat.BcVec[j].nl, Eigen::Matrix<T, 1, Eigen::Dynamic>(n));
      std::vector<Eigen::Matrix<T, 1, Eigen::Dynamic> > rbc_con_vecs(
          Lmat.BcVec[j].nr, Eigen::Matrix<T, 1, Eigen::Dynamic>(n));
      Eigen::Matrix<T, 1, Eigen::Dynamic> ones(N + 1 + 2 * n);
      Eigen::Matrix<T, 1, Eigen::Dynamic> onem(N + 1 + 2 * n);
      Eigen::Matrix<T, 1, Eigen::Dynamic> onesn(n);
      Eigen::Matrix<T, 1, Eigen::Dynamic> onemn(n);
      // This follows as T_n(1) = 1.0.
      ones.setConstant(1.0);
      onesn.setConstant(1.0);
      onem.setConstant(1.0);
      onemn.setConstant(1.0);
      // This follows as T_n(-1) = (-1)^n.
      for (int k = 1; k < N + 1 + 2 * n; k = k + 2) {
        onem[k] = -1.0;
      }
      for (int k = 1; k < n; k = k + 2) {
        onemn[k] = -1.0;
      }
      if (n > 0) {
        ones[0] = 0.5;
        onesn[0] = 0.5;
        onem[0] = 0.5;
        onemn[0] = 0.5;
      }
      // Next just multiply matrices based on the order.
      for (int k = 0; k < Lmat.BcVec[j].nl; k++) {
        lbc_vecs[k].resize(N + 1 + 2 * n);
        lbc_vecs[k].setConstant(0.0);
        for (int l = 0; l < Lmat.BcVec[j].ord + 1; l++) {
          lbc_vecs[k] += Lmat.BcVec[j].coefl(k, l) *
                         (onem * Mat.mats[l + (n - Lmat.BcVec[j].ord)]);
        }
      }

      for (int k = 0; k < Lmat.BcVec[j].nr; k++) {
        rbc_vecs[k].resize(N + 1 + 2 * n);
        rbc_vecs[k].setConstant(0.0);
        for (int l = 0; l < Lmat.BcVec[j].ord + 1; l++) {
          rbc_vecs[k] += Lmat.BcVec[j].coefr(k, l) *
                         (ones * Mat.mats[l + (n - Lmat.BcVec[j].ord)]);
        }
      }
      for (int k = 0; k < Lmat.BcVec[j].nl; k++) {
        lbc_con_vecs[k].setConstant(0.0);
        for (int l = 0; l < Lmat.BcVec[j].ord + 1; l++) {
          if (l + n - Lmat.BcVec[j].ord - 1 > -1) {
            lbc_con_vecs[k] +=
                Lmat.BcVec[j].coefl(k, l) *
                (onemn * Mat.con_mats[l + n - Lmat.BcVec[j].ord - 1]);
          }
        }
      }
      for (int k = 0; k < Lmat.BcVec[j].nr; k++) {
        rbc_con_vecs[k].setConstant(0.0);
        for (int l = 0; l < Lmat.BcVec[j].ord + 1; l++) {
          if (l + n - Lmat.BcVec[j].ord - 1 > -1) {
            rbc_con_vecs[k] +=
                Lmat.BcVec[j].coefr(k, l) *
                (onesn * Mat.con_mats[l + n - Lmat.BcVec[j].ord - 1]);
          }
        }
      }

      for (int k = 0; k < Lmat.BcVec[j].nl; k++) {
        mat_temp.block(row_counter, col_counter, 1,
                       Lmat.BcVec[j].nl + Lmat.BcVec[j].nr) =
            lbc_vecs[k].head(Lmat.BcVec[j].nl + Lmat.BcVec[j].nr);
        row_counter++;
      }
      for (int k = Lmat.BcVec[j].nl; k < Lmat.BcVec[j].nl + Lmat.BcVec[j].nr;
           k++) {
        mat_temp.block(row_counter, col_counter, 1,
                       Lmat.BcVec[j].nl + Lmat.BcVec[j].nr) =
            rbc_vecs[k - Lmat.BcVec[j].nl].head(Lmat.BcVec[j].nl +
                                                Lmat.BcVec[j].nr);
        row_counter++;
      }
      col_counter = col_counter + Lmat.BcVec[j].nbc();
      for (int k = 0; k < Lmat.BcVec[j].nl; k++) {
        subs_mat.block(subs_mat_row_counter, subs_mat_col_counter, 1,
                       (N + 1 - Lmat.BcVec[j].nbc())) =
            lbc_vecs[k].block(0, Lmat.BcVec[j].nbc(), 1,
                              N + 1 - Lmat.BcVec[j].nbc());
        subs_mat.block(subs_mat_row_counter,
                       subs_mat_col_counter + (N + 1 - Lmat.BcVec[j].nbc()), 1,
                       n) = lbc_con_vecs[k];
        subs_mat_row_counter++;
      }

      for (int k = Lmat.BcVec[j].nl; k < Lmat.BcVec[j].nbc(); k++) {
        subs_mat.block(subs_mat_row_counter, subs_mat_col_counter, 1,
                       (N + 1 - Lmat.BcVec[j].nbc())) =
            rbc_vecs[k - Lmat.BcVec[j].nl].block(0, Lmat.BcVec[j].nbc(), 1,
                                                 N + 1 - Lmat.BcVec[j].nbc());
        subs_mat.block(subs_mat_row_counter,
                       subs_mat_col_counter + (N + 1 - Lmat.BcVec[j].nbc()), 1,
                       n) = rbc_con_vecs[k - Lmat.BcVec[j].nl];
        subs_mat_row_counter++;
      }
      subs_mat_col_counter += (N + 1 + n - Lmat.BcVec[j].nbc());
    }
    if (mat_temp.rows() != 0 || mat_temp.cols() != 0) {
      Eigen::ColPivHouseholderQR<
          Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >
          consSolver(mat_temp);
      if (!consSolver.isInvertible()) {
        std::cout << "the matrix is not invertible." << '\n';
        exit(1);
      }
      subs_mat = -consSolver.inverse() * subs_mat.eval();
    }
    row_counter = 0;
    col_counter = 0;
    int master_row_counter = 0;
    int master_col_counter = 0;
    subs_mat_row_counter = 0;
    subs_mat_col_counter = 0;
    for (int j = 0; j < c; j++) {
      int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                  ? highest_each_columnL[j]
                  : highest_each_columnM[j];
      Mat.compute(n);
      for (int i = 0; i < r; i++) {

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> solver_matL(
            N + 1, N + 1 + n - Lmat.BcVec[j].nbc());
        solver_matL.setConstant(0.0);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> solver_matM(
            N + 1, N + 1 + n - Lmat.BcVec[j].nbc());
        solver_matM.setConstant(0.0);
        int diffn = n - Lmat(i, j).n;
        if (Lmat(i, j).NCC == 0) {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            solver_matL.block(0, 0, N + 1, N + 1 + n - Lmat.BcVec[j].nbc()) +=
                Lmat(i, j).coef[k] *
                (Mat.mats2[k + diffn].block(0, Lmat.BcVec[j].nbc(), N + 1,
                                            N + 1 + n - Lmat.BcVec[j].nbc()));
            subs_mat_col_mul.block(i * (N + 1), subs_mat_col_counter, N + 1,
                                   Lmat.BcVec[j].nbc()) +=
                Lmat(i, j).coef[k] *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, Lmat.BcVec[j].nbc()));
          }
        } else {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            solver_matL.block(0, 0, N + 1, N + 1 + n - Lmat.BcVec[j].nbc()) +=
                Lmat(i, j).coefFun[k].MultMat().block(0, 0, N + 1, N + 1) *
                (Mat.mats2[k + diffn].block(0, Lmat.BcVec[j].nbc(), N + 1,
                                            N + 1 + n - Lmat.BcVec[j].nbc()));
            subs_mat_col_mul.block(i * (N + 1), subs_mat_col_counter, N + 1,
                                   Lmat.BcVec[j].nbc()) +=
                Lmat(i, j).coefFun[k].MultMat().block(0, 0, N + 1, N + 1) *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, Lmat.BcVec[j].nbc()));
          }
        }
        diffn = n - Mmat(i, j).n;
        if (Mmat(i, j).NCC == 0) {
          for (int k = 0; k < Mmat(i, j).n + 1; k++) {
            solver_matM.block(0, 0, N + 1, N + 1 + n - Lmat.BcVec[j].nbc()) +=
                Mmat(i, j).coef[k] *
                (Mat.mats2[k + diffn].block(0, Lmat.BcVec[j].nbc(), N + 1,
                                            N + 1 + n - Lmat.BcVec[j].nbc()));
            subs_mat_col_mul2.block(i * (N + 1), subs_mat_col_counter, N + 1,
                                    Lmat.BcVec[j].nbc()) +=
                Mmat(i, j).coef[k] *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, Lmat.BcVec[j].nbc()));
          }
        } else {
          for (int k = 0; k < Mmat(i, j).n + 1; k++) {
            solver_matM.block(0, 0, N + 1, N + 1 + n - Lmat.BcVec[j].nbc()) +=
                Mmat(i, j).coefFun[k].MultMat().block(0, 0, N + 1, N + 1) *
                (Mat.mats2[k + diffn].block(0, Lmat.BcVec[j].nbc(), N + 1,
                                            N + 1 + n - Lmat.BcVec[j].nbc()));
            subs_mat_col_mul2.block(i * (N + 1), subs_mat_col_counter, N + 1,
                                    Lmat.BcVec[j].nbc()) +=
                Mmat(i, j).coefFun[k].MultMat().block(0, 0, N + 1, N + 1) *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, Lmat.BcVec[j].nbc()));
          }
        }
        subs_mat_row_counter += N + 1 + n - Lmat.BcVec[j].nbc();
        masterL.block(i * (N + 1), master_col_counter, N + 1,
                      N + 1 + n - Lmat.BcVec[j].nbc()) = solver_matL;
        masterM.block(i * (N + 1), master_col_counter, N + 1,
                      N + 1 + n - Lmat.BcVec[j].nbc()) = solver_matM;
      }
      subs_mat_col_counter += Lmat.BcVec[j].nbc();
      subs_mat_row_counter = 0;
      master_row_counter = 0;
      master_col_counter += N + 1 + n - Lmat.BcVec[j].nbc();
    }
    if (mat_temp.rows() != 0 || mat_temp.cols() != 0) {
      masterL.block(0, 0, c * (N + 1), c * (N + 1) + Lbc.m) +=
          masterL.block(0, 0, c * (N + 1), c * (N + 1) + Lbc.m).eval() +
          (subs_mat_col_mul * subs_mat);
      masterM.block(0, 0, c * (N + 1), c * (N + 1) + Lbc.m) +=
          masterM.block(0, 0, c * (N + 1), c * (N + 1) + Lbc.m).eval() +
          (subs_mat_col_mul2 * subs_mat);
    }
    // Now append the constraints:

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> constraints_bigL(
        Lbc.m, c * (N + 1) + Lbc.m),
        constraints_smallL(Lbc.m, total_boundary_conditions),
        constraints_bigM(Lbc.m, c * (N + 1) + Lbc.m),
        constraints_smallM(Lbc.m, total_boundary_conditions);
    constraints_bigL.setConstant(0.0);
    constraints_smallL.setConstant(0.0);
    constraints_bigM.setConstant(0.0);
    constraints_smallM.setConstant(0.0);
    col_counter = 0;
    row_counter = 0;
    int col_counter2 = 0;
    int row_counter2 = 0;
    for (int i = 0; i < Lbc.m; i++) {
      for (int j = 0; j < Lbc.n; j++) {
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> temp =
            Lbc(i, j, highest_each_column[j]);

        constraints_bigL.block(i, col_counter, 1,
                               N + 1 + highest_each_column[j] -
                                   num_bc_each_var[j]) =
            temp.block(0, num_bc_each_var[j], 1,
                       N + 1 + highest_each_column[j] - num_bc_each_var[j]);

        constraints_smallL.block(i, col_counter2, 1, num_bc_each_var[j]) =
            temp.block(0, 0, 1, num_bc_each_var[j]);
        //  std::cout << "i = " << i << ", j = " << j << ",
        //  highest_each_column["
        //            << j << "]= " << highest_each_column[j] << '\n';
        temp = Mbc(i, j, highest_each_column[j]);

        constraints_bigM.block(i, col_counter, 1,
                               N + 1 + highest_each_column[j] -
                                   num_bc_each_var[j]) =
            temp.block(0, num_bc_each_var[j], 1,
                       N + 1 + highest_each_column[j] - num_bc_each_var[j]);

        constraints_smallM.block(i, col_counter2, 1, num_bc_each_var[j]) =
            temp.block(0, 0, 1, num_bc_each_var[j]);
        col_counter += N + 1 + highest_each_column[j] - num_bc_each_var[j];
        col_counter2 += num_bc_each_var[j];
      }
      col_counter = 0;
      col_counter2 = 0;
    }
    if (mat_temp.rows() != 0 || mat_temp.cols() != 0) {
      masterL.block(c * (N + 1), 0, Lbc.m, c * (N + 1) + Lbc.m) =
          constraints_bigL + constraints_smallL * subs_mat;
      masterM.block(c * (N + 1), 0, Lbc.m, c * (N + 1) + Lbc.m) =
          constraints_bigM + constraints_smallM * subs_mat;
    } else {
      masterL.block(c * (N + 1), 0, Lbc.m, c * (N + 1) + Lbc.m) =
          constraints_bigL;
      masterM.block(c * (N + 1), 0, Lbc.m, c * (N + 1) + Lbc.m) =
          constraints_bigM;
    }

#ifndef SIS_USE_LAPACK
    Eigen::ColPivHouseholderQR<
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >
        solver(masterM);
    Eigen::ComplexEigenSolver<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >
        eigs;
    bool Is_M_Invertible;
    if (solver.isInvertible()) {
      Is_M_Invertible = true;
      eigs.compute(solver.inverse() * masterL);
      //  std::cout << "Eigenvalues :\n" << eigs.eigenvalues() << '\n';
    } else {
      Is_M_Invertible = false;
      solver.compute(masterL);
      eigs.compute(solver.inverse() * masterM);
      //  std::cout << "Eigenvalues :\n" << 1 / eigs.eigenvalues().array() <<
      //  '\n';
    }
    // std::cin >> bre;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
        evecMat_master(c * (N + 1) + Lbc.m, c * (N + 1) + Lbc.m);
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> temp_vec;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> fir_n_mat,
        full_mat;
    if (subs_mat.size() > 0) {
      fir_n_mat = subs_mat * eigs.eigenvectors();
    }
    full_mat.resize(c * (N + 1) + Lbc.m + fir_n_mat.rows(),
                    c * (N + 1) + Lbc.m);

    row_counter = 0;
    row_counter2 = 0;
    col_counter = 0;
    // Repack eigenvectors using fir_n_mat:
    for (int i = 0; i < c; i++) {
      int n = (highest_each_columnL[i] > highest_each_columnM[i])
                  ? highest_each_columnL[i]
                  : highest_each_columnM[i];
      if (subs_mat.size() > 0) {
        full_mat.block(row_counter, 0, num_bc_each_var[i],
                       c * (N + 1) + Lbc.m) =
            fir_n_mat.block(col_counter, 0, num_bc_each_var[i],
                            c * (N + 1) + Lbc.m);
      }
      row_counter += num_bc_each_var[i];
      col_counter += num_bc_each_var[i];

      full_mat.block(row_counter, 0, N + 1 + n - num_bc_each_var[i],
                     c * (N + 1) + Lbc.m) =
          eigs.eigenvectors().block(row_counter2, 0,
                                    N + 1 + n - num_bc_each_var[i],
                                    c * (N + 1) + Lbc.m);
      row_counter += N + 1 + n - num_bc_each_var[i];
      row_counter2 += N + 1 + n - num_bc_each_var[i];
    }

    row_counter = 0;
    for (int i = 0; i < c; i++) {
      int n = (highest_each_columnL[i] > highest_each_columnM[i])
                  ? highest_each_columnL[i]
                  : highest_each_columnM[i];
      Mat.compute(n);
      evecMat_master.block(i * (N + 1), 0, N + 1, c * (N + 1) + Lbc.m) =
          Mat.mats2[n] *
          full_mat.block(row_counter, 0, N + 1 + n, c * (N + 1) + Lbc.m);
      row_counter += N + 1 + n;
    }

    std::vector<std::vector<T> > eigenval_trunc_sorter;

    eigenval_trunc_sorter.resize(c * (N + 1) + Lbc.m);
    for (int i = 0; i < c * (N + 1) + Lbc.m; i++) {
      eigenval_trunc_sorter[i].resize(2);
    }

    Eigen::Matrix<int, Eigen::Dynamic, 1> MPcount(c * (N + 1) + Lbc.m);
    MPcount.setConstant(0);
    std::vector<EigenSorter<T> > sorter;
    sorter.resize(c * (N + 1) + Lbc.m);
    for (int i = 0; i < c * (N + 1) + Lbc.m; i++) {
      sorter[i].compute(evecMat_master.block(0, i, c * (N + 1), 1), c);
      if (!Is_M_Invertible) {
        // If M is not invertible, there will be zero eigenvalues, so that 1/0 =
        // inf eigenvalue, remove them.
        if (abs(eigs.eigenvalues()[i]) < 1e-11) {
          sorter[i].numMp = N + 1;
          sorter[i].Mp = false;
          sorter[i].minVal = 1e11;
        }
      }
      // sorter[i].compute(evecMat_master.block(0, i, (N + 1), 1));
      if (sorter[i].Mp == true) {
        MPcount[i] = 1;
      }
    }
    int MPcount_sum = MPcount.sum();
    eigenvectorsMat.resize(num_vals);
    for (int i = 0; i < num_vals; i++) {
      eigenvectorsMat[i].resize(c, 1);
    }
    eigenvalues.resize(num_vals);
    if (MPcount_sum >= num_vals) {
      for (int i = 0; i < c * (N + 1) + Lbc.m; i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
        // std::cout << eigs.eigenvalues()[i] << " " << sorter[i].numMp <<
        // "\n"; std::cin > > bre;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < num_vals; i++) {
        if (Is_M_Invertible) {
          eigenvalues[i] = eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
        } else {
          eigenvalues[i] =
              1.0 / eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
        }
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];

        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < num_vals; i++) {
          eigenvectorsMat[i][j] = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }
    } else {
      std::cout << "Last " << num_vals - MPcount_sum
                << " eigenvectors are not resolved to machine precision."
                << '\n';
      for (int i = 0; i < c * (N + 1) + Lbc.m; i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < MPcount_sum; i++) {
        if (Is_M_Invertible) {
          eigenvalues[i] = eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
        } else {
          eigenvalues[i] =
              1.0 / eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
        }
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];
        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < MPcount_sum; i++) {
          eigenvectorsMat[i][j] = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }

      std::vector<std::vector<T> > eigenval_trunc_sorter2;
      eigenval_trunc_sorter2.resize(c * (N + 1) + Lbc.m - MPcount_sum);
      for (int i = MPcount_sum; i < c * (N + 1) + Lbc.m; i++) {
        eigenval_trunc_sorter2[i - MPcount_sum].resize(2);
        eigenval_trunc_sorter2[i - MPcount_sum][0] =
            sorter[int(eigenval_trunc_sorter[i][1])].minVal;
        eigenval_trunc_sorter2[i - MPcount_sum][1] =
            eigenval_trunc_sorter[i][1];
      }
      std::sort(eigenval_trunc_sorter2.begin(), eigenval_trunc_sorter2.end());
      for (int i = MPcount_sum; i < num_vals; i++) {
        if (Is_M_Invertible) {
          eigenvalues[i] = eigs.eigenvalues()[int(
              eigenval_trunc_sorter2[i - MPcount_sum][1])];
        } else {
          eigenvalues[i] =
              1.0 / eigs.eigenvalues()[int(
                        eigenval_trunc_sorter2[i - MPcount_sum][1])];
        }
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];
        Mat.compute(n);
        temp_vec.resize(N + 1 + n, 1);
        for (int i = MPcount_sum; i < num_vals; i++) {
          eigenvectorsMat[i][j] = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter2[i - MPcount_sum][1]),
              N + 1, 1);
          eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }
    }

#else
    char jobvl = 'N';                // Don't compute left evecs
    char jobvr = 'V';                // Compute right evecs
    std::complex<double> wkopt;      // Eistimate optimum workspace
    std::complex<double> *work;      // allocate optimum workspace
    alpha.resize(masterL.rows(), 1); // alpha for gen. eig. prob.
    beta.resize(masterL.rows(), 1);  // beta for gen. eig. prob.

    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> masterL_,
        masterM_, vl(masterL.rows(), masterL.rows()),
        vr(masterL.rows(), masterL.rows()), eigenvalues_temp(masterL.rows(), 1),
        alpha_temp(masterL.rows(), 1), beta_temp(masterL.rows(), 1);
    masterL_ = masterL;
    masterM_ = masterM;
    // vl : left evecs, vr: right evecs.
    int ldL = masterL.outerStride(); // ld for leading dimension
    int ldM = masterM.outerStride();
    int ldvl = vl.outerStride();
    int ldvr = vr.outerStride();
    int sizeL = masterL.rows();
    int lwork = -1; // set lwork to -1 to estimate workspace.
    double rwork[8 * sizeL];

    // call this to estimate workspace
    zggev_(&jobvl, &jobvr, &sizeL, masterL_.data(), &ldL, masterM_.data(), &ldM,
           alpha_temp.data(), beta_temp.data(), vl.data(), &ldvl, vr.data(),
           &ldvr, &wkopt, &lwork, rwork, &info);

    // Now allocate workspace:
    lwork = (int)real(wkopt);
    work = (std::complex<double> *)malloc(lwork * sizeof(std::complex<double>));

    // Solve eigenvalue problem:
    zggev_(&jobvl, &jobvr, &sizeL, masterL_.data(), &ldL, masterM_.data(), &ldM,
           alpha_temp.data(), beta_temp.data(), vl.data(), &ldvl, vr.data(),
           &ldvr, work, &lwork, rwork, &info);

    // Free workspace.
    free((void *)work);

    eigenvalues_temp = alpha_temp.array() / beta_temp.array();

    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
        evecMat_master(c * (N + 1), c * (N + 1) + Lbc.m);
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> temp_vec;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> fir_n_mat =
        subs_mat * vr;
    row_counter = 0;
    col_counter = 0;
    for (int i = 0; i < c; i++) {
      int n = (highest_each_columnL[i] > highest_each_columnM[i])
                  ? highest_each_columnL[i]
                  : highest_each_columnM[i];
      Mat.compute(n);
      // std::cout << "n = " << n << '\n';
      // std::cin > > bre;

      Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> evecMat(
          N + 1 + n, c * (N + 1) + Lbc.m);
      evecMat.setConstant(0.0);
      evecMat.block(num_bc_each_var[i], 0, N + 1 + n - num_bc_each_var[i],
                    c * (N + 1) + Lbc.m) =
          vr.block(row_counter, 0, N + 1 + n - num_bc_each_var[i],
                   c * (N + 1) + Lbc.m);
      row_counter += N + 1 + n - num_bc_each_var[i];
      evecMat.block(0, 0, num_bc_each_var[i], c * (N + 1) + Lbc.m) =
          fir_n_mat.block(col_counter, 0, num_bc_each_var[i],
                          c * (N + 1) + Lbc.m);
      col_counter += num_bc_each_var[i];
      evecMat_master.block(i * (N + 1), 0, N + 1, c * (N + 1) + Lbc.m) =
          Mat.mats2[n] * evecMat;
    }
    std::vector<std::vector<T> > eigenval_trunc_sorter;

    eigenval_trunc_sorter.resize(c * (N + 1) + Lbc.m);
    for (int i = 0; i < c * (N + 1) + Lbc.m; i++) {
      eigenval_trunc_sorter[i].resize(2);
    }

    Eigen::Matrix<int, Eigen::Dynamic, 1> MPcount(c * (N + 1) + Lbc.m);
    MPcount.setConstant(0);
    std::vector<EigenSorter<T> > sorter;
    sorter.resize(c * (N + 1) + Lbc.m);
    for (int i = 0; i < c * (N + 1) + Lbc.m; i++) {
      sorter[i].compute(evecMat_master.block(0, i, c * (N + 1), 1), c);
      // sorter[i].compute(evecMat_master.block(0, i, (N + 1), 1));
      if (sorter[i].Mp == true) {
        MPcount[i] = 1;
      }
    }
    int MPcount_sum = MPcount.sum();

    if (num_vals > masterL.rows()) {
      std::cout << "Only " << masterL.rows()
                << " eigenvalues can be calculated."
                << "Storing only that many." << '\n';
      num_vals = masterL.rows();
    }
    eigenvectorsMat.resize(num_vals);
    for (int i = 0; i < num_vals; i++) {
      eigenvectorsMat[i].resize(c, 1);
    }
    eigenvalues.resize(num_vals);
    if (MPcount_sum >= num_vals) {
      for (int i = 0; i < c * (N + 1); i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
        // std::cout << eigs.eigenvalues()[i] << " " << sorter[i].numMp <<
        // "\n"; std::cin > > bre;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < num_vals; i++) {
        eigenvalues[i] = eigenvalues_temp(int(eigenval_trunc_sorter[i][1]), 0);
        alpha[i] = alpha_temp(int(eigenval_trunc_sorter[i][1]), 0);
        beta[i] = beta_temp(int(eigenval_trunc_sorter[i][1]), 0);
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];

        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < num_vals; i++) {
          eigenvectorsMat[i][j] = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }
    } else {
      std::cout << "Last " << num_vals - MPcount_sum
                << " eigenvectors are not resolved to machine precision."
                << '\n';
      for (int i = 0; i < c * (N + 1); i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < MPcount_sum; i++) {
        eigenvalues[i] = eigenvalues_temp(int(eigenval_trunc_sorter[i][1]), 0);
        alpha[i] = alpha_temp(int(eigenval_trunc_sorter[i][1]), 0);
        beta[i] = beta_temp(int(eigenval_trunc_sorter[i][1]), 0);
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];
        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < MPcount_sum; i++) {
          eigenvectorsMat[i][j] = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }

      std::vector<std::vector<T> > eigenval_trunc_sorter2;
      eigenval_trunc_sorter2.resize(c * (N + 1) + Lbc.m - MPcount_sum);
      for (int i = MPcount_sum; i < c * (N + 1) + Lbc.m; i++) {
        eigenval_trunc_sorter2[i - MPcount_sum].resize(2);
        eigenval_trunc_sorter2[i - MPcount_sum][0] =
            sorter[int(eigenval_trunc_sorter[i][1])].minVal;
        eigenval_trunc_sorter2[i - MPcount_sum][1] =
            eigenval_trunc_sorter[i][1];
      }
      std::sort(eigenval_trunc_sorter2.begin(), eigenval_trunc_sorter2.end());
      for (int i = MPcount_sum; i < num_vals; i++) {
        eigenvalues[i] = eigenvalues_temp(
            int(eigenval_trunc_sorter2[i - MPcount_sum][1]), 0);
        alpha[i] =
            alpha_temp(int(eigenval_trunc_sorter2[i - MPcount_sum][1]), 0);
        beta[i] = beta_temp(int(eigenval_trunc_sorter2[i - MPcount_sum][1]), 0);
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];
        Mat.compute(n);
        temp_vec.resize(N + 1 + n, 1);
        for (int i = MPcount_sum; i < num_vals; i++) {
          eigenvectorsMat[i][j] = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter2[i - MPcount_sum][1]),
              N + 1, 1);
          eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }
    }

#endif
  }

  /// This is the most generic form of the eigenvalue solver, in which boundary
  /// conditions that don't have an associated eigenvalue are to be embedded (
  /// using Lbc_embed_), and boundary conditions that have an associated
  /// eigevalue associated are to be appended.
  /// This means that you can go about solving problems where you would
  /// append boundary conditions anyway, instead of embedding them, see the
  /// example in our paper, and example/Ex_19.cpp.
  void compute_with_constraints(const LinopMat<T> &Lmat_,
                                const LinopMat<T> &Mmat_, int num_vals,
                                const BcMat<T> &Lbc_embed_,
                                const BcMat<T> &Lbc_append_,
                                const BcMat<T> &Mbc_append_) {
    LinopMat<T> Lmat = Lmat_;
    LinopMat<T> Mmat = Mmat_;
    BcMat<T> Lbc_embed = Lbc_embed_;
    BcMat<T> Lbc_append = Lbc_append_;
    BcMat<T> Mbc_append = Mbc_append_;
    int bre;
    int total_of_all_orders = 0;
    int total_boundary_conditions = 0;
    int total_constraints = 0;
    if (Lmat.r != Lmat.c) {
      std::cout << "Solution only possible with square LinopMats. Exiting ..."
                << '\n';
      exit(1);
    }
    if (Mmat.r != Mmat.c) {
      std::cout << "Solution only possible with square LinopMats. Exiting ..."
                << '\n';
      exit(1);
    }
    if (Mmat.c != Lmat.c) {
      std::cout << "Both matrices have to be of same size. Exiting ..." << '\n';
      exit(1);
    }

    if (Lmat.c != Lbc_embed.n) {
      std::cout
          << "Number of columns in Lbc_embed have to be equal to number of "
          << "columns in Lmat . Exiting .... In line " << __LINE__ << '\n';
      exit(1);
    }
    if (Lmat.c != Lbc_append.n) {
      std::cout
          << "Number of columns in Lbc_append have to be equal to number of "
          << "columns in Lmat . Exiting .... In line " << __LINE__ << '\n';
      exit(1);
    }
    if (Lmat.c != Mbc_append.n) {
      std::cout
          << "Number of columns in Mbc_append have to be equal to number of "
          << "columns in Lmat . Exiting .... In line " << __LINE__ << '\n';
      exit(1);
    }
    if (Lbc_append.m != Mbc_append.m) {
      std::cout << "Number of rows in Mbc_append have to be equal to number of "
                << "rows in Lbc_append. Exiting .... In line " << __LINE__
                << '\n';
      exit(1);
    }

    int r = Lmat.r, c = Lmat.c;
    // Find the highest derivative in each column. To do this create a vector
    // highest_each_column, and a temp_int_vec that will hold all values of a
    // given column to and the maximum will be stored in highest_each_column
    std::vector<int> highest_each_columnL, highest_each_columnM,
        highest_each_column;
    highest_each_columnL.resize(c);
    highest_each_columnM.resize(c);
    highest_each_column.resize(c);

    std::vector<int> temp_vec_int;
    temp_vec_int.resize(r);
    for (int j = 0; j < c; j++) {
      for (int i = 0; i < r; i++) {
        temp_vec_int[i] = Lmat(i, j).n;
      }
      highest_each_columnL[j] =
          *std::max_element(temp_vec_int.begin(), temp_vec_int.end());
    }
    for (int j = 0; j < c; j++) {
      for (int i = 0; i < r; i++) {
        temp_vec_int[i] = Mmat(i, j).n;
      }
      highest_each_columnM[j] =
          *std::max_element(temp_vec_int.begin(), temp_vec_int.end());
    }
    for (int i = 0; i < c; i++) {
      total_of_all_orders += (highest_each_columnL[i] > highest_each_columnM[i])
                                 ? highest_each_columnL[i]
                                 : highest_each_columnM[i];
      highest_each_column[i] =
          (highest_each_columnL[i] > highest_each_columnM[i])
              ? highest_each_columnL[i]
              : highest_each_columnM[i];
    }
    total_boundary_conditions = Lbc_embed.m + Lbc_append.m;

    // total_of_all_orders has to be equal to total number of boundary
    // conditions, else the problem is ill-posed, if ill-posed, cout the same
    // and exit.
    if (total_of_all_orders != total_boundary_conditions) {
      std::cout << "The problem is ill-posed, the total of the highest "
                   "orders of all "
                   "dependent variables has to be equal to the total number of "
                   "boundary conditions specified."
                   "\n "
                   "Total no. of boundary conditions: "
                << total_boundary_conditions
                << "\n"
                   "Total of all orders: "
                << total_of_all_orders << "\n Exiting ...\n";
      exit(1);
    }
    // Declare the master matrix L:
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> masterL(
        r * (N + 1), c * (N + 1) + total_boundary_conditions);
    // Declare the master matrix M:
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> masterM(
        r * (N + 1), c * (N + 1) + total_boundary_conditions);

    masterL.setConstant(0.0);
    masterM.setConstant(0.0);
    MatGen<T> Mat;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> subs_mat;
    subs_mat.resize(total_boundary_conditions, (N + 1) * c);
    subs_mat.setConstant(0.0);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat_temp(
        total_boundary_conditions, total_boundary_conditions),
        constraints_embed(Lbc_embed.m, c * (N + 1) + total_of_all_orders);

    mat_temp.setConstant(0.0);
    int row_counter = 0, col_counter = 0;

    for (int i = 0; i < Lbc_embed.m; i++) {
      for (int j = 0; j < Lbc_embed.n; j++) {
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> temp =
            Lbc_embed(i, j, highest_each_column[j]);
        constraints_embed.block(i, col_counter, 1, temp.cols()) = temp;
        col_counter += temp.cols();
      }
      col_counter = 0;
    }

    Eigen::ColPivHouseholderQR<
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >
        qr;
    qr.compute(constraints_embed);
    if (qr.rank() != constraints_embed.rows()) {
      std::cout << "The boundary conditions supplied are not "
                << "  linearly independent." << '\n';
      std::cout << "qr.rank = " << qr.rank()
                << ", no. bcs: " << total_boundary_conditions << ". Exiting ..."
                << '\n';
      exit(1);
    }

    // Permutation matrix
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P;
    P = qr.colsPermutation();
    // Permute constraints
    constraints_embed = constraints_embed * P;
    mat_temp = constraints_embed.block(0, 0, constraints_embed.rows(),
                                       constraints_embed.rows());

    Eigen::ColPivHouseholderQR<
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >
        consSolver(mat_temp);
    subs_mat = -consSolver.inverse() *
               constraints_embed.block(
                   0, constraints_embed.rows(), constraints_embed.rows(),
                   constraints_embed.cols() - constraints_embed.rows());

    row_counter = 0;
    col_counter = 0;
    int master_row_counter = 0;
    int master_col_counter = 0;
    for (int j = 0; j < c; j++) {
      int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                  ? highest_each_columnL[j]
                  : highest_each_columnM[j];
      Mat.compute(n);
      for (int i = 0; i < r; i++) {
        int diffn = n - Lmat(i, j).n;
        if (Lmat(i, j).NCC == 0) {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            masterL.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Lmat(i, j).coef[k] *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, N + 1 + n));
          }
        } else {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            masterL.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Lmat(i, j).coefFun[k].MultMat().block(0, 0, N + 1, N + 1) *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, N + 1 + n));
          }
        }
        diffn = n - Mmat(i, j).n;
        if (Mmat(i, j).NCC == 0) {
          for (int k = 0; k < Mmat(i, j).n + 1; k++) {
            masterM.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Mmat(i, j).coef[k] *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, N + 1 + n));
          }
        } else {
          for (int k = 0; k < Mmat(i, j).n + 1; k++) {
            masterM.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Mmat(i, j).coefFun[k].MultMat().block(0, 0, N + 1, N + 1) *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, N + 1 + n));
          }
        }
        master_row_counter += N + 1;
      }
      master_row_counter = 0;
      master_col_counter += N + 1 + n;
    }
    // Permute columns of M and L:
    masterL = masterL * P;
    masterM = masterM * P;

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> masterL2(
        r * (N + 1) + Lbc_append.m, c * (N + 1) + Lbc_append.m),
        masterM2(r * (N + 1) + Lbc_append.m, c * (N + 1) + Lbc_append.m);

    masterL2.block(0, 0, r * (N + 1), r * (N + 1) + Lbc_append.m) =
        masterL.block(0, constraints_embed.rows(), c * (N + 1),
                      c * (N + 1) + Lbc_append.m) +
        (masterL.block(0, 0, c * (N + 1), constraints_embed.rows()) * subs_mat);
    masterM2.block(0, 0, r * (N + 1), r * (N + 1) + Lbc_append.m) =
        masterM.block(0, constraints_embed.rows(), c * (N + 1),
                      c * (N + 1) + Lbc_append.m) +
        (masterM.block(0, 0, c * (N + 1), constraints_embed.rows()) * subs_mat);

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> constraints_appendL(
        Lbc_append.m, c * (N + 1) + total_of_all_orders),
        constraints_appendM(Lbc_append.m, c * (N + 1) + total_of_all_orders);
    col_counter = 0;
    for (int i = 0; i < Lbc_append.m; i++) {
      for (int j = 0; j < Lbc_append.n; j++) {
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> temp =
            Lbc_append(i, j, highest_each_column[j]);
        constraints_appendL.block(i, col_counter, 1, temp.cols()) = temp;
        col_counter += temp.cols();
      }
      col_counter = 0;
    }
    col_counter = 0;
    for (int i = 0; i < Lbc_append.m; i++) {
      for (int j = 0; j < Lbc_append.n; j++) {
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> temp =
            Mbc_append(i, j, highest_each_column[j]);
        constraints_appendM.block(i, col_counter, 1, temp.cols()) = temp;
        col_counter += temp.cols();
      }
      col_counter = 0;
    }
    constraints_appendM = constraints_appendM * P;
    constraints_appendL = constraints_appendL * P;

    masterL2.block(r * (N + 1), 0, Lbc_append.m, r * (N + 1) + Lbc_append.m) =
        constraints_appendL.block(0, Lbc_embed.m, Lbc_append.m,
                                  r * (N + 1) + Lbc_append.m) +
        (constraints_appendL.block(0, 0, Lbc_append.m, Lbc_embed.m) * subs_mat);

    masterM2.block(r * (N + 1), 0, Lbc_append.m, r * (N + 1) + Lbc_append.m) =
        constraints_appendM.block(0, Lbc_embed.m, Lbc_append.m,
                                  r * (N + 1) + Lbc_append.m) +
        (constraints_appendM.block(0, 0, Lbc_append.m, Lbc_embed.m) * subs_mat);

#ifndef SIS_USE_LAPACK
    Eigen::ColPivHouseholderQR<
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >
        solver(masterM2);
    Eigen::ComplexEigenSolver<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >
        eigs;
    bool Is_M_Invertible;
    if (solver.isInvertible()) {
      //  std::cout << "M is invertible" << '\n';
      //  std::cin >> bre;
      Is_M_Invertible = true;
      eigs.compute(solver.inverse() * masterL2);
      // std::cout << "Eigenvalues :\n" << eigs.eigenvalues() << '\n';
    } else {
      // std::cout << "M is not invertible." << '\n';
      // std::cin >> bre;
      Is_M_Invertible = false;
      solver.compute(masterL2);
      eigs.compute(solver.inverse() * masterM2);
      // std::cout << "Eigenvalues :\n" << 1 / eigs.eigenvalues().array() <<
      // '\n';
    }
    // std::cin >> bre;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> evecMat_master(
        c * (N + 1), c * (N + 1) + Lbc_append.m);
    Eigen::Matrix<T, Eigen::Dynamic, 1> temp_vec;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> fir_n_mat =
        subs_mat * eigs.eigenvectors();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> evecMat(
        c * (N + 1) + total_boundary_conditions, c * (N + 1) + Lbc_append.m);
    evecMat << fir_n_mat, eigs.eigenvectors();

    // Unpermute eveMat:
    evecMat = P * evecMat;
    row_counter = 0;
    col_counter = 0;
    for (int i = 0; i < c; i++) {
      int n = (highest_each_columnL[i] > highest_each_columnM[i])
                  ? highest_each_columnL[i]
                  : highest_each_columnM[i];
      Mat.compute(n);
      evecMat_master.block(i * (N + 1), 0, N + 1, c * (N + 1) + Lbc_append.m) =
          Mat.mats2[n] *
          evecMat.block(row_counter, 0, N + 1 + n, c * (N + 1) + Lbc_append.m);
      row_counter += N + 1 + n;
    }
    std::vector<std::vector<T> > eigenval_trunc_sorter;

    eigenval_trunc_sorter.resize(c * (N + 1) + Lbc_append.m);
    for (int i = 0; i < c * (N + 1) + Lbc_append.m; i++) {
      eigenval_trunc_sorter[i].resize(2);
    }

    Eigen::Matrix<int, Eigen::Dynamic, 1> MPcount(c * (N + 1) + Lbc_append.m);
    MPcount.setConstant(0);
    std::vector<EigenSorter<T> > sorter;
    sorter.resize(c * (N + 1) + Lbc_append.m);
    for (int i = 0; i < c * (N + 1) + Lbc_append.m; i++) {
      sorter[i].compute(evecMat_master.block(0, i, c * (N + 1), 1), c);
      if (!Is_M_Invertible) {
        // If M is not invertible, there will be zero eigenvalues, so that 1/0 =
        // inf eigenvalue, remove them.
        if (abs(eigs.eigenvalues()[i]) < 1e-11) {
          sorter[i].numMp = N + 1;
          sorter[i].Mp = false;
          sorter[i].minVal = 1e11;
        }
      }
      // sorter[i].compute(evecMat_master.block(0, i, (N + 1), 1));
      if (sorter[i].Mp == true) {
        MPcount[i] = 1;
      }
    }
    if (num_vals > masterL.rows()) {
      std::cout << "Only " << masterL.rows()
                << " eigenvalues can be calculated."
                << "Storing only that many." << '\n';
      num_vals = masterL.rows();
    }
    int MPcount_sum = MPcount.sum();
    eigenvectorsMat.resize(num_vals);
    for (int i = 0; i < num_vals; i++) {
      eigenvectorsMat[i].resize(c, 1);
    }
    eigenvalues.resize(num_vals);
    if (MPcount_sum >= num_vals) {
      for (int i = 0; i < c * (N + 1) + Lbc_append.m; i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
        // std::cout << eigs.eigenvalues()[i] << " " << sorter[i].numMp <<
        // "\n"; std::cin > > bre;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < num_vals; i++) {
        if (Is_M_Invertible) {
          eigenvalues[i] = eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
        } else {
          eigenvalues[i] =
              1.0 / eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
        }
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];

        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < num_vals; i++) {
          eigenvectorsMat[i][j] = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
        }
      }
    } else {
      std::cout << "Last " << num_vals - MPcount_sum
                << " eigenvectors are not resolved to machine precision."
                << '\n';
      for (int i = 0; i < c * (N + 1) + Lbc_append.m; i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < MPcount_sum; i++) {
        if (Is_M_Invertible) {
          eigenvalues[i] = eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
        } else {
          eigenvalues[i] =
              1.0 / eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
        }
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];
        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < MPcount_sum; i++) {
          eigenvectorsMat[i][j] = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
        }
      }

      std::vector<std::vector<T> > eigenval_trunc_sorter2;
      eigenval_trunc_sorter2.resize(c * (N + 1) + Lbc_append.m - MPcount_sum);
      for (int i = MPcount_sum; i < c * (N + 1) + Lbc_append.m; i++) {
        eigenval_trunc_sorter2[i - MPcount_sum].resize(2);
        eigenval_trunc_sorter2[i - MPcount_sum][0] =
            sorter[int(eigenval_trunc_sorter[i][1])].minVal;
        eigenval_trunc_sorter2[i - MPcount_sum][1] =
            eigenval_trunc_sorter[i][1];
      }
      std::sort(eigenval_trunc_sorter2.begin(), eigenval_trunc_sorter2.end());
      for (int i = MPcount_sum; i < num_vals; i++) {
        if (Is_M_Invertible) {
          eigenvalues[i] = eigs.eigenvalues()[int(
              eigenval_trunc_sorter2[i - MPcount_sum][1])];
        } else {
          eigenvalues[i] =
              1.0 / eigs.eigenvalues()[int(
                        eigenval_trunc_sorter2[i - MPcount_sum][1])];
        }
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];
        Mat.compute(n);
        temp_vec.resize(N + 1 + n, 1);
        for (int i = MPcount_sum; i < num_vals; i++) {
          eigenvectorsMat[i][j] = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter2[i - MPcount_sum][1]),
              N + 1, 1);
        }
      }
    }
#else
    std::cout << "Using lapack routine..." << '\n';
    char jobvl = 'N';                 // Don't compute left evecs
    char jobvr = 'V';                 // Compute right evecs
    std::complex<double> wkopt;       // Eistimate optimum workspace
    std::complex<double> *work;       // allocate optimum workspace
    alpha.resize(masterL2.rows(), 1); // alpha for gen. eig. prob.
    beta.resize(masterL2.rows(), 1);  // beta for gen. eig. prob.

    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> masterL2_,
        masterM2_, vl(masterL2.rows(), masterL2.rows()),
        vr(masterL2.rows(), masterL2.rows()),
        eigenvalues_temp(masterL2.rows(), 1), alpha_temp(masterL2.rows(), 1),
        beta_temp(masterL2.rows(), 1);
    masterL2_ = masterL2;
    masterM2_ = masterM2;
    // vl : left evecs, vr: right evecs.
    int ldL = masterL2.outerStride(); // ld for leading dimension
    int ldM = masterM2.outerStride();
    int ldvl = vl.outerStride();
    int ldvr = vr.outerStride();
    int sizeL = masterL2.rows();
    int lwork = -1; // set lwork to -1 to estimate workspace.
    double rwork[8 * sizeL];

    // call this to estimate workspace
    zggev_(&jobvl, &jobvr, &sizeL, masterL2_.data(), &ldL, masterM2_.data(),
           &ldM, alpha_temp.data(), beta_temp.data(), vl.data(), &ldvl,
           vr.data(), &ldvr, &wkopt, &lwork, rwork, &info);

    // Now allocate workspace:
    lwork = (int)real(wkopt);
    work = (std::complex<double> *)malloc(lwork * sizeof(std::complex<double>));

    // Solve eigenvalue problem:
    zggev_(&jobvl, &jobvr, &sizeL, masterL2_.data(), &ldL, masterM2_.data(),
           &ldM, alpha_temp.data(), beta_temp.data(), vl.data(), &ldvl,
           vr.data(), &ldvr, work, &lwork, rwork, &info);

    // Free workspace.
    free((void *)work);

    eigenvalues_temp = alpha_temp.array() / beta_temp.array();

    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
        evecMat_master(c * (N + 1), c * (N + 1) + Lbc_append.m);
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> temp_vec;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> fir_n_mat =
        subs_mat * vr;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> evecMat(
        c * (N + 1) + total_boundary_conditions, c * (N + 1) + Lbc_append.m);
    evecMat << fir_n_mat, vr;
    // Unpermute eveMat:
    evecMat = P * evecMat;
    row_counter = 0;
    col_counter = 0;
    for (int i = 0; i < c; i++) {
      int n = (highest_each_columnL[i] > highest_each_columnM[i])
                  ? highest_each_columnL[i]
                  : highest_each_columnM[i];
      Mat.compute(n);
      evecMat_master.block(i * (N + 1), 0, N + 1, c * (N + 1) + Lbc_append.m) =
          Mat.mats2[n] *
          evecMat.block(row_counter, 0, N + 1 + n, c * (N + 1) + Lbc_append.m);
      row_counter += N + 1 + n;
    }
    std::vector<std::vector<T> > eigenval_trunc_sorter;

    eigenval_trunc_sorter.resize(c * (N + 1) + Lbc_append.m);
    for (int i = 0; i < c * (N + 1) + Lbc_append.m; i++) {
      eigenval_trunc_sorter[i].resize(2);
    }

    Eigen::Matrix<int, Eigen::Dynamic, 1> MPcount(c * (N + 1) + Lbc_append.m);
    MPcount.setConstant(0);
    std::vector<EigenSorter<T> > sorter;
    sorter.resize(c * (N + 1) + Lbc_append.m);
    for (int i = 0; i < c * (N + 1) + Lbc_append.m; i++) {
      sorter[i].compute(evecMat_master.block(0, i, c * (N + 1), 1), c);
      // sorter[i].compute(evecMat_master.block(0, i, (N + 1), 1));
      if (sorter[i].Mp == true) {
        MPcount[i] = 1;
      }
    }
    if (num_vals > masterL2.rows()) {
      std::cout << "Only " << masterL2.rows()
                << " eigenvalues can be calculated."
                << "Storing only that many." << '\n';
      num_vals = masterL2.rows();
    }

    int MPcount_sum = MPcount.sum();
    eigenvectorsMat.resize(num_vals);
    for (int i = 0; i < num_vals; i++) {
      eigenvectorsMat[i].resize(c, 1);
    }
    eigenvalues.resize(num_vals);
    beta.resize(num_vals);
    alpha.resize(num_vals);
    if (MPcount_sum >= num_vals) {
      for (int i = 0; i < c * (N + 1) + Lbc_append.m; i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
        // std::cout << eigs.eigenvalues()[i] << " " << sorter[i].numMp <<
        // "\n"; std::cin > > bre;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < num_vals; i++) {
        eigenvalues[i] = eigenvalues_temp(int(eigenval_trunc_sorter[i][1]), 0);
        alpha[i] = alpha_temp(int(eigenval_trunc_sorter[i][1]), 0);
        beta[i] = beta_temp(int(eigenval_trunc_sorter[i][1]), 0);
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];

        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < num_vals; i++) {
          eigenvectorsMat[i][j] = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }
    } else {
      std::cout << "Last " << num_vals - MPcount_sum
                << " eigenvectors are not resolved to machine precision."
                << '\n';
      for (int i = 0; i < c * (N + 1) + Lbc_append.m; i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < MPcount_sum; i++) {
        eigenvalues[i] = eigenvalues_temp(int(eigenval_trunc_sorter[i][1]), 0);
        alpha[i] = alpha_temp(int(eigenval_trunc_sorter[i][1]), 0);
        beta[i] = beta_temp(int(eigenval_trunc_sorter[i][1]), 0);
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];
        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < MPcount_sum; i++) {
          eigenvectorsMat[i][j] = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }

      std::vector<std::vector<T> > eigenval_trunc_sorter2;
      eigenval_trunc_sorter2.resize(c * (N + 1) + Lbc_append.m - MPcount_sum);
      for (int i = MPcount_sum; i < c * (N + 1) + Lbc_append.m; i++) {
        eigenval_trunc_sorter2[i - MPcount_sum].resize(2);
        eigenval_trunc_sorter2[i - MPcount_sum][0] =
            sorter[int(eigenval_trunc_sorter[i][1])].minVal;
        eigenval_trunc_sorter2[i - MPcount_sum][1] =
            eigenval_trunc_sorter[i][1];
      }
      std::sort(eigenval_trunc_sorter2.begin(), eigenval_trunc_sorter2.end());
      for (int i = MPcount_sum; i < num_vals; i++) {
        eigenvalues[i] = eigenvalues_temp(
            int(eigenval_trunc_sorter2[i - MPcount_sum][1]), 0);
        alpha[i] =
            alpha_temp(int(eigenval_trunc_sorter2[i - MPcount_sum][1]), 0);
        beta[i] = beta_temp(int(eigenval_trunc_sorter2[i - MPcount_sum][1]), 0);
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];
        Mat.compute(n);
        temp_vec.resize(N + 1 + n, 1);
        for (int i = MPcount_sum; i < num_vals; i++) {
          eigenvectorsMat[i][j] = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter2[i - MPcount_sum][1]),
              N + 1, 1);
          eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }
    }

#endif
  }
};

template <class T> class GeneralizedEigenSolver<std::complex<T> > {
private:
public:
  //std::valarray<ChebfunMat<std::complex<T> > > eigenvectorsMat;
  ChebfunMat<std::complex<T> > eigenvectors;
  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> eigenvalues;
  Eigen::Matrix<int, Eigen::Dynamic, 1> MPorNot;
  /// \brief Number of eigenvalues that have converged to machine precision.
  int converged;
  /// \brief Null constructor
  GeneralizedEigenSolver() : converged(0){};

#ifdef SIS_USE_LAPACK
  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> alpha, beta;
  int info; // gives info on lapack solver.
#endif

/// \brief Call this with an input Linear operator to solve for eigenvalues
/// and vectors. The number of Eigen values/vectors is num_vals, num_vals
/// has to be less than N*r.
  void compute(const LinopMat<std::complex<T> > &Lmat_,
               int num_vals, const BcMat<std::complex<T> > &bc_) {
    LinopMat<std::complex<T> > Mmat(Lmat_.r, Lmat_.c);
    Mmat.setIdentity();
    compute(Lmat_, Mmat, num_vals, bc_);
  }

  /// \brief Call this with an input Linear operator to solve for eigenvalues
  /// and vectors.
    void compute(const LinopMat<std::complex<T> > &Lmat_,
                 const BcMat<std::complex<T> > &bc_) {
      LinopMat<std::complex<T> > Mmat(Lmat_.r, Lmat_.c);
      int num_vals = Lmat_.r * (N + 1);
      Mmat.setIdentity();
      compute(Lmat_, Mmat, num_vals, bc_);
    }
  /// \brief The main solver for LinopMat. Read about class BcMat to see how
  /// this works, also see examples
  /// example/Ex_16.cpp and test/Visco_3D_pipe.cpp of how this is applied.
  /// This function is useful when boundary conditions are mixed between
  /// variables. Also, if you have boundary conditions that have an associated
  /// eigenvalue, see compute_with_constraints.
  void compute(const LinopMat<std::complex<T> > &Lmat_,
               const LinopMat<std::complex<T> > &Mmat_, int num_vals,
               const BcMat<std::complex<T> > &Lbc_) {
    LinopMat<std::complex<T> > Lmat = Lmat_;
    LinopMat<std::complex<T> > Mmat = Mmat_;
    BcMat<std::complex<T> > Lbc = Lbc_;
    int bre;
    int total_of_all_orders = 0;
    int total_boundary_conditions = 0;
    if (Lmat.r != Lmat.c) {
      std::cout << "Solution only possible with square LinopMats. Exiting ..."
                << '\n';
      exit(1);
    }
    if (Mmat.r != Mmat.c) {
      std::cout << "Solution only possible with square LinopMats. Exiting ..."
                << '\n';
      exit(1);
    }
    if (Mmat.c != Lmat.c) {
      std::cout << "Both matrices have to be of same size. Exiting ..." << '\n';
      exit(1);
    }
    int r = Lmat.r, c = Lmat.c;
    // Find the highest derivative in each column. To do this create a vector
    // highest_each_column, and a temp_int_vec that will hold all values of a
    // given column to and the maximum will be stored in highest_each_column
    std::vector<int> highest_each_columnL, highest_each_columnM,
        highest_each_column;
    highest_each_columnL.resize(c);
    highest_each_columnM.resize(c);
    highest_each_column.resize(c);

    std::vector<int> temp_vec_int;
    temp_vec_int.resize(r);
    for (int j = 0; j < c; j++) {
      for (int i = 0; i < r; i++) {
        temp_vec_int[i] = Lmat(i, j).n;
      }
      highest_each_columnL[j] =
          *std::max_element(temp_vec_int.begin(), temp_vec_int.end());
    }
    for (int j = 0; j < c; j++) {
      for (int i = 0; i < r; i++) {
        temp_vec_int[i] = Mmat(i, j).n;
      }
      highest_each_columnM[j] =
          *std::max_element(temp_vec_int.begin(), temp_vec_int.end());
    }
    for (int i = 0; i < c; i++) {
      total_of_all_orders += (highest_each_columnL[i] > highest_each_columnM[i])
                                 ? highest_each_columnL[i]
                                 : highest_each_columnM[i];
      highest_each_column[i] =
          (highest_each_columnL[i] > highest_each_columnM[i])
              ? highest_each_columnL[i]
              : highest_each_columnM[i];
      //  std::cout << "highest_each_column["<< i << "]: " <<
      //  highest_each_column[i] << '\n';
    }
    total_boundary_conditions = Lbc.m;

    // total_of_all_orders has to be equal to total number of boundary
    // conditions, else the problem is ill-posed, if ill-posed, cout the same
    // and exit.
    if (total_of_all_orders != total_boundary_conditions) {
      std::cout << "The problem is ill-posed, the total of the highest "
                   "orders of all "
                   "dependent variables has to be equal to the total number of "
                   "boundary conditions specified."
                   "\n "
                   "Total no. of boundary conditions: "
                << total_boundary_conditions
                << "\n"
                   "Total of all orders: "
                << total_of_all_orders << "\n Exiting ...\n";
      //  exit(1);
    }
    // Declare the master matrix L:
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> masterL(
        r * (N + 1), c * (N + 1) + total_boundary_conditions);
    // Declare the master matrix M:
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> masterM(
        r * (N + 1), c * (N + 1) + total_boundary_conditions);

    masterL.setConstant(0.0);
    masterM.setConstant(0.0);
    MatGen<T> Mat;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> subs_mat;
    subs_mat.resize(total_boundary_conditions, (N + 1) * c);
    subs_mat.setConstant(0.0);
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> mat_temp(
        total_boundary_conditions, total_boundary_conditions),
        constraints(total_boundary_conditions,
                    c * (N + 1) + total_of_all_orders);

    mat_temp.setConstant(0.0);
    int row_counter = 0, col_counter = 0;

    for (int i = 0; i < Lbc.m; i++) {
      for (int j = 0; j < Lbc.n; j++) {
        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> temp =
            Lbc(i, j, highest_each_column[j]);
        constraints.block(i, col_counter, 1, temp.cols()) = temp;
        col_counter += temp.cols();
      }
      col_counter = 0;
    }

    Eigen::ColPivHouseholderQR<
        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
        qr;
    qr.compute(constraints);
    if (qr.rank() != constraints.rows()) {
      std::cout << "The boundary conditions supplied are not "
                << "  linearly independent." << '\n';
      std::cout << "qr.rank = " << qr.rank()
                << ", no. bcs: " << total_boundary_conditions << ". Exiting ..."
                << '\n';
      exit(1);
    }

    // Permutation matrix
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P;
    P = qr.colsPermutation();
    // Permute constraints
    constraints = constraints * P;
    mat_temp = constraints.block(0, 0, constraints.rows(), constraints.rows());

    Eigen::ColPivHouseholderQR<
        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
        consSolver(mat_temp);
    subs_mat = -consSolver.inverse() *
               constraints.block(0, constraints.rows(), constraints.rows(),
                                 constraints.cols() - constraints.rows());

    row_counter = 0;
    col_counter = 0;
    int master_row_counter = 0;
    int master_col_counter = 0;
    for (int j = 0; j < c; j++) {
      int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                  ? highest_each_columnL[j]
                  : highest_each_columnM[j];
      Mat.compute(n);
      for (int i = 0; i < r; i++) {
        int diffn = n - Lmat(i, j).n;
        if (Lmat(i, j).NCC == 0) {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            masterL.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Lmat(i, j).coef[k] *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, N + 1 + n));
          }
        } else {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            masterL.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Lmat(i, j).coefFun[k].MultMat().block(0, 0, N + 1, N + 1) *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, N + 1 + n));
          }
        }
        diffn = n - Mmat(i, j).n;
        if (Mmat(i, j).NCC == 0) {
          for (int k = 0; k < Mmat(i, j).n + 1; k++) {
            masterM.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Mmat(i, j).coef[k] *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, N + 1 + n));
          }
        } else {
          for (int k = 0; k < Mmat(i, j).n + 1; k++) {
            masterM.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Mmat(i, j).coefFun[k].MultMat().block(0, 0, N + 1, N + 1) *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, N + 1 + n));
          }
        }
        master_row_counter += N + 1;
      }
      master_row_counter = 0;
      master_col_counter += N + 1 + n;
    }
    // Permute columns of M and L:
    masterL = masterL * P;
    masterM = masterM * P;

    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> masterL2,
        masterM2;

    masterL2 =
        masterL.block(0, constraints.rows(), c * (N + 1), c * (N + 1)) +
        (masterL.block(0, 0, c * (N + 1), constraints.rows()) * subs_mat);

    masterM2 =
        masterM.block(0, constraints.rows(), c * (N + 1), c * (N + 1)) +
        (masterM.block(0, 0, c * (N + 1), constraints.rows()) * subs_mat);

#if defined SIS_USE_LAPACK
    char jobvl = 'N';                 // Don't compute left evecs
    char jobvr = 'V';                 // Compute right evecs
    std::complex<double> wkopt;       // Eistimate optimum workspace
    std::complex<double> *work;       // allocate optimum workspace
    alpha.resize(masterL2.rows(), 1); // alpha for gen. eig. prob.
    beta.resize(masterL2.rows(), 1);  // beta for gen. eig. prob.

    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> vl(
        masterL2.rows(), masterL2.rows()),
        vr(masterL2.rows(), masterL2.rows()),
        eigenvalues_temp(masterL2.rows(), 1), alpha_temp(masterL2.rows(), 1),
        beta_temp(masterL2.rows(), 1);
    // vl : left evecs, vr: right evecs.
    int ldL = masterL2.outerStride(); // ld for leading dimension
    int ldM = masterM2.outerStride();
    int ldvl = vl.outerStride();
    int ldvr = vr.outerStride();
    int sizeL = masterL2.rows();
    int lwork = -1; // set lwork to -1 to estimate workspace.
    double rwork[8 * sizeL];

    // call this to estimate workspace
    zggev_(&jobvl, &jobvr, &sizeL, masterL2.data(), &ldL, masterM2.data(), &ldM,
           alpha_temp.data(), beta_temp.data(), vl.data(), &ldvl, vr.data(),
           &ldvr, &wkopt, &lwork, rwork, &info);

    // Now allocate workspace:
    lwork = (int)real(wkopt);
    work = (std::complex<double> *)malloc(lwork * sizeof(std::complex<double>));

    // Solve eigenvalue problem:
    zggev_(&jobvl, &jobvr, &sizeL, masterL2.data(), &ldL, masterM2.data(), &ldM,
           alpha_temp.data(), beta_temp.data(), vl.data(), &ldvl, vr.data(),
           &ldvr, work, &lwork, rwork, &info);

    // Free workspace.
    free((void *)work);

    eigenvalues_temp = alpha_temp.array() / beta_temp.array();

    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
        evecMat_master(c * (N + 1), c * (N + 1));
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> temp_vec;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> fir_n_mat =
        subs_mat * vr;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> evecMat(
        c * (N + 1) + total_boundary_conditions, c * (N + 1));
    evecMat << fir_n_mat, vr;
    // Unpermute eveMat:
    evecMat = P * evecMat;
    row_counter = 0;
    col_counter = 0;
    for (int i = 0; i < c; i++) {
      int n = (highest_each_columnL[i] > highest_each_columnM[i])
                  ? highest_each_columnL[i]
                  : highest_each_columnM[i];
      Mat.compute(n);
      evecMat_master.block(i * (N + 1), 0, N + 1, c * (N + 1)) =
          Mat.mats2[n] * evecMat.block(row_counter, 0, N + 1 + n, c * (N + 1));
      row_counter += N + 1 + n;
    }
#ifndef SIS_DONT_SORT
    std::vector<std::vector<T> > eigenval_trunc_sorter;

    eigenval_trunc_sorter.resize(c * (N + 1));
    for (int i = 0; i < c * (N + 1); i++) {
      eigenval_trunc_sorter[i].resize(2);
    }

    Eigen::Matrix<int, Eigen::Dynamic, 1> MPcount(c * (N + 1));
    MPcount.setConstant(0);
    std::vector<EigenSorter<T> > sorter;
    sorter.resize(c * (N + 1));
    for (int i = 0; i < c * (N + 1); i++) {
      sorter[i].compute(evecMat_master.block(0, i, c * (N + 1), 1), c);
      // sorter[i].compute(evecMat_master.block(0, i, (N + 1), 1));
      if (std::abs(beta_temp(i, 0)) < 1e-11) {
        sorter[i].numMp = N + 1;
        sorter[i].Mp = false;
        sorter[i].minVal = 1e11;
      }
      if (sorter[i].Mp == true) {
        MPcount[i] = 1;
      }
    }
    if (num_vals > masterL2.rows()) {
      std::cout << "Only " << masterL2.rows()
                << " eigenvalues can be calculated."
                << "Storing only that many." << '\n';
      num_vals = masterL2.rows();
    }

    int MPcount_sum = MPcount.sum();
    converged = MPcount_sum;

    eigenvectors.resize(c, num_vals);
    //eigenvectorsMat.resize(num_vals);
    //for (int i = 0; i < num_vals; i++) {
    //  eigenvectorsMat[i].resize(c, 1);
    //}
    eigenvalues.resize(num_vals);
    MPorNot.resize(num_vals);
    beta.resize(num_vals);
    alpha.resize(num_vals);
    converged = MPcount_sum;
    if (MPcount_sum >= num_vals) {
      for (int i = 0; i < c * (N + 1); i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
        // std::cout << eigs.eigenvalues()[i] << " " << sorter[i].numMp <<
        // "\n"; std::cin > > bre;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < num_vals; i++) {
        eigenvalues[i] = eigenvalues_temp(int(eigenval_trunc_sorter[i][1]), 0);
        alpha[i] = alpha_temp(int(eigenval_trunc_sorter[i][1]), 0);
        beta[i] = beta_temp(int(eigenval_trunc_sorter[i][1]), 0);
        MPorNot[i] = 1;
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];

        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < num_vals; i++) {
          eigenvectors(j,i) = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          eigenvectors(j,i).dct_flag = SIS_CHEB_SPACE;
          //eigenvectorsMat[i][j] = evecMat_master.block(
          //    j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          //eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }
    } else {
      std::cout << "Last " << num_vals - MPcount_sum
                << " eigenvectors are not resolved to machine precision."
                << '\n';
      for (int i = 0; i < c * (N + 1); i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < MPcount_sum; i++) {
        eigenvalues[i] = eigenvalues_temp(int(eigenval_trunc_sorter[i][1]), 0);
        alpha[i] = alpha_temp(int(eigenval_trunc_sorter[i][1]), 0);
        beta[i] = beta_temp(int(eigenval_trunc_sorter[i][1]), 0);
        MPorNot[i] = 1;
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];
        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < MPcount_sum; i++) {
          eigenvectors(j,i) = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          eigenvectors(j,i).dct_flag = SIS_CHEB_SPACE;
        //  eigenvectorsMat[i][j] = evecMat_master.block(
        //      j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
        //  eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }

      std::vector<std::vector<T> > eigenval_trunc_sorter2;
      eigenval_trunc_sorter2.resize(c * (N + 1) - MPcount_sum);
      for (int i = MPcount_sum; i < c * (N + 1); i++) {
        eigenval_trunc_sorter2[i - MPcount_sum].resize(2);
        eigenval_trunc_sorter2[i - MPcount_sum][0] =
            sorter[int(eigenval_trunc_sorter[i][1])].minVal;
        eigenval_trunc_sorter2[i - MPcount_sum][1] =
            eigenval_trunc_sorter[i][1];
      }
      std::sort(eigenval_trunc_sorter2.begin(), eigenval_trunc_sorter2.end());
      for (int i = MPcount_sum; i < num_vals; i++) {
        eigenvalues[i] = eigenvalues_temp(
            int(eigenval_trunc_sorter2[i - MPcount_sum][1]), 0);
        alpha[i] =
            alpha_temp(int(eigenval_trunc_sorter2[i - MPcount_sum][1]), 0);
        beta[i] = beta_temp(int(eigenval_trunc_sorter2[i - MPcount_sum][1]), 0);
        MPorNot[i] = 0;
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];
        Mat.compute(n);
        temp_vec.resize(N + 1 + n, 1);
        for (int i = MPcount_sum; i < num_vals; i++) {
          eigenvectors(j,i) = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter2[i - MPcount_sum][1]),
              N + 1, 1);
          eigenvectors(j,i).dct_flag = SIS_CHEB_SPACE;
      //    eigenvectorsMat[i][j] = evecMat_master.block(
      //        j * (N + 1), int(eigenval_trunc_sorter2[i - MPcount_sum][1]),
      //        N + 1, 1);
      //    eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }
    }
#else
    eigenvectors.resize(c, c * (N + 1));
    //eigenvectorsMat.resize(c * (N + 1));
    //for (int i = 0; i < c * (N + 1); i++) {
    //  eigenvectorsMat[i].resize(c, 1);
    //}
    eigenvalues.resize(c * (N + 1));
    MPorNot.resize(c * (N + 1));
    beta.resize(c * (N + 1));
    alpha.resize(c * (N + 1));
    converged = -2;
    for (int i = 0; i < c * (N + 1); i++) {
      eigenvalues[i] = eigenvalues_temp(i, 0);
      alpha[i] = alpha_temp(i, 0);
      beta[i] = beta_temp(i, 0);
      MPorNot[i] = -2;
    }
    for (int j = 0; j < c; j++) {
      int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                  ? highest_each_columnL[j]
                  : highest_each_columnM[j];

      temp_vec.resize(N + 1 + n, 1);
      for (int i = 0; i < c * (N + 1); i++) {
        eigenvectors(j,i) = evecMat_master.block(j * (N + 1), i, N + 1, 1);
        eigenvectors(j,i).dct_flag = SIS_CHEB_SPACE;
    //    eigenvectorsMat[i][j] = evecMat_master.block(j * (N + 1), i, N + 1, 1);
    //    eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
      }
    }
#endif
#else

    Eigen::ColPivHouseholderQR<
        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
        solver(masterM2);
    Eigen::ComplexEigenSolver<
        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
        eigs;
    bool Is_M_Invertible;
    if (solver.isInvertible()) {
      //  std::cout << "M is invertible" << '\n';
      //  std::cin >> bre;
      Is_M_Invertible = true;
      eigs.compute(solver.inverse() * masterL2);
      // std::cout << "Eigenvalues :\n" << eigs.eigenvalues() << '\n';
    } else {
      // std::cout << "M is not invertible." << '\n';
      // std::cin >> bre;
      Is_M_Invertible = false;
      solver.compute(masterL2);
      eigs.compute(solver.inverse() * masterM2);
      // std::cout << "Eigenvalues :\n" << 1 / eigs.eigenvalues().array() <<
      // '\n';
    }
    // std::cin >> bre;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
        evecMat_master(c * (N + 1), c * (N + 1));
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> temp_vec;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> fir_n_mat =
        subs_mat * eigs.eigenvectors();
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> evecMat(
        c * (N + 1) + total_boundary_conditions, c * (N + 1));
    evecMat << fir_n_mat, eigs.eigenvectors();

    // Unpermute eveMat:
    evecMat = P * evecMat;
    row_counter = 0;
    col_counter = 0;
    for (int i = 0; i < c; i++) {
      int n = (highest_each_columnL[i] > highest_each_columnM[i])
                  ? highest_each_columnL[i]
                  : highest_each_columnM[i];
      Mat.compute(n);
      evecMat_master.block(i * (N + 1), 0, N + 1, c * (N + 1)) =
          Mat.mats2[n] * evecMat.block(row_counter, 0, N + 1 + n, c * (N + 1));
      row_counter += N + 1 + n;
    }

#ifndef SIS_DONT_SORT
    std::vector<std::vector<T> > eigenval_trunc_sorter;

    eigenval_trunc_sorter.resize(c * (N + 1));
    for (int i = 0; i < c * (N + 1); i++) {
      eigenval_trunc_sorter[i].resize(2);
    }

    Eigen::Matrix<int, Eigen::Dynamic, 1> MPcount(c * (N + 1));
    MPcount.setConstant(0);
    std::vector<EigenSorter<T> > sorter;
    sorter.resize(c * (N + 1));
    for (int i = 0; i < c * (N + 1); i++) {
      sorter[i].compute(evecMat_master.block(0, i, c * (N + 1), 1), c);
      if (!Is_M_Invertible) {
        // If M is not invertible, there will be zero eigenvalues, so that 1/0 =
        // inf eigenvalue, remove them.
        if (abs(eigs.eigenvalues()[i]) < 1e-11) {
          sorter[i].numMp = N + 1;
          sorter[i].Mp = false;
          sorter[i].minVal = 1e11;
        }
      }
      // sorter[i].compute(evecMat_master.block(0, i, (N + 1), 1));
      if (sorter[i].Mp == true) {
        MPcount[i] = 1;
      }
    }
    if (num_vals > masterL2.rows()) {
      std::cout << "Only " << masterL2.rows()
                << " eigenvalues can be calculated."
                << "Storing only that many." << '\n';
      num_vals = masterL2.rows();
    }

    int MPcount_sum = MPcount.sum();
    converged = MPcount_sum;
    eigenvectors.resize(c,num_vals);
    //eigenvectorsMat.resize(num_vals);
    //for (int i = 0; i < num_vals; i++) {
    //  eigenvectorsMat[i].resize(c, 1);
    //}
    eigenvalues.resize(num_vals);
    MPorNot.resize(num_vals);
    if (MPcount_sum >= num_vals) {
      for (int i = 0; i < c * (N + 1); i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
        // std::cout << eigs.eigenvalues()[i] << " " << sorter[i].numMp <<
        // "\n"; std::cin > > bre;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < num_vals; i++) {
        if (Is_M_Invertible) {
          eigenvalues[i] = eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
          MPorNot[i] = 1;
        } else {
          eigenvalues[i] =
              1.0 / eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
          MPorNot[i] = 1;
        }
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];

        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < num_vals; i++) {
          eigenvectors(j,i) = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          eigenvectors(j,i).dct_flag = SIS_CHEB_SPACE;
          //eigenvectorsMat[i][j] = evecMat_master.block(
          //    j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          //eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }
    } else {
      std::cout << "Last " << num_vals - MPcount_sum
                << " eigenvectors are not resolved to machine precision."
                << '\n';
      for (int i = 0; i < c * (N + 1); i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < MPcount_sum; i++) {
        if (Is_M_Invertible) {
          eigenvalues[i] = eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
          MPorNot[i] = 1;

        } else {
          eigenvalues[i] =
              1.0 / eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
          MPorNot[i] = 1;
        }
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];
        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < MPcount_sum; i++) {
          eigenvectors(j,i) = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          eigenvectors(j,i).dct_flag = SIS_CHEB_SPACE;
          //eigenvectorsMat[i][j] = evecMat_master.block(
          //    j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          //eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }

      std::vector<std::vector<T> > eigenval_trunc_sorter2;
      eigenval_trunc_sorter2.resize(c * (N + 1) - MPcount_sum);
      for (int i = MPcount_sum; i < c * (N + 1); i++) {
        eigenval_trunc_sorter2[i - MPcount_sum].resize(2);
        eigenval_trunc_sorter2[i - MPcount_sum][0] =
            sorter[int(eigenval_trunc_sorter[i][1])].minVal;
        eigenval_trunc_sorter2[i - MPcount_sum][1] =
            eigenval_trunc_sorter[i][1];
      }
      std::sort(eigenval_trunc_sorter2.begin(), eigenval_trunc_sorter2.end());
      for (int i = MPcount_sum; i < num_vals; i++) {
        if (Is_M_Invertible) {
          eigenvalues[i] = eigs.eigenvalues()[int(
              eigenval_trunc_sorter2[i - MPcount_sum][1])];
          MPorNot[i] = 0;

        } else {
          eigenvalues[i] =
              1.0 / eigs.eigenvalues()[int(
                        eigenval_trunc_sorter2[i - MPcount_sum][1])];
          MPorNot[i] = 0;
        }
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];
        Mat.compute(n);
        temp_vec.resize(N + 1 + n, 1);
        for (int i = MPcount_sum; i < num_vals; i++) {
          eigenvectors(j,i) = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter2[i - MPcount_sum][1]),
              N + 1, 1);
          eigenvectors(j,i).dct_flag = SIS_CHEB_SPACE;
          //eigenvectorsMat[i][j] = evecMat_master.block(
          //    j * (N + 1), int(eigenval_trunc_sorter2[i - MPcount_sum][1]),
          //    N + 1, 1);
          //eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }
    }

#else
    eigenvectors.resize(c, c * (N + 1))
    //eigenvectorsMat.resize(c * (N + 1));
    //for (int i = 0; i < c * (N + 1); i++) {
    //  eigenvectorsMat[i].resize(c, 1);
    //}
    eigenvalues.resize(c * (N + 1));
    MPorNot.resize(c * (N + 1));
    converged = -2;
    for (int i = 0; i < c * (N + 1); i++) {
      if (Is_M_Invertible) {
        eigenvalues[i] = eigs.eigenvalues()[i];
        MPorNot[i] = -2;
      } else {
        eigenvalues[i] = 1.0 / eigs.eigenvalues()[i];
        MPorNot[i] = -2;
      }
    }
    for (int j = 0; j < c; j++) {
      int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                  ? highest_each_columnL[j]
                  : highest_each_columnM[j];

      temp_vec.resize(N + 1 + n, 1);
      for (int i = 0; i < c * (N + 1); i++) {
        eigenvectors(j,i) = evecMat_master.block(j * (N + 1), i, N + 1, 1);
        eigenvectors(j,i).dct_flag = SIS_CHEB_SPACE;
      //  eigenvectorsMat[i][j] = evecMat_master.block(j * (N + 1), i, N + 1, 1);
      //  eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
      }
    }
#endif
#endif
  // Change everything to physical space:
  eigenvectors.c2p();
  }

#ifdef SIS_USE_FEAST
  void feast_compute(const LinopMat<std::complex<T> > &Lmat_,
                     const LinopMat<std::complex<T> > &Mmat_, int num_vals,
                     const BcMat<std::complex<T> > &Lbc_) {
    LinopMat<std::complex<T> > Lmat = Lmat_;
    LinopMat<std::complex<T> > Mmat = Mmat_;
    BcMat<std::complex<T> > Lbc = Lbc_;
    int bre;
    int total_of_all_orders = 0;
    int total_boundary_conditions = 0;
    if (Lmat.r != Lmat.c) {
      std::cout << "Solution only possible with square LinopMats. Exiting ..."
                << '\n';
      exit(1);
    }
    if (Mmat.r != Mmat.c) {
      std::cout << "Solution only possible with square LinopMats. Exiting ..."
                << '\n';
      exit(1);
    }
    if (Mmat.c != Lmat.c) {
      std::cout << "Both matrices have to be of same size. Exiting ..." << '\n';
      exit(1);
    }
    int r = Lmat.r, c = Lmat.c;
    // Find the highest derivative in each column. To do this create a vector
    // highest_each_column, and a temp_int_vec that will hold all values of a
    // given column to and the maximum will be stored in highest_each_column
    std::vector<int> highest_each_columnL, highest_each_columnM,
        highest_each_column;
    highest_each_columnL.resize(c);
    highest_each_columnM.resize(c);
    highest_each_column.resize(c);

    std::vector<int> temp_vec_int;
    temp_vec_int.resize(r);
    for (int j = 0; j < c; j++) {
      for (int i = 0; i < r; i++) {
        temp_vec_int[i] = Lmat(i, j).n;
      }
      highest_each_columnL[j] =
          *std::max_element(temp_vec_int.begin(), temp_vec_int.end());
    }
    for (int j = 0; j < c; j++) {
      for (int i = 0; i < r; i++) {
        temp_vec_int[i] = Mmat(i, j).n;
      }
      highest_each_columnM[j] =
          *std::max_element(temp_vec_int.begin(), temp_vec_int.end());
    }
    for (int i = 0; i < c; i++) {
      total_of_all_orders += (highest_each_columnL[i] > highest_each_columnM[i])
                                 ? highest_each_columnL[i]
                                 : highest_each_columnM[i];
      highest_each_column[i] =
          (highest_each_columnL[i] > highest_each_columnM[i])
              ? highest_each_columnL[i]
              : highest_each_columnM[i];
      //  std::cout << "highest_each_column["<< i << "]: " <<
      //  highest_each_column[i] << '\n';
    }
    total_boundary_conditions = Lbc.m;

    // total_of_all_orders has to be equal to total number of boundary
    // conditions, else the problem is ill-posed, if ill-posed, cout the same
    // and exit.
    if (total_of_all_orders != total_boundary_conditions) {
      std::cout << "The problem is ill-posed, the total of the highest "
                   "orders of all "
                   "dependent variables has to be equal to the total number of "
                   "boundary conditions specified."
                   "\n "
                   "Total no. of boundary conditions: "
                << total_boundary_conditions
                << "\n"
                   "Total of all orders: "
                << total_of_all_orders << "\n Exiting ...\n";
      exit(1);
    }
    // Declare the master matrix L:
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> masterL(
        r * (N + 1), c * (N + 1) + total_boundary_conditions);
    // Declare the master matrix M:
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> masterM(
        r * (N + 1), c * (N + 1) + total_boundary_conditions);

    masterL.setConstant(0.0);
    masterM.setConstant(0.0);
    MatGen<T> Mat;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> subs_mat;
    subs_mat.resize(total_boundary_conditions, (N + 1) * c);
    subs_mat.setConstant(0.0);
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> mat_temp(
        total_boundary_conditions, total_boundary_conditions),
        constraints(total_boundary_conditions,
                    c * (N + 1) + total_of_all_orders);

    mat_temp.setConstant(0.0);
    int row_counter = 0, col_counter = 0;

    for (int i = 0; i < Lbc.m; i++) {
      for (int j = 0; j < Lbc.n; j++) {
        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> temp =
            Lbc(i, j, highest_each_column[j]);
        constraints.block(i, col_counter, 1, temp.cols()) = temp;
        col_counter += temp.cols();
      }
      col_counter = 0;
    }

    Eigen::ColPivHouseholderQR<
        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
        qr;
    qr.compute(constraints);
    if (qr.rank() != constraints.rows()) {
      std::cout << "The boundary conditions supplied are not "
                << "  linearly independent." << '\n';
      std::cout << "qr.rank = " << qr.rank()
                << ", no. bcs: " << total_boundary_conditions << ". Exiting ..."
                << '\n';
      exit(1);
    }

    // Permutation matrix
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P;
    P = qr.colsPermutation();
    // Permute constraints
    constraints = constraints * P;
    mat_temp = constraints.block(0, 0, constraints.rows(), constraints.rows());

    Eigen::ColPivHouseholderQR<
        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
        consSolver(mat_temp);
    subs_mat = -consSolver.inverse() *
               constraints.block(0, constraints.rows(), constraints.rows(),
                                 constraints.cols() - constraints.rows());

    row_counter = 0;
    col_counter = 0;
    int master_row_counter = 0;
    int master_col_counter = 0;
    for (int j = 0; j < c; j++) {
      int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                  ? highest_each_columnL[j]
                  : highest_each_columnM[j];
      Mat.compute(n);
      for (int i = 0; i < r; i++) {
        int diffn = n - Lmat(i, j).n;
        if (Lmat(i, j).NCC == 0) {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            masterL.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Lmat(i, j).coef[k] *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, N + 1 + n));
          }
        } else {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            masterL.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Lmat(i, j).coefFun[k].MultMat().block(0, 0, N + 1, N + 1) *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, N + 1 + n));
          }
        }
        diffn = n - Mmat(i, j).n;
        if (Mmat(i, j).NCC == 0) {
          for (int k = 0; k < Mmat(i, j).n + 1; k++) {
            masterM.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Mmat(i, j).coef[k] *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, N + 1 + n));
          }
        } else {
          for (int k = 0; k < Mmat(i, j).n + 1; k++) {
            masterM.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Mmat(i, j).coefFun[k].MultMat().block(0, 0, N + 1, N + 1) *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, N + 1 + n));
          }
        }
        master_row_counter += N + 1;
      }
      master_row_counter = 0;
      master_col_counter += N + 1 + n;
    }
    // Permute columns of M and L:
    masterL = masterL * P;
    masterM = masterM * P;

    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> masterL2,
        masterM2;

    masterL2 =
        masterL.block(0, constraints.rows(), c * (N + 1), c * (N + 1)) +
        (masterL.block(0, 0, c * (N + 1), constraints.rows()) * subs_mat);

    masterM2 =
        masterM.block(0, constraints.rows(), c * (N + 1), c * (N + 1)) +
        (masterM.block(0, 0, c * (N + 1), constraints.rows()) * subs_mat);

    int ldL = masterL2.outerStride(); // ld for leading dimension
    int ldM = masterM2.outerStride();
    double epsout;
    double Emid[2];
    int loop;
    int M;
    // feast::M0 = 2;
    double res[2 * feast::M0];
    feast::feast_init();
    // feast::feastparam[6] = 8;
    feast::feastparam[0] = 1;
    feast::feastparam[7] = 56;
    // feast::feastparam[13] = 2;
    // Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> lambda(
    //  masterL2.rows(), 1), vecs(masterL2.rows(), 2 * masterL2.rows()),
    //  vr;

    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> lambda(
        feast::M0, 1),
        vecs(masterM2.rows(), 2 * feast::M0), vr;
    lambda.setConstant(0.0);
    vecs.setConstant(0.0);
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
        eigenvalues_temp;
    Emid[0] = std::real(feast::center);
    Emid[1] = std::imag(feast::center);

    std::cout << "in " << __LINE__ << '\n' << std::flush;
    zfeast_gegv_(&N, masterL2.data(), &ldL, masterM2.data(), &ldM,
                 feast::feastparam, &epsout, &loop, Emid, &feast::radius,
                 &feast::M0, lambda.data(), vecs.data(), &M, res, &feast::info);
    feast::display();
    // feast::M0 = M + 10;
    // lambda.resize(feast::M0, 1);
    // vecs.resize(feast::M0, 2 * feast::M0);
    // feast::feast_init();
    // zfeast_gegv_(&N, masterL2.data(), &ldL, masterM2.data(), &ldM,
    // feast::feastparam, &epsout, &loop, Emid, &feast::radius, &feast::M0,
    //lambda.data(), vecs.data(), &M, res, &feast::info);
    std::cout << "M = " << M << '\n';
    // std::cout << "lambda = " << lambda.block(0,0,M,1) << '\n';
    eigenvalues_temp = lambda.block(0, 0, M, 1);
    vr.resize(masterL2.rows(), M);
    vr = vecs.block(0, 0, vecs.rows(), M);
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
        evecMat_master(c * (N + 1), M);
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> temp_vec;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> fir_n_mat =
        subs_mat * vr;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> evecMat(
        c * (N + 1) + total_boundary_conditions, M);
    evecMat << fir_n_mat, vr;
    // Unpermute eveMat:
    evecMat = P * evecMat;
    row_counter = 0;
    col_counter = 0;
    for (int i = 0; i < c; i++) {
      int n = (highest_each_columnL[i] > highest_each_columnM[i])
                  ? highest_each_columnL[i]
                  : highest_each_columnM[i];
      Mat.compute(n);
      evecMat_master.block(i * (N + 1), 0, N + 1, M) =
          Mat.mats2[n] * evecMat.block(row_counter, 0, N + 1 + n, M);
      row_counter += N + 1 + n;
    }
    std::vector<std::vector<T> > eigenval_trunc_sorter;

    eigenval_trunc_sorter.resize(M);
    for (int i = 0; i < M; i++) {
      eigenval_trunc_sorter[i].resize(2);
    }

    Eigen::Matrix<int, Eigen::Dynamic, 1> MPcount(M);
    MPcount.setConstant(0);
    std::vector<EigenSorter<T> > sorter;
    sorter.resize(M);
    for (int i = 0; i < M; i++) {
      sorter[i].compute(evecMat_master.block(0, i, c * (N + 1), 1), c);
      // sorter[i].compute(evecMat_master.block(0, i, (N + 1), 1));
      if (sorter[i].Mp == true) {
        MPcount[i] = 1;
      }
    }
    if (num_vals > M) {
      std::cout << "Only " << M << " eigenvalues can be calculated."
                << "Storing only that many." << '\n';
      num_vals = M;
    }

    int MPcount_sum = MPcount.sum();
    converged = MPcount_sum;

    eigenvectorsMat.resize(num_vals);
    for (int i = 0; i < num_vals; i++) {
      eigenvectorsMat[i].resize(c, 1);
    }
    eigenvalues.resize(num_vals);
    converged = MPcount_sum;
    if (MPcount_sum >= num_vals) {
      for (int i = 0; i < M; i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
        // std::cout << eigs.eigenvalues()[i] << " " << sorter[i].numMp <<
        // "\n"; std::cin > > bre;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < num_vals; i++) {
        eigenvalues[i] = eigenvalues_temp(int(eigenval_trunc_sorter[i][1]), 0);
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];

        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < num_vals; i++) {
          eigenvectorsMat[i][j] = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }
    } else {
      std::cout << "Last " << num_vals - MPcount_sum
                << " eigenvectors are not resolved to machine precision."
                << '\n';
      for (int i = 0; i < M; i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < MPcount_sum; i++) {
        eigenvalues[i] = eigenvalues_temp(int(eigenval_trunc_sorter[i][1]), 0);
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];
        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < MPcount_sum; i++) {
          eigenvectorsMat[i][j] = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }

      std::vector<std::vector<T> > eigenval_trunc_sorter2;
      eigenval_trunc_sorter2.resize(M - MPcount_sum);
      for (int i = MPcount_sum; i < M; i++) {
        eigenval_trunc_sorter2[i - MPcount_sum].resize(2);
        eigenval_trunc_sorter2[i - MPcount_sum][0] =
            sorter[int(eigenval_trunc_sorter[i][1])].minVal;
        eigenval_trunc_sorter2[i - MPcount_sum][1] =
            eigenval_trunc_sorter[i][1];
      }
      std::sort(eigenval_trunc_sorter2.begin(), eigenval_trunc_sorter2.end());
      for (int i = MPcount_sum; i < num_vals; i++) {
        eigenvalues[i] = eigenvalues_temp(
            int(eigenval_trunc_sorter2[i - MPcount_sum][1]), 0);
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];
        Mat.compute(n);
        temp_vec.resize(N + 1 + n, 1);
        for (int i = MPcount_sum; i < num_vals; i++) {
          eigenvectorsMat[i][j] = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter2[i - MPcount_sum][1]),
              N + 1, 1);
          eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }
    }
  }
#endif
  /// Use this to solve eigenvalue problems where eigenvalues appears in the
  /// boundary condition, example, in eigenvalue problems with fluid-fluid
  /// interfaces. Again
  /// solvability must hold. That total of all the highest orders of every
  /// independent variable has to be equal to the total number of boundary
  /// conditions specified + the number of constraints.
  /// Further, please ensure that the boundary conditions and constraints are
  /// linearly independent, else this will throw an error.
  void compute_with_constraints(const LinopMat<std::complex<T> > &Lmat_,
                                const LinopMat<std::complex<T> > &Mmat_,
                                int num_vals,
                                const BcMat<std::complex<T> > &Lbc_,
                                const BcMat<std::complex<T> > &Mbc_) {

    LinopMat<std::complex<T> > Lmat = Lmat_;
    LinopMat<std::complex<T> > Mmat = Mmat_;
    BcMat<std::complex<T> > Lbc = Lbc_;
    BcMat<std::complex<T> > Mbc = Mbc_;
    int bre;
    int total_of_all_orders = 0;
    int total_boundary_conditions = 0;
    if (Lmat.r != Lmat.c) {
      std::cout << "Solution only possible with square LinopMats. Exiting ..."
                << '\n';
      exit(1);
    }
    if (Mmat.r != Mmat.c) {
      std::cout << "Solution only possible with square LinopMats. Exiting ..."
                << '\n';
      exit(1);
    }
    if (Mmat.c != Lmat.c) {
      std::cout << "Both matrices have to be of same size. Exiting ..." << '\n';
      exit(1);
    }
    int r = Lmat.r, c = Lmat.c;
    // Find the highest derivative in each column. To do this create a vector
    // highest_each_column, and a temp_int_vec that will hold all values of a
    // given column to and the maximum will be stored in highest_each_column
    std::vector<int> highest_each_columnL, highest_each_columnM,
        highest_each_column, num_bc_each_var;
    highest_each_columnL.resize(c);
    highest_each_columnM.resize(c);
    highest_each_column.resize(c);
    num_bc_each_var.resize(c);
    std::vector<int> temp_vec_int;
    temp_vec_int.resize(r);
    for (int j = 0; j < c; j++) {
      for (int i = 0; i < r; i++) {
        temp_vec_int[i] = Lmat(i, j).n;
      }
      highest_each_columnL[j] =
          *std::max_element(temp_vec_int.begin(), temp_vec_int.end());
    }
    for (int j = 0; j < c; j++) {
      for (int i = 0; i < r; i++) {
        temp_vec_int[i] = Mmat(i, j).n;
      }
      highest_each_columnM[j] =
          *std::max_element(temp_vec_int.begin(), temp_vec_int.end());
    }
    for (int i = 0; i < c; i++) {
      total_of_all_orders += (highest_each_columnL[i] > highest_each_columnM[i])
                                 ? highest_each_columnL[i]
                                 : highest_each_columnM[i];
      highest_each_column[i] =
          (highest_each_columnL[i] > highest_each_columnM[i])
              ? highest_each_columnL[i]
              : highest_each_columnM[i];
    }
    for (int i = 0; i < c; i++) {
      total_boundary_conditions += Lmat.BcVec[i].nbc();
      num_bc_each_var[i] = Lmat.BcVec[i].nbc();
    }

    if (Lbc.m != Mbc.m) {
      std::cout << "The Lbc and Mbc have to be of same dimensions" << '\n'
                << "Exiting in " << __LINE__ << "...\n";

      exit(1);
    }
    // total_of_all_orders has to be equal to total number of boundary
    // conditions, else the problem is ill-posed, if ill-posed, cout the same
    // and exit.
    if (total_of_all_orders != (total_boundary_conditions + Lbc.m)) {
      std::cout << "The problem is ill-posed, the total of the highest "
                   "orders of all "
                   "dependent variables has to be equal to the total number of "
                   "boundary conditions plus number of constraints specified."
                   "\n "
                   "Total no. of boundary conditions: "
                << total_boundary_conditions
                << "\n"
                   "Total of all orders: "
                << total_of_all_orders << "\n Number of Constraints : " << Lbc.m
                << "\n Exiting in " << __LINE__ << " ...\n ";
      exit(1);
    }
    // Declare the master matrix L:
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> masterL(
        r * (N + 1) + Lbc.m, c * (N + 1) + Lbc.m);
    // Declare the master matrix M:
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> masterM(
        r * (N + 1) + Lbc.m, c * (N + 1) + Lbc.m);

    masterL.setConstant(0.0);
    masterM.setConstant(0.0);
    MatGen<T> Mat;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> subs_mat;

    subs_mat.resize(total_boundary_conditions, (N + 1) * c + Lbc.m);

    subs_mat.setConstant(0.0);
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> mat_temp(
        total_boundary_conditions, total_boundary_conditions);
    mat_temp.setConstant(0.0);
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
        subs_mat_col_mul(c * (N + 1), total_boundary_conditions);
    subs_mat_col_mul.setConstant(0.0);
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
        subs_mat_col_mul2(c * (N + 1), total_boundary_conditions);
    subs_mat_col_mul2.setConstant(0.0);
    int subs_mat_row_counter = 0;
    int subs_mat_col_counter = 0;
    int row_counter = 0, col_counter = 0;
    for (int j = 0; j < c; j++) {
      // overall order is the one which is greater in each column:
      int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                  ? highest_each_columnL[j]
                  : highest_each_columnM[j];
      // Store the LinopMat with the larger degree to setup boundary conditions
      // LinopMat<std::complex<T> > Tmat =
      //    (highest_each_columnL[j] > highest_each_columnM[j]) ? Lmat : Mmat;

      Mat.compute(n);
      std::vector<Eigen::Matrix<std::complex<T>, 1, Eigen::Dynamic> > lbc_vecs(
          Lmat.BcVec[j].nl,
          Eigen::Matrix<std::complex<T>, 1, Eigen::Dynamic>(N + 1 + 2 * n));
      std::vector<Eigen::Matrix<std::complex<T>, 1, Eigen::Dynamic> > rbc_vecs(
          Lmat.BcVec[j].nr,
          Eigen::Matrix<std::complex<T>, 1, Eigen::Dynamic>(N + 1 + 2 * n));
      std::vector<Eigen::Matrix<std::complex<T>, 1, Eigen::Dynamic> >
          lbc_con_vecs(Lmat.BcVec[j].nl,
                       Eigen::Matrix<std::complex<T>, 1, Eigen::Dynamic>(n));
      std::vector<Eigen::Matrix<std::complex<T>, 1, Eigen::Dynamic> >
          rbc_con_vecs(Lmat.BcVec[j].nr,
                       Eigen::Matrix<std::complex<T>, 1, Eigen::Dynamic>(n));
      Eigen::Matrix<T, 1, Eigen::Dynamic> ones(N + 1 + 2 * n);
      Eigen::Matrix<T, 1, Eigen::Dynamic> onem(N + 1 + 2 * n);
      Eigen::Matrix<T, 1, Eigen::Dynamic> onesn(n);
      Eigen::Matrix<T, 1, Eigen::Dynamic> onemn(n);
      // This follows as T_n(1) = 1.0.
      ones.setConstant(1.0);
      onesn.setConstant(1.0);
      onem.setConstant(1.0);
      onemn.setConstant(1.0);
      // This follows as T_n(-1) = (-1)^n.
      for (int k = 1; k < N + 1 + 2 * n; k = k + 2) {
        onem[k] = -1.0;
      }
      for (int k = 1; k < n; k = k + 2) {
        onemn[k] = -1.0;
      }
      if (n > 0) {
        ones[0] = 0.5;
        onesn[0] = 0.5;
        onem[0] = 0.5;
        onemn[0] = 0.5;
      }
      // Next just multiply matrices based on the order.
      for (int k = 0; k < Lmat.BcVec[j].nl; k++) {
        lbc_vecs[k].resize(N + 1 + 2 * n);
        lbc_vecs[k].setConstant(std::complex<T>(0.0, 0.0));
        for (int l = 0; l < Lmat.BcVec[j].ord + 1; l++) {
          lbc_vecs[k] += Lmat.BcVec[j].coefl(k, l) *
                         (onem * Mat.mats[l + (n - Lmat.BcVec[j].ord)]);
        }
      }

      for (int k = 0; k < Lmat.BcVec[j].nr; k++) {
        rbc_vecs[k].resize(N + 1 + 2 * n);
        rbc_vecs[k].setConstant(std::complex<T>(0.0, 0.0));
        for (int l = 0; l < Lmat.BcVec[j].ord + 1; l++) {
          rbc_vecs[k] += Lmat.BcVec[j].coefr(k, l) *
                         (ones * Mat.mats[l + (n - Lmat.BcVec[j].ord)]);
        }
      }
      for (int k = 0; k < Lmat.BcVec[j].nl; k++) {
        lbc_con_vecs[k].setConstant(0.0);
        for (int l = 0; l < Lmat.BcVec[j].ord + 1; l++) {
          if (l + n - Lmat.BcVec[j].ord - 1 > -1) {
            lbc_con_vecs[k] +=
                Lmat.BcVec[j].coefl(k, l) *
                (onemn * Mat.con_mats[l + n - Lmat.BcVec[j].ord - 1]);
          }
        }
      }
      for (int k = 0; k < Lmat.BcVec[j].nr; k++) {
        rbc_con_vecs[k].setConstant(0.0);
        for (int l = 0; l < Lmat.BcVec[j].ord + 1; l++) {
          if (l + n - Lmat.BcVec[j].ord - 1 > -1) {
            rbc_con_vecs[k] +=
                Lmat.BcVec[j].coefr(k, l) *
                (onesn * Mat.con_mats[l + n - Lmat.BcVec[j].ord - 1]);
          }
        }
      }

      for (int k = 0; k < Lmat.BcVec[j].nl; k++) {
        mat_temp.block(row_counter, col_counter, 1,
                       Lmat.BcVec[j].nl + Lmat.BcVec[j].nr) =
            lbc_vecs[k].head(Lmat.BcVec[j].nl + Lmat.BcVec[j].nr);
        row_counter++;
      }
      for (int k = Lmat.BcVec[j].nl; k < Lmat.BcVec[j].nl + Lmat.BcVec[j].nr;
           k++) {
        mat_temp.block(row_counter, col_counter, 1,
                       Lmat.BcVec[j].nl + Lmat.BcVec[j].nr) =
            rbc_vecs[k - Lmat.BcVec[j].nl].head(Lmat.BcVec[j].nl +
                                                Lmat.BcVec[j].nr);
        row_counter++;
      }
      col_counter = col_counter + Lmat.BcVec[j].nbc();
      for (int k = 0; k < Lmat.BcVec[j].nl; k++) {
        subs_mat.block(subs_mat_row_counter, subs_mat_col_counter, 1,
                       (N + 1 - Lmat.BcVec[j].nbc())) =
            lbc_vecs[k].block(0, Lmat.BcVec[j].nbc(), 1,
                              N + 1 - Lmat.BcVec[j].nbc());
        subs_mat.block(subs_mat_row_counter,
                       subs_mat_col_counter + (N + 1 - Lmat.BcVec[j].nbc()), 1,
                       n) = lbc_con_vecs[k];
        subs_mat_row_counter++;
      }

      for (int k = Lmat.BcVec[j].nl; k < Lmat.BcVec[j].nbc(); k++) {
        subs_mat.block(subs_mat_row_counter, subs_mat_col_counter, 1,
                       (N + 1 - Lmat.BcVec[j].nbc())) =
            rbc_vecs[k - Lmat.BcVec[j].nl].block(0, Lmat.BcVec[j].nbc(), 1,
                                                 N + 1 - Lmat.BcVec[j].nbc());
        subs_mat.block(subs_mat_row_counter,
                       subs_mat_col_counter + (N + 1 - Lmat.BcVec[j].nbc()), 1,
                       n) = rbc_con_vecs[k - Lmat.BcVec[j].nl];
        subs_mat_row_counter++;
      }
      subs_mat_col_counter += (N + 1 + n - Lmat.BcVec[j].nbc());
    }
    if (mat_temp.rows() != 0 || mat_temp.cols() != 0) {
      Eigen::ColPivHouseholderQR<
          Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
          consSolver(mat_temp);
      if (!consSolver.isInvertible()) {
        std::cout << "the matrix is not invertible." << '\n';
        exit(1);
      }
      subs_mat = -consSolver.inverse() * subs_mat.eval();
    }
    row_counter = 0;
    col_counter = 0;
    int master_row_counter = 0;
    int master_col_counter = 0;
    subs_mat_row_counter = 0;
    subs_mat_col_counter = 0;
    for (int j = 0; j < c; j++) {
      int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                  ? highest_each_columnL[j]
                  : highest_each_columnM[j];
      Mat.compute(n);
      for (int i = 0; i < r; i++) {

        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
            solver_matL(N + 1, N + 1 + n - Lmat.BcVec[j].nbc());
        solver_matL.setConstant(0.0);
        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
            solver_matM(N + 1, N + 1 + n - Lmat.BcVec[j].nbc());
        solver_matM.setConstant(0.0);
        int diffn = n - Lmat(i, j).n;
        if (Lmat(i, j).NCC == 0) {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            solver_matL.block(0, 0, N + 1, N + 1 + n - Lmat.BcVec[j].nbc()) +=
                Lmat(i, j).coef[k] *
                (Mat.mats2[k + diffn].block(0, Lmat.BcVec[j].nbc(), N + 1,
                                            N + 1 + n - Lmat.BcVec[j].nbc()));
            subs_mat_col_mul.block(i * (N + 1), subs_mat_col_counter, N + 1,
                                   Lmat.BcVec[j].nbc()) +=
                Lmat(i, j).coef[k] *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, Lmat.BcVec[j].nbc()));
          }
        } else {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            solver_matL.block(0, 0, N + 1, N + 1 + n - Lmat.BcVec[j].nbc()) +=
                Lmat(i, j).coefFun[k].MultMat().block(0, 0, N + 1, N + 1) *
                (Mat.mats2[k + diffn].block(0, Lmat.BcVec[j].nbc(), N + 1,
                                            N + 1 + n - Lmat.BcVec[j].nbc()));
            subs_mat_col_mul.block(i * (N + 1), subs_mat_col_counter, N + 1,
                                   Lmat.BcVec[j].nbc()) +=
                Lmat(i, j).coefFun[k].MultMat().block(0, 0, N + 1, N + 1) *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, Lmat.BcVec[j].nbc()));
          }
        }
        diffn = n - Mmat(i, j).n;
        if (Mmat(i, j).NCC == 0) {
          for (int k = 0; k < Mmat(i, j).n + 1; k++) {
            solver_matM.block(0, 0, N + 1, N + 1 + n - Lmat.BcVec[j].nbc()) +=
                Mmat(i, j).coef[k] *
                (Mat.mats2[k + diffn].block(0, Lmat.BcVec[j].nbc(), N + 1,
                                            N + 1 + n - Lmat.BcVec[j].nbc()));
            subs_mat_col_mul2.block(i * (N + 1), subs_mat_col_counter, N + 1,
                                    Lmat.BcVec[j].nbc()) +=
                Mmat(i, j).coef[k] *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, Lmat.BcVec[j].nbc()));
          }
        } else {
          for (int k = 0; k < Mmat(i, j).n + 1; k++) {
            solver_matM.block(0, 0, N + 1, N + 1 + n - Lmat.BcVec[j].nbc()) +=
                Mmat(i, j).coefFun[k].MultMat().block(0, 0, N + 1, N + 1) *
                (Mat.mats2[k + diffn].block(0, Lmat.BcVec[j].nbc(), N + 1,
                                            N + 1 + n - Lmat.BcVec[j].nbc()));
            subs_mat_col_mul2.block(i * (N + 1), subs_mat_col_counter, N + 1,
                                    Lmat.BcVec[j].nbc()) +=
                Mmat(i, j).coefFun[k].MultMat().block(0, 0, N + 1, N + 1) *
                (Mat.mats2[k + diffn].block(0, 0, N + 1, Lmat.BcVec[j].nbc()));
          }
        }
        subs_mat_row_counter += N + 1 + n - Lmat.BcVec[j].nbc();
        masterL.block(i * (N + 1), master_col_counter, N + 1,
                      N + 1 + n - Lmat.BcVec[j].nbc()) = solver_matL;
        masterM.block(i * (N + 1), master_col_counter, N + 1,
                      N + 1 + n - Lmat.BcVec[j].nbc()) = solver_matM;
      }
      subs_mat_col_counter += Lmat.BcVec[j].nbc();
      subs_mat_row_counter = 0;
      master_row_counter = 0;
      master_col_counter += N + 1 + n - Lmat.BcVec[j].nbc();
    }
    if (mat_temp.rows() != 0 || mat_temp.cols() != 0) {
      masterL.block(0, 0, c * (N + 1), c * (N + 1) + Lbc.m) +=
          masterL.block(0, 0, c * (N + 1), c * (N + 1) + Lbc.m).eval() +
          (subs_mat_col_mul * subs_mat);
      masterM.block(0, 0, c * (N + 1), c * (N + 1) + Lbc.m) +=
          masterM.block(0, 0, c * (N + 1), c * (N + 1) + Lbc.m).eval() +
          (subs_mat_col_mul2 * subs_mat);
    }
    // Now append the constraints:

    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
        constraints_bigL(Lbc.m, c * (N + 1) + Lbc.m),
        constraints_smallL(Lbc.m, total_boundary_conditions),
        constraints_bigM(Lbc.m, c * (N + 1) + Lbc.m),
        constraints_smallM(Lbc.m, total_boundary_conditions);
    constraints_bigL.setConstant(0.0);
    constraints_smallL.setConstant(0.0);
    constraints_bigM.setConstant(0.0);
    constraints_smallM.setConstant(0.0);
    col_counter = 0;
    row_counter = 0;
    int col_counter2 = 0;
    int row_counter2 = 0;
    for (int i = 0; i < Lbc.m; i++) {
      for (int j = 0; j < Lbc.n; j++) {
        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> temp =
            Lbc(i, j, highest_each_column[j]);

        constraints_bigL.block(i, col_counter, 1,
                               N + 1 + highest_each_column[j] -
                                   num_bc_each_var[j]) =
            temp.block(0, num_bc_each_var[j], 1,
                       N + 1 + highest_each_column[j] - num_bc_each_var[j]);

        constraints_smallL.block(i, col_counter2, 1, num_bc_each_var[j]) =
            temp.block(0, 0, 1, num_bc_each_var[j]);
        //  std::cout << "i = " << i << ", j = " << j << ",
        //  highest_each_column["
        //            << j << "]= " << highest_each_column[j] << '\n';
        temp = Mbc(i, j, highest_each_column[j]);

        constraints_bigM.block(i, col_counter, 1,
                               N + 1 + highest_each_column[j] -
                                   num_bc_each_var[j]) =
            temp.block(0, num_bc_each_var[j], 1,
                       N + 1 + highest_each_column[j] - num_bc_each_var[j]);

        constraints_smallM.block(i, col_counter2, 1, num_bc_each_var[j]) =
            temp.block(0, 0, 1, num_bc_each_var[j]);
        col_counter += N + 1 + highest_each_column[j] - num_bc_each_var[j];
        col_counter2 += num_bc_each_var[j];
      }
      col_counter = 0;
      col_counter2 = 0;
    }
    if (mat_temp.rows() != 0 || mat_temp.cols() != 0) {
      masterL.block(c * (N + 1), 0, Lbc.m, c * (N + 1) + Lbc.m) =
          constraints_bigL + constraints_smallL * subs_mat;
      masterM.block(c * (N + 1), 0, Lbc.m, c * (N + 1) + Lbc.m) =
          constraints_bigM + constraints_smallM * subs_mat;
    } else {
      masterL.block(c * (N + 1), 0, Lbc.m, c * (N + 1) + Lbc.m) =
          constraints_bigL;
      masterM.block(c * (N + 1), 0, Lbc.m, c * (N + 1) + Lbc.m) =
          constraints_bigM;
    }

#ifndef SIS_USE_LAPACK
    Eigen::ColPivHouseholderQR<
        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
        solver(masterM);
    Eigen::ComplexEigenSolver<
        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
        eigs;
    bool Is_M_Invertible;
    if (solver.isInvertible()) {
      Is_M_Invertible = true;
      eigs.compute(solver.inverse() * masterL);
      //  std::cout << "Eigenvalues :\n" << eigs.eigenvalues() << '\n';
    } else {
      Is_M_Invertible = false;
      solver.compute(masterL);
      eigs.compute(solver.inverse() * masterM);
      //  std::cout << "Eigenvalues :\n" << 1 / eigs.eigenvalues().array() <<
      //  '\n';
    }
    // std::cin >> bre;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
        evecMat_master(c * (N + 1) + Lbc.m, c * (N + 1) + Lbc.m);
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> temp_vec;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> fir_n_mat,
        full_mat;
    if (subs_mat.size() > 0) {
      fir_n_mat = subs_mat * eigs.eigenvectors();
    }
    full_mat.resize(c * (N + 1) + Lbc.m + fir_n_mat.rows(),
                    c * (N + 1) + Lbc.m);

    row_counter = 0;
    row_counter2 = 0;
    col_counter = 0;
    // Repack eigenvectors using fir_n_mat:
    for (int i = 0; i < c; i++) {
      int n = (highest_each_columnL[i] > highest_each_columnM[i])
                  ? highest_each_columnL[i]
                  : highest_each_columnM[i];
      if (subs_mat.size() > 0) {
        full_mat.block(row_counter, 0, num_bc_each_var[i],
                       c * (N + 1) + Lbc.m) =
            fir_n_mat.block(col_counter, 0, num_bc_each_var[i],
                            c * (N + 1) + Lbc.m);
      }
      row_counter += num_bc_each_var[i];
      col_counter += num_bc_each_var[i];

      full_mat.block(row_counter, 0, N + 1 + n - num_bc_each_var[i],
                     c * (N + 1) + Lbc.m) =
          eigs.eigenvectors().block(row_counter2, 0,
                                    N + 1 + n - num_bc_each_var[i],
                                    c * (N + 1) + Lbc.m);
      row_counter += N + 1 + n - num_bc_each_var[i];
      row_counter2 += N + 1 + n - num_bc_each_var[i];
    }

    row_counter = 0;
    for (int i = 0; i < c; i++) {
      int n = (highest_each_columnL[i] > highest_each_columnM[i])
                  ? highest_each_columnL[i]
                  : highest_each_columnM[i];
      Mat.compute(n);
      evecMat_master.block(i * (N + 1), 0, N + 1, c * (N + 1) + Lbc.m) =
          Mat.mats2[n] *
          full_mat.block(row_counter, 0, N + 1 + n, c * (N + 1) + Lbc.m);
      row_counter += N + 1 + n;
    }

    std::vector<std::vector<T> > eigenval_trunc_sorter;

    eigenval_trunc_sorter.resize(c * (N + 1) + Lbc.m);
    for (int i = 0; i < c * (N + 1) + Lbc.m; i++) {
      eigenval_trunc_sorter[i].resize(2);
    }

    Eigen::Matrix<int, Eigen::Dynamic, 1> MPcount(c * (N + 1) + Lbc.m);
    MPcount.setConstant(0);
    std::vector<EigenSorter<T> > sorter;
    sorter.resize(c * (N + 1) + Lbc.m);
    for (int i = 0; i < c * (N + 1) + Lbc.m; i++) {
      sorter[i].compute(evecMat_master.block(0, i, c * (N + 1), 1), c);
      if (!Is_M_Invertible) {
        // If M is not invertible, there will be zero eigenvalues, so that 1/0 =
        // inf eigenvalue, remove them.
        if (abs(eigs.eigenvalues()[i]) < 1e-11) {
          sorter[i].numMp = N + 1;
          sorter[i].Mp = false;
          sorter[i].minVal = 1e11;
        }
      }
      // sorter[i].compute(evecMat_master.block(0, i, (N + 1), 1));
      if (sorter[i].Mp == true) {
        MPcount[i] = 1;
      }
    }
    if (num_vals > masterL.rows()) {
      std::cout << "Only " << masterL.rows()
                << " eigenvalues can be calculated."
                << "Storing only that many." << '\n';
      num_vals = masterL.rows();
    }
    int MPcount_sum = MPcount.sum();
    converged = MPcount_sum;

    eigenvectors.resize(c, num_vals);
    //eigenvectorsMat.resize(num_vals);
    //for (int i = 0; i < num_vals; i++) {
    //  eigenvectorsMat[i].resize(c, 1);
    //}
    eigenvalues.resize(num_vals);
    if (MPcount_sum >= num_vals) {
      for (int i = 0; i < c * (N + 1) + Lbc.m; i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
        // std::cout << eigs.eigenvalues()[i] << " " << sorter[i].numMp <<
        // "\n"; std::cin > > bre;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < num_vals; i++) {
        if (Is_M_Invertible) {
          eigenvalues[i] = eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
        } else {
          eigenvalues[i] =
              1.0 / eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
        }
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];

        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < num_vals; i++) {
          eigenvectors(j,i) = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          eigenvectors(j,i).dct_flag = SIS_CHEB_SPACE;
          //eigenvectorsMat[i][j] = evecMat_master.block(
          //    j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          //eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }
    } else {
      std::cout << "Last " << num_vals - MPcount_sum
                << " eigenvectors are not resolved to machine precision."
                << '\n';
      for (int i = 0; i < c * (N + 1) + Lbc.m; i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < MPcount_sum; i++) {
        if (Is_M_Invertible) {
          eigenvalues[i] = eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
        } else {
          eigenvalues[i] =
              1.0 / eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
        }
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];
        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < MPcount_sum; i++) {
          eigenvectors(j,i)  = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          eigenvectors(j,i).dct_flag = SIS_CHEB_SPACE;
          //eigenvectorsMat[i][j] = evecMat_master.block(
          //    j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          //eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }

      std::vector<std::vector<T> > eigenval_trunc_sorter2;
      eigenval_trunc_sorter2.resize(c * (N + 1) + Lbc.m - MPcount_sum);
      for (int i = MPcount_sum; i < c * (N + 1) + Lbc.m; i++) {
        eigenval_trunc_sorter2[i - MPcount_sum].resize(2);
        eigenval_trunc_sorter2[i - MPcount_sum][0] =
            sorter[int(eigenval_trunc_sorter[i][1])].minVal;
        eigenval_trunc_sorter2[i - MPcount_sum][1] =
            eigenval_trunc_sorter[i][1];
      }
      std::sort(eigenval_trunc_sorter2.begin(), eigenval_trunc_sorter2.end());
      for (int i = MPcount_sum; i < num_vals; i++) {
        if (Is_M_Invertible) {
          eigenvalues[i] = eigs.eigenvalues()[int(
              eigenval_trunc_sorter2[i - MPcount_sum][1])];
        } else {
          eigenvalues[i] =
              1.0 / eigs.eigenvalues()[int(
                        eigenval_trunc_sorter2[i - MPcount_sum][1])];
        }
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];
        Mat.compute(n);
        temp_vec.resize(N + 1 + n, 1);
        for (int i = MPcount_sum; i < num_vals; i++) {
          eigenvectors(j,i) = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter2[i - MPcount_sum][1]),
              N + 1, 1);
          eigenvectors(j,i).dct_flag = SIS_CHEB_SPACE;
          //eigenvectorsMat[i][j] = evecMat_master.block(
          //    j * (N + 1), int(eigenval_trunc_sorter2[i - MPcount_sum][1]),
          //    N + 1, 1);
          //eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }
    }

#else
    char jobvl = 'N';                // Don't compute left evecs
    char jobvr = 'V';                // Compute right evecs
    std::complex<double> wkopt;      // Eistimate optimum workspace
    std::complex<double> *work;      // allocate optimum workspace
    alpha.resize(masterL.rows(), 1); // alpha for gen. eig. prob.
    beta.resize(masterL.rows(), 1);  // beta for gen. eig. prob.

    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> vl(
        masterL.rows(), masterL.rows()),
        vr(masterL.rows(), masterL.rows()), eigenvalues_temp(masterL.rows(), 1),
        alpha_temp(masterL.rows(), 1), beta_temp(masterL.rows(), 1);
    // vl : left evecs, vr: right evecs.
    int ldL = masterL.outerStride(); // ld for leading dimension
    int ldM = masterM.outerStride();
    int ldvl = vl.outerStride();
    int ldvr = vr.outerStride();
    int sizeL = masterL.rows();
    int lwork = -1; // set lwork to -1 to estimate workspace.
    double rwork[8 * sizeL];

    // call this to estimate workspace
    zggev_(&jobvl, &jobvr, &sizeL, masterL.data(), &ldL, masterM.data(), &ldM,
           alpha_temp.data(), beta_temp.data(), vl.data(), &ldvl, vr.data(),
           &ldvr, &wkopt, &lwork, rwork, &info);

    // Now allocate workspace:
    lwork = (int)real(wkopt);
    work = (std::complex<double> *)malloc(lwork * sizeof(std::complex<double>));

    // Solve eigenvalue problem:
    zggev_(&jobvl, &jobvr, &sizeL, masterL.data(), &ldL, masterM.data(), &ldM,
           alpha_temp.data(), beta_temp.data(), vl.data(), &ldvl, vr.data(),
           &ldvr, work, &lwork, rwork, &info);

    // Free workspace.
    free((void *)work);

    eigenvalues_temp = alpha_temp.array() / beta_temp.array();

    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
        evecMat_master(c * (N + 1), c * (N + 1) + Lbc.m);
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> temp_vec;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> fir_n_mat =
        subs_mat * vr;
    row_counter = 0;
    col_counter = 0;
    for (int i = 0; i < c; i++) {
      int n = (highest_each_columnL[i] > highest_each_columnM[i])
                  ? highest_each_columnL[i]
                  : highest_each_columnM[i];
      Mat.compute(n);
      // std::cout << "n = " << n << '\n';
      // std::cin > > bre;

      Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> evecMat(
          N + 1 + n, c * (N + 1) + Lbc.m);
      evecMat.setConstant(0.0);
      evecMat.block(num_bc_each_var[i], 0, N + 1 + n - num_bc_each_var[i],
                    c * (N + 1) + Lbc.m) =
          vr.block(row_counter, 0, N + 1 + n - num_bc_each_var[i],
                   c * (N + 1) + Lbc.m);
      row_counter += N + 1 + n - num_bc_each_var[i];
      evecMat.block(0, 0, num_bc_each_var[i], c * (N + 1) + Lbc.m) =
          fir_n_mat.block(col_counter, 0, num_bc_each_var[i],
                          c * (N + 1) + Lbc.m);
      col_counter += num_bc_each_var[i];
      evecMat_master.block(i * (N + 1), 0, N + 1, c * (N + 1) + Lbc.m) =
          Mat.mats2[n] * evecMat;
    }
    std::vector<std::vector<T> > eigenval_trunc_sorter;

    eigenval_trunc_sorter.resize(c * (N + 1) + Lbc.m);
    for (int i = 0; i < c * (N + 1) + Lbc.m; i++) {
      eigenval_trunc_sorter[i].resize(2);
    }

    Eigen::Matrix<int, Eigen::Dynamic, 1> MPcount(c * (N + 1) + Lbc.m);
    MPcount.setConstant(0);
    std::vector<EigenSorter<T> > sorter;
    sorter.resize(c * (N + 1) + Lbc.m);
    for (int i = 0; i < c * (N + 1) + Lbc.m; i++) {
      sorter[i].compute(evecMat_master.block(0, i, c * (N + 1), 1), c);
      // sorter[i].compute(evecMat_master.block(0, i, (N + 1), 1));
      if (sorter[i].Mp == true) {
        MPcount[i] = 1;
      }
    }
    int MPcount_sum = MPcount.sum();
    converged = MPcount_sum;

    if (num_vals > masterL.rows()) {
      std::cout << "Only " << masterL.rows()
                << " eigenvalues can be calculated."
                << "Storing only that many." << '\n';
      num_vals = masterL.rows();
    }

    eigenvectors.resize(c,num_vals);
    //eigenvectorsMat.resize(num_vals);
    //for (int i = 0; i < num_vals; i++) {
    //  eigenvectorsMat[i].resize(c, 1);
  //  }
    eigenvalues.resize(num_vals);
    if (MPcount_sum >= num_vals) {
      for (int i = 0; i < c * (N + 1); i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
        // std::cout << eigs.eigenvalues()[i] << " " << sorter[i].numMp <<
        // "\n"; std::cin > > bre;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < num_vals; i++) {
        eigenvalues[i] = eigenvalues_temp(int(eigenval_trunc_sorter[i][1]), 0);
        alpha[i] = alpha_temp(int(eigenval_trunc_sorter[i][1]), 0);
        beta[i] = beta_temp(int(eigenval_trunc_sorter[i][1]), 0);
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];

        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < num_vals; i++) {
          eigenvectors(j,i) = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          eigenvectors(j,i).dct_flag = SIS_CHEB_SPACE;
          //eigenvectorsMat[i][j] = evecMat_master.block(
          //    j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          //eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }
    } else {
      std::cout << "Last " << num_vals - MPcount_sum
                << " eigenvectors are not resolved to machine precision."
                << '\n';
      for (int i = 0; i < c * (N + 1); i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < MPcount_sum; i++) {
        eigenvalues[i] = eigenvalues_temp(int(eigenval_trunc_sorter[i][1]), 0);
        alpha[i] = alpha_temp(int(eigenval_trunc_sorter[i][1]), 0);
        beta[i] = beta_temp(int(eigenval_trunc_sorter[i][1]), 0);
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];
        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < MPcount_sum; i++) {
          eigenvectors(j,i) = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          eigenvectors(j,i).dct_flag = SIS_CHEB_SPACE;
          //eigenvectorsMat[i][j] = evecMat_master.block(
          //    j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          //eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }

      std::vector<std::vector<T> > eigenval_trunc_sorter2;
      eigenval_trunc_sorter2.resize(c * (N + 1) + Lbc.m - MPcount_sum);
      for (int i = MPcount_sum; i < c * (N + 1) + Lbc.m; i++) {
        eigenval_trunc_sorter2[i - MPcount_sum].resize(2);
        eigenval_trunc_sorter2[i - MPcount_sum][0] =
            sorter[int(eigenval_trunc_sorter[i][1])].minVal;
        eigenval_trunc_sorter2[i - MPcount_sum][1] =
            eigenval_trunc_sorter[i][1];
      }
      std::sort(eigenval_trunc_sorter2.begin(), eigenval_trunc_sorter2.end());
      for (int i = MPcount_sum; i < num_vals; i++) {
        eigenvalues[i] = eigenvalues_temp(
            int(eigenval_trunc_sorter2[i - MPcount_sum][1]), 0);
        alpha[i] =
            alpha_temp(int(eigenval_trunc_sorter2[i - MPcount_sum][1]), 0);
        beta[i] = beta_temp(int(eigenval_trunc_sorter2[i - MPcount_sum][1]), 0);
      }
      for (int j = 0; j < c; j++) {
        int n = (highest_each_columnL[j] >= highest_each_columnM[j])
                    ? highest_each_columnL[j]
                    : highest_each_columnM[j];
        Mat.compute(n);
        temp_vec.resize(N + 1 + n, 1);
        for (int i = MPcount_sum; i < num_vals; i++) {
          eigenvectors(j,i) = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter2[i - MPcount_sum][1]),
              N + 1, 1);
          eigenvectors(j,i).dct_flag = SIS_CHEB_SPACE;
          //eigenvectorsMat[i][j] = evecMat_master.block(
          //    j * (N + 1), int(eigenval_trunc_sorter2[i - MPcount_sum][1]),
          //    N + 1, 1);
          //eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }
    }

#endif
  // Convert everything to physical space
  eigenvalues.c2p();
  }


  /// \brief This function will sort the eigenvalues and vectors by the largest real
  /// part.
  void sortByLargestReal() {

    int n = eigenvalues.size();
    std::vector<std::vector<T> > sorter;
    sorter.resize(n);
    for (int i = 0; i < n; i++) {
      sorter[i].resize(2);
    }

    for (int i = 0; i < n; i++) {
      sorter[i][0] = eigenvalues[i].real();
      sorter[i][1] = i;
    }

    // Sort in ascending order
    std::sort(sorter.begin(), sorter.end());
    // Then reverse it.
    std::reverse(sorter.begin(), sorter.end());

    ChebfunMat<std::complex<T> > tempMat(eigenvectors.r, eigenvectors.c);
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> tempvals;
    Eigen::Matrix<int, Eigen::Dynamic, 1> tempMP;
    tempvals = eigenvalues;
    tempMat = eigenvectors;
    tempMP = MPorNot;
#ifdef SIS_USE_LAPACK
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> tempAlpha, tempBeta;
    tempAlpha = alpha;
    tempBeta = beta;
#endif
    for (int i = 0; i < n; i++) {
      eigenvalues[i] = tempvals[sorter[i][1]];
      for (int j = 0; j < tempMat.r; j++)
        eigenvectors(j,i) = tempMat(j,sorter[i][1]);
      MPorNot[i] = tempMP[sorter[i][1]];

#ifdef SIS_USE_LAPACK
      alpha[i] = tempAlpha[sorter[i][1]];
      beta[i] = tempBeta[sorter[i][1]];
#endif
    }
  }

  /// \brief This will remove all unconverged  and infinite eigenvalues from the list.
  void keepConverged() {
    if (converged < eigenvalues.size()) {
      int c = eigenvectors.c;
      int r = eigenvectors.r;
      ChebfunMat<std::complex<T> > eigenvectors_temp;
      Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> eigenvalues_temp;
#ifdef SIS_USE_LAPACK
      Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> alpha_temp = alpha;
      Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> beta_temp = beta;
      alpha.resize(converged);
      beta.resize(converged);
#endif
      eigenvectors_temp = eigenvectors;
      eigenvalues_temp = eigenvalues;
      eigenvectors.resize(eigenvectors_temp.r,converged);

      eigenvalues.resize(converged);
      int count = 0;
      for (int i = 0; i < eigenvalues_temp.size(); i++) {
        if (MPorNot[i] == 1) {
          for (int j = 0; j < eigenvectors_temp.r; j++)
            eigenvectors(j,count) = eigenvectors_temp(j,i);
          eigenvalues[count] = eigenvalues_temp[i];
#ifdef SIS_USE_LAPACK
          alpha[count] = alpha_temp[i];
          beta[count] = beta_temp[i];
#endif
          count++;
        }
      }
      MPorNot.resize(converged);
      MPorNot.setConstant(1);
    }
  }

  /// \brief This will keep only n and remove the rest.
  /// Use this carefully. Only first n will be kept, irrespective of convergence.
  void keep(int n) {
    if (n < eigenvalues.size()) {
      int c = eigenvectors.c;
      int r = eigenvectors.r;
      ChebfunMat<std::complex<T> >  eigenvectors_temp;
      Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> eigenvalues_temp;
      eigenvectors_temp = eigenvectors;
      eigenvalues_temp = eigenvalues;
#ifdef SIS_USE_LAPACK
      Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> alpha_temp = alpha;
      Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> beta_temp = beta;
#endif
      eigenvectors.resize(eigenvectors_temp.r,n);
      //eigenvectorsMat.resize(n);
      //for (int i = 0; i < n; i++) {
      //  eigenvectorsMat[i].resize(c, 1);
      //}

      eigenvalues.resize(n);
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < eigenvectors_temp.r; j++)
          eigenvectors(j,i) = eigenvectors_temp(j,i);
        eigenvalues[i] = eigenvalues_temp[i];
#ifdef SIS_USE_LAPACK
        alpha[i] = alpha_temp[i];
        beta[i] = beta_temp[i];
#endif
      }
    }
  }

  /// \brief This will remove all infinite eigenvalues.
  /// Useful when solving generalized eigenvalue problems
  void removeInf() {
    std::vector<T> noninf_indices;
    for (int i = 0; i < eigenvalues.size(); i++) {
#ifndef SIS_USE_LAPACK
      if (std::abs(eigenvalues[i]) > 1.0e8 ||
          std::isnan(eigenvalues[i].real()) ||
          std::isnan(eigenvalues[i].imag())) {
      } else {
        noninf_indices.push_back(i);
      }
#else
      if (std::abs(beta[i]) < 1.0e-8) {
      } else {
        noninf_indices.push_back(i);
      }
#endif
    }

    int n = noninf_indices.size();
    int c = eigenvectors.c;
    int r = eigenvectors.r;
    ChebfunMat<std::complex<T> >  eigenvectors_temp;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> eigenvalues_temp;
    Eigen::Matrix<int, Eigen::Dynamic, 1> MPorNot_temp;
#ifdef SIS_USE_LAPACK
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> alpha_temp = alpha;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> beta_temp = beta;
#endif
    eigenvectors_temp = eigenvectors;
    eigenvalues_temp = eigenvalues;
    MPorNot_temp = MPorNot;
    eigenvectors.resize(eigenvectors_temp.r,n);
    //for (int i = 0; i < n; i++) {
    //  eigenvectors[i].resize(c, 1);
    //}

    eigenvalues.resize(n);
    MPorNot.resize(n);
#ifdef SIS_USE_LAPACK
    alpha.resize(n);
    beta.resize(n);
#endif
    int counter = 0;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < eigenvectors_temp.r; j++)
        eigenvectors(j,i) = eigenvectors_temp(j,noninf_indices[i]);
      eigenvalues[i] = eigenvalues_temp[noninf_indices[i]];
      MPorNot[i] = MPorNot_temp[noninf_indices[i]];
#ifdef SIS_USE_LAPACK
      alpha[i] = alpha_temp[noninf_indices[i]];
      beta[i] = beta_temp[noninf_indices[i]];
#endif
    }
  }
  /// \brief This is used to use a discretization to compute eigenvalues
  void compute(
      const Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> &L_,
      const Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> &M_,
      Discretize<std::complex<T> > Dis) {
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> masterM2 =
                                                                       M_,
                                                                   masterL2 =
                                                                       L_;
    int num_vals = masterM2.rows();
    int c = masterM2.rows() / (N + 1);
    int total_boundary_conditions = Dis.subs_mat.rows();
    MatGen<std::complex<T> > Mat;
    int row_counter = 0;
    int col_counter = 0;
#ifndef SIS_USE_LAPACK
    Eigen::ColPivHouseholderQR<
        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
        solver(masterM2);
    Eigen::ComplexEigenSolver<
        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
        eigs;
    bool Is_M_Invertible;
    if (solver.isInvertible()) {
      //  std::cout << "M is invertible" << '\n';
      //  std::cin >> bre;
      Is_M_Invertible = true;
      eigs.compute(solver.inverse() * masterL2);
      // std::cout << "Eigenvalues :\n" << eigs.eigenvalues() << '\n';
    } else {
      // std::cout << "M is not invertible." << '\n';
      // std::cin >> bre;
      Is_M_Invertible = false;
      solver.compute(masterL2);
      eigs.compute(solver.inverse() * masterM2);
      // std::cout << "Eigenvalues :\n" << 1 / eigs.eigenvalues().array() <<
      // '\n';
    }
    // std::cin >> bre;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
        evecMat_master(c * (N + 1), c * (N + 1));
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> temp_vec;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> fir_n_mat =
        Dis.subs_mat * eigs.eigenvectors();
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> evecMat(
        c * (N + 1) + total_boundary_conditions, c * (N + 1));
    evecMat << fir_n_mat, eigs.eigenvectors();

    // Unpermute eveMat:

    evecMat = Dis.P * evecMat;

    for (int i = 0; i < c; i++) {
      int n = Dis.highest_each_column[i];
      Mat.compute(n);
      evecMat_master.block(i * (N + 1), 0, N + 1, c * (N + 1)) =
          Mat.mats2[n] * evecMat.block(row_counter, 0, N + 1 + n, c * (N + 1));
      row_counter += N + 1 + n;
    }
    std::vector<std::vector<T> > eigenval_trunc_sorter;

    eigenval_trunc_sorter.resize(c * (N + 1));
    for (int i = 0; i < c * (N + 1); i++) {
      eigenval_trunc_sorter[i].resize(2);
    }

    Eigen::Matrix<int, Eigen::Dynamic, 1> MPcount(c * (N + 1));
    MPcount.setConstant(0);
    std::vector<EigenSorter<T> > sorter;
    sorter.resize(c * (N + 1));
    for (int i = 0; i < c * (N + 1); i++) {
      sorter[i].compute(evecMat_master.block(0, i, c * (N + 1), 1), c);
      if (!Is_M_Invertible) {
        // If M is not invertible, there will be zero eigenvalues, so that 1/0 =
        // inf eigenvalue, remove them.
        if (abs(eigs.eigenvalues()[i]) < 1e-11) {
          sorter[i].numMp = N + 1;
          sorter[i].Mp = false;
          sorter[i].minVal = 1e11;
        }
      }
      // sorter[i].compute(evecMat_master.block(0, i, (N + 1), 1));
      if (sorter[i].Mp == true) {
        MPcount[i] = 1;
      }
    }
    if (num_vals > masterL2.rows()) {
      std::cout << "Only " << masterL2.rows()
                << " eigenvalues can be calculated."
                << "Storing only that many." << '\n';
      num_vals = masterL2.rows();
    }

    int MPcount_sum = MPcount.sum();
    converged = MPcount_sum;
    eigenvectors.resize(c,num_vals);
    //eigenvectorsMat.resize(num_vals);
    //for (int i = 0; i < num_vals; i++) {
    //  eigenvectorsMat[i].resize(c, 1);
    //}
    eigenvalues.resize(num_vals);
    if (MPcount_sum >= num_vals) {
      for (int i = 0; i < c * (N + 1); i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
        // std::cout << eigs.eigenvalues()[i] << " " << sorter[i].numMp <<
        // "\n"; std::cin > > bre;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < num_vals; i++) {
        if (Is_M_Invertible) {
          eigenvalues[i] = eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
        } else {
          eigenvalues[i] =
              1.0 / eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
        }
      }
      for (int j = 0; j < c; j++) {
        int n = Dis.highest_each_column[j];

        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < num_vals; i++) {
          eigenvectors(j,i) = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          eigenvectors(j,i).dct_flag = SIS_CHEB_SPACE;
          //eigenvectorsMat[i][j] = evecMat_master.block(
          //    j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          //eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }
    } else {
      std::cout << "Last " << num_vals - MPcount_sum
                << " eigenvectors are not resolved to machine precision."
                << '\n';
      for (int i = 0; i < c * (N + 1); i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < MPcount_sum; i++) {
        if (Is_M_Invertible) {
          eigenvalues[i] = eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
        } else {
          eigenvalues[i] =
              1.0 / eigs.eigenvalues()[int(eigenval_trunc_sorter[i][1])];
        }
      }
      for (int j = 0; j < c; j++) {
        int n = Dis.highest_each_column[j];
        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < MPcount_sum; i++) {
          eigenvectors(j,i) = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          eigenvectors(j,i).dct_flag = SIS_CHEB_SPACE;
          //eigenvectorsMat[i][j] = evecMat_master.block(
          //    j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          //eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }

      std::vector<std::vector<T> > eigenval_trunc_sorter2;
      eigenval_trunc_sorter2.resize(c * (N + 1) - MPcount_sum);
      for (int i = MPcount_sum; i < c * (N + 1); i++) {
        eigenval_trunc_sorter2[i - MPcount_sum].resize(2);
        eigenval_trunc_sorter2[i - MPcount_sum][0] =
            sorter[int(eigenval_trunc_sorter[i][1])].minVal;
        eigenval_trunc_sorter2[i - MPcount_sum][1] =
            eigenval_trunc_sorter[i][1];
      }
      std::sort(eigenval_trunc_sorter2.begin(), eigenval_trunc_sorter2.end());
      for (int i = MPcount_sum; i < num_vals; i++) {
        if (Is_M_Invertible) {
          eigenvalues[i] = eigs.eigenvalues()[int(
              eigenval_trunc_sorter2[i - MPcount_sum][1])];
        } else {
          eigenvalues[i] =
              1.0 / eigs.eigenvalues()[int(
                        eigenval_trunc_sorter2[i - MPcount_sum][1])];
        }
      }
      for (int j = 0; j < c; j++) {
        int n = Dis.highest_each_column[j];
        Mat.compute(n);
        temp_vec.resize(N + 1 + n, 1);
        for (int i = MPcount_sum; i < num_vals; i++) {
          eigenvectors(j,i) = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter2[i - MPcount_sum][1]),
              N + 1, 1);
          eigenvectors(j,i).dct_flag = SIS_CHEB_SPACE;
        //  eigenvectorsMat[i][j] = evecMat_master.block(
        //      j * (N + 1), int(eigenval_trunc_sorter2[i - MPcount_sum][1]),
        //      N + 1, 1);
        //  eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }
    }
#else
    char jobvl = 'N';                 // Don't compute left evecs
    char jobvr = 'V';                 // Compute right evecs
    std::complex<double> wkopt;       // Eistimate optimum workspace
    std::complex<double> *work;       // allocate optimum workspace
    alpha.resize(masterL2.rows(), 1); // alpha for gen. eig. prob.
    beta.resize(masterL2.rows(), 1);  // beta for gen. eig. prob.

    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> vl(
        masterL2.rows(), masterL2.rows()),
        vr(masterL2.rows(), masterL2.rows()),
        eigenvalues_temp(masterL2.rows(), 1), alpha_temp(masterL2.rows(), 1),
        beta_temp(masterL2.rows(), 1);
    // vl : left evecs, vr: right evecs.
    int ldL = masterL2.outerStride(); // ld for leading dimension
    int ldM = masterM2.outerStride();
    int ldvl = vl.outerStride();
    int ldvr = vr.outerStride();
    int sizeL = masterL2.rows();
    int lwork = -1; // set lwork to -1 to estimate workspace.
    double rwork[8 * sizeL];

    // call this to estimate workspace
    zggev_(&jobvl, &jobvr, &sizeL, masterL2.data(), &ldL, masterM2.data(), &ldM,
           alpha_temp.data(), beta_temp.data(), vl.data(), &ldvl, vr.data(),
           &ldvr, &wkopt, &lwork, rwork, &info);

    // Now allocate workspace:
    lwork = (int)real(wkopt);
    work = (std::complex<double> *)malloc(lwork * sizeof(std::complex<double>));

    // Solve eigenvalue problem:
    zggev_(&jobvl, &jobvr, &sizeL, masterL2.data(), &ldL, masterM2.data(), &ldM,
           alpha_temp.data(), beta_temp.data(), vl.data(), &ldvl, vr.data(),
           &ldvr, work, &lwork, rwork, &info);

    // Free workspace.
    free((void *)work);

    eigenvalues_temp = alpha_temp.array() / beta_temp.array();

    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
        evecMat_master(c * (N + 1), c * (N + 1));
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> temp_vec;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> fir_n_mat =
        Dis.subs_mat * vr;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> evecMat(
        c * (N + 1) + total_boundary_conditions, c * (N + 1));
    evecMat << fir_n_mat, vr;
    // Unpermute eveMat:
    evecMat = Dis.P * evecMat;
    for (int i = 0; i < c; i++) {
      int n = Dis.highest_each_column[i];
      Mat.compute(n);
      evecMat_master.block(i * (N + 1), 0, N + 1, c * (N + 1)) =
          Mat.mats2[n] * evecMat.block(row_counter, 0, N + 1 + n, c * (N + 1));
      row_counter += N + 1 + n;
    }
    std::vector<std::vector<T> > eigenval_trunc_sorter;

    eigenval_trunc_sorter.resize(c * (N + 1));
    for (int i = 0; i < c * (N + 1); i++) {
      eigenval_trunc_sorter[i].resize(2);
    }

    Eigen::Matrix<int, Eigen::Dynamic, 1> MPcount(c * (N + 1));
    MPcount.setConstant(0);
    std::vector<EigenSorter<T> > sorter;
    sorter.resize(c * (N + 1));
    for (int i = 0; i < c * (N + 1); i++) {
      sorter[i].compute(evecMat_master.block(0, i, c * (N + 1), 1), c);
      // sorter[i].compute(evecMat_master.block(0, i, (N + 1), 1));
      if (sorter[i].Mp == true) {
        MPcount[i] = 1;
      }
    }
    if (num_vals > masterL2.rows()) {
      std::cout << "Only " << masterL2.rows()
                << " eigenvalues can be calculated."
                << "Storing only that many." << '\n';
      num_vals = masterL2.rows();
    }

    int MPcount_sum = MPcount.sum();
    converged = MPcount_sum;

    eigenvectors.resize(c,num_vals);
    //eigenvectorsMat.resize(num_vals);
    //for (int i = 0; i < num_vals; i++) {
    //  eigenvectorsMat[i].resize(c, 1);
    //}
    eigenvalues.resize(num_vals);
    beta.resize(num_vals);
    alpha.resize(num_vals);
    converged = MPcount_sum;
    if (MPcount_sum >= num_vals) {
      for (int i = 0; i < c * (N + 1); i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
        // std::cout << eigs.eigenvalues()[i] << " " << sorter[i].numMp <<
        // "\n"; std::cin > > bre;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < num_vals; i++) {
        eigenvalues[i] = eigenvalues_temp(int(eigenval_trunc_sorter[i][1]), 0);
        alpha[i] = alpha_temp(int(eigenval_trunc_sorter[i][1]), 0);
        beta[i] = beta_temp(int(eigenval_trunc_sorter[i][1]), 0);
      }
      for (int j = 0; j < c; j++) {
        int n = Dis.highest_each_column[j];

        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < num_vals; i++) {
          eigenvectors(j,i) = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          eigenvectors(j,i).dct_flag = SIS_CHEB_SPACE;
          //eigenvectorsMat[i][j] = evecMat_master.block(
          //    j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          //eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }
    } else {
      std::cout << "Last " << num_vals - MPcount_sum
                << " eigenvectors are not resolved to machine precision."
                << '\n';
      for (int i = 0; i < c * (N + 1); i++) {
        eigenval_trunc_sorter[i][0] = T(sorter[i].numMp);
        eigenval_trunc_sorter[i][1] = i;
      }
      std::sort(eigenval_trunc_sorter.begin(), eigenval_trunc_sorter.end());
      for (int i = 0; i < MPcount_sum; i++) {
        eigenvalues[i] = eigenvalues_temp(int(eigenval_trunc_sorter[i][1]), 0);
        alpha[i] = alpha_temp(int(eigenval_trunc_sorter[i][1]), 0);
        beta[i] = beta_temp(int(eigenval_trunc_sorter[i][1]), 0);
      }
      for (int j = 0; j < c; j++) {
        int n = Dis.highest_each_column[j];
        temp_vec.resize(N + 1 + n, 1);
        for (int i = 0; i < MPcount_sum; i++) {
          eigenvectors(j,i) = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          eigenvectors(j,i).dct_flag = SIS_CHEB_SPACE;
          //eigenvectorsMat[i][j] = evecMat_master.block(
          //    j * (N + 1), int(eigenval_trunc_sorter[i][1]), N + 1, 1);
          //eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }

      std::vector<std::vector<T> > eigenval_trunc_sorter2;
      eigenval_trunc_sorter2.resize(c * (N + 1) - MPcount_sum);
      for (int i = MPcount_sum; i < c * (N + 1); i++) {
        eigenval_trunc_sorter2[i - MPcount_sum].resize(2);
        eigenval_trunc_sorter2[i - MPcount_sum][0] =
            sorter[int(eigenval_trunc_sorter[i][1])].minVal;
        eigenval_trunc_sorter2[i - MPcount_sum][1] =
            eigenval_trunc_sorter[i][1];
      }
      std::sort(eigenval_trunc_sorter2.begin(), eigenval_trunc_sorter2.end());
      for (int i = MPcount_sum; i < num_vals; i++) {
        eigenvalues[i] = eigenvalues_temp(
            int(eigenval_trunc_sorter2[i - MPcount_sum][1]), 0);
        alpha[i] =
            alpha_temp(int(eigenval_trunc_sorter2[i - MPcount_sum][1]), 0);
        beta[i] = beta_temp(int(eigenval_trunc_sorter2[i - MPcount_sum][1]), 0);
      }
      for (int j = 0; j < c; j++) {
        int n = Dis.highest_each_column[j];
        Mat.compute(n);
        temp_vec.resize(N + 1 + n, 1);
        for (int i = MPcount_sum; i < num_vals; i++) {
          eigenvectors(j,i) = evecMat_master.block(
              j * (N + 1), int(eigenval_trunc_sorter2[i - MPcount_sum][1]),
              N + 1, 1);
          eigenvectors(j,i).dct_flag = SIS_CHEB_SPACE;
          //eigenvectorsMat[i][j] = evecMat_master.block(
          //    j * (N + 1), int(eigenval_trunc_sorter2[i - MPcount_sum][1]),
          //    N + 1, 1);
          //eigenvectorsMat[i][j].dct_flag = SIS_CHEB_SPACE;
        }
      }
    }

#endif
  //Convert everything to physical space
  eigenvectors.c2p();
  }
};

/// \brief This class computes the eigenvalues and eigenvectors (functions)
/// of a Linear operator Linop. See documentation of Eigen on how to access
/// eigenvalues and eigen vectors.
template <class T> class EigenSolver : public GeneralizedEigenSolver<T> {};

/// \brief This class computes the eigenvalues and eigenvectors (functions)
/// of a Linear operator Linop, overloaded to a complex Linop.
template <class T>
class EigenSolver<std::complex<T> >
    : public GeneralizedEigenSolver<std::complex<T> > {};

/// \brief This class sets up integration Matrices. This class must be
/// intiated by the highest order of based on which integration matrices will
/// be made.
template <class T> class MatGen {
private:
  int n;

public:
  std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > mats;
  std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > mats2;
  std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > con_mats;

  /// \brief Null Constructor
  MatGen() {}
  /// \brief Constructor
  MatGen(int n_) { compute(n_); }

  /// \brief Call this to generate integration matrices for highest order n.
  void compute(int n_) {
    clear();
    n = n_;

    mats.resize(n + 1);
    for (int i = 0; i < n + 1; i++) {
      mats[i].resize(N + 1 + 2 * n, N + 1 + 2 * n);
      mats[i].setConstant(0.0);
    }
    mats[0].block(0, 0, N + 1, N + 1) =
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(N + 1,
                                                                   N + 1);
    for (int i = 1; i < n + 1; i++) {
      mats[i].row(0) = 0.5 * mats[i - 1].row(1);
      for (int j = 1; j < N + 2 * n; j++) {
        mats[i].row(j) =
            (mats[i - 1].row(j - 1) - mats[i - 1].row(j + 1)) / (2.0 * j);
      }
    }

    if (n > 0) {
      con_mats.resize(n);
      for (int i = 0; i < n; i++) {
        con_mats[i].resize(n, n);
        con_mats[i].setConstant(0.0);
      }
      // Next we set the matrices for the constants
      con_mats[n - 1].setIdentity(n, n);
      for (int i = n - 2; i > -1; i--) {
        for (int j = 0; j < n; j++) {
          con_mats[i].col(j) = diff<T>(con_mats[i + 1].col(j));
        }
      }
    }
    mats2.resize(n + 1);
    for (int i = 0; i < n + 1; i++) {
      mats2[i].resize(N + 1, N + 1 + n);
      mats2[i].setConstant(0.0);
    }
    for (int i = 0; i < n + 1; i++) {
      mats2[i].block(0, 0, N + 1, N + 1) = mats[i].block(0, 0, N + 1, N + 1);
    }
    for (int i = 0; i < n; i++) {
      mats2[i + 1].block(0, N + 1, n, n) = con_mats[i];
    }
  }
  void clear() {
    mats.clear();
    mats2.clear();
    con_mats.clear();
  }
};

/// \brief This class sets up integration Matrices. This class must be
/// intiated by the highest order of based on which integration matrices will
/// be made.
template <class T> class MatGen<std::complex<T> > {
private:
  int n;

public:
  std::vector<Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
      mats;
  std::vector<Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
      mats2;
  std::vector<Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
      con_mats;

  /// \brief Null Constructor
  MatGen() {}
  /// \brief Constructor
  MatGen(int n_) { compute(n_); }

  /// \brief Call this to generate integration matrices for highest order n.
  void compute(int n_) {
    clear();
    n = n_;

    mats.resize(n + 1);
    for (int i = 0; i < n + 1; i++) {
      mats[i].resize(N + 1 + 2 * n, N + 1 + 2 * n);
      mats[i].setConstant(std::complex<T>(0.0, 0.0));
    }
    mats[0].block(0, 0, N + 1, N + 1) =
        Eigen::Matrix<std::complex<T>, Eigen::Dynamic,
                      Eigen::Dynamic>::Identity(N + 1, N + 1);
    for (int i = 1; i < n + 1; i++) {
      mats[i].row(0) = 0.5 * mats[i - 1].row(1);
      for (int j = 1; j < N + 2 * n; j++) {
        mats[i].row(j) =
            (mats[i - 1].row(j - 1) - mats[i - 1].row(j + 1)) / (2.0 * j);
      }
    }
    con_mats.resize(n);
    for (int i = 0; i < n; i++) {
      con_mats[i].resize(n, n);
      con_mats[i].setConstant(std::complex<T>(0.0, 0.0));
    }
    // Next we set the matrices for the constants
    con_mats[n - 1].setIdentity(n, n);
    for (int i = n - 2; i > -1; i--) {
      for (int j = 0; j < n; j++) {
        con_mats[i].col(j) = diff<std::complex<T> >(con_mats[i + 1].col(j));
      }
    }
    mats2.resize(n + 1);
    for (int i = 0; i < n + 1; i++) {
      mats2[i].resize(N + 1, N + 1 + n);
      mats2[i].setConstant(std::complex<T>(0.0, 0.0));
    }
    for (int i = 0; i < n + 1; i++) {
      mats2[i].block(0, 0, N + 1, N + 1) = mats[i].block(0, 0, N + 1, N + 1);
    }
    for (int i = 0; i < n; i++) {
      mats2[i + 1].block(0, N + 1, n, n) = con_mats[i];
    }
  }
  void clear() {
    mats.clear();
    mats2.clear();
    con_mats.clear();
  }
};

/// \brief This class represents a block matrix operator. It is a matrix of
/// operators.
template <class T> class LinopMat {
  friend class EigenSolver<T>;
  friend class GeneralizedEigenSolver<T>;
  friend class GeneralizedEigenSolver<std::complex<T> >;
  friend class LinopMat<std::complex<T> >;

private:
  int count; ///< Rows, columns and counter of number of inputs.

public:
  /// \brief Number of rows
  int r;
  /// \brief Number of columns
  int c;
  /// \brief This 1D vector of Linops holds all Linops of the LinopMat in
  /// the row-major format.
  std::valarray<Linop<T> > LinopVec;
  int countc;
  int countr;
  bool go_to_next_row;
  /// \brief This holds a vector of boundary conditions. There should be as
  /// many boundary conditions as there are columns, for example consider a
  /// (2,3) LinopMat: \f[ \left[ \begin{array}{ccc}
  /// L & M & O\\
  /// Q & R & S
  /// \end{array}
  /// \right]\left[
  /// \begin{array}{c}
  /// u\\
  /// v\\
  /// w
  /// \end{array}
  /// \right]
  /// \f]
  /// One needs to specify the boundary conditions for three dependent
  /// variables, i.e., \f$u\f$, \f$v\f$ and \f$w\f$

  /// \brief Null constructor
  LinopMat() : LinopVec(), r(0), c(0), count(0) {}

  /// \brief Initializes a Linop Matrix of r_ rows and c_ columns.
  LinopMat(int r_, int c_) {
    r = r_;
    c = c_;
    count = 0;
    LinopVec.resize(r * c);
  };

  /// \brief Copy constructor
  LinopMat(const LinopMat<T> &in) {
    r = in.r;
    c = in.c;
    resize(r, c);
    count = in.count;
    for (int i = 0; i < r * c; i++) {
      LinopVec[i] = in.LinopVec[i];
    }
  }
  /// \brief Assignment operator
  void operator=(const LinopMat &in) {
    //  std::cout << "in copier" << '\n' << std::flush;
    r = in.r;
    c = in.c;
    resize(r, c);
    count = in.count;
    for (int i = 0; i < r * c; i++) {
      LinopVec[i] = in.LinopVec[i];
    }
  }
  /// \brief Use this to input multiple Linops to the LinopMat using comma
  /// separators. Input type of constant.
  LinopMat<T> &operator<<(T b) {
    int bre;
    LinopVec[count] = b;
    count++;
    return *this;
  };

  /// \brief Use this to to refer to Linops. For 1D Matrices, we use the Row
  /// major format when refering to [i]. See wiki for row-major format.
  Linop<T> &operator[](int i) {
    Linop<T> &temp = LinopVec[i];
    return temp;
  };

  /// \brief Use this to to refer to Linop. For 2D Matrices, refer to matrix
  /// element Linop by (i,j) with the row and column index. Indices start from
  /// 0.
  Linop<T> &operator()(int i, int j) { return LinopVec[j + c * i]; };

  /// \brief Use this to input multiple Linops to the LinopMat using comma
  /// separators. Input type Eigen array.
  LinopMat<T> &operator<<(Eigen::Array<T, Eigen::Dynamic, 1> b) {
    int bre;
    LinopVec[count] = b;
    count++;
    return *this;
  };

  /// \brief Use this to input multiple Linops to the LinopMat using comma
  LinopMat<T> &operator<<(Linop<T> b) {
    int bre;
    LinopVec[count] = b;
    count++;
    return *this;
  };

  /// \brief This clears all contents in the LinopMat, and then creates a
  /// fresh LinopMat of size r_ x c_.
  void resize(int r_, int c_) {
    r = r_;
    c = c_;
    LinopVec.resize(r * c);
    count = 0;
  }

  /// \brief Overloads comma separator to input Chebfuns into a ChebfunMat.
  /// Input type constant.
  LinopMat<T> &operator,(T b) {
    LinopVec[count] = b;
    count++;
    return *this;
  }
  /// \brief Overloads comma separator to input Chebfuns into a ChebfunMat.
  /// Input type Eigen array.
  LinopMat<T> &operator,(Eigen::Array<T, Eigen::Dynamic, 1> b) {
    LinopVec[count] = b;
    count++;
    return *this;
  }
  /// \brief Overloads comma separator to input Chebfuns into a ChebfunMat.
  /// Input type Chebfun.
  LinopMat<T> &operator,(Linop<T> b) {
    LinopVec[count] = b;
    count++;
    return *this;
  }

  /// \brief Solves an input ChebfunMat, Dimensions of ChebfunMat have to be
  /// (c,1)
  ChebfunMat<T> solve(const ChebfunMat<T> &in) {}

  /// \brief Use this to input multiple LinopMats to the LinopMat using comma
  /// separators. Call at after initializing a LinopMat, or after calling
  /// resize, , similar to Eigen's << operator. CAREFUL: If dimensions are not
  /// compatible, it will not throw an error, and could end up in segv.
  LinopMat<T> &operator<<(LinopMat<T> b) {
    setConstant(0.0);
    go_to_next_row = false;
    int bre;
    if (b.c > c || b.r > r) {
      std::cout << "Incompatible shape to input LinopMat."
                << " In" << __LINE__ << '\n';
      exit(1);
    }
    for (int i = 0; i < b.r; i++) {
      for (int j = 0; j < b.c; j++) {
        operator()(i, j) = b(i, j);
      }
    }
    countr = 0;
    countc = b.c;
    if (countc == c) {
      go_to_next_row = true;
      countr += b.r;
    }
    // if (ind == 2) {
    //  std::cout << "countc: " << countc << '\n';
    //  std::cout << "countr: " << countr << '\n';
    //  std::cout << "count: " << count << '\n';
    //  std::cout << "go_to_next_row: " << go_to_next_row << '\n';
    //  std::cin >> bre;
    //}
    count++;
    return *this;
  };

  LinopMat<T> &operator,(LinopMat<T> b) {
    int bre;
    if (!go_to_next_row) {
      // std::cout << "in if" << '\n';
      for (int i = 0; i < b.r; i++) {
        for (int j = 0; j < b.c; j++) {
          operator()(countr + i, countc + j) = b(i, j);
          //  std::cout << "(" << countr + i<< "," << countc + j << ")" << '\n';
        }
      }
      countc += b.c;
      if (countc == c) {
        go_to_next_row = true;
        countr += b.r;
        countc = 0;
      }
    } else {
      // std::cout << "in else" << '\n';
      go_to_next_row = false;
      for (int i = 0; i < b.r; i++) {
        for (int j = 0; j < b.c; j++) {
          operator()(countr + i, countc + j) = b(i, j);
          //  std::cout << "(" << countr + i<< "," << countc + j << ")" << '\n';
        }
      }
      countc += b.c;
      if (countc == c) {
        go_to_next_row = true;
        countr += b.r;
        countc = 0;
      }
    }
    //  if (ind == 2) {
    //    std::cout << "countc: " << countc << '\n';
    //    std::cout << "countr: " << countr << '\n';
    //    std::cout << "count: " << count << '\n';
    //    std::cout << "go_to_next_row: " << go_to_next_row << '\n';
    //    std::cin >> bre;
    //  }
    count++;
    return *this;
  };

  LinopMat<T> cTranspose() {
    LinopMat<T> temp;
    temp.resize(r, c);
    for (int i = 0; i < r; i++) {
      for (int j = 0; j < c; j++) {
        temp(i, j) = conj(operator()(j, i));
      }
    }
    return temp;
  }

  /// \brief Sets every member in LinopMat to constant.
  void setConstant(T in) {
    for (int i = 0; i < r; i++) {
      for (int j = 0; j < c; j++) {
        operator()(i, j) = in;
      }
    }
  };

  /// \brief Sets LinopMat to identity.
  void setIdentity() {
    if (r != c) {
      std::cout << "setIdentity() is only for square LinopMats" << '\n';
      exit(1);
    }
    setConstant(0.0);
    for (int i = 0; i < r; i++) {
      operator()(i, i) = 1.0;
    }
  };

  //  LinopMat<T> operator-(LinopMat<T> in) {
  //    LinopMat<T> temp(in.r, in.c);
  //    for (int i = 0; i < in.r; i++) {
  //      for (int j = 0; j < in.c; j++) {
  //        temp(i, j) = -in(i, j);
  //      }
  //    }
  //   return temp;//
  //  }
  //~LinopMat(){
  //  LinopVec.~valarray();
  //}
};

/// \brief This class represents a block matrix operator. It is a matrix of
/// operators.
template <class T> class LinopMat<std::complex<T> > {
  friend class EigenSolver<std::complex<T> >;
  friend class GeneralizedEigenSolver<std::complex<T> >;
  friend class BcMat<std::complex<T> >;

private:
  int count;

public:
  /// \brief Number of rows
  int r;
  /// \brief Number of columns
  int c;
  /// \brief This 1D vector of Linops holds all Linops of the LinopMat in
  /// the row-major format.
  std::valarray<Linop<std::complex<T> > > LinopVec;

  int countc;
  int countr;
  bool go_to_next_row;

  /// \brief Null constructor
  LinopMat() : LinopVec(), r(0), c(0), count(0) {}

  /// \brief Initializes a Linop Matrix of r_ rows and c_ columns.
  LinopMat(int r_, int c_) {
    r = r_;
    c = c_;
    count = 0;
    LinopVec.resize(r * c);
  };

  /// \brief Copy constructor
  LinopMat(const LinopMat<std::complex<T> > &in) {
    r = in.r;
    c = in.c;
    resize(r, c);
    count = in.count;
    for (int i = 0; i < r * c; i++) {
      LinopVec[i] = in.LinopVec[i];
    }
  }

  /// \brief Copy constructor
  LinopMat(const LinopMat<T> &in) {
    r = in.r;
    c = in.c;
    resize(r, c);
    count = in.count;
    for (int i = 0; i < r * c; i++) {
      LinopVec[i] = in.LinopVec[i];
    }
  }
  /// \brief Assignment operator
  void operator=(const LinopMat<std::complex<T> > &in) {
    r = in.r;
    c = in.c;
    resize(r, c);
    count = in.count;
    for (int i = 0; i < r * c; i++) {
      LinopVec[i] = in.LinopVec[i];
    }
  }

  /// \brief Assignment operator
  void operator=(const LinopMat<T> &in) {
    r = in.r;
    c = in.c;
    resize(r, c);
    count = in.count;
    for (int i = 0; i < r * c; i++) {
      LinopVec[i] = in.LinopVec[i];
    }
  }
  /// \brief Use this to input multiple Linops to the LinopMat using comma
  /// separators. Input type of constant.
  LinopMat<std::complex<T> > &operator<<(std::complex<T> b) {
    int bre;
    LinopVec[count] = b;
    count++;
    return *this;
  };
  /// \brief Use this to input multiple Linops to the LinopMat using comma
  /// separators. Input type of constant.
  LinopMat<std::complex<T> > &operator<<(T b) {
    int bre;
    LinopVec[count] = b;
    count++;
    return *this;
  };

  /// \brief Use this to to refer to Linops. For 1D Matrices, we use the Row
  /// major format when refering to [i]. See wiki for row-major format.
  Linop<std::complex<T> > &operator[](int i) { return LinopVec[i]; };

  /// \brief Use this to to refer to Linop. For 2D Matrices, refer to matrix
  /// element Linop by (i,j) with the row and column index. Indices start from
  /// 0.
  Linop<std::complex<T> > &operator()(int i, int j) {
    // std::cout << ind++ << '\n';
    return LinopVec[j + c * i];
  };

  /// \brief Use this to apply a LinopMat on a ChebfunMat.
  ChebfunMat<std::complex<T> > operator()(ChebfunMat<std::complex<T> > in) {
    ChebfunMat<std::complex<T> > out(r, in.c);
    if (c != in.r){
      std::cout << "Incompatible shapes in applying LinopMat to ChebfunMat."
                << "In line "<< __LINE__ <<". Exiting "<< '\n';
      exit(1);
    }
    for (int i = 0 ; i < out.r; i++){
      for (int j = 0; j < out.c; j++){
        out(i,j).v = std::complex<T>(0.0,0.0);
      }
    }
    for (int i = 0; i < r; i++) {
      for (int j = 0; j < in.c; j++) {
        for (int k = 0; k < c; k++) {
          out(i, j) = out(i, j) + (operator()(i, k) (in(k, j)));
        }
      }
    }
    return out;
  };

  /// \brief Use this to input multiple Linops to the LinopMat using comma
  /// separators. Input type Eigen array.
  LinopMat<std::complex<T> > &
  operator<<(Eigen::Array<std::complex<T>, Eigen::Dynamic, 1> b) {
    int bre;
    LinopVec[count] = b;
    count++;
    return *this;
  };

  /// \brief Use this to input multiple Linops to the LinopMat using comma
  /// separators. Input type valarray.
  LinopMat<std::complex<T> > &operator<<(std::valarray<std::complex<T> > b) {
    int bre;
    LinopVec[count] = b;
    count++;
    return *this;
  };
  /// \brief Use this to input multiple Linops to the LinopMat using comma
  /// separators. Input type valarray.
  LinopMat<std::complex<T> > &operator<<(std::valarray<T> b) {
    int bre;
    LinopVec[count] = b;
    count++;
    return *this;
  };

  /// \brief Use this to input multiple Linops to the LinopMat using comma
  LinopMat<std::complex<T> > &operator<<(Linop<std::complex<T> > b) {
    int bre;
    LinopVec[count] = b;
    count++;
    return *this;
  };

  /// \brief Use this to input multiple Linops to the LinopMat using comma
  LinopMat<std::complex<T> > &operator<<(Linop<T> b) {
    int bre;
    LinopVec[count] = b;
    count++;
    return *this;
  };

  /// \brief This clears all contents in the LinopMat, and then creates a
  /// fresh LinopMat of size r_ x c_.
  void resize(int r_, int c_) {
    r = r_;
    c = c_;
    LinopVec.resize(r * c);
    count = 0;
  }

  /// \brief Overloads comma separator to input Chebfuns into a ChebfunMat.
  /// Input type constant.
  LinopMat<std::complex<T> > &operator,(std::complex<T> b) {
    LinopVec[count] = b;
    count++;
    return *this;
  }

  /// \brief Overloads comma separator to input Chebfuns into a ChebfunMat.
  /// Input type constant.
  LinopMat<std::complex<T> > &operator,(T b) {
    LinopVec[count] = b;
    count++;
    return *this;
  }

  /// \brief Overloads comma separator to input Chebfuns into a ChebfunMat.
  /// Input type Eigen array.
  LinopMat<std::complex<T> > &operator,(
      Eigen::Array<std::complex<T>, Eigen::Dynamic, 1> b) {
    LinopVec[count] = b;
    count++;
    return *this;
  }

  /// \brief Overloads comma separator to input Chebfuns into a ChebfunMat.
  /// Input type valarray.
  LinopMat<std::complex<T> > &operator,(std::valarray<std::complex<T> > b) {
    LinopVec[count] = b;
    count++;
    return *this;
  }
  /// \brief Overloads comma separator to input Chebfuns into a ChebfunMat.
  /// Input type valarray.
  LinopMat<std::complex<T> > &operator,(std::valarray<T> b) {
    LinopVec[count] = b;
    count++;
    return *this;
  }

  /// \brief Overloads comma separator to input Chebfuns into a ChebfunMat.
  /// Input type Chebfun.
  LinopMat<std::complex<T> > &operator,(Linop<std::complex<T> > b) {
    LinopVec[count] = b;
    count++;
    return *this;
  }

  /// \brief Overloads comma separator to input Chebfuns into a ChebfunMat.
  /// Input type Chebfun.
  LinopMat<std::complex<T> > &operator,(Linop<T> b) {
    LinopVec[count] = b;
    count++;
    return *this;
  }

  /// \brief Sets every member in LinopMat to constant.
  void setConstant(T in) {
    for (int i = 0; i < r; i++) {
      for (int j = 0; j < c; j++) {
        operator()(i, j) = in;
      }
    }
  };

  /// \brief Sets every member in LinopMat to complex constant.
  void setConstant(std::complex<T> in) {
    for (int i = 0; i < r; i++) {
      for (int j = 0; j < c; j++) {
        operator()(i, j) = in;
      }
    }
  };

  /// \brief Sets LinopMat to identity.
  void setIdentity() {
    if (r != c) {
      std::cout << "setIdentity() is only for square LinopMats" << '\n';
      exit(1);
    }
    setConstant(0.0);
    for (int i = 0; i < r; i++) {
      operator()(i, i) = 1.0;
    }
  };
  /// \brief Use this to input multiple LinopMats to the LinopMat using comma
  /// separators. Call at after initializing a LinopMat, or after calling
  /// resize, , similar to Eigen's << operator. CAREFUL: If dimensions are not
  /// compatible, it will not throw an error, and could end up in segv.
  LinopMat<std::complex<T> > &operator<<(LinopMat<std::complex<T> > b) {
    setConstant(0.0);
    go_to_next_row = false;
    int bre;
    if (b.c > c || b.r > r) {
      std::cout << "Incompatible shape to input LinopMat."
                << " In" << __LINE__ << '\n';
      exit(1);
    }
    for (int i = 0; i < b.r; i++) {
      for (int j = 0; j < b.c; j++) {
        operator()(i, j) = b(i, j);
      }
    }
    countr = 0;
    countc = b.c;
    if (countc == c) {
      go_to_next_row = true;
      countr += b.r;
      countc = 0;
    }
    // if (ind == 2) {
    //  std::cout << "countc: " << countc << '\n';
    //  std::cout << "countr: " << countr << '\n';
    //  std::cout << "count: " << count << '\n';
    //  std::cout << "go_to_next_row: " << go_to_next_row << '\n';
    //  std::cin >> bre;
    //}
    count++;
    return *this;
  };

  LinopMat<std::complex<T> > &operator,(LinopMat<std::complex<T> > b) {
    int bre;
    if (!go_to_next_row) {
      // std::cout << "in if" << '\n';
      for (int i = 0; i < b.r; i++) {
        for (int j = 0; j < b.c; j++) {
          operator()(countr + i, countc + j) = b(i, j);
          //  std::cout << "(" << countr + i<< "," << countc + j << ")" << '\n';
        }
      }
      countc += b.c;
      if (countc == c) {
        go_to_next_row = true;
        countr += b.r;
        countc = 0;
      }
    } else {
      // std::cout << "in else" << '\n';
      go_to_next_row = false;
      for (int i = 0; i < b.r; i++) {
        for (int j = 0; j < b.c; j++) {
          operator()(countr + i, countc + j) = b(i, j);
          //    std::cout << "(" << countr + i<< "," << countc + j << ")" <<
          //    '\n';
        }
      }
      countc += b.c;
      if (countc == c) {
        go_to_next_row = true;
        countr += b.r;
        countc = 0;
      }
    }
    // if (ind == 2) {
    //  std::cout << "countc: " << countc << '\n';
    //  std::cout << "countr: " << countr << '\n';
    //  std::cout << "count: " << count << '\n';
    //  std::cout << "go_to_next_row: " << go_to_next_row << '\n';
    //  std::cin >> bre;
    //}
    count++;
    return *this;
  };

  LinopMat<std::complex<T> > &operator<<(LinopMat<T> b) {
    setConstant(0.0);
    go_to_next_row = false;
    int bre;
    if (b.c > c || b.r > r) {
      std::cout << "Incompatible shape to input LinopMat."
                << " In" << __LINE__ << '\n';
      exit(1);
    }
    for (int i = 0; i < b.r; i++) {
      for (int j = 0; j < b.c; j++) {
        operator()(i, j) = b(i, j);
      }
    }
    countr = 0;
    countc = b.c;
    if (countc == c) {
      go_to_next_row = true;
      countr += b.r;
    }
    // if (ind == 2) {
    //  std::cout << "countc: " << countc << '\n';
    //  std::cout << "countr: " << countr << '\n';
    //  std::cout << "count: " << count << '\n';
    //  std::cout << "go_to_next_row: " << go_to_next_row << '\n';
    //  std::cin >> bre;
    //}
    count++;
    return *this;
  };

  LinopMat<std::complex<T> > &operator,(LinopMat<T> b) {
    int bre;
    if (!go_to_next_row) {
      // std::cout << "in if" << '\n';
      for (int i = 0; i < b.r; i++) {
        for (int j = 0; j < b.c; j++) {
          operator()(countr + i, countc + j) = b(i, j);
          //  std::cout << "(" << countr + i<< "," << countc + j << ")" << '\n';
        }
      }
      countc += b.c;
      if (countc == c) {
        go_to_next_row = true;
        countr += b.r;
        countc = 0;
      }
    } else {
      // std::cout << "in else" << '\n';
      go_to_next_row = false;
      for (int i = 0; i < b.r; i++) {
        for (int j = 0; j < b.c; j++) {
          operator()(countr + i, countc + j) = b(i, j);
          //  std::cout << "(" << countr + i<< "," << countc + j << ")" << '\n';
        }
      }
      countc += b.c;
      if (countc == c) {
        go_to_next_row = true;
        countr += b.r;
        countc = 0;
      }
    }
    //  if (ind == 2) {
    //    std::cout << "countc: " << countc << '\n';
    //    std::cout << "countr: " << countr << '\n';
    //    std::cout << "count: " << count << '\n';
    //    std::cout << "go_to_next_row: " << go_to_next_row << '\n';
    //    std::cin >> bre;
    //  }
    count++;
    return *this;
  };

  LinopMat<std::complex<T> > cTranspose() {
    LinopMat<std::complex<T> > temp;
    temp.resize(r, c);
    for (int i = 0; i < r; i++) {
      for (int j = 0; j < c; j++) {
        temp(i, j) = conj(operator()(j, i));
      }
    }
    return temp;
  }

  //  LinopMat<std::complex<T> > operator-(LinopMat<std::complex<T> > in) {
  //    LinopMat<std::complex<T> > temp(in.r, in.c);
  //    for (int i = 0; i < in.r; i++) {
  //      for (int j = 0; j < in.c; j++) {
  //        temp(i, j) = -in(i, j);
  //      }
  //    }
  //    return temp;
  //  }
  //~LinopMat(){
  //  LinopVec.~valarray();
  //}
};

/// \brief BcMat will hold general Boundary conditions as LinopMats at
/// evealuation points, as given by operator L and evaluation points, eval.
template <class T> class BcMat {
  friend class EigenSolver<T>;
  friend class GeneralizedEigenSolver<T>;
  friend class GeneralizedEigenSolver<std::complex<T> >;
  friend class LinopMat<std::complex<T> >;
  friend class SingularValueDecomposition<T>;
  friend class SingularValueDecomposition<std::complex<T> >;
  friend class Discretize<T>;
  friend class Discretize<std::complex<T> >;

private:
  int m, n;

public:
  LinopMat<T> L;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> eval;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> vals;

  BcMat() {
    m = 0;
    n = 0;
  }
  BcMat(int m_, int n_) {
    m = m_;
    n = n_;
    eval.resize(m_, n_);
    vals.resize(m_, 1);
    L.resize(m, n);
    eval.setConstant(0.0);
    vals.setConstant(0.0);
  }
  void resize(int m_, int n_) {
    m = m_;
    n = n_;
    eval.resize(m_, n_);
    vals.resize(m_, 1);
    eval.setConstant(0.0);
    vals.setConstant(0.0);
    L.resize(m, n);
  }
  void operator=(const BcMat<T> &in) {
    m = in.m;
    n = in.n;
    L = in.L;
    eval = in.eval;
    vals = in.vals;
  }

  int rows() {
      return m;
    }

    int cols() {
      return n;
    }

  /// \brief Calling BcMat(i,j, ord) will produce row vector representing
  /// [------][a0 a1 ... an C0 C1]^T, assuming ord = 2. i and j refers to an
  /// element in the LinopMat L. Suppose L(1,1) is Dyy, and eval(1,1) = -1, then
  /// the row matrix will represent the evaluation of that condition in the
  /// basis of a0 to C1.
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> operator()(int i, int j,
                                                              int ord) {
    int n = ord;
    Eigen::Matrix<T, 1, Eigen::Dynamic> vals(1, N + 1 + 2 * n), valcons(n);
    for (int k = 0; k < N + 1 + 2 * n; k++) {
      vals[k] = cos(k * acos(eval(i, j)));
    }
    for (int k = 0; k < n; k++) {
      valcons[k] = cos(k * acos(eval(i, j)));
    }
    vals[0] *= 0.5;
    valcons[0] *= 0.5;
    MatGen<T> Mat;
    Mat.compute(n);
    Eigen::Matrix<T, 1, Eigen::Dynamic> out(1, N + 1 + n);
    out.setConstant(0.0);
    int diffn = n - L(i, j).n;
    for (int k = 0; k < L(i, j).n + 1; k++) {
      if (L(i, j).NCC == 0) {

        out.block(0, 0, 1, N + 1) +=
            (L(i, j).coef[k] * vals * Mat.mats[k + diffn]).head(N + 1);
        if (k + diffn - 1 > -1) {
          out.block(0, N + 1, 1, n) +=
              (L(i, j).coef[k] * valcons * Mat.con_mats[k + diffn - 1]);
        }
      } else {
        out.block(0, 0, 1, N + 1) +=
            (L(i, j).coefFun[k](eval(i, j)) * vals * Mat.mats[k + diffn])
                .head(N + 1);
        if (k + diffn - 1 > -1) {
          out.block(0, N + 1, 1, n) += (L(i, j).coefFun[k](eval(i, j)) *
                                        valcons * Mat.con_mats[k + diffn - 1]);
        }
      }
    }
    return out;
  }
};

template <class T> class BcMat<std::complex<T> > {
  friend class EigenSolver<T>;
  friend class GeneralizedEigenSolver<T>;
  friend class GeneralizedEigenSolver<std::complex<T> >;
  friend class LinopMat<std::complex<T> >;
  friend class SingularValueDecomposition<T>;
  friend class SingularValueDecomposition<std::complex<T> >;
  friend class Discretize<T>;
  friend class Discretize<std::complex<T> >;

private:
  int m, n;

public:
  LinopMat<std::complex<T> > L;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> eval;
  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> vals;
  BcMat() {
    m = 0;
    n = 0;
  }
  BcMat(int m_, int n_) {
    m = m_;
    n = n_;
    eval.resize(m_, n_);
    vals.resize(m_, 1);
    L.resize(m, n);
  //  eval.setConstant(0.0);
  //  vals.setConstant(0.0);
  }
  void resize(int m_, int n_) {
    m = m_;
    n = n_;
    eval.resize(m_, n_);
    vals.resize(m_, 1);
//    eval.setConstant(0.0);
//    vals.setConstant(0.0);
    L.resize(m, n);
  }
  void operator=(const BcMat<std::complex<T> > &in) {
    m = in.m;
    n = in.n;
    L = in.L;
    eval = in.eval;
    vals = in.vals;
  }

  int rows() {
    return m;
  }

  int cols() {
    return n;
  }

  /// \brief Calling BcMat(i,j, ord) will produce row vector representing
  /// [------][a0 a1 ... an C0 C1]^T, assuming ord = 2. i and j refers to an
  /// element in the LinopMat L. Suppose L(1,1) is Dyy, and eval(1,1) = -1, then
  /// the row matrix will represent the evaluation of that condition in the
  /// basis of a0 to C1.
  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
  operator()(int i, int j, int ord) {
    int n = ord;
    int bre;
    Eigen::Matrix<T, 1, Eigen::Dynamic> vals(1, N + 1 + 2 * n), valcons(n);
    for (int k = 0; k < N + 1 + 2 * n; k++) {
      vals[k] = cos(T(k) * acos(eval(i, j)));
    }
    for (int k = 0; k < n; k++) {
      valcons[k] = cos(T(k) * acos(eval(i, j)));
    }

    vals[0] *= 0.5;
    if (n > 0) {
      valcons[0] *= 0.5;
    }
    MatGen<T> Mat;
    Mat.compute(n);
    Eigen::Matrix<std::complex<T>, 1, Eigen::Dynamic> out(1, N + 1 + n);
    out.setConstant(0.0);
    int diffn = n - L(i, j).n;
    for (int k = 0; k < L(i, j).n + 1; k++) {
      if (L(i, j).NCC == 0) {
        // std::cout << "Im here;;;" << '\n';
        // if (ind == 2){
        // std::cout << "L("<<i<<", "<<j<<").coef["<<k<<"]:"<<L(i, j).coef[k] <<
        // '\n'; std::cout << "size of vals: " << vals.size() << '\n'; std::cout
        // << "Mat.mats["<< k + diffn << "]: "<< size(Mat.mats[k + diffn])  <<
        // '\n';
        //}
        out.block(0, 0, 1, N + 1) +=
            (L(i, j).coef[k] * vals * Mat.mats[k + diffn]).head(N + 1);
        if (k + diffn - 1 > -1) {
          out.block(0, N + 1, 1, n) +=
              (L(i, j).coef[k] * valcons * Mat.con_mats[k + diffn - 1]);
        }
      } else {
        out.block(0, 0, 1, N + 1) +=
            (L(i, j).coefFun[k](eval(i, j)) * vals * Mat.mats[k + diffn])
                .head(N + 1);
        if (k + diffn - 1 > -1) {
          out.block(0, N + 1, 1, n) += (L(i, j).coefFun[k](eval(i, j)) *
                                        valcons * Mat.con_mats[k + diffn - 1]);
        }
      }
    }
    return out;
  }
};

/// \brief Class to compute binomial coefficients \f$ \binom{n}{k} \f$ returns
/// a double output.
template <class T> class nChoosek {
public:
  nChoosek(){};

  T operator()(int n, int k) { return T(nCk_(n, k)); }

  int nCk_(int n, int k) {
    if (k == 0 || k == n) {
      return 1;
    } else if (n < k) {
      return 0;
    }
    return nCk_(n - 1, k - 1) + nCk_(n - 1, k);
  }
};

/// \brief Class to compute binomial coefficients \f$ \binom{n}{k} \f$, returns
/// complex<double> output.
template <class T> class nChoosek<std::complex<T> > {
public:
  nChoosek(){};

  T operator()(int n, int k) { return T(nCk_(n, k)); }

  int nCk_(int n, int k) {
    if (k == 0 || k == n) {
      return 1;
    } else if (n < k) {
      return 0;
    }
    return nCk_(n - 1, k - 1) + nCk_(n - 1, k);
  }
};

/// \]brief Evaluates all the Chebfuns in a  2D vector of ChebfunMats in the
/// row major format at the point "a" in the domain. The size of the 2D vector
/// specified through r and c.
template <class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
feval2D(std::valarray<ChebfunMat<T> > Amat, int r, int c, T a) {

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> outmat(r * Amat[0].r,
                                                          r * Amat[0].r);
  outmat.setConstant(0.0);
  int bre;
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < c; j++) {
      for (int l = 0; l < Amat[i * c + j].r; l++) {
        for (int m = 0; m < Amat[i * c + j].c; m++) {
          outmat(i * Amat[i * c + j].r + l, j * Amat[i * c + j].c + m) =
              Amat[i * c + j](l, m)(a);
          // std::cout << "(" << i << "," << j << "," << l << "," << m << ")"
          //          << '\n';
          // std::cout << "outmat: \n" << outmat << '\n';
          // std::cin >> bre;
        }
      }
    }
  }

  return outmat;
};

/// \]brief Evaluates all the Chebfuns in a  2D vector of ChebfunMats in the
/// row major format at the point "a" in the domain. The size of the 2D vector
/// specified through r and c.
template <class T>
Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
feval2D(std::valarray<ChebfunMat<std::complex<T> > > Amat, int r, int c, T a) {

  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> outmat(
      r * Amat[0].r, r * Amat[0].r);
  outmat.setConstant(0.0);
  int bre;
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < c; j++) {
      for (int l = 0; l < Amat[i * c + j].r; l++) {
        for (int m = 0; m < Amat[i * c + j].c; m++) {
          outmat(i * Amat[i * c + j].r + l, j * Amat[i * c + j].c + m) =
              Amat[i * c + j](l, m)(a);
          // std::cout << "(" << i << "," << j << "," << l << "," << m << ")"
          //        << '\n';
          // std::cout << "outmat: \n" << outmat << '\n';
          // std::cin >> bre;
        }
      }
    }
  }
  return outmat;
};
/// \brief This class computes various SingularValues of a differential
/// block matrix
/// operator using using it's adjoint.
/// Class has various utilities, like computing the adjoint, adjoint boundary
/// conditions, and also computing singular values of the frequency response
/// operator.
template <class T>
class SingularValueDecomposition : public GeneralizedEigenSolver<T> {
private:
public:
  /// \brief Computes the singular values/functions of a Linear block matrix
  /// operator.
  void compute(const LinopMat<T> &A_, const BcMat<T> &Lbc_,
               const BcMat<T> &Rbc_, int num_vals) {
    int bre;
    LinopMat<T> A = A_;
    BcMat<T> Lbc = Lbc_;
    BcMat<T> Rbc = Rbc_;
    LinopMat<T> Astar = Adjoint(A);
    BcMat<T> Bc_adjoint = AdjointBc_analytical(A, Lbc_, Rbc_);

    // Define the composite operator:
    LinopMat<T> Acomp(2 * A.r, 2 * A.r);
    LinopMat<T> Z(A.r, A.r);
    Z.setConstant(0.0);
    Acomp << A, Z, //
        Z, Astar;
    //  for (int i = 0; i < Acomp.r; i ++){
    //    for (int j = 0; j < Acomp.c;j ++){
    //      if(Acomp(i,j).NCC == 0){
    //        std::cout << "("<< i <<"," << j << "):" << '\n';
    //        std::cout << Acomp(i,j).coef << '\n';
    //        std::cin >> bre;
    //      } else {
    //        for (int k = 0; k < Acomp(i,j).n + 1; k ++){
    //          std::cout << "("<< i <<"," << j << "," << k << "):" << '\n';
    //          std::cout << Acomp(i,j).coefFun[k] << '\n';
    //          std::cin >> bre;
    //        }
    //      }
    //    }
    //  }

    //  std::cout << "Bc_adjoint: \n" << '\n';
    //  for (int i = 0; i < Bc_adjoint.m; i++) {
    //    for (int j = 0; j < Bc_adjoint.n; j++) {
    //      std::cout << "(" << i << "," << j << ")" << '\n';
    //      std::cout << Bc_adjoint.L(i, j).coef << '\n';
    //    }
    //}
    //    std::cout << "Bc_adjoint.eval: \n" << Bc_adjoint.eval << '\n';
    //    std::cin >> bre;
    // Form composite boundary conditions:
    BcMat<T> BcComp(Lbc.m + Rbc.m + Bc_adjoint.m, 2 * A.r);
    LinopMat<T> Z1(Lbc.m, Astar.c), Z2(Rbc.m, Astar.c), Z3(Bc_adjoint.m, A.c);
    Z1.setConstant(0.0);
    Z2.setConstant(0.0);
    Z3.setConstant(0.0);
    ind = 2;

    for (int i = 0; i < Lbc.m; i++) {
      for (int j = 0; j < Lbc.n; j++) {
        Lbc.eval(i, j) = -1.0;
      }
    }
    for (int i = 0; i < Rbc.m; i++) {
      for (int j = 0; j < Rbc.n; j++) {
        Rbc.eval(i, j) = 1.0;
      }
    }
    // std::cout << "("<< Lbc.L.r << ","   ")"<< << '\n';
    // std::cout << "Entering ::::::::::::::" << '\n';
    // std::cout << "BcComp: " << "(" << BcComp.m <<","<< BcComp.n<< ")" <<
    // '\n'; std::cout << "Lbc: " << "(" << Lbc.m <<","<< Lbc.n<< ")" << '\n';
    // std::cout << "Z1: " << "(" << Z1.r <<","<< Z1.c<< ")" << '\n';
    // std::cout << "Rbc: " << "(" << Rbc.m <<","<< Rbc.n<< ")" << '\n';
    // std::cout << "Z2: " << "(" << Z2.r <<","<< Z2.c<< ")" << '\n';
    // std::cout << "Z3: " << "(" << Z3.r <<","<< Z3.c<< ")" << '\n';
    // std::cout << "Bc_adjoint: " << "(" << Bc_adjoint.m <<","<< Bc_adjoint.n<<
    // ")" << '\n';
    BcComp.L << Lbc.L, Z1, //
        Rbc.L, Z2,         //
        Z3, Bc_adjoint.L;

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> z1(Lbc.m, Astar.c),
        z2(Rbc.m, Astar.c), z3(Bc_adjoint.m, A.c);
    z1.setConstant(0.0);
    z2.setConstant(0.0);
    z3.setConstant(0.0);
    BcComp.eval << Lbc.eval, z1, //
        Rbc.eval, z2,            //
        z3, Bc_adjoint.eval;
    // std::cout << "done >>>" << std::flush << '\n';
    LinopMat<T> Mmat(Acomp.r, Acomp.c), I(A.r, A.r);
    I.setIdentity();
    Mmat << Z, I, //
        I, Z;
    //        std::cout << "BcComp: \n" << '\n';
    //        for (int i = 0; i < BcComp.m; i++) {
    //          for (int j = 0; j < BcComp.n; j++) {
    //            std::cout << "(" << i << "," << j << ")" << '\n';
    //            std::cout << BcComp.L(i, j).coef << '\n';
    //          }
    //        }
    //        std::cin >> bre;
    // std::cout << "Bccomp.eval: " << BcComp.eval << std::endl;
    // std::cin >> bre;
    GeneralizedEigenSolver<T>::compute(Acomp, Mmat, num_vals, BcComp);
  };

  /// \brief Computes the Singular value(s) of a system in the input-output
  /// form,
  /// \f{align}
  /// [\mathcal{A}\,\psi(\cdot)](y)  \;&=\;
  /// [\mathcal{B}\,\boldmath{d}(\cdot)](y) \\
  /// \phi(y)  \;&=\;
  /// \mathcal{C}\,\psi(y), \f} This is used to the find the frequency response
  /// of block matrix operators.
  void compute(const LinopMat<T> &A_, const LinopMat<T> &B_,
               const LinopMat<T> &C_, const BcMat<T> &Lbc_,
               const BcMat<T> &Rbc_, int num_vals){

  };

  /// \brief Computes the Singular value(s) of a system in the input-output
  /// form,
  /// \f{align}
  /// [\mathcal{A}\,\psi(\cdot)](y)  \;&=\;  B(y)\,d
  /// \phi(y)  \;&=\; \mathcal{C}\,\psi(y),
  /// \f}
  /// here, B(y) is a function and not an operator.
  void compute(const LinopMat<T> &A_, const ChebfunMat<T> &B_,
               const LinopMat<T> &C_, const BcMat<T> &Lbc_,
               const BcMat<T> &Rbc_, int num_vals){

  };

  /// \brief Computes the adjoint of a linear operator. This is a differential
  /// adjoint, and not the conjugate transpose of the operators in the LinopMat.
  /// To get the latter, use cTranspose();
  /// Adjoints are for generic discriptor systems, like
  /// \f$M\partial_{t}\phi(y) = L\phi(y)\f$.
  LinopMat<T> Adjoint(const LinopMat<T> &Lmat_) {
    int bre;
    LinopMat<T> Lmat = Lmat_;
    int r = Lmat.r;
    int c = Lmat.c;
    // std::cin >> bre;
    LinopMat<T> outLmat;
    std::vector<int> highest_each_columnL;
    highest_each_columnL.resize(c);
    std::vector<int> temp_vec_int;
    temp_vec_int.resize(r);
    std::vector<int> alphaSize;
    int total_of_all_orders = 0;

    int counter = 0;
    for (int j = 0; j < c; j++) {
      for (int i = 0; i < r; i++) {
        temp_vec_int[i] = Lmat(i, j).n;
      }
      highest_each_columnL[j] =
          *std::max_element(temp_vec_int.begin(), temp_vec_int.end());
    }

    for (int i = 0; i < c; i++) {
      total_of_all_orders += highest_each_columnL[i];
    }
    int n = *std::max_element(highest_each_columnL.begin(),
                              highest_each_columnL.end());
    // n is the overall order.
    // std::cout << "Total of all orders = " << total_of_all_orders << '\n';
    //  std::cin >> bre;
    std::valarray<ChebfunMat<T> > alphais;
    std::valarray<ChebfunMat<T> > betais;
    outLmat.resize(Lmat.c, Lmat.r);
    alphais.resize(n + 1);
    betais.resize(n + 1);
    int outer_counter = 0;
    // Assign alphais
    for (int k = 0; k < n + 1; k++) {
      alphais[k].resize(Lmat.r, Lmat.c);
      alphais[k].setConstant(0.0);
      for (int i = 0; i < Lmat.r; i++) {
        for (int j = 0; j < Lmat.c; j++) {
          if (Lmat(i, j).n - k >= 0) {
            if (Lmat(i, j).NCC == 0) {
              alphais[k](i, j) = Lmat(i, j).coef[Lmat(i, j).n - k];

            } else {
              alphais[k](i, j) = Lmat(i, j).coefFun[Lmat(i, j).n - k];
            }
          }
        }
      }
    }
    /*std::cout << "alphais[0](0,0).coef" << '\n';
    std::cout << alphais[0](0, 0) << '\n';
    std::cout << alphais[1](0, 0) << '\n';
    std::cout << alphais[2](0, 0) << '\n';
    */
    if (ind == 2) {
      std::cout << "in " << __LINE__ << '\n';
      for (int k = 0; k < alphais.size(); k++) {
        for (int i = 0; i < 2; i++) {
          for (int j = 0; j < 2; j++) {
            std::cout << alphais[k](i, j) << '\n';
            std::cout << "(" << i << "," << j << ")" << '\n';
            std::cout << "k: " << k << '\n';
            // std::cin >> bre;
          }
        }
      }
    }
    // std::cout << "DONE alphai ,..." << '\n';
    //  std::cin >> bre;
    nChoosek<T> nck;
    // Then Multiply binomial coefficients and assign to betais.
    for (int i = 0; i < n + 1; i++) {
      betais[i].resize(Lmat.r, Lmat.c);
      betais[i].setConstant(0.0);
      for (int j = i; j < n + 1; j++) {
        betais[i] =
            betais[i] + (pow(-1, j) * nck(j, i) * diff(alphais[i], j - i));
      }
    }

    /*std::cout << "Im in " << __LINE__ << '\n';
    std::cout << "betais[0](0,0).coef" << '\n';
    std::cout << betais[0](0, 0) << '\n';
    std::cout << betais[1](0, 0) << '\n';
    std::cout << betais[2](0, 0) << '\n';
    */
    // take the cTranspose, we will represent betais same as betais_star.
    for (int i = 0; i < n + 1; i++) {
      betais[i] = betais[i].cTranspose();
    }

    if (ind == 2) {
      std::cout << "in " << __LINE__ << '\n';
      for (int k = 0; k < betais.size(); k++) {
        for (int i = 0; i < 2; i++) {
          for (int j = 0; j < 2; j++) {
            std::cout << betais[k](i, j) << '\n';
            std::cout << "(" << i << "," << j << ")" << '\n';
            std::cout << "k: " << k << '\n';
            //          std::cin >> bre;
          }
        }
      }
    }
    // Now, assemble the adjoint operator from betais_star
    // Assign alphais.
    // std::cout << "Lmatr,c : " << r << "," << c << '\n';
    for (int i = 0; i < Lmat.c; i++) {
      for (int j = 0; j < Lmat.r; j++) {
        outLmat(i, j).ncc(n);
        for (int k = 0; k < n + 1; k++) {
          outLmat(i, j).coefFun[n - k] = betais[k](i, j);
        }
      }
    }

    // std::cout << "DONE!!" << std::flush << '\n';

    return outLmat;
  };

  /// \brief Computes adjoint boundary conditions analytically, needs boundary
  /// conditions all right boundary conditions in rbc, and all left boundary
  /// conditions in lbc. All dependent variables need not have boundary
  /// conditions specified, but it is essential that the number of boundary
  /// conditions be equal to the sum of derivatives of all independent
  /// variables, and all the boundary conditions be linearly independent. The
  /// adjoint boundary conditions will be returned in the output, as a BcMat,
  /// will all boundary conditions for left and right in it, in the orders left
  /// followed by right.
  ///

  // \todo It is possible to have incosistent equations and yet find an Adjoint
  // bc that compensates for the rank deficiency.
  BcMat<T> AdjointBc_analytical(const LinopMat<T> &Lmat_, const BcMat<T> &Lbc_,
                                const BcMat<T> &Rbc_) {
    int bre;
    LinopMat<T> Lmat = Lmat_;
    BcMat<T> Lbc = Lbc_;
    BcMat<T> Rbc = Rbc_;
    int r = Lmat.r;
    int c = Lmat.c;
    std::vector<int> highest_each_columnL;
    highest_each_columnL.resize(c);
    std::vector<int> temp_vec_int;
    temp_vec_int.resize(r);
    std::vector<int> alphaSize;
    int total_of_all_orders = 0;
    // std::cout << "Entered ABCAN" << '\n';

    int counter = 0;
    for (int j = 0; j < c; j++) {
      for (int i = 0; i < r; i++) {
        temp_vec_int[i] = Lmat(i, j).n;
      }
      highest_each_columnL[j] =
          *std::max_element(temp_vec_int.begin(), temp_vec_int.end());
    }

    for (int i = 0; i < c; i++) {
      total_of_all_orders += highest_each_columnL[i];
    }
    int n = *std::max_element(highest_each_columnL.begin(),
                              highest_each_columnL.end());
    // n is the overall order.
    //  std::cin >> bre;
    std::valarray<ChebfunMat<T> > alphais;
    std::valarray<ChebfunMat<T> > betais;
    alphais.resize(n + 1);
    betais.resize(n + 1);
    int outer_counter = 0;
    // Assign alphais
    for (int k = 0; k < n + 1; k++) {
      alphais[k].resize(Lmat.r, Lmat.c);
      alphais[k].setConstant(0.0);
      for (int i = 0; i < Lmat.r; i++) {
        for (int j = 0; j < Lmat.c; j++) {
          if (Lmat(i, j).n - k >= 0) {
            if (Lmat(i, j).NCC == 0) {
              alphais[k](i, j) = Lmat(i, j).coef[Lmat(i, j).n - k];
            } else {
              alphais[k](i, j) = Lmat(i, j).coefFun[Lmat(i, j).n - k];
            }
          }
        }
      }
    }
    // std::cout << "in " << __LINE__ << '\n' << std::flush;
    // std::cin >> bre;
    // Need a 2D vector of ChebfunMats, storing in row major format.
    //
    std::valarray<ChebfunMat<T> > Amat;
    Amat.resize(n * n);
    for (int i = 0; i < n * n; i++) {
      Amat[i].resize(Lmat.r, Lmat.c);
      Amat[i].setConstant(0.0);
    }
    // std::cout << "in " << __LINE__ << '\n' << std::flush;
    nChoosek<T> nck;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n - i; j++) {
        for (int k = i; k < n - j; k++) {
          Amat[i * n + j] = Amat[i * n + j] + (pow(-1, k) * nck(k, i) *
                                               diff(alphais[k + j + 1], k - i));
        }
      }
    }
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Aplus(n * r, n * r),
        Aminus(n * r, n * r);

    Aplus = feval2D(Amat, n, n, 1.0);
    Aminus = feval2D(Amat, n, n, -1.0);
    // std::cout << "Aplus : \n" << Aplus << '\n';
    // std::cout << "Aminus : \n" << Aminus << '\n';
    // std::cin >> bre;
    Amat.resize(0);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Blbc(Lbc.m, n * Lbc.n),
        Brbc(Rbc.m, n * Rbc.n);
    Blbc.setConstant(0.0);
    Brbc.setConstant(0.0);
    // std::cout << "Blbc: " << size(Blbc) << '\n';
    for (int k = 0; k < n + 1; k++) {
      for (int i = 0; i < Lbc.m; i++) {
        for (int j = 0; j < Lbc.n; j++) {
          if (Lbc.L(i, j).n - k >= 0) {
            if (Lbc.L(i, j).NCC == 0) {
              //        std::cout << "(" << i << "," << j << "," << k << ")" <<
              //        '\n';
              Blbc(i, Lbc.n * k + j) = Lbc.L(i, j).coef[Lbc.L(i, j).n - k];
              //      std::cout << "Blbc: \n" << Blbc << '\n';
              //      std::cin >> bre;
            } else {
              std::cout << "boundary conditions cannot have non-constant"
                        << "coefficients" << '\n';
              exit(1);
            }
          }
        }
      }
    }
    for (int k = 0; k < n + 1; k++) {
      for (int i = 0; i < Rbc.m; i++) {
        for (int j = 0; j < Rbc.n; j++) {
          if (Rbc.L(i, j).n - k >= 0) {
            if (Rbc.L(i, j).NCC == 0) {
              Brbc(i, Rbc.n * k + j) = Rbc.L(i, j).coef[Rbc.L(i, j).n - k];
              //    std::cout << "(" << i << "," << j << "," << k << ")" <<
              //    '\n'; std::cout << "Brbc: \n" << Brbc << '\n'; std::cin >>
              //    bre;
            } else {
              std::cout << "boundary conditions cannot have non-constant"
                        << "coefficients" << '\n';
              exit(1);
            }
          }
        }
      }
    }

    Eigen::CompleteOrthogonalDecomposition<
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >
        cod;
    cod.compute(Blbc);
    if (cod.rank() != Blbc.rows()) {
      std::cout << "the bounndary conditions are not linearly independent."
                << "Exiting ..." << '\n';
      exit(1);
    }
    // Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> forcL(lbc_val.size(),
    // 1),
    //    forcR(rbc_val.size(), 1);
    // for (int i = 0; i < lbc_val.size(); i++) {
    //  forcL(i, 0) = lbc_val[i];
    //  }
    //  for (int i = 0; i < rbc_val.size(); i++) {
    //    forcR(i, 0) = rbc_val[i];
    //  }

    // Find the particular solutions:
    //  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> minNormSol_lbc =
    //      cod.solve(forcL);
    // std::cout << "minNormSol_lbc: " << minNormSol_lbc << std::flush << '\n';
    // Find the null-space:
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> V =
        cod.matrixZ().adjoint();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Null_space =
        V.block(0, cod.rank(), V.rows(), V.cols() - cod.rank());
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P = cod.colsPermutation();
    Null_space = P * Null_space;

    // Determine as many LI solutions using the Null-space.
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> solutionL(
        Null_space.rows(), Null_space.cols());
    //  for (int i = 0; i < Null_space.cols(); i++) {
    //    solutionL.col(i) = Null_space.col(i) + minNormSol_lbc;
    //  }
    solutionL = Null_space;
    // Do the same for RBC:
    cod.compute(Brbc);
    // std::cout << "in " << __LINE__ << '\n';
    if (cod.rank() != Brbc.rows()) {
      std::cout << "the bounndary conditions are not linearly independent."
                << "Exiting ..." << '\n';
      exit(1);
    }
    // Find the particular solutions:
    // Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> minNormSol_rbc =
    //    cod.solve(forcR);

    // Find the null-space:
    V = cod.matrixZ().adjoint();

    Null_space = V.block(0, cod.rank(), V.rows(), V.cols() - cod.rank());

    P = cod.colsPermutation();

    Null_space = P * Null_space;

    // Determine as many LI solutions using the Null-space.
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> solutionR(
        Null_space.rows(), Null_space.cols());
    // for (int i = 0; i < Null_space.cols(); i++) {
    //  solutionR.col(i) = Null_space.col(i) + minNormSol_rbc;
    //}
    solutionR = Null_space;
    // Now determine the boundary conditions in lbc:
    solutionL.adjointInPlace();
    Aminus.adjointInPlace();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> lbc = solutionL * Aminus;

    solutionR.adjointInPlace();
    Aplus.adjointInPlace();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rbc = solutionR * Aplus;

    // The lbc and rbc are now in the basis u1, u2, express in the normal basis:
    BcMat<T> out;
    out.resize(lbc.rows() + rbc.rows(), Lmat.r);
    //    std::cout << "lbc: \n " << lbc << '\n';
    //    std::cout << "rbc: \n " << rbc << '\n';
    //    std::cin >> bre;
    for (int i = 0; i < lbc.rows(); i++) {
      for (int j = 0; j < Lmat.c; j++) {
        out.L(i, j).set(n - 1);
        out.eval(i, j) = -1.0;
        for (int k = 0; k < n; k++) {
          out.L(i, j).coef[n - k - 1] = lbc(i, Lmat.c * k + j);
        }
      }
    }

    for (int i = 0; i < rbc.rows(); i++) {
      for (int j = 0; j < Lmat.c; j++) {
        out.L(i + lbc.rows(), j).set(n - 1);
        out.eval(i + lbc.rows(), j) = 1.0;
        for (int k = 0; k < n; k++) {
          //          std::cout << "(" << i << "," << j << "," << k << ")" <<
          //          '\n';
          out.L(i + lbc.rows(), j).coef[n - k - 1] = rbc(i, Lmat.c * k + j);
        }
      }
    }
    //  std::cout << "out: \n" << '\n';
    //  for (int i = 0; i < out.m; i++) {
    //    for (int j = 0; j < out.n; j++) {
    //      std::cout << "(" << i << "," << j << ")" << '\n';
    //      std::cout << out.L(i, j).coef << '\n';
    //    }
    //  }
    //  std::cout << "out.eval: \n" << out.eval << '\n';
    //  std::cout << "in " << __LINE__ << '\n';

    // return the expression for boundary conditions:
    return out;
  };
};

/// \brief Given a linear block matrix operator and appropriate boundary
/// conditions, this class will produce an Eigen matrix representing the
/// discretized version. The implementation will naturally involve column
/// pivoting, and the pivot matrix is also stored.
template <class T> class Discretize : public MatGen<T> {
public:
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> subs_mat;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat_temp;
  std::vector<int> highest_each_column;
  int num_bc;
  /// Generates the Eigen matrix. Note that columns are pivoted according to P.
  /// Lmat is the operator, bc is the boundary condition and n_ is the highest
  /// order.
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
  operator()(const LinopMat<T> &Lmat_, const BcMat<T> &bc_) {

    LinopMat<T> Lmat = Lmat_;
    BcMat<T> bc = bc_;
    int r = Lmat.r, c = Lmat.c;
    if (Lmat_.c != bc_.n) {
      std::cout << "The number of comlumns have to be the same in " << __LINE__
                << '\n';
      std::cout << ". Exiting..." << '\n';
      exit(1);
    }
    num_bc = bc.m;
    highest_each_column.resize(c);
    int bre;
    int total_of_all_orders = 0;
    int total_boundary_conditions = 0;
    std::vector<int> temp_vec_int;
    temp_vec_int.resize(r);
    for (int j = 0; j < c; j++) {
      for (int i = 0; i < r; i++) {
        temp_vec_int[i] = Lmat(i, j).n;
      }
      highest_each_column[j] =
          *std::max_element(temp_vec_int.begin(), temp_vec_int.end());
      total_of_all_orders += highest_each_column[j];
    }

    total_boundary_conditions = bc.m;

    // total_of_all_orders has to be equal to total number of boundary
    // conditions, else the problem is ill-posed, if ill-posed, cout the same
    // and exit.
    if (total_of_all_orders != total_boundary_conditions) {
      std::cout << "The problem is ill-posed, the total of the highest "
                   "orders of all "
                   "dependent variables has to be equal to the total number of "
                   "boundary conditions specified."
                   "\n "
                   "Total no. of boundary conditions: "
                << total_boundary_conditions
                << "\n"
                   "Total of all orders: "
                << total_of_all_orders << "\n Exiting ...\n";
      exit(1);
    }

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> masterL(
        r * (N + 1), c * (N + 1) + total_boundary_conditions);

    masterL.setConstant(0.0);

    subs_mat.resize(total_boundary_conditions, (N + 1) * c);
    subs_mat.setConstant(0.0);

    mat_temp.resize(total_boundary_conditions, total_boundary_conditions);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> constraints(
        total_boundary_conditions, c * (N + 1) + total_of_all_orders);
    mat_temp.setConstant(0.0);
    int row_counter = 0, col_counter = 0;
    for (int i = 0; i < bc.m; i++) {
      for (int j = 0; j < bc.n; j++) {
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> temp =
            bc(i, j, highest_each_column[j]);
        constraints.block(i, col_counter, 1, temp.cols()) = temp;
        col_counter += temp.cols();
      }
      col_counter = 0;
    }

    Eigen::ColPivHouseholderQR<
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >
        qr;
    qr.compute(constraints);
    if (qr.rank() != constraints.rows()) {
      std::cout << "The boundary conditions supplied are not "
                << "  linearly independent." << '\n';
      std::cout << "qr.rank = " << qr.rank()
                << ", no. bcs: " << total_boundary_conditions << ". Exiting ..."
                << '\n';
      exit(1);
    }

    P = qr.colsPermutation();
    // Permute constraints
    constraints = constraints * P;
    mat_temp = constraints.block(0, 0, constraints.rows(), constraints.rows());

    Eigen::ColPivHouseholderQR<
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >
        consSolver(mat_temp);
    subs_mat = -consSolver.inverse() *
               constraints.block(0, constraints.rows(), constraints.rows(),
                                 constraints.cols() - constraints.rows());

    row_counter = 0;
    col_counter = 0;
    int master_row_counter = 0;
    int master_col_counter = 0;
    for (int j = 0; j < c; j++) {
      int n = highest_each_column[j];
      MatGen<T>::compute(n);
      for (int i = 0; i < r; i++) {
        int diffn = n - Lmat(i, j).n;
        if (Lmat(i, j).NCC == 0) {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            masterL.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Lmat(i, j).coef[k] *
                (MatGen<T>::mats2[k + diffn].block(0, 0, N + 1, N + 1 + n));
          }
        } else {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            masterL.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Lmat(i, j).coefFun[k].MultMat().block(0, 0, N + 1, N + 1) *
                (MatGen<T>::mats2[k + diffn].block(0, 0, N + 1, N + 1 + n));
          }
        }
        master_row_counter += N + 1;
      }
      master_row_counter = 0;
      master_col_counter += N + 1 + n;
    }

    // Permute columns of M and L:
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> masterL2;
    masterL = masterL * P;
    masterL2 =
        masterL.block(0, constraints.rows(), c * (N + 1), c * (N + 1)) +
        (masterL.block(0, 0, c * (N + 1), constraints.rows()) * subs_mat);
    return masterL2;
  }

  /// \brief Discretization based on a previous call with Boundary conditions.
  /// Then other LinopMats can be discretized based on a previous Lbc called
  /// with operator()(const LinopMat<T> &Lmat, const BcMat<T> &bc).
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
  operator()(const LinopMat<T> &Lmat) {
    int r = Lmat.r, c = Lmat.c;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> masterL(
        r * (N + 1), c * (N + 1) + num_bc);

    masterL.setConstant(0.0);

    int row_counter = 0;
    int col_counter = 0;
    int master_row_counter = 0;
    int master_col_counter = 0;

    for (int j = 0; j < c; j++) {
      int n = highest_each_column[j];
      MatGen<T>::MatGen::compute(n);
      for (int i = 0; i < r; i++) {
        int diffn = n - Lmat(i, j).n;
        if (Lmat(i, j).NCC == 0) {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            masterL.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Lmat(i, j).coef[k] *
                (MatGen<T>::mats2[k + diffn].block(0, 0, N + 1, N + 1 + n));
          }
        } else {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            masterL.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Lmat(i, j).coefFun[k].MultMat().block(0, 0, N + 1, N + 1) *
                (MatGen<T>::mats2[k + diffn].block(0, 0, N + 1, N + 1 + n));
          }
        }
        master_row_counter += N + 1;
      }
      master_row_counter = 0;
      master_col_counter += N + 1 + n;
    }

    // Permute columns of M and L:
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> masterL2;
    masterL = masterL * P;
    masterL2 = masterL.block(0, num_bc, c * (N + 1), c * (N + 1)) +
               (masterL.block(0, 0, c * (N + 1), num_bc) * subs_mat);

    return masterL2;
  }
};

/// \brief Given a linear block matrix operator and appropriate boundary
/// conditions, this class will produce an Eigen matrix representing the
/// discretized version. The implementation will naturally involve column
/// pivoting, and the pivot matrix is also stored.
template <class T>
class Discretize<std::complex<T> > : public MatGen<std::complex<T> > {
public:
  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> P;
  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> subs_mat;
  // Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> mat_temp;
  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> A1;
  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> invmat_temp;
  std::vector<int> highest_each_column;

  int num_bc;

  /// \brief Null Constructor
  Discretize() : num_bc(0) {}
  /// Generates the Eigen matrix. Note that columns are pivoted according to P.
  /// Lmat is the operator, bc is the boundary condition and n_ is the highest
  /// order.
  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
  operator()(const LinopMat<std::complex<T> > &Lmat_,
             const BcMat<std::complex<T> > &bc_) {
    int bre;
    LinopMat<std::complex<T> > Lmat = Lmat_;
    BcMat<std::complex<T> > bc = bc_;
    int r = Lmat.r, c = Lmat.c;
    if (Lmat_.c != bc_.n) {
      std::cout << "The number of comlumns have to be the same in " << __LINE__
                << '\n';
      std::cout << ". Exiting..." << '\n';
      exit(1);
    }
    num_bc = bc.m;
    highest_each_column.resize(c);
    int total_of_all_orders = 0;
    int total_boundary_conditions = 0;
    std::vector<int> temp_vec_int;
    temp_vec_int.resize(r);
    for (int j = 0; j < c; j++) {
      for (int i = 0; i < r; i++) {
        temp_vec_int[i] = Lmat(i, j).n;
      }
      highest_each_column[j] =
          *std::max_element(temp_vec_int.begin(), temp_vec_int.end());
      //  std::cout << "highest_each_column[" << j << "]" <<
      //  highest_each_column[j]
      //          << '\n';
      total_of_all_orders += highest_each_column[j];
    }
    // std::cin >> bre;
    total_boundary_conditions = bc.m;

    // total_of_all_orders has to be equal to total number of boundary
    // conditions, else the problem is ill-posed, if ill-posed, cout the same
    // and exit.
    if (total_of_all_orders != total_boundary_conditions) {
      std::cout << "The problem is ill-posed, the total of the highest "
                   "orders of all "
                   "dependent variables has to be equal to the total number of "
                   "boundary conditions specified."
                   "\n "
                   "Total no. of boundary conditions: "
                << total_boundary_conditions
                << "\n"
                   "Total of all orders: "
                << total_of_all_orders << "\n Exiting ...\n";
      exit(1);
    }

    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> masterL(
        r * (N + 1), c * (N + 1) + total_boundary_conditions);

    masterL.setConstant(0.0);

    subs_mat.resize(total_boundary_conditions, (N + 1) * c);
    subs_mat.setConstant(0.0);
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> mat_temp;
    mat_temp.resize(total_boundary_conditions, total_boundary_conditions);
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> constraints(
        total_boundary_conditions, c * (N + 1) + total_of_all_orders);
    mat_temp.setConstant(0.0);
    int row_counter = 0, col_counter = 0;
    for (int i = 0; i < bc.m; i++) {
      for (int j = 0; j < bc.n; j++) {
        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> temp =
            bc(i, j, highest_each_column[j]);
        constraints.block(i, col_counter, 1, temp.cols()) = temp;
        col_counter += temp.cols();
      }
      col_counter = 0;
    }

    Eigen::ColPivHouseholderQR<
        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
        qr;
    qr.compute(constraints);
    if (qr.rank() != constraints.rows()) {
      std::cout << "The boundary conditions supplied are not "
                << "  linearly independent." << '\n';
      std::cout << "qr.rank = " << qr.rank()
                << ", no. bcs: " << total_boundary_conditions << ". Exiting ..."
                << '\n';
      //  exit(1);
    }
    P = qr.colsPermutation();

    // Permute constraints
    constraints = constraints * P;
    mat_temp = constraints.block(0, 0, constraints.rows(), constraints.rows());

    Eigen::ColPivHouseholderQR<
        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
        consSolver(mat_temp);
    subs_mat = -consSolver.inverse() *
               constraints.block(0, constraints.rows(), constraints.rows(),
                                 constraints.cols() - constraints.rows());
    invmat_temp = consSolver.inverse();
    row_counter = 0;
    col_counter = 0;
    int master_row_counter = 0;
    int master_col_counter = 0;
    for (int j = 0; j < c; j++) {
      int n = highest_each_column[j];
      MatGen<std::complex<T> >::compute(n);
      for (int i = 0; i < r; i++) {
        int diffn = n - Lmat(i, j).n;
        if (Lmat(i, j).NCC == 0) {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            masterL.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Lmat(i, j).coef[k] *
                (MatGen<std::complex<T> >::mats2[k + diffn].block(0, 0, N + 1,
                                                                  N + 1 + n));
          }
        } else {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            masterL.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Lmat(i, j).coefFun[k].MultMat().block(0, 0, N + 1, N + 1) *
                (MatGen<std::complex<T> >::mats2[k + diffn].block(0, 0, N + 1,
                                                                  N + 1 + n));
          }
        }
        master_row_counter += N + 1;
      }
      master_row_counter = 0;
      master_col_counter += N + 1 + n;
    }

    // Permute columns of M and L:
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> masterL2;
    masterL = masterL * P;
    masterL2 =
        masterL.block(0, constraints.rows(), c * (N + 1), c * (N + 1)) +
        (masterL.block(0, 0, c * (N + 1), constraints.rows()) * subs_mat);
    // std::cout << "masterL2: \n" << masterL2 << '\n';
    // std::cin >> bre;
    A1 = masterL.block(0, 0, c * (N + 1), constraints.rows());
    return masterL2;
  }

  /// \brief Discretization based on a previous call with Boundary conditions.
  /// Then other LinopMats can be discretized based on a previous Lbc called
  /// with operator()(const LinopMat<T> &Lmat, const BcMat<T> &bc).
  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
  operator()(const LinopMat<std::complex<T> > &Lmat_) {
    LinopMat<std::complex<T> > Lmat = Lmat_;
    int r = Lmat.r, c = Lmat.c;

    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> masterL(
        r * (N + 1), c * (N + 1) + num_bc);
    masterL.setConstant(0.0);
    int row_counter = 0;
    int col_counter = 0;
    int master_row_counter = 0;
    int master_col_counter = 0;
    for (int j = 0; j < c; j++) {
      int n = highest_each_column[j];
      MatGen<std::complex<T> >::compute(n);
      for (int i = 0; i < r; i++) {
        int diffn = n - Lmat(i, j).n;
        if (Lmat(i, j).NCC == 0) {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            masterL.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Lmat(i, j).coef[k] *
                (MatGen<std::complex<T> >::mats2[k + diffn].block(0, 0, N + 1,
                                                                  N + 1 + n));
          }
        } else {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            masterL.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Lmat(i, j).coefFun[k].MultMat().block(0, 0, N + 1, N + 1) *
                (MatGen<std::complex<T> >::mats2[k + diffn].block(0, 0, N + 1,
                                                                  N + 1 + n));
          }
        }
        master_row_counter += N + 1;
      }
      master_row_counter = 0;
      master_col_counter += N + 1 + n;
    }

    // Permute columns of M and L:
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> masterL2;
    masterL = masterL * P;
    masterL2 = masterL.block(0, num_bc, r * (N + 1), c * (N + 1)) +
               (masterL.block(0, 0, r * (N + 1), num_bc) * subs_mat);
    return masterL2;
  }

  /// Chebyshev differentiation operator for LinopMat
  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>

  ChebDiff(const LinopMat<std::complex<T> > &Lmat_) {
    LinopMat<std::complex<T> > Lmat = Lmat_;
    int r = Lmat.r, c = Lmat.c;
    highest_each_column.resize(c);
    int bre;
    int total_of_all_orders = 0;
    int total_boundary_conditions = 0;
    std::vector<int> temp_vec_int;
    temp_vec_int.resize(r);
    for (int j = 0; j < c; j++) {
      for (int i = 0; i < r; i++) {
        temp_vec_int[i] = Lmat(i, j).n;
      }
      highest_each_column[j] =
          *std::max_element(temp_vec_int.begin(), temp_vec_int.end());
      std::cout << "highest_each_column[" << j << "]" << highest_each_column[j]
                << '\n';
      total_of_all_orders += highest_each_column[j];
    }
    int n = *std::max_element(highest_each_column.begin(),
                              highest_each_column.end());
    std::vector<Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
        diffMats;
    diffMats.resize(n + 1);

    for (int i = 0; i < n + 1; i++) {
      diffMats[i].resize(N + 1, N + 1);
    }

    diffMats[0] = Eigen::Matrix<std::complex<T>, Eigen::Dynamic,
                                Eigen::Dynamic>::Identity(N + 1, N + 1);
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> Dmat(N + 1,
                                                                        N + 1),
        masterL(r * (N + 1), c * (N + 1)), temp;

    Dmat.setConstant(0.0);
    masterL.setConstant(0.0);
    for (int k = 0; k < N + 1; k++) {
      for (int p = k + 1; p < N + 1; p++) {
        if ((p + k) % 2) {
          Dmat(k, p) = p;
        }
      }
    }
    std::valarray<T> y(N + 1), y1(N + 1);
    setChebPts(y);
    Chebfun<std::complex<T> > fun;
    fun.vr = y * y;
    fun.vi = 0.0;
    fun.p2c();
    temp = fun.evc();

    for (int i = 1; i < n + 1; i++) {
      diffMats[i] = Dmat * diffMats[i - 1];
    }
    // temp = diffMats[1] * temp;
    // std::cout << "temp = \n" << temp << '\n';
    // std::cin >> bre;

    int master_row_counter = 0;
    int master_col_counter = 0;
    for (int j = 0; j < c; j++) {
      int n = highest_each_column[j];
      for (int i = 0; i < r; i++) {
        int diffn = n - Lmat(i, j).n;
        if (Lmat(i, j).NCC == 0) {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            masterL.block(master_row_counter, master_col_counter, N + 1,
                          N + 1) +=
                Lmat(i, j).coef[Lmat(i, j).n - k] * (diffMats[k]);
          }
        } else {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            masterL.block(master_row_counter, master_col_counter, N + 1,
                          N + 1) +=
                Lmat(i, j).coefFun[Lmat(i, j).n - k].MultMat().block(
                    0, 0, N + 1, N + 1) *
                (diffMats[k]);
          }
        }
        master_row_counter += N + 1;
      }
      master_row_counter = 0;
      master_col_counter += N + 1;
    }

    return masterL;
  }

  /// Generates the Eigen matrix of discretization as full matrix, without
  /// embedding constraints, by appending them instead. Constraints are at the
  /// bottom.
  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
  MatAppend(const LinopMat<std::complex<T> > &Lmat_,
            const BcMat<std::complex<T> > &bc_) {
    int bre;
    LinopMat<std::complex<T> > Lmat = Lmat_;
    BcMat<std::complex<T> > bc = bc_;
    int r = Lmat.r, c = Lmat.c;
    if (Lmat_.c != bc_.n) {
      std::cout << "The number of comlumns have to be the same in " << __LINE__
                << '\n';
      std::cout << ". Exiting..." << '\n';
      exit(1);
    }
    num_bc = bc.m;
    highest_each_column.resize(c);
    int total_of_all_orders = 0;
    int total_boundary_conditions = 0;
    std::vector<int> temp_vec_int;
    temp_vec_int.resize(r);
    for (int j = 0; j < c; j++) {
      for (int i = 0; i < r; i++) {
        temp_vec_int[i] = Lmat(i, j).n;
      }
      highest_each_column[j] =
          *std::max_element(temp_vec_int.begin(), temp_vec_int.end());
      //  std::cout << "highest_each_column[" << j << "]" <<
      //  highest_each_column[j]
      //          << '\n';
      total_of_all_orders += highest_each_column[j];
    }
    // std::cin >> bre;
    total_boundary_conditions = bc.m;

    // total_of_all_orders has to be equal to total number of boundary
    // conditions, else the problem is ill-posed, if ill-posed, cout the same
    // and exit.
    if (total_of_all_orders != total_boundary_conditions) {
      std::cout << "The problem is ill-posed, the total of the highest "
                   "orders of all "
                   "dependent variables has to be equal to the total number of "
                   "boundary conditions specified."
                   "\n "
                   "Total no. of boundary conditions: "
                << total_boundary_conditions
                << "\n"
                   "Total of all orders: "
                << total_of_all_orders << "\n Exiting ...\n";
      exit(1);
    }

    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> masterL(
        r * (N + 1), c * (N + 1) + total_boundary_conditions);

    masterL.setConstant(0.0);

    subs_mat.resize(total_boundary_conditions, (N + 1) * c);
    subs_mat.setConstant(0.0);
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> mat_temp;
    mat_temp.resize(total_boundary_conditions, total_boundary_conditions);
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> constraints(
        total_boundary_conditions, c * (N + 1) + total_of_all_orders);
    mat_temp.setConstant(0.0);
    int row_counter = 0, col_counter = 0;
    for (int i = 0; i < bc.m; i++) {
      for (int j = 0; j < bc.n; j++) {
        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> temp =
            bc(i, j, highest_each_column[j]);
        constraints.block(i, col_counter, 1, temp.cols()) = temp;
        col_counter += temp.cols();
      }
      col_counter = 0;
    }

    row_counter = 0;
    col_counter = 0;
    int master_row_counter = 0;
    int master_col_counter = 0;
    for (int j = 0; j < c; j++) {
      int n = highest_each_column[j];
      MatGen<std::complex<T> >::compute(n);
      for (int i = 0; i < r; i++) {
        int diffn = n - Lmat(i, j).n;
        if (Lmat(i, j).NCC == 0) {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            masterL.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Lmat(i, j).coef[k] *
                (MatGen<std::complex<T> >::mats2[k + diffn].block(0, 0, N + 1,
                                                                  N + 1 + n));
          }
        } else {
          for (int k = 0; k < Lmat(i, j).n + 1; k++) {
            masterL.block(master_row_counter, master_col_counter, N + 1,
                          N + 1 + n) +=
                Lmat(i, j).coefFun[k].MultMat().block(0, 0, N + 1, N + 1) *
                (MatGen<std::complex<T> >::mats2[k + diffn].block(0, 0, N + 1,
                                                                  N + 1 + n));
          }
        }
        master_row_counter += N + 1;
      }
      master_row_counter = 0;
      master_col_counter += N + 1 + n;
    }

    // Permute columns of M and L:
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> masterL2(
        r * (N + 1) + total_boundary_conditions,
        c * (N + 1) + total_boundary_conditions);

    masterL2 << masterL, constraints;
    return masterL2;
  }
};

/// \brief This class computes various SingularValues of a differential
/// block matrix
/// operator using using it's adjoint.
/// Class has various utilities, like computing the adjoint, adjoint boundary
/// conditions, and also computing singular values of the frequency response
/// operator.
template <class T>
class SingularValueDecomposition<std::complex<T> >
    : public GeneralizedEigenSolver<std::complex<T> > {
private:
public:
  // Options: SIS_SVD_RIGHT: compute only right singular eigenvectors
  //          SIS_SVD_LEFT: compute only left singular eigenvectors.
  //          SIS_SVD: compute both. The former two are twice as fast.
  int svd_flag;
  int isSingular;

  /// \brief This is the H-infinity norm.
  double gamma_opt;
  /// \brief This is the value of \f$\omega\f$ at which H-infinity norm occurs.
  double omega_opt;
  SingularValueDecomposition() : svd_flag(SIS_SVD), isSingular(0) {}
  /// \brief Computes the singular values/functions of a Linear block matrix
  /// operator.
  void compute(const LinopMat<std::complex<T> > &A_,
               const BcMat<std::complex<T> > &Lbc_,
               const BcMat<std::complex<T> > &Rbc_, int num_vals) {
    int bre;
    LinopMat<std::complex<T> > A = A_;
    BcMat<std::complex<T> > Lbc = Lbc_;
    BcMat<std::complex<T> > Rbc = Rbc_;
    LinopMat<std::complex<T> > Astar = Adjoint(A);
    BcMat<std::complex<T> > Bc_adjoint = AdjointBc_analytical(A, Lbc_, Rbc_);

    // Define the composite operator:
    LinopMat<std::complex<T> > Acomp(2 * A.r, 2 * A.r);
    LinopMat<std::complex<T> > Z(A.r, A.r);
    Z.setConstant(0.0);
    Acomp << A, Z, //
        Z, Astar;
    //    for (int i = 0; i < Astar.r; i ++){
    //      for (int j = 0; j < Astar.c;j ++){
    //        if(Astar(i,j).NCC == 0){
    //          std::cout << "ij = ("<< i <<"," << j << "):" << '\n';
    //          std::cout << Astar(i,j).coef << '\n';
    //          std::cin >> bre;
    //        } else {
    //          for (int k = 0; k < Astar(i,j).n + 1; k ++){
    //            std::cout << "("<< i <<"," << j << "," << k << "):" << '\n';
    //            std::cout << Astar(i,j).coefFun[k] << '\n';
    //            std::cin >> bre;
    //          }
    //        }
    //      }
    //    }

    //  std::cout << "Bc_adjoint: \n" << '\n';
    //  for (int i = 0; i < Bc_adjoint.m; i++) {
    //    for (int j = 0; j < Bc_adjoint.n; j++) {
    //      std::cout << "(" << i << "," << j << ")" << '\n';
    //      std::cout << Bc_adjoint.L(i, j).coef << '\n';
    //    }
    //}
    //    std::cout << "Bc_adjoint.eval: \n" << Bc_adjoint.eval << '\n';
    //    std::cin >> bre;
    // Form composite boundary conditions:
    BcMat<std::complex<T> > BcComp(Lbc.m + Rbc.m + Bc_adjoint.m, 2 * A.r);
    LinopMat<std::complex<T> > Z1(Lbc.m, Astar.c), Z2(Rbc.m, Astar.c),
        Z3(Bc_adjoint.m, A.c);
    Z1.setConstant(0.0);
    Z2.setConstant(0.0);
    Z3.setConstant(0.0);

    for (int i = 0; i < Lbc.m; i++) {
      for (int j = 0; j < Lbc.n; j++) {
        Lbc.eval(i, j) = -1.0;
      }
    }
    for (int i = 0; i < Rbc.m; i++) {
      for (int j = 0; j < Rbc.n; j++) {
        Rbc.eval(i, j) = 1.0;
      }
    }
    // std::cout << "("<< Lbc.L.r << ","   ")"<< << '\n';
    // std::cout << "Entering ::::::::::::::" << '\n';
    // std::cout << "BcComp: " << "(" << BcComp.m <<","<< BcComp.n<< ")" <<
    // '\n'; std::cout << "Lbc: " << "(" << Lbc.m <<","<< Lbc.n<< ")" << '\n';
    // std::cout << "Z1: " << "(" << Z1.r <<","<< Z1.c<< ")" << '\n';
    // std::cout << "Rbc: " << "(" << Rbc.m <<","<< Rbc.n<< ")" << '\n';
    // std::cout << "Z2: " << "(" << Z2.r <<","<< Z2.c<< ")" << '\n';
    // std::cout << "Z3: " << "(" << Z3.r <<","<< Z3.c<< ")" << '\n';
    // std::cout << "Bc_adjoint: " << "(" << Bc_adjoint.m <<","<< Bc_adjoint.n<<
    // ")" << '\n';
    BcComp.L << Lbc.L, Z1, //
        Rbc.L, Z2,         //
        Z3, Bc_adjoint.L;

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> z1(Lbc.m, Astar.c),
        z2(Rbc.m, Astar.c), z3(Bc_adjoint.m, A.c);
    z1.setConstant(0.0);
    z2.setConstant(0.0);
    z3.setConstant(0.0);

    BcComp.eval << Lbc.eval, z1, //
        Rbc.eval, z2,            //
        z3, Bc_adjoint.eval;
    // std::cout << "done >>>" << std::flush << '\n';

    LinopMat<std::complex<T> > Mmat(Acomp.r, Acomp.c), I(A.r, A.r);
    I.setIdentity();
    Mmat << Z, I, //
        I, Z;
    //        std::cout << "BcComp: \n" << '\n';
    //        for (int i = 0; i < BcComp.m; i++) {
    //          for (int j = 0; j < BcComp.n; j++) {
    //            std::cout << "(" << i << "," << j << ")" << '\n';
    //            std::cout << BcComp.L(i, j).coef << '\n';
    //          }
    //        }
    //        std::cin >> bre;
    // std::cout << "Bccomp.eval: " << BcComp.eval << std::endl;
    // std::cin >> bre;
    GeneralizedEigenSolver<std::complex<T> >::compute(Acomp, Mmat, num_vals,
                                                      BcComp);
  };

  /// \brief Computes singular values of the frequency response operator
  /// of a system in the input-output
  /// form,
  /// \f{align}
  /// [\mathcal{A}\,\psi(\cdot)](y)  \;&=\;
  /// [\mathcal{B}\,\boldmath{d}(\cdot)](y)\\ \phi(y)  \;&=\;
  /// \mathcal{C}\,\psi(y), \f}.
  void compute(const LinopMat<std::complex<T> > &A_,
               const LinopMat<std::complex<T> > &B_,
               const LinopMat<std::complex<T> > &C_,
               const BcMat<std::complex<T> > &Lbc_,
               const BcMat<std::complex<T> > &Rbc_, int num_vals) {
    int bre;
    LinopMat<std::complex<T> > A = A_;
    LinopMat<std::complex<T> > B = B_;
    LinopMat<std::complex<T> > C = C_;
    BcMat<std::complex<T> > Lbc = Lbc_;
    BcMat<std::complex<T> > Rbc = Rbc_;
    LinopMat<std::complex<T> > Astar = Adjoint(A);
    LinopMat<std::complex<T> > Bstar = Adjoint(B);
    LinopMat<std::complex<T> > Cstar = Adjoint(C);
    BcMat<std::complex<T> > Astar_bc = AdjointBc_analytical(A, Lbc, Rbc);

    // std::cout << "Astar_bc: \n" << '\n';
    // for (int i = 0; i < Astar_bc.m; i++) {
    // for (int j = 0; j < Astar_bc.n; j++) {
    //   std::cout << "ij = (" << i << "," << j << ")" << '\n';
    //   std::cout << Astar_bc.L(i, j).coef << '\n';
    // }
    // }
    // std::cin >> bre;
    // std::cout << "Astar_bc.eval: " << Astar_bc.eval << std::endl;
    // std::cin >> bre;
    // std::valarray<std::complex<T> > rr(N + 1);
    // setChebPts(rr);
    // for (int i = 0; i < Astar.r; i ++){
    //   for (int j = 0; j < Astar.c;j ++){
    //     if(Astar(i,j).NCC == 0){
    //       std::cout << "ij = ("<< i <<"," << j << "):" << '\n';
    //       std::cout << Astar(i,j).coef << '\n';
    //       std::cin >> bre;
    //     } else {
    //       for (int k = 0; k < Astar(i,j).n + 1; k ++){
    //         if (Astar(i,j).coefFun[k].dct_flag == SIS_PHYS_SPACE){
    //           Astar(i,j).coefFun[k].v = Astar(i,j).coefFun[k].v / (rr);
    //         } else {
    //           Astar(i,j).coefFun[k].c2p();
    //           Astar(i,j).coefFun[k].v = Astar(i,j).coefFun[k].v / (rr);
    //           Astar(i,j).coefFun[k].p2c();
    //         }
    //         std::cout << "("<< i <<"," << j << "," << k << "):" << '\n';
    //         std::cout << Astar(i,j).coefFun[k](-1) << '\n';
    //         std::cin >> bre;
    //       }
    //     }
    //   }
    //}

    BcMat<std::complex<T> > A_bc(Lbc.m + Rbc.m, Lbc.n);
    Lbc.eval.setConstant(-1.0);
    Rbc.eval.setConstant(1.0);

    A_bc.L << Lbc.L, //
        Rbc.L;

    A_bc.eval << Lbc.eval, //
        Rbc.eval;          //

    LinopMat<std::complex<T> > BBstar, CstarC;
    BBstar = B * Bstar;
    CstarC = Cstar * C;

    if (svd_flag == SIS_SVD) {
      LinopMat<std::complex<T> > Lmat(A.r + Astar.r, A.c + Astar.c),
          Mmat(A.r + Astar.r, A.c + Astar.c);
      LinopMat<std::complex<T> > Z1(A.r, A.c), Z2(A_bc.m, Astar_bc.n),
          Z3(Astar_bc.m, A_bc.n);
      Z1.setConstant(0.0);
      Z2.setConstant(0.0);
      Z3.setConstant(0.0);
      Lmat << A, Z1, //
          Z1, Astar;
      Mmat << Z1, BBstar, //
          CstarC, Z1;

      BcMat<std::complex<T> > bcf(A_bc.m + Astar_bc.m, A_bc.n + A_bc.n);
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> z2(A_bc.eval.rows(),
                                                          Astar_bc.eval.cols()),
          z3(Astar_bc.eval.rows(), A_bc.eval.cols());
      z2.setConstant(0.0);
      z3.setConstant(0.0);
      bcf.L << A_bc.L, Z2, //
          Z3, Astar_bc.L;
      bcf.eval << A_bc.eval, z2, //
          z3, Astar_bc.eval;
      GeneralizedEigenSolver<std::complex<T> >::compute(Mmat, Lmat,
                                                        Lmat.r * (N + 1), bcf);
      GeneralizedEigenSolver<std::complex<T> >::removeInf();
      GeneralizedEigenSolver<std::complex<T> >::sortByLargestReal();
      GeneralizedEigenSolver<std::complex<T> >::keep(num_vals);
      // GeneralizedEigenSolver<std::complex<T> >::keepConverged();
      // GeneralizedEigenSolver<std::complex<T> >::sortByLargestReal();
    } else if (svd_flag == SIS_SVD_RIGHT) {
      Discretize<std::complex<T> > ADis, AstarDis;
      Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> A0, A0star,
          B0B0star, C0starC0, L, M, invA0, invA0star;
      A0 = ADis(A, A_bc);
      A0star = AstarDis(Astar, Astar_bc);
      B0B0star = AstarDis((BBstar));
      C0starC0 = ADis((CstarC));
      Eigen::ColPivHouseholderQR<
          Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
          qr(A0star);
      if (!qr.isInvertible()) {
        std::cout << "Cannot compute the right singular vectors alone."
                  << "Change svd_flag to SIS_SVD. Exiting ... " << '\n';
        exit(1);
      }
      invA0star = qr.inverse();
      qr.compute(A0);
      if (!qr.isInvertible()) {
        std::cout << "Cannot compute the right singular vectors alone."
                  << "Change svd_flag to SIS_SVD. Exiting ... " << '\n';
        exit(1);
      }
      invA0 = qr.inverse();
      L = B0B0star * invA0star * C0starC0;
      M = A0;
      Eigen::ComplexEigenSolver<
          Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
          eigs;

      GeneralizedEigenSolver<std::complex<T> >::compute(L, A0, ADis);

      GeneralizedEigenSolver<std::complex<T> >::removeInf();
      GeneralizedEigenSolver<std::complex<T> >::sortByLargestReal();
      GeneralizedEigenSolver<std::complex<T> >::keep(num_vals);
    } else if (svd_flag == SIS_SVD_LEFT) {
      Discretize<std::complex<T> > ADis, AstarDis;
      Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> A0, A0star,
          B0B0star, C0starC0, L, M, invA0, invA0star;
      A0 = ADis(A, A_bc);
      A0star = AstarDis(Astar, Astar_bc);
      B0B0star = AstarDis((BBstar));
      C0starC0 = ADis((CstarC));
      Eigen::ColPivHouseholderQR<
          Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
          qr(A0star);
      if (!qr.isInvertible()) {
        std::cout << "Cannot compute the right singular vectors alone."
                  << "Change svd_flag to SIS_SVD. Exiting ... " << '\n';
        exit(1);
      }
      invA0star = qr.inverse();
      qr.compute(A0);
      if (!qr.isInvertible()) {
        std::cout << "Cannot compute the right singular vectors alone."
                  << "Change svd_flag to SIS_SVD. Exiting ... " << '\n';
        exit(1);
      }
      invA0 = qr.inverse();
      L = C0starC0 * invA0 * B0B0star;

      Eigen::ComplexEigenSolver<
          Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
          eigs;

      GeneralizedEigenSolver<std::complex<T> >::compute(L, A0star, AstarDis);

      GeneralizedEigenSolver<std::complex<T> >::removeInf();
      GeneralizedEigenSolver<std::complex<T> >::sortByLargestReal();
      GeneralizedEigenSolver<std::complex<T> >::keep(num_vals);
    }
  };

  /// \brief Computes singular values of the frequency response operator
  /// of a system in the input-output
  /// form,
  /// \f{align}
  /// [\mathcal{A}\,\psi(\cdot)](y)  \;&=\;
  /// [\mathcal{B}\,\boldmath{d}(\cdot)](y)\\ \phi(y)  \;&=\;
  /// \mathcal{C}\,\psi(y), \f}.
  /// When using this, you need to explicitly specify the adjoint boundary
  /// conditions. This is useful in certain discriptor systems.
  void compute(const LinopMat<std::complex<T> > &A_,
               const LinopMat<std::complex<T> > &B_,
               const LinopMat<std::complex<T> > &C_,
               const BcMat<std::complex<T> > &Lbc_,
               const BcMat<std::complex<T> > &Rbc_,
               const BcMat<std::complex<T> > &Lbc_adjoint_,
               const BcMat<std::complex<T> > &Rbc_adjoint_, int num_vals) {
    int bre;
    LinopMat<std::complex<T> > A = A_;
    LinopMat<std::complex<T> > B = B_;
    LinopMat<std::complex<T> > C = C_;
    BcMat<std::complex<T> > Lbc = Lbc_;
    BcMat<std::complex<T> > Rbc = Rbc_;
    BcMat<std::complex<T> > Lbc_adjoint = Lbc_adjoint_;
    BcMat<std::complex<T> > Rbc_adjoint = Rbc_adjoint_;
    LinopMat<std::complex<T> > Astar = Adjoint(A);
    LinopMat<std::complex<T> > Bstar = Adjoint(B);
    LinopMat<std::complex<T> > Cstar = Adjoint(C);

    // std::cout << "Astar_bc: \n" << '\n';
    // for (int i = 0; i < Astar_bc.m; i++) {
    // for (int j = 0; j < Astar_bc.n; j++) {
    //   std::cout << "ij = (" << i << "," << j << ")" << '\n';
    //   std::cout << Astar_bc.L(i, j).coef << '\n';
    // }
    // }
    // std::cin >> bre;
    // std::cout << "Astar_bc.eval: " << Astar_bc.eval << std::endl;
    // std::cin >> bre;
    // std::valarray<std::complex<T> > rr(N + 1);
    // setChebPts(rr);
    //   for (int i = 0; i < Astar.r; i ++){
    //     for (int j = 0; j < Astar.c;j ++){
    //       if(Astar(i,j).NCC == 0){
    //         std::cout << "ij = ("<< i <<"," << j << "):" << '\n';
    //         std::cout << Astar(i,j).coef << '\n';
    //         std::cin >> bre;
    //       } else {
    //         for (int k = 0; k < Astar(i,j).n + 1; k ++){
    //           if (Astar(i,j).coefFun[k].dct_flag == SIS_PHYS_SPACE){
    //             Astar(i,j).coefFun[k].v = Astar(i,j).coefFun[k].v;
    //           } else {
    //             Astar(i,j).coefFun[k].c2p();
    //             Astar(i,j).coefFun[k].v = Astar(i,j).coefFun[k].v;
    //             Astar(i,j).coefFun[k].p2c();
    //           }
    //           std::cout << "("<< i <<"," << j << "," << k << "):" << '\n';
    //           std::cout << Astar(i,j).coefFun[k] << '\n';
    //           std::cin >> bre;
    //         }
    //       }
    //     }
    //  }

    BcMat<std::complex<T> > A_bc(Lbc.m + Rbc.m, Lbc.n);
    BcMat<std::complex<T> > Astar_bc(Lbc_adjoint.m + Rbc_adjoint.m,
                                     Lbc_adjoint.n);
    Lbc.eval.setConstant(-1.0);
    Rbc.eval.setConstant(1.0);

    Lbc_adjoint.eval.setConstant(-1.0);
    Rbc_adjoint.eval.setConstant(1.0);

    A_bc.L << Lbc.L, //
        Rbc.L;
    Astar_bc.L << Lbc_adjoint.L, //
        Rbc_adjoint.L;

    A_bc.eval << Lbc.eval,             //
        Rbc.eval;                      //
    Astar_bc.eval << Lbc_adjoint.eval, //
        Rbc_adjoint.eval;

    LinopMat<std::complex<T> > BBstar, CstarC;
    BBstar = B * Bstar;
    CstarC = Cstar * C;

    if (svd_flag == SIS_SVD) {
      LinopMat<std::complex<T> > Lmat(A.r + Astar.r, A.c + Astar.c),
          Mmat(A.r + Astar.r, A.c + Astar.c);
      LinopMat<std::complex<T> > Z1(A.r, A.c), Z2(A_bc.m, Astar_bc.n),
          Z3(Astar_bc.m, A_bc.n);
      Z1.setConstant(0.0);
      Z2.setConstant(0.0);
      Z3.setConstant(0.0);
      Lmat << A, Z1, //
          Z1, Astar;
      Mmat << Z1, BBstar, //
          CstarC, Z1;

      BcMat<std::complex<T> > bcf(A_bc.m + Astar_bc.m, A_bc.n + A_bc.n);
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> z2(A_bc.eval.rows(),
                                                          Astar_bc.eval.cols()),
          z3(Astar_bc.eval.rows(), A_bc.eval.cols());
      z2.setConstant(0.0);
      z3.setConstant(0.0);
      bcf.L << A_bc.L, Z2, //
          Z3, Astar_bc.L;
      bcf.eval << A_bc.eval, z2, //
          z3, Astar_bc.eval;
      GeneralizedEigenSolver<std::complex<T> >::compute(Mmat, Lmat,
                                                        Lmat.r * (N + 1), bcf);
      // std::cout << "evals:\n"
      //<< GeneralizedEigenSolver<std::complex<T> >::eigenvalues << '\n';
      // std::cin >> bre;
      // std::cout << "evals:\n"
      //<< GeneralizedEigenSolver<std::complex<T> >::eigenvectorsMat[2][0] <<
      //'\n'; std::cin >> bre; std::cout << "evals:\n"
      //<< GeneralizedEigenSolver<std::complex<T> >::eigenvectorsMat[2][1] <<
      //'\n'; std::cin >> bre; std::cout << "evals:\n"
      //<< GeneralizedEigenSolver<std::complex<T> >::eigenvectorsMat[2][2] <<
      //'\n'; std::cin >> bre; std::cout << "evals:\n"
      //<< GeneralizedEigenSolver<std::complex<T> >::eigenvectorsMat[2][3] <<
      //'\n'; std::cin >> bre;

      GeneralizedEigenSolver<std::complex<T> >::removeInf();
      // GeneralizedEigenSolver<std::complex<T> >::keepConverged();
      GeneralizedEigenSolver<std::complex<T> >::sortByLargestReal();
      GeneralizedEigenSolver<std::complex<T> >::keep(num_vals);
    } else if (svd_flag == SIS_SVD_RIGHT) {
      std::cout << "computing only right singular vectors" << '\n';
      Discretize<std::complex<T> > ADis, AstarDis;
      Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> A0, A0star,
          B0B0star, C0starC0, L, M, invA0, invA0star;
      A0 = ADis(A, A_bc);
      A0star = AstarDis(Astar, Astar_bc);
      B0B0star = AstarDis((BBstar));
      C0starC0 = ADis((CstarC));
      Eigen::ColPivHouseholderQR<
          Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
          qr(A0star);
      if (!qr.isInvertible()) {
        std::cout << "Cannot compute the right singular vectors alone."
                  << "Change svd_flag to SIS_SVD. Exiting ... " << '\n';
        exit(1);
      }
      invA0star = qr.inverse();

      L = B0B0star * invA0star * C0starC0;

      GeneralizedEigenSolver<std::complex<T> >::compute(L, A0, ADis);

      GeneralizedEigenSolver<std::complex<T> >::removeInf();
      GeneralizedEigenSolver<std::complex<T> >::sortByLargestReal();
      GeneralizedEigenSolver<std::complex<T> >::keep(num_vals);
    } else if (svd_flag == SIS_SVD_LEFT) {
      Discretize<std::complex<T> > ADis, AstarDis;
      Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> A0, A0star,
          B0B0star, C0starC0, L, M, invA0, invA0star;
      A0 = ADis(A, A_bc);
      A0star = AstarDis(Astar, Astar_bc);
      B0B0star = AstarDis((BBstar));
      C0starC0 = ADis((CstarC));
      Eigen::ColPivHouseholderQR<
          Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
          qr;
      qr.compute(A0);
      if (!qr.isInvertible()) {
        std::cout << "Cannot compute the right singular vectors alone."
                  << "Change svd_flag to SIS_SVD. Exiting ... " << '\n';
        exit(1);
      }
      invA0 = qr.inverse();
      L = C0starC0 * invA0 * B0B0star;

      GeneralizedEigenSolver<std::complex<T> >::compute(L, A0star, AstarDis);

      GeneralizedEigenSolver<std::complex<T> >::removeInf();
      GeneralizedEigenSolver<std::complex<T> >::sortByLargestReal();
      GeneralizedEigenSolver<std::complex<T> >::keep(num_vals);
    }
  };

#ifdef SIS_USE_FEAST
  /// \brief This is for using feast to compute eigenvalues. Still in trial.
  void feast_compute(const LinopMat<std::complex<T> > &A_,
                     const LinopMat<std::complex<T> > &B_,
                     const LinopMat<std::complex<T> > &C_,
                     const BcMat<std::complex<T> > &Lbc_,
                     const BcMat<std::complex<T> > &Rbc_,
                     const BcMat<std::complex<T> > &Lbc_adjoint_,
                     const BcMat<std::complex<T> > &Rbc_adjoint_) {
    int bre;
    LinopMat<std::complex<T> > A = A_;
    LinopMat<std::complex<T> > B = B_;
    LinopMat<std::complex<T> > C = C_;
    BcMat<std::complex<T> > Lbc = Lbc_;
    BcMat<std::complex<T> > Rbc = Rbc_;
    BcMat<std::complex<T> > Lbc_adjoint = Lbc_adjoint_;
    BcMat<std::complex<T> > Rbc_adjoint = Rbc_adjoint_;
    LinopMat<std::complex<T> > Astar = Adjoint(A);
    LinopMat<std::complex<T> > Bstar = Adjoint(B);
    LinopMat<std::complex<T> > Cstar = Adjoint(C);

    BcMat<std::complex<T> > A_bc(Lbc.m + Rbc.m, Lbc.n);
    BcMat<std::complex<T> > Astar_bc(Lbc_adjoint.m + Rbc_adjoint.m,
                                     Lbc_adjoint.n);
    Lbc.eval.setConstant(-1.0);
    Rbc.eval.setConstant(1.0);

    Lbc_adjoint.eval.setConstant(-1.0);
    Rbc_adjoint.eval.setConstant(1.0);

    A_bc.L << Lbc.L, //
        Rbc.L;
    Astar_bc.L << Lbc_adjoint.L, //
        Rbc_adjoint.L;

    A_bc.eval << Lbc.eval,             //
        Rbc.eval;                      //
    Astar_bc.eval << Lbc_adjoint.eval, //
        Rbc_adjoint.eval;

    LinopMat<std::complex<T> > BBstar, CstarC;
    BBstar = B * Bstar;
    CstarC = Cstar * C;

    LinopMat<std::complex<T> > Lmat(A.r + Astar.r, A.c + Astar.c),
        Mmat(A.r + Astar.r, A.c + Astar.c);
    LinopMat<std::complex<T> > Z1(A.r, A.c), Z2(A_bc.m, Astar_bc.n),
        Z3(Astar_bc.m, A_bc.n);
    Z1.setConstant(0.0);
    Z2.setConstant(0.0);
    Z3.setConstant(0.0);
    Lmat << A, Z1, //
        Z1, Astar;
    Mmat << Z1, BBstar, //
        CstarC, Z1;

    BcMat<std::complex<T> > bcf(A_bc.m + Astar_bc.m, A_bc.n + A_bc.n);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> z2(A_bc.eval.rows(),
                                                        Astar_bc.eval.cols()),
        z3(Astar_bc.eval.rows(), A_bc.eval.cols());
    z2.setConstant(0.0);
    z3.setConstant(0.0);
    bcf.L << A_bc.L, Z2, //
        Z3, Astar_bc.L;
    bcf.eval << A_bc.eval, z2, //
        z3, Astar_bc.eval;

    GeneralizedEigenSolver<std::complex<T> >::feast_compute(Mmat, Lmat, 10,
                                                            bcf);
    // GeneralizedEigenSolver<std::complex<T> >::removeInf();
    GeneralizedEigenSolver<std::complex<T> >::sortByLargestReal();

    // std::cout << "evals:\n"
    //<< GeneralizedEigenSolver<std::complex<T> >::eigenvalues << '\n';
    // std::cin >> bre;
    // std::cout << "evals:\n"
    //<< GeneralizedEigenSolver<std::complex<T> >::eigenvectorsMat[2][0] <<
    //'\n'; std::cin >> bre; std::cout << "evals:\n"
    //<< GeneralizedEigenSolver<std::complex<T> >::eigenvectorsMat[2][1] <<
    //'\n'; std::cin >> bre; std::cout << "evals:\n"
    //<< GeneralizedEigenSolver<std::complex<T> >::eigenvectorsMat[2][2] <<
    //'\n'; std::cin >> bre; std::cout << "evals:\n"
    //<< GeneralizedEigenSolver<std::complex<T> >::eigenvectorsMat[2][3] <<
    //'\n'; std::cin >> bre;
  };
#endif
  /// \brief Computes the H-infinity norm of a system in the input-output
  /// form,
  /// \f{align}
  /// [\mathcal{I\omega E - A}\,\psi(\cdot)](y)  \;&=\;
  /// [\mathcal{B}\,\boldmath{d}(\cdot)](y) \phi(y)  \;&=\;
  /// \mathcal{C}\,\psi(y), \f}.
  /// Still needs validation and improvement. Works only for the evolution form.
  void HinfNorm(const LinopMat<std::complex<T> > &E_,
                const LinopMat<std::complex<T> > &A_,
                const LinopMat<std::complex<T> > &B_,
                const LinopMat<std::complex<T> > &C_,
                const BcMat<std::complex<T> > &Lbc_,
                const BcMat<std::complex<T> > &Rbc_,
                const BcMat<std::complex<T> > &Lbc_adjoint_,
                const BcMat<std::complex<T> > &Rbc_adjoint_) {
    // An initial guess to omega.
    double seed = 10.0;
    omega_opt = seed;
    int bre;

    std::complex<T> ii(0.0, 1.0);
    LinopMat<std::complex<T> > E = E_;
    LinopMat<std::complex<T> > A = A_;
    LinopMat<std::complex<T> > B = B_;
    LinopMat<std::complex<T> > C = C_;
    BcMat<std::complex<T> > Lbc = Lbc_;
    BcMat<std::complex<T> > Rbc = Rbc_;
    BcMat<std::complex<T> > Lbc_adjoint = Lbc_adjoint_;
    BcMat<std::complex<T> > Rbc_adjoint = Rbc_adjoint_;
    LinopMat<std::complex<T> > Estar = Adjoint(E);
    LinopMat<std::complex<T> > Astar = Adjoint(A);
    LinopMat<std::complex<T> > Bstar = Adjoint(B);
    LinopMat<std::complex<T> > Cstar = Adjoint(C);
    LinopMat<std::complex<T> > BBstar = B * Bstar;
    LinopMat<std::complex<T> > CstarC = Cstar * C;
    LinopMat<std::complex<T> > A2, A2star, HamiltonL, HamiltonM;
    BcMat<std::complex<T> > A_bc(Lbc.m + Rbc.m, Lbc.n);
    BcMat<std::complex<T> > Astar_bc(Lbc_adjoint.m + Rbc_adjoint.m,
                                     Lbc_adjoint.n);

    A_bc.L << Lbc.L, //
        Rbc.L;
    Astar_bc.L << Lbc_adjoint.L, //
        Rbc_adjoint.L;

    A_bc.eval << Lbc.eval,             //
        Rbc.eval;                      //
    Astar_bc.eval << Lbc_adjoint.eval, //
        Rbc_adjoint.eval;

    BcMat<std::complex<T> > bcf(A_bc.m + Astar_bc.m, A_bc.n + A_bc.n);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> z2(A_bc.eval.rows(),
                                                        Astar_bc.eval.cols()),
        z3(Astar_bc.eval.rows(), A_bc.eval.cols());
    z2.setConstant(0.0);
    z3.setConstant(0.0);
    LinopMat<std::complex<T> > Z1(A.r, A.c), Z2(A_bc.m, Astar_bc.n),
        Z3(Astar_bc.m, A_bc.n);
    Z1.setConstant(0.0);
    Z2.setConstant(0.0);
    Z3.setConstant(0.0);

    bcf.L << A_bc.L, Z2, //
        Z3, Astar_bc.L;
    bcf.eval << A_bc.eval, z2, //
        z3, Astar_bc.eval;

    double gamma_lb, gamma_ub, gamma;
    double Hinf_eps = 1e-4;

    A2 = (ii * omega_opt * E) - A;
    compute(A2, B, C, Lbc, Rbc, Lbc_adjoint, Rbc_adjoint, A2.r * (N + 1));
    gamma_lb = GeneralizedEigenSolver<std::complex<T> >::eigenvalues[0].real();
    // gamma_lb= std::sqrt(gamma_lb);
    while (true) {
      gamma = (1.0 + 2.0 * Hinf_eps) * gamma_lb;
      std::cout << "gamma : " << gamma << '\n';
      //  std::cin >> bre;
      //  std::cout << "A: ("<< A.r <<"," << A.c << ")" << '\n';
      //  std::cout << "BBstar: ("<< BBstar.r <<"," << BBstar.c << ")" << '\n';
      //  std::cout << "CstarC: ("<< CstarC.r <<"," << CstarC.c << ")" << '\n';
      //  std::cout << "Astar: ("<< Astar.r <<"," << Astar.c << ")" << '\n';
      //  std::cin >> bre;

      HamiltonL.resize(2 * A2.r, 2 * A2.c);
      HamiltonM.resize(2 * A2.r, 2 * A2.c);
      HamiltonL << A, (BBstar / (-gamma * gamma)), //
          CstarC, Astar;
      HamiltonM << E, Z1, //
          Z1, Estar;
      GeneralizedEigenSolver<std::complex<T> >::compute(HamiltonL, HamiltonM,
                                                        4 * A.r * (N + 1), bcf);
      // GeneralizedEigenSolver<std::complex<T> >::keepConverged();
      GeneralizedEigenSolver<std::complex<T> >::sortByLargestReal();
      std::ofstream outf("eval_saver.txt");

      Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> mattemp(
          GeneralizedEigenSolver<std::complex<T> >::eigenvalues.rows(), 2);
      mattemp << GeneralizedEigenSolver<std::complex<T> >::alpha,
          GeneralizedEigenSolver<std::complex<T> >::beta;

      std::cout << "eigenvalues from Hamiltonian: \n"
                << GeneralizedEigenSolver<std::complex<T> >::eigenvalues
                << '\n';
      outf << "eigenvalues from Hamiltonian: \n" << mattemp << '\n';

      outf.close();
      std::cin >> bre;
      std::vector<T> omegas;
      for (int i = 0;
           i < GeneralizedEigenSolver<std::complex<T> >::eigenvalues.size();
           i = i + 1) {
        if (std::abs(GeneralizedEigenSolver<std::complex<T> >::eigenvalues[i]
                         .imag()) < 1e-9) {
          omegas.push_back(std::abs(
              GeneralizedEigenSolver<std::complex<T> >::eigenvalues[i].real()));
        }
      }
      if (omegas.size() == 0) {
        A2 = (ii * omega_opt * E) - A;
        int num_vals = A2.r * 2 * (N + 1);
        compute(A2, B, C, Lbc, Rbc, Lbc_adjoint, Rbc_adjoint, num_vals);
        gamma_ub = gamma;
        break;
      }

      std::sort(omegas.begin(), omegas.end());
      std::valarray<T> ms;
      ms.resize(omegas.size() - 1);
      for (int i = 0; i < omegas.size() - 1; i++) {
        ms[i] = (omegas[i] + omegas[i + 1]) / 2.0;
        std::cout << "ms[" << i << "] : " << ms[i] << '\n';
      }
      omegas.clear();
      std::vector<std::vector<T> > gammas;
      gammas.resize(ms.size());
      for (int i = 0; i < ms.size(); i++) {
        // std::cout << "ms["<< i << "] : " << ms[i]  << '\n';
        A2 = (ii * ms[i] * E) - A;
        int num_vals = A2.r * 2 * (N + 1);
        compute(A2, B, C, Lbc, Rbc, Lbc_adjoint, Rbc_adjoint, num_vals);
        // std::cout << "in " << __LINE__ << '\n';
        gammas[i].resize(2);
        // std::cout << "size: " << GeneralizedEigenSolver<std::complex<T>
        // >::eigenvalues.size() << '\n';
        gammas[i][0] =
            GeneralizedEigenSolver<std::complex<T> >::eigenvalues[0].real();
        gammas[i][1] = i;
      }

      std::sort(gammas.begin(), gammas.end());
      std::cout << "gamma end: " << gammas[gammas.size() - 1][0] << '\n';
      std::cout << "gamma 0: " << gammas[0][0] << '\n';
      std::cin >> bre;
      gamma_lb = gammas[gammas.size() - 1][0];
      omega_opt = ms[gammas[gammas.size() - 1][1]];
      std::cout << "omega_opt = " << omega_opt << '\n';
    }
    gamma_opt = (gamma_lb + gamma_ub) / 2.0;
  }

  /// \brief Computes the formal adjoint of a Linear operator.
  /// This is a differential
  /// adjoint, and not the conjugate transpose of the operators in the LinopMat.
  /// To get the latter, use cTranspose();
  /// Adjoints are for generic discriptor systems, like
  /// \f$M\partial_{t}\phi(y) = L\phi(y)\f$.
  LinopMat<std::complex<T> > Adjoint(const LinopMat<std::complex<T> > &Lmat_) {
    int bre;
    LinopMat<std::complex<T> > Lmat = Lmat_;
    int r = Lmat.r;
    int c = Lmat.c;
    //  std::cin >> bre;
    LinopMat<std::complex<T> > outLmat;
    std::vector<int> highest_each_columnL;
    highest_each_columnL.resize(c);
    std::vector<int> temp_vec_int;
    temp_vec_int.resize(r);
    std::vector<int> alphaSize;
    int total_of_all_orders = 0;

    int counter = 0;
    for (int j = 0; j < c; j++) {
      for (int i = 0; i < r; i++) {
        temp_vec_int[i] = Lmat(i, j).n;
      }
      highest_each_columnL[j] =
          *std::max_element(temp_vec_int.begin(), temp_vec_int.end());
    }

    for (int i = 0; i < c; i++) {
      total_of_all_orders += highest_each_columnL[i];
    }
    int n = *std::max_element(highest_each_columnL.begin(),
                              highest_each_columnL.end());
    // n is the overall order.
    // std::cout << "Total of all orders = " << total_of_all_orders << '\n';
    //  std::cin >> bre;
    std::valarray<ChebfunMat<std::complex<T> > > alphais;
    std::valarray<ChebfunMat<std::complex<T> > > betais;
    outLmat.resize(Lmat.c, Lmat.r);
    alphais.resize(n + 1);
    betais.resize(n + 1);
    int outer_counter = 0;
    // Assign alphais
    for (int k = 0; k < n + 1; k++) {
      alphais[k].resize(Lmat.r, Lmat.c);
      alphais[k].setConstant(0.0);
      for (int i = 0; i < Lmat.r; i++) {
        for (int j = 0; j < Lmat.c; j++) {
          if (Lmat(i, j).n - k >= 0) {
            if (Lmat(i, j).NCC == 0) {
              alphais[k](i, j) = Lmat(i, j).coef[Lmat(i, j).n - k];

            } else {
              alphais[k](i, j) = Lmat(i, j).coefFun[Lmat(i, j).n - k];
            }
          }
        }
      }
    }
    /*std::cout << "alphais[0](0,0).coef" << '\n';
    std::cout << alphais[0](0, 0) << '\n';
    std::cout << alphais[1](0, 0) << '\n';
    std::cout << alphais[2](0, 0) << '\n';
    */
    //    if (ind == 2) {
    //      std::cout << "in " << __LINE__ << '\n';
    //      for (int k = 0; k < alphais.size(); k++) {
    //        for (int i = 0; i < 2; i++) {
    //          for (int j = 0; j < 2; j++) {
    //            std::cout << alphais[k](i, j) << '\n';
    //            std::cout << "(" << i << "," << j << ")" << '\n';
    //            std::cout << "k: " << k << '\n';
    //            //          std::cin >> bre;
    //          }
    //        }
    //      }
    //    }
    /// std::cout << "DONE alphai ,..." << '\n';
    //    std::cin >> bre;
    nChoosek<std::complex<T> > nck;
    // Then Multiply binomial coefficients and assign to betais.
    for (int i = 0; i < n + 1; i++) {
      betais[i].resize(Lmat.r, Lmat.c);
      betais[i].setConstant(0.0);
      for (int j = i; j < n + 1; j++) {
        betais[i] =
            betais[i] + (pow(-1, j) * nck(j, i) * diff(alphais[j], j - i));
      }
    }

    /*std::cout << "Im in " << __LINE__ << '\n';
    std::cout << "betais[0](0,0).coef" << '\n';
    std::cout << betais[0](0, 0) << '\n';
    std::cout << betais[1](0, 0) << '\n';
    std::cout << betais[2](0, 0) << '\n';
    */
    // take the cTranspose, we will represent betais same as betais_star.
    for (int i = 0; i < n + 1; i++) {
      betais[i] = betais[i].cTranspose();
    }

    //    if (ind == 2) {
    //      std::cout << "in " << __LINE__ << '\n';
    //      for (int k = 0; k < betais.size(); k++) {
    //        for (int i = 0; i < 2; i++) {
    //          for (int j = 0; j < 2; j++) {
    //            std::cout << betais[k](i, j) << '\n';
    //            std::cout << "(" << i << "," << j << ")" << '\n';
    //            std::cout << "k: " << k << '\n';
    //            //        std::cin >> bre;
    //          }
    //        }
    //      }
    //    }
    // Now, assemble the adjoint operator from betais_star
    // Assign alphais.
    /// std::cout << "Lmatr,c : " << r << "," << c << '\n';
    for (int i = 0; i < Lmat.c; i++) {
      for (int j = 0; j < Lmat.r; j++) {
        outLmat(i, j).ncc(n);
        for (int k = 0; k < n + 1; k++) {
          outLmat(i, j).coefFun[n - k] = betais[k](i, j);
        }
      }
    }
    //  std::cout << "in Adjoint::::::" << '\n';
    //  for (int i = 0; i < outLmat.r; i ++){
    //    for (int j = 0; j < outLmat.c; j ++){
    //      for (int k = 0; k < outLmat(i,j).n + 1; k ++){
    //        std::cout << "ijk: ("<< i <<","<<j<<","<<k<< ")" << '\n';
    //        std::cout << outLmat(i,j).coefFun[k] << '\n';
    //        std::cin >> bre;
    //      }
    //    }
    //  }
    //    // Reduce the orders so that orders of Linops are reduced according to
    //    a
    // transpose of orders.
    LinopMat<std::complex<T> > tempMat(outLmat.r, outLmat.c);
    for (int i = 0; i < outLmat.r; i++) {
      for (int j = 0; j < outLmat.c; j++) {
        tempMat(i, j).ncc(Lmat(j, i).n);
        for (int k = 0; k < Lmat(j, i).n + 1; k++) {
          tempMat(i, j).coefFun[Lmat(j, i).n - k] =
              outLmat(i, j).coefFun[outLmat(i, j).n - k];
        }
      }
    }
    //  total_of_all_orders = 0;

    //  counter = 0;
    //  for (int j = 0; j < c; j++) {
    //    for (int i = 0; i < r; i++) {
    //      temp_vec_int[i] = Lmat(i, j).n;
    //    }
    //    highest_each_columnL[j] =
    //        *std::max_element(temp_vec_int.begin(), temp_vec_int.end());
    //    std::cout << "hecl["<< j << "]: " << highest_each_columnL[j] << '\n';
    //  }

    //  for (int i = 0; i < c; i++) {
    //    total_of_all_orders += highest_each_columnL[i];
    //  }
    //  n = *std::max_element(highest_each_columnL.begin(),
    //                            highest_each_columnL.end());

    return tempMat;
  };

  /// \brief Computes adjoint boundary conditions analytically, needs boundary
  /// conditions all right boundary conditions in rbc, and all left boundary
  /// conditions in lbc. All dependent variables need not have boundary
  /// conditions specified, but it is essential that the number of boundary
  /// conditions be equal to the sum of derivatives of all independent
  /// variables, and all the boundary conditions be linearly independent. The
  /// adjoint boundary conditions will be returned in the output, as BcMat,
  /// will all boundary conditions for left and right in it, in the orders left
  /// followed by right.
  ///

  // \todo It is possible to have incosistent equations and yet find an Adjoint
  // bc that compensates for the rank deficiency.
  BcMat<std::complex<T> >
  AdjointBc_analytical(const LinopMat<std::complex<T> > &Lmat_,
                       const BcMat<std::complex<T> > &Lbc_,
                       const BcMat<std::complex<T> > &Rbc_) {
    int bre;
    LinopMat<std::complex<T> > Lmat = Lmat_;
    BcMat<std::complex<T> > Lbc = Lbc_;
    BcMat<std::complex<T> > Rbc = Rbc_;
    int r = Lmat.r;
    int c = Lmat.c;
    std::vector<int> highest_each_columnL;
    highest_each_columnL.resize(c);
    std::vector<int> temp_vec_int;
    temp_vec_int.resize(r);
    std::vector<int> alphaSize;
    int total_of_all_orders = 0;

    int counter = 0;
    for (int j = 0; j < c; j++) {
      for (int i = 0; i < r; i++) {
        temp_vec_int[i] = Lmat(i, j).n;
      }
      highest_each_columnL[j] =
          *std::max_element(temp_vec_int.begin(), temp_vec_int.end());
    }

    for (int i = 0; i < c; i++) {
      total_of_all_orders += highest_each_columnL[i];
    }
    int n = *std::max_element(highest_each_columnL.begin(),
                              highest_each_columnL.end());
    // n is the overall order.
    //  std::cin >> bre;
    std::valarray<ChebfunMat<std::complex<T> > > alphais;
    std::valarray<ChebfunMat<std::complex<T> > > betais;
    alphais.resize(n + 1);
    betais.resize(n + 1);
    int outer_counter = 0;
    // Assign alphais
    for (int k = 0; k < n + 1; k++) {
      alphais[k].resize(Lmat.r, Lmat.c);
      alphais[k].setConstant(0.0);
      for (int i = 0; i < Lmat.r; i++) {
        for (int j = 0; j < Lmat.c; j++) {
          if (Lmat(i, j).n - k >= 0) {
            if (Lmat(i, j).NCC == 0) {
              alphais[k](i, j) = Lmat(i, j).coef[Lmat(i, j).n - k];
            } else {
              alphais[k](i, j) = Lmat(i, j).coefFun[Lmat(i, j).n - k];
            }
          }
        }
      }
    }
    // Need a 2D vector of ChebfunMats, storing in row major format.
    //
    std::valarray<ChebfunMat<std::complex<T> > > Amat;
    Amat.resize(n * n);
    for (int i = 0; i < n * n; i++) {
      Amat[i].resize(Lmat.r, Lmat.c);
      Amat[i].setConstant(0.0);
    }
    nChoosek<std::complex<T> > nck;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n - i; j++) {
        for (int k = i; k < n - j; k++) {
          Amat[i * n + j] = Amat[i * n + j] + (pow(-1, k) * nck(k, i) *
                                               diff(alphais[k + j + 1], k - i));
        }
      }
    }

    if (isSingular) {
      std::valarray<std::complex<double> > rr(N + 1);
      setChebPts(rr);
      r = 0.5 * r + 0.5;
      for (int i = 0; i < n * n; i++) {
        for (int l = 0; l < Amat[i].r; l++) {
          for (int m = 0; m < Amat[i].c; m++) {
            if (Amat[i](l, m).dct_flag == SIS_PHYS_SPACE) {
              Amat[i](l, m).v = Amat[i](l, m).v / (rr * rr);
            } else {
              Amat[i](l, m).c2p();
              Amat[i](l, m).v = Amat[i](l, m).v / (rr * rr);
              Amat[i](l, m).p2c();
            }
          }
        }
      }
    }
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> Aplus(n * r,
                                                                         n * r),
        Aminus(n * r, n * r);

    Aplus = feval2D(Amat, n, n, 1.0);
    Aminus = feval2D(Amat, n, n, -1.0);
    //  std::cout << "Aplus : \n" << Aplus << '\n';
    //  std::cout << "Aminus : \n" << Aminus << '\n';
    //  std::cin >> bre;
    Amat.resize(0);
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> Blbc(
        Lbc.m, n * Lbc.n),
        Brbc(Rbc.m, n * Rbc.n);
    Blbc.setConstant(0.0);
    Brbc.setConstant(0.0);
    ///   std::cout << "Blbc: " << size(Blbc) << '\n';
    for (int k = 0; k < n + 1; k++) {
      for (int i = 0; i < Lbc.m; i++) {
        for (int j = 0; j < Lbc.n; j++) {
          if (Lbc.L(i, j).n - k >= 0) {
            if (Lbc.L(i, j).NCC == 0) {
              //          std::cout << "(" << i << "," << j << "," << k << ")"
              //          << '\n';
              Blbc(i, Lbc.n * k + j) = Lbc.L(i, j).coef[Lbc.L(i, j).n - k];
              //          std::cout << "Blbc: \n" << Blbc << '\n';
              //          std::cin >> bre;
            } else {
              std::cout << "boundary conditions cannot have non-constant"
                        << "coefficients" << '\n';
              exit(1);
            }
          }
        }
      }
    }
    for (int k = 0; k < n + 1; k++) {
      for (int i = 0; i < Rbc.m; i++) {
        for (int j = 0; j < Rbc.n; j++) {
          if (Rbc.L(i, j).n - k >= 0) {
            if (Rbc.L(i, j).NCC == 0) {
              Brbc(i, Rbc.n * k + j) = Rbc.L(i, j).coef[Rbc.L(i, j).n - k];
              ///         std::cout << "(" << i << "," << j << "," << k << ")"
              ///         << '\n'; std::cout << "Brbc: \n" << Brbc << '\n';
              ///         std::cin >> bre;
            } else {
              std::cout << "boundary conditions cannot have non-constant"
                        << "coefficients" << '\n';
              exit(1);
            }
          }
        }
      }
    }
    // std::cout << "blbc:" << size(Blbc) << '\n';
    // std::cout << "blbc: \n" << Blbc << '\n';
    Eigen::CompleteOrthogonalDecomposition<
        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
        cod;
    cod.compute(Blbc);
    if (cod.rank() != Blbc.rows()) {
      std::cout << "the bounndary conditions are not linearly independent."
                << "Exiting ..." << '\n';
      exit(1);
    }
    //  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> forcL(
    //      lbc_val.size(), 1),
    //      forcR(rbc_val.size(), 1);
    //  for (int i = 0; i < lbc_val.size(); i++) {
    //    forcL(i, 0) = lbc_val[i];
    //  }
    //  for (int i = 0; i < rbc_val.size(); i++) {
    //    forcR(i, 0) = rbc_val[i];
    //  }

    // Find the particular solutions:
    // Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
    //    minNormSol_lbc = cod.solve(forcL);
    /// std::cout << "minNormSol_lbc: " << minNormSol_lbc << std::flush << '\n';
    // Find the null-space:
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> V =
        cod.matrixZ().adjoint();
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> Null_space =
        V.block(0, cod.rank(), V.rows(), V.cols() - cod.rank());
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P = cod.colsPermutation();
    Null_space = P * Null_space;
    //  std::cout << "Null_space:" << size(Null_space) << '\n';
    cod.compute(Null_space);
    //    std::cout << "rank of Null: " << cod.rank() << '\n';

    // Determine as many LI solutions using the Null-space.
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> solutionL(
        Null_space.rows(), Null_space.cols());
    // for (int i = 0; i < Null_space.cols(); i++) {
    //  solutionL.col(i) = Null_space.col(i) + minNormSol_lbc;
    //}
    solutionL = Null_space;
    // Do the same for RBC:
    cod.compute(Brbc);
    if (cod.rank() != Brbc.rows()) {
      std::cout << "the bounndary conditions are not linearly independent."
                << "Exiting ..." << '\n';
      exit(1);
    }
    // Find the particular solutions:
    // Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
    //    minNormSol_rbc = cod.solve(forcR);

    // Find the null-space:
    V = cod.matrixZ().adjoint();

    Null_space = V.block(0, cod.rank(), V.rows(), V.cols() - cod.rank());

    P = cod.colsPermutation();

    Null_space = P * Null_space;

    // Determine as many LI solutions using the Null-space.
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> solutionR(
        Null_space.rows(), Null_space.cols());
    // for (int i = 0; i < Null_space.cols(); i++) {
    //  solutionR.col(i) = Null_space.col(i) + minNormSol_rbc;
    //}
    solutionR = Null_space;

    // Now determine the boundary conditions in lbc:
    solutionL.adjointInPlace();
    Aminus.adjointInPlace();
    cod.compute(Aminus);

    //  std::cout << " size of Aminus:" << size(Aminus) << " rank of Aminus: "
    //  << cod.rank() << '\n'; std::cin >> bre;
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
        lbc1 = solutionL * Aminus,
        Ptemp, lbc;
    //      std::cout << "lbc1: \n" << lbc1  << '\n';
    //      std::cout << "in :::::; :::::" << '\n';
    //      std::cin >> bre;
    cod.compute(lbc1.transpose());
    Ptemp = cod.colsPermutation();
    lbc1 = Ptemp * lbc1;
    lbc = lbc1.block(0, 0, cod.rank(), lbc1.cols());

    solutionR.adjointInPlace();
    Aplus.adjointInPlace();
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>
        rbc1 = solutionR * Aplus,
        Ptemp2, rbc;
    cod.compute(rbc1.transpose());
    Ptemp2 = cod.colsPermutation();
    rbc1 = Ptemp2 * rbc1;
    rbc = rbc1.block(0, 0, cod.rank(), rbc1.cols());

    // The lbc and rbc are now in the basis u1, u2, express in the normal basis:
    BcMat<std::complex<T> > out;
    out.resize(lbc.rows() + rbc.rows(), Lmat.r);
    ///     std::cout << "lbc: \n " << lbc << '\n';
    ///     std::cout << "rbc: \n " << rbc << '\n';
    ///     std::cin >> bre;
    for (int i = 0; i < lbc.rows(); i++) {
      for (int j = 0; j < Lmat.c; j++) {
        out.L(i, j).set(n - 1);
        out.eval(i, j) = -1.0;
        for (int k = 0; k < n; k++) {
          out.L(i, j).coef[n - k - 1] = lbc(i, Lmat.c * k + j);
        }
      }
    }

    for (int i = 0; i < rbc.rows(); i++) {
      for (int j = 0; j < Lmat.c; j++) {
        out.L(i + lbc.rows(), j).set(n - 1);
        out.eval(i + lbc.rows(), j) = 1.0;
        for (int k = 0; k < n; k++) {
          /// std::cout << "(" << i << "," << j << "," << k << ")" << '\n';
          out.L(i + lbc.rows(), j).coef[n - k - 1] = rbc(i, Lmat.c * k + j);
        }
      }
    }

    // This is to filter out higher derivatives in the boundary conditions,
    BcMat<std::complex<T> > tempMat(out.m, out.n);
    for (int i = 0; i < out.m; i++) {
      for (int j = 0; j < out.n; j++) {
        tempMat.L(i, j).set(highest_each_columnL[j] - 1);
        for (int k = 0; k < highest_each_columnL[j]; k++) {
          tempMat.L(i, j).coef[highest_each_columnL[j] - k - 1] =
              out.L(i, j).coef[out.L(i, j).n - k];
        }
      }
    }
    tempMat.eval = out.eval;

    //    std::cout << "out: \n" << '\n';
    //    for (int i = 0; i < out.m; i++) {
    //      for (int j = 0; j < out.n; j++) {
    //        std::cout << "(" << i << "," << j << ")" << '\n';
    //        std::cout << out.L(i, j).coef << '\n';
    //      }
    //    }
    //    std::cout << "out.eval: \n" << out.eval << '\n';
    //    std::cout << "in " << __LINE__ << '\n';

    // return the expression for boundary conditions:
    return tempMat;
  };
};

/// \brief Exports all data of 2D3C to a file. There are two dimensions, the
/// wall-normal and the spanwise, while there are three components,
/// u, v, w, and also pressure. vtk file can be directly exported into
/// paraview, or into Python using meshio
/// Supply filename without .vtk extension, that will be automatically added.
void vtkExportCartesian2D3C(const std::string &flnm,
                            const std::valarray<double> &y,
                            const std::valarray<double> &z,
                            const std::valarray<double> &u,
                            const std::valarray<double> &v,
                            const std::valarray<double> &w,
                            const std::valarray<double> &p) {
  std::string filename;
  filename = flnm + ".vtk";
  std::ofstream outf(filename);
  int Nz = z.size();
  int Ny = y.size();
  outf << "# vtk DataFile Version 2.0\n"
       << "2D3C velocity vector data \n"
       << "ASCII" << '\n'
       << "DATASET RECTILINEAR_GRID\n"
       << "DIMENSIONS " << 1 << " " << Ny << " " << Nz << "\n"
       << "X_COORDINATES " << 1 << " double\n";
  outf << 0.0 << " "
       << "\n"
       << "Y_COORDINATES " << Ny << " double\n";
  for (int i = 0; i < Ny; i++) {
    outf << y[i] << " ";
  }
  outf << "\n"
       << "Z_COORDINATES " << Nz << " double\n";
  for (int i = 0; i < Nz; i++) {
    outf << z[i] << " ";
  }
  outf << "\n\n"
       << "POINT_DATA " << Nz * (Ny) << "\n"
       << "SCALARS pressure double 1\n"
       << "LOOKUP_TABLE default\n";
  for (int j = 0; j < Nz; j++) {
    for (int i = 0; i < Ny; i++) {
      outf << p[(Ny)*i + j] << "\t";
    }
    outf << "\n";
  }
  outf << "\n\n"
       << "VECTORS Velocity double\n";
  for (int j = 0; j < Nz; j++) {
    for (int i = 0; i < Ny; i++) {
      outf << u[Nz * i + j] << " " << v[Nz * i + j] << " " << w[Nz * i + j]
           << "\n";
    }
    outf << "\n";
  }

  outf.close();
}

/// \brief Exports all data of 3D velocity to a file.
/// VTK file can be directly exported into
/// paraview, or into Python using meshio.
/// Supply filename without .vtk extension, that will be automatically added.
void vtkExportCartesian3D(
    const std::string &flnm, const std::valarray<double> &x,
    const std::valarray<double> &y, const std::valarray<double> &z,
    const std::valarray<double> &u, const std::valarray<double> &v,
    const std::valarray<double> &w, const std::valarray<double> &p) {
  std::string filename;
  filename = flnm + ".vtk";
  std::ofstream outf(filename);
  int Nz = z.size();
  int Nx = x.size();
  int Ny = y.size();
  outf << "# vtk DataFile Version 2.0\n"
       << "3D velocity vector data \n"
       << "ASCII" << '\n'
       << "DATASET RECTILINEAR_GRID\n"
       << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n"
       << "X_COORDINATES " << Nx << " double\n";
  for (int i = 0; i < Nx; i++) {
    outf << x[i] << " ";
  }
  outf << "\n"
       << "Y_COORDINATES " << Ny << " double\n";
  for (int i = 0; i < Ny; i++) {
    outf << y[i] << " ";
  }
  outf << "\n"
       << "Z_COORDINATES " << Nz << " double\n";
  for (int i = 0; i < Nz; i++) {
    outf << z[i] << " ";
  }
  outf << "\n\n"
       << "POINT_DATA " << Nx * Nz * (Ny) << "\n"
       << "SCALARS pressure double 1\n"
       << "LOOKUP_TABLE default\n";
  for (int k = 0; k < Nz; k++) {
    for (int i = 0; i < Ny; i++) {
      for (int j = 0; j < Nx; j++) {
        outf << p[i * Nx * Nz + j * Nz + k] << "\t";
      }
    }
    outf << "\n";
  }
  outf << "\n\n"
       << "VECTORS Velocity double\n";
  for (int k = 0; k < Nz; k++) {
    for (int i = 0; i < Ny; i++) {
      for (int j = 0; j < Nx; j++) {
        outf << u[i * Nx * Nz + j * Nz + k] << " "
             << v[i * Nx * Nz + j * Nz + k] << " "
             << w[i * Nx * Nz + j * Nz + k] << "\n";
      }
    }
    outf << "\n";
  }

  outf.close();
}

/// \brief Exports stress data to a file. VTK file can be directly exported into
/// paraview, or into Python using meshio.
/// Supply filename without .vtk extension, that will be automatically added.
void vtkExportCartesianStress3D(
    const std::string &flnm, const std::valarray<double> &x,
    const std::valarray<double> &y, const std::valarray<double> &z,
    const std::valarray<double> &t11, const std::valarray<double> &t12,
    const std::valarray<double> &t13, const std::valarray<double> &t22,
    const std::valarray<double> &t23, const std::valarray<double> &t33) {
  std::string filename;
  filename = flnm + ".vtk";
  std::ofstream outf(filename);
  int Nz = z.size();
  int Nx = x.size();
  int Ny = y.size();
  outf << "# vtk DataFile Version 2.0\n"
       << "3D velocity vector data \n"
       << "ASCII" << '\n'
       << "DATASET RECTILINEAR_GRID\n"
       << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n"
       << "X_COORDINATES " << Nx << " double\n";
  for (int i = 0; i < Nx; i++) {
    outf << x[i] << " ";
  }
  outf << "\n"
       << "Y_COORDINATES " << Ny << " double\n";
  for (int i = 0; i < Ny; i++) {
    outf << y[i] << " ";
  }
  outf << "\n"
       << "Z_COORDINATES " << Nz << " double\n";
  for (int i = 0; i < Nz; i++) {
    outf << z[i] << " ";
  }
  outf << "\n\n"
       << "POINT_DATA " << Nx * Nz * (Ny) << "\n"
       << "TENSORS Stress double\n";
  for (int k = 0; k < Nz; k++) {
    for (int i = 0; i < Ny; i++) {
      for (int j = 0; j < Nx; j++) {
        outf << t11[i * Nx * Nz + j * Nz + k] << " "
             << t12[i * Nx * Nz + j * Nz + k] << " "
             << t13[i * Nx * Nz + j * Nz + k] << "\n"
             << t12[i * Nx * Nz + j * Nz + k] << " "
             << t22[i * Nx * Nz + j * Nz + k] << " "
             << t23[i * Nx * Nz + j * Nz + k] << "\n"
             << t13[i * Nx * Nz + j * Nz + k] << " "
             << t23[i * Nx * Nz + j * Nz + k] << " "
             << t33[i * Nx * Nz + j * Nz + k] << "\n\n";
      }
    }
  }

  outf.close();
}

template <class T> Linop<T> conj(Linop<T> in) {
  if (in.NCC == 0) {
    in.coef = in.coef.conjugate().eval();
  } else {
    for (int i = 0; i < in.n + 1; i++) {
      in.coefFun[i] = conj(in.coefFun[i]);
    }
  }
  return in;
};
/// \brief Defining addition of ChebfunMats.
template <class T> ChebfunMat<T> operator+(ChebfunMat<T> a, ChebfunMat<T> b) {
  ChebfunMat<T> out;
  int bre;
  /// std::cout << "arc: " << b.r << "," << b.c << '\n';
  if (a.r != b.r || a.c != b.c) {
    std::cout << "addition of vector of ChebfunMats not possible"
              << ". In " << __LINE__ << '\n';
    exit(1);
  } else {
    out.resize(a.r, a.c);
    for (int i = 0; i < out.r; i++) {
      for (int j = 0; j < out.c; j++) {
        // std::cout << "a("<< i << ","<<j << ")" << '\n';
        // std::cout << a(i,j) << '\n';
        // std::cout << "b("<< i << ","<<j << ")" << '\n';
        // std::cout << b(i,j) << '\n';
        // std::cin >> bre;
        out(i, j) = a(i, j) + b(i, j);
      }
    }
  }
  return out;
};

/// \brief Defining addition of vector of ChebfunMats. Used for eigenvectors,
/// Singular values, adjoints etc.
template <class T>
std::vector<ChebfunMat<T> > operator+(std::vector<sis::ChebfunMat<T> > a,
                                      std::vector<sis::ChebfunMat<T> > b) {
  if (a.size() != b.size()) {
    std::cout << "addition of vector of ChebfunMats not possible"
              << ". In " << __LINE__ << '\n';
    exit(1);
  }
  std::vector<sis::ChebfunMat<T> > out;
  out.resize(a.size());
  for (int i = 0; i < a.size(); i++) {
    if (a[i].r != b[i].r || a[i].c != b[i].c) {
      std::cout << "addition of vector of ChebfunMats not possible"
                << ". In " << __LINE__ << '\n';
      exit(1);
    } else {
      out[i] = a[i] + b[i];
    }
  }
  return out;
};

/// \brief This function overloads std::cout<< to valarrays.
template <class T>
std::ostream &operator<<(std::ostream &stream, std::valarray<T> a) {
  stream << "\n";
  for (int i = 0; i < a.size(); i++) {
    stream << a[i] << "\n";
  };
  return stream;
};

template <class T>
std::ostream &operator<<(std::ostream &stream, std::vector<T> a) {
  stream << "\n";
  for (int i = 0; i < a.size(); i++) {
    stream << a[i] << "\n";
  };
  return stream;
};

// Minus of a Linop:
template <class T> Linop<T> operator-(const Linop<T> &in) {
  Linop<T> out;
  out = in;
  if (out.NCC == 0) {
    out.coef = -out.coef;
  } else {
    for (int i = 0; i < out.n + 1; i++) {
      out.coefFun[i].v = -out.coefFun[i].v;
    }
  }
  return out;
};

// Minus of a Chebfun:
template <class T> Chebfun<T> operator-(const Chebfun<T> &in) {
  Chebfun<T> out;
  out = in;
  out.v = -out.v;
  return out;
};

// Addition of two Linop:
template <class T> Linop<T> operator+(Linop<T> l_, Linop<T> r_) {
  Linop<T> out;
  if (l_.NCC == 0) {
    if (r_.NCC == 0) {
      if (l_.n > r_.n) {
        int diffn = l_.n - r_.n;
        out.set(l_.n);
        for (int i = r_.n; i > -1; i--) {
          out.coef[i + diffn] = l_.coef[i + diffn] + r_.coef[i];
        }
        for (int i = 0; i < diffn; i++) {
          out.coef[i] = l_.coef[i];
        }
      } else {
        int diffn = r_.n - l_.n;
        out.set(r_.n);
        for (int i = l_.n; i > -1; i--) {
          out.coef[i + diffn] = r_.coef[i + diffn] + l_.coef[i];
        }
        for (int i = 0; i < diffn; i++) {
          out.coef[i] = r_.coef[i];
        }
      }
    } else {
      if (l_.n > r_.n) {
        //  std::cout << "l.n ,r.n" << l_.n << "," << r_.n << '\n';
        int diffn = l_.n - r_.n;
        out.ncc(l_.n);
        for (int i = r_.n; i > -1; i--) {
          out.coefFun[i + diffn].v = l_.coef[i + diffn] + r_.coefFun[i].v;
        }
        for (int i = 0; i < diffn; i++) {
          out.coefFun[i].v = l_.coef[i];
        }
      } else {
        int diffn = r_.n - l_.n;
        out.ncc(r_.n);
        for (int i = l_.n; i > -1; i--) {
          out.coefFun[i + diffn].v = r_.coefFun[i + diffn].v + l_.coef[i];
        }
        for (int i = 0; i < diffn; i++) {
          out.coefFun[i].v = r_.coefFun[i].v;
        }
      }
    }
  } else {
    if (r_.NCC == 0) {
      if (l_.n > r_.n) {
        int diffn = l_.n - r_.n;
        out.ncc(l_.n);
        for (int i = r_.n; i > -1; i--) {
          out.coefFun[i + diffn].v = l_.coefFun[i + diffn].v + r_.coef[i];
        }
        for (int i = 0; i < diffn; i++) {
          out.coefFun[i].v = l_.coefFun[i].v;
        }
      } else {
        int diffn = r_.n - l_.n;
        out.ncc(r_.n);
        for (int i = l_.n; i > -1; i--) {
          out.coefFun[i + diffn].v = r_.coef[i + diffn] + l_.coefFun[i].v;
        }
        for (int i = 0; i < diffn; i++) {
          out.coefFun[i].v = r_.coef[i];
        }
      }
    } else {
      if (l_.n > r_.n) {
        int diffn = l_.n - r_.n;
        out.ncc(l_.n);
        for (int i = r_.n; i > -1; i--) {
          out.coefFun[i + diffn].v = l_.coefFun[i + diffn].v + r_.coefFun[i].v;
        }
        for (int i = 0; i < diffn; i++) {
          out.coefFun[i].v = l_.coefFun[i].v;
        }
      } else {
        int diffn = r_.n - l_.n;
        out.ncc(r_.n);
        for (int i = l_.n; i > -1; i--) {
          out.coefFun[i + diffn].v = r_.coefFun[i + diffn].v + l_.coefFun[i].v;
        }
        for (int i = 0; i < diffn; i++) {
          out.coefFun[i].v = r_.coefFun[i].v;
        }
      }
    }
  }
  return out;
};

/// \brief Converts real and imaginary parts of a Linop to a complex Linop.
template <class T>
Linop<std::complex<T> > dou2com(Linop<T> real, Linop<T> imag) {
  Linop<std::complex<T> > out;
  if (real.NCC == 0) {
    if (imag.NCC == 0) {
      if (real.n > imag.n) {
        int diffn = real.n - imag.n;
        out.set(real.n);
        for (int i = imag.n; i > -1; i--) {
          out.coef[i + diffn] =
              std::complex<T>(real.coef[i + diffn], imag.coef[i]);
        }
        for (int i = 0; i < diffn; i++) {
          out.coef[i] = real.coef[i];
        }
      } else {
        int diffn = imag.n - real.n;
        out.set(imag.n);
        for (int i = real.n; i > -1; i--) {
          out.coef[i + diffn] =
              std::complex<T>(real.coef[i], imag.coef[i + diffn]);
        }
        for (int i = 0; i < diffn; i++) {
          out.coef[i] = std::complex<T>(0.0, imag.coef[i]);
        }
      }
    } else {
      if (real.n > imag.n) {
        int diffn = real.n - imag.n;
        out.ncc(real.n);
        for (int i = imag.n; i > -1; i--) {
          std::valarray<T> temp(N + 1);
          temp = real.coef[i + diffn];
          out.coefFun[i + diffn].v = dou2com(temp, imag.coefFun[i].v);
        }
        for (int i = 0; i < diffn; i++) {
          out.coefFun[i].v = std::complex<T>(real.coef[i], 0.0);
        }
      } else {
        int diffn = imag.n - real.n;
        out.ncc(imag.n);
        for (int i = real.n; i > -1; i--) {
          std::valarray<T> temp(N + 1);
          temp = real.coef[i];
          out.coefFun[i + diffn].v = dou2com(temp, imag.coefFun[i + diffn].v);
        }
        for (int i = 0; i < diffn; i++) {
          std::valarray<T> temp(N + 1);
          temp = 0.0;
          out.coefFun[i].v = dou2com(temp, imag.coefFun[i].v);
        }
      }
    }
  } else {
    if (imag.NCC == 0) {
      if (real.n > imag.n) {
        int diffn = real.n - imag.n;
        out.ncc(real.n);
        for (int i = imag.n; i > -1; i--) {
          std::valarray<T> temp(N + 1);
          temp = imag.coef[i];
          out.coefFun[i + diffn].v = dou2com(real.coefFun[i + diffn].v, temp);
        }
        for (int i = 0; i < diffn; i++) {
          std::valarray<T> temp(N + 1);
          temp = 0.0;
          out.coefFun[i].v = dou2com(real.coefFun[i].v, temp);
        }
      } else {
        int diffn = imag.n - real.n;
        out.ncc(imag.n);
        for (int i = real.n; i > -1; i--) {
          std::valarray<T> temp(N + 1);
          temp = imag.coef[i + diffn];
          out.coefFun[i + diffn].v = dou2com(real.coefFun[i].v, temp);
        }
        for (int i = 0; i < diffn; i++) {
          out.coefFun[i].v = std::complex<T>(0.0, imag.coef[i]);
        }
      }
    } else {
      if (real.n > imag.n) {
        int diffn = real.n - imag.n;
        out.ncc(real.n);
        for (int i = imag.n; i > -1; i--) {
          out.coefFun[i + diffn].v =
              dou2com(real.coefFun[i + diffn].v, imag.coefFun[i].v);
        }
        for (int i = 0; i < diffn; i++) {
          std::valarray<T> temp(N + 1);
          temp = 0.0;
          out.coefFun[i].v = dou2com(real.coefFun[i].v, temp);
        }
      } else {
        int diffn = imag.n - real.n;
        out.ncc(imag.n);
        for (int i = real.n; i > -1; i--) {
          out.coefFun[i + diffn].v =
              dou2com(real.coefFun[i].v, imag.coefFun[i + diffn].v);
        }
        for (int i = 0; i < diffn; i++) {
          std::valarray<T> temp(N + 1);
          temp = 0.0;
          out.coefFun[i].v = dou2com(temp, imag.coefFun[i].v);
        }
      }
    }
  }
  return out;
};

/// \brief Addition of two Linops, complex with real:
template <class T>
Linop<std::complex<T> > operator+(Linop<std::complex<T> > l_, Linop<T> r_) {
  Linop<T> outr, outi;
  Linop<std::complex<T> > out;
  outr = l_.real() + r_;
  outi = l_.imag();
  out = dou2com(outr, outi);
  return out;
};

/// \brief Addition of two Linops, real with complex:
template <class T>
Linop<std::complex<T> > operator+(Linop<T> l_, Linop<std::complex<T> > r_) {
  Linop<T> outr, outi;
  Linop<std::complex<T> > out;
  outr = l_ + r_.real();
  outi = r_.imag();
  out = dou2com(outr, outi);
  return out;
};

/// \brief Addition of Linop to Chebfun:
template <class T> Linop<T> operator+(Linop<T> l_, Chebfun<T> r_) {
  Linop<T> temp;
  temp = r_;
  return l_ + temp;
};

/// \brief Addition of two Linop with Chebfun:
template <class T>
Linop<std::complex<T> > operator+(Linop<std::complex<T> > l_,
                                  Chebfun<std::complex<T> > r_) {
  Linop<std::complex<T> > temp;
  temp = r_;
  return l_ + temp;
};

/// \brief Addition of two Linop with Chebfun:
template <class T>
Linop<std::complex<T> > operator+(Linop<T> l_, Chebfun<std::complex<T> > r_) {

  Linop<T> outr, outi;
  outr = l_ + r_.real();
  outi = r_.imag();
  return dou2com(outr, outi);
};

/// \brief Addition of Linop with Chebfun:
template <class T>
Linop<std::complex<T> > operator+(Linop<std::complex<T> > l_, Chebfun<T> r_) {
  Linop<T> outr, outi;
  outr = l_.real() + r_;
  outi = l_.imag();
  return dou2com(l_, r_);
};

/// \brief Addition of Chebfun to Linop:
template <class T> Linop<T> operator+(Chebfun<T> l_, Linop<T> r_) {
  Linop<T> temp;
  temp = l_;
  return r_ + temp;
};

/// \brief Addition of Chebfun to Linop:
template <class T>
Linop<std::complex<T> > operator+(Chebfun<std::complex<T> > l_,
                                  Linop<std::complex<T> > r_) {
  Linop<std::complex<T> > temp;
  temp = l_;
  return temp + r_;
};

/// \brief Addition of Chebfun to Linop:
template <class T>
Linop<std::complex<T> > operator+(Chebfun<T> l_, Linop<std::complex<T> > r_) {
  Linop<T> outr, outi;
  outr = l_ + r_.real();
  outi = r_.imag();
  return dou2com(outr, outi);
};
/// \brief Addition of Chebfun to Linop:
template <class T>
Linop<std::complex<T> > operator+(Chebfun<std::complex<T> > l_, Linop<T> r_) {
  Linop<T> outr, outi;
  outr = l_.real() + r_;
  outi = l_.imag();
  return dou2com(outr, outi);
};

/// \brief Addition of Linop to Chebfun:
template <class T> Linop<T> operator+(Linop<T> l_, std::valarray<T> r_) {
  Linop<T> temp;
  temp = r_;
  return l_ + temp;
};

/// \brief Addition of two Linop with Chebfun:
template <class T>
Linop<std::complex<T> > operator+(Linop<std::complex<T> > l_,
                                  std::valarray<std::complex<T> > r_) {
  Linop<std::complex<T> > temp;
  temp = r_;
  return l_ + temp;
};

/// \brief Addition of two Linop with Chebfun:
template <class T>
Linop<std::complex<T> > operator+(Linop<T> l_,
                                  std::valarray<std::complex<T> > r_) {

  Linop<std::complex<T> > temp1;
  Linop<std::complex<T> > temp2;
  temp1 = l_;
  temp2 = r_;
  return temp1 + temp2;
};

/// \brief Addition of Linop with Chebfun:
template <class T>
Linop<std::complex<T> > operator+(Linop<std::complex<T> > l_,
                                  std::valarray<T> r_) {
  Linop<T> outr, outi;
  outr = l_.real() + r_;
  outi = l_.imag();
  return dou2com(outr, outi);
};

/// \brief Addition of Chebfun to Linop:
template <class T> Linop<T> operator+(std::valarray<T> l_, const Linop<T> &r_) {
  Linop<T> temp;
  temp = l_;
  return r_ + temp;
};

/// \brief Addition of Chebfun to Linop:
template <class T>
Linop<std::complex<T> > operator+(const std::valarray<std::complex<T> > &l_,
                                  const Linop<std::complex<T> > &r_) {
  Linop<std::complex<T> > temp;
  temp = l_;
  return temp + r_;
};

/// \brief Addition of Chebfun to Linop:
template <class T>
Linop<std::complex<T> > operator+(std::valarray<T> l_,
                                  Linop<std::complex<T> > r_) {
  return dou2com(l_, r_);
};
/// \brief Addition of Chebfun to Linop:
template <class T>
Linop<std::complex<T> > operator+(std::valarray<std::complex<T> > l_,
                                  Linop<T> r_) {
  Linop<T> outr, outi;
  outr = real(l_) + r_;
  outi = imag(l_);
  return dou2com(outr, outi);
};

/// \brief Subtraction for two Linop:
template <class T> Linop<T> operator-(Linop<T> l_, Linop<T> r_) {
  return l_ + (-r_);
};

/// \brief Subtraction of two Linops, complex with real:
template <class T>
Linop<std::complex<T> > operator-(Linop<std::complex<T> > l_, Linop<T> r_) {
  return l_ + (-r_);
};

/// \brief Subtraction of two Linops, real with complex:
template <class T>
Linop<std::complex<T> > operator-(Linop<T> l_, Linop<std::complex<T> > r_) {
  return l_ + (-r_);
};

/// \brief Subtraction of Linop to Chebfun:
template <class T> Linop<T> operator-(Linop<T> l_, Chebfun<T> r_) {
  return l_ + (-r_);
};

/// \brief Subtraction of two Linop with Chebfun:
template <class T>
Linop<std::complex<T> > operator-(Linop<std::complex<T> > l_,
                                  Chebfun<std::complex<T> > r_) {
  return l_ + (-r_);
};

/// \brief Subtraction of two Linop with Chebfun:
template <class T>
Linop<std::complex<T> > operator-(Linop<T> l_, Chebfun<std::complex<T> > r_) {
  return l_ + (-r_);
};

/// \brief Subtraction of Linop with Chebfun:
template <class T>
Linop<std::complex<T> > operator-(Linop<std::complex<T> > l_, Chebfun<T> r_) {
  return l_ + (-r_);
};

/// \brief Subtraction of Chebfun to Linop:
template <class T> Linop<T> operator-(Chebfun<T> l_, Linop<T> r_) {
  return l_ + (-r_);
};

/// \brief Subtraction of Chebfun to Linop:
template <class T>
Linop<std::complex<T> > operator-(Chebfun<std::complex<T> > l_,
                                  Linop<std::complex<T> > r_) {
  return l_ + (-r_);
};

/// \brief Subtraction of Chebfun to Linop:
template <class T>
Linop<std::complex<T> > operator-(Chebfun<T> l_, Linop<std::complex<T> > r_) {
  return l_ + (-r_);
};
/// \brief Subtraction of Chebfun to Linop:
template <class T>
Linop<std::complex<T> > operator-(Chebfun<std::complex<T> > l_, Linop<T> r_) {
  return l_ + (-r_);
};

/// \brief Subtraction of Linop to Chebfun:
template <class T> Linop<T> operator-(Linop<T> l_, std::valarray<T> r_) {
  return l_ + std::valarray<T>(-r_);
};

/// \brief Subtraction of two Linop with Chebfun:
template <class T>
Linop<std::complex<T> > operator-(Linop<T> l_,
                                  std::valarray<std::complex<T> > r_) {
  return l_ + std::valarray<std::complex<T> >(-r_);
};

/// \brief Subtraction of Linop with Chebfun:
template <class T>
Linop<std::complex<T> > operator-(Linop<std::complex<T> > l_,
                                  std::valarray<T> r_) {
  return l_ + std::valarray<T>(-r_);
};

/// \brief Subtraction of Chebfun to Linop:
template <class T>
Linop<T> operator-(const std::valarray<T> &l_, const Linop<T> &r_) {
  return l_ + (-r_);
};

/// \brief Subtraction of Chebfun to Linop:
template <class T>
Linop<std::complex<T> > operator-(std::valarray<T> l_,
                                  Linop<std::complex<T> > r_) {
  return l_ + (-r_);
};
/// \brief Subtraction of Chebfun to Linop:
template <class T>
Linop<std::complex<T> > operator-(std::valarray<std::complex<T> > l_,
                                  Linop<T> r_) {
  return l_ + (-r_);
};

template <class T>
Eigen::Array<T, Eigen::Dynamic, 1>
operator*(std::valarray<T> left,
          const Eigen::Array<T, Eigen::Dynamic, 1> &right) {
  Eigen::Array<T, Eigen::Dynamic, 1> out;
  out.resize(left.size());
  for (int i = 0; i < left.size(); i++) {
    out[i] = left[i] * right[i];
  }
  return out;
};

template <class T> Linop<T> operator*(std::valarray<T> left, Linop<T> right) {
  Linop<T> temp;
  temp = right;
  if (right.NCC == 0) {
    temp.ncc();
    for (int i = 0; i < right.n + 1; i++) {
      temp.coefFun[i].v = left * right.coef[i];
    }
  } else {
    temp.ncc();
    for (int i = 0; i < right.n + 1; i++) {
      temp.coefFun[i].v = left * right.coefFun[i].v;
    }
  }
  return temp;
};

template <class T>
Linop<std::complex<T> > operator*(std::valarray<std::complex<T> > left,
                                  Linop<T> right) {
  Linop<std::complex<T> > temp;
  temp = right;
  return left * temp;
};

template <class T>
Linop<std::complex<T> > operator*(std::valarray<T> left,
                                  Linop<std::complex<T> > right) {
  std::valarray<std::complex<T> > temp;
  std::valarray<T> zvec;
  zvec = 0.0;
  temp = dou2com(temp, zvec);
  return temp * right;
};

/////////
template <class T> Linop<T> operator*(Chebfun<T> left, Linop<T> right) {
  Linop<T> temp;
  temp = right;
  if (right.NCC == 0) {
    temp.ncc();
    for (int i = 0; i < right.n + 1; i++) {
      temp.coefFun[i] = left * right.coef[i];
    }
  } else {
    temp.ncc();
    for (int i = 0; i < right.n + 1; i++) {
      temp.coefFun[i] = left * right.coefFun[i];
    }
  }
  return temp;
};

template <class T>
Linop<std::complex<T> > operator*(Chebfun<std::complex<T> > left,
                                  Linop<T> right) {
  Linop<std::complex<T> > temp;
  temp = right;
  return left * temp;
};

template <class T>
Linop<std::complex<T> > operator*(Chebfun<T> left,
                                  Linop<std::complex<T> > right) {
  Chebfun<std::complex<T> > temp;
  temp = left;
  return temp * right;
};

////////
template <class T> Linop<T> operator*(T left, Linop<T> right) {
  Linop<T> temp;
  temp = right;
  if (right.NCC == 0) {
    for (int i = 0; i < right.n + 1; i++) {
      temp.coef[i] = left * right.coef[i];
    }
  } else {
    temp.ncc();
    for (int i = 0; i < right.n + 1; i++) {
      temp.coefFun[i].v = left * right.coefFun[i].v;
    }
  }
  return temp;
};

template <class T>
Linop<std::complex<T> > operator*(std::complex<T> left, Linop<T> right) {
  Linop<std::complex<T> > temp;
  temp = right;
  return left * temp;
};

template <class T>
Linop<std::complex<T> > operator*(T left, Linop<std::complex<T> > right) {
  std::complex<T> temp(left, 0.0);
  return temp * right;
};

/// \brief Adding scalar to Linop.
template <class T> Linop<T> operator+(const Linop<T> &left, T right) {
  Linop<T> _right;
  _right = right;
  return _right + left;
};

/// \brief Adding Linop to scalar.
template <class T> Linop<T> operator+(T left, const Linop<T> &right) {
  Linop<T> _left;
  _left = left;
  return _left + right;
};

/// \brief Adding scalar to Linop complex.
template <class T>
Linop<std::complex<T> > operator+(const Linop<std::complex<T> > &left,
                                  T right) {
  Linop<std::complex<T> > _right;
  _right = right;
  return _right + left;
};

/// \brief Adding Linop to scalar.
template <class T>
Linop<std::complex<T> > operator+(std::complex<T> left, const Linop<T> &right) {
  Linop<std::complex<T> > _left;
  _left = left;
  return _left + right;
};

/// \brief Adding scalar to Linop complex.
template <class T>
Linop<std::complex<T> > operator+(const Linop<T> &left, std::complex<T> right) {
  Linop<std::complex<T> > _right;
  _right = right;
  return _right + left;
};

/// \brief Adding Linop to scalar.
template <class T>
Linop<std::complex<T> > operator+(T left,
                                  const Linop<std::complex<T> > &right) {
  Linop<std::complex<T> > _left;
  _left = left;
  return _left + right;
};
/// \brief Subracting scalar from Linop.
template <class T> Linop<T> operator-(const Linop<T> &left, T right) {
  Linop<T> _right;
  _right = right;
  return left - _right;
};

/// \brief Adding Linop to scalar.
template <class T> Linop<T> operator-(T left, const Linop<T> &right) {
  Linop<T> _left;
  ;
  _left = left;
  return _left - right;
};

/// \brief Adding Linop to scallar
template <class T>
Linop<std::complex<T> > operator-(T left,
                                  const Linop<std::complex<T> > &right) {
  Linop<std::complex<T> > _left;
  _left = left;
  return _left - right;
};

template <class T>
Linop<std::complex<T> > operator-(std::complex<T> left, const Linop<T> &right) {
  Linop<std::complex<T> > _left;
  _left = left;
  return _left - right;
};

template <class T>
Linop<std::complex<T> > operator-(const Linop<std::complex<T> > &left,
                                  T right) {
  Linop<std::complex<T> > _right;
  _right = right;
  return left - _right;
};

template <class T>
Linop<std::complex<T> > operator-(const Linop<T> &left, std::complex<T> right) {
  Linop<std::complex<T> > _right;
  _right = right;
  return left - _right;
};

/// \brief Multiplying Chebfun to a constant:
template <class T> Chebfun<T> operator*(T a, Chebfun<T> b) {
  Chebfun<T> out;
  out.v = a * b.v;
  out.dct_flag = b.dct_flag;
  return out;
};

/// \brief Multiplying Chebfun to a constant:
template <class T> Chebfun<T> operator*(Chebfun<T> b, T a) {
  Chebfun<T> out;
  out.v = a * b.v;
  out.dct_flag = b.dct_flag;
  return out;
};

/// \brief Multiplying Chebfun to a constant:
template <class T>
Chebfun<std::complex<T> > operator*(Chebfun<std::complex<T> > b, T a) {
  Chebfun<std::complex<T> > out;
  out.fftw_init();
  out.v = a * b.v;
  out.dct_flag = b.dct_flag;
  return out;
};

/// \brief Multiplying Chebfun to a constant:
template <class T>
Chebfun<std::complex<T> > operator*(T a, Chebfun<std::complex<T> > b) {
  Chebfun<std::complex<T> > out;
  out.v = a * b.v;
  out.dct_flag = b.dct_flag;
  return out;
};

/// \brief Multiplying ChebfunMat to a constant, all Chebfuns in the ChebfunMat
/// is multiplied by the constant.
template <class T> ChebfunMat<T> operator*(T a, ChebfunMat<T> in) {
  ChebfunMat<T> out;
  out.resize(in.r, in.c);
  int bre;

  for (int i = 0; i < in.r; i++) {
    for (int j = 0; j < in.c; j++) {
      out(i, j) = a * in(i, j);
    }
  }

  return out;
};

/// \brief Multiplying ChebfunMat to a constant, all Chebfuns in the ChebfunMat
/// is multiplied by the constant.
template <class T>
ChebfunMat<std::complex<T> > operator*(T a, ChebfunMat<std::complex<T> > in) {
  ChebfunMat<std::complex<T> > out;
  out.resize(in.r, in.c);
  for (int i = 0; i < in.r; i++) {
    for (int j = 0; j < in.c; j++) {
      out(i, j) = a * in(i, j);
    }
  }
  return out;
};

/// \brief Multiplying ChebfunMat to a constant, all Chebfuns in the ChebfunMat
/// is multiplied by the constant.
template <class T>
ChebfunMat<std::complex<T> > operator*(std::complex<T> a, ChebfunMat<T> in) {
  ChebfunMat<std::complex<T> > out;
  out.resize(in.r, in.c);
  for (int i = 0; i < in.r; i++) {
    for (int j = 0; j < in.c; j++) {
      out(i, j) = a * in(i, j);
    }
  }
};

//////////////
/// \brief Multiplying a vector of ChebfunMats to a constant, all Chebfuns in
/// all ChebfunMats
/// are multiplied by the constant.
template <class T>
std::valarray<ChebfunMat<T> > operator*(T a, std::valarray<ChebfunMat<T> > in) {
  std::valarray<ChebfunMat<T> > out;
  out.resize(in.size());
  for (int i = 0; i < in.size(); i++) {
    out[i] = a * in[i];
  }
};

/// \brief Multiplying a vector of ChebfunMats to a constant, all Chebfuns in
/// all ChebfunMats
/// are multiplied by the constant.
template <class T>
std::valarray<ChebfunMat<std::complex<T> > >
operator*(std::complex<T> a, std::valarray<ChebfunMat<T> > in) {
  std::valarray<ChebfunMat<std::complex<T> > > out;
  out.resize(in.size());
  for (int i = 0; i < in.size(); i++) {
    out[i] = a * in[i];
  }
};

/// \brief Multiplying a vector of ChebfunMats to a constant, all Chebfuns in
/// all ChebfunMats
/// are multiplied by the constant.
template <class T>
std::valarray<ChebfunMat<T> >
operator*(T a, std::valarray<ChebfunMat<std::complex<T> > > in) {
  std::valarray<ChebfunMat<std::complex<T> > > out;
  out.resize(in.size());
  for (int i = 0; i < in.size(); i++) {
    out[i] = a * in[i];
  }
};

/////////////

/// \brief Multiplying two Chebfuns. If in same space, returns in same space,
/// else defauts to SIS_PHYS_SPACE.
template <class T> Chebfun<T> operator*(Chebfun<T> a, Chebfun<T> b) {
  Chebfun<T> out;
  if (a.dct_flag == b.dct_flag) {
    if (a.dct_flag == SIS_PHYS_SPACE) {
      out.v = a.v * b.v;
      out.dct_flag = SIS_PHYS_SPACE;
      return out;
    } else {
      a.c2p();
      b.c2p();
      out.v = a.v + b.v;
      out.dct_flag = a.dct_flag;
      a.p2c();
      b.p2c();
      out.p2c();
      return out;
    }
  } else {
    if (a.dct_flag == SIS_PHYS_SPACE) {
      b.c2p();
      out.v = a.v + b.v;
      out.dct_flag = SIS_PHYS_SPACE;
      b.p2c();
      return out;
    } else {
      a.c2p();
      out.v = a.v + b.v;
      out.dct_flag = SIS_PHYS_SPACE;
      a.p2c();
      return out;
    }
  }
};

template <class T>
Chebfun<std::complex<T> > operator*(Chebfun<T> a, Chebfun<std::complex<T> > b) {

  Chebfun<std::complex<T> > out, a_;
  a_ = a;
  return (a_ * b);
};

template <class T>
Chebfun<std::complex<T> > operator*(Chebfun<std::complex<T> > b, Chebfun<T> a) {
  Chebfun<std::complex<T> > out, a_;
  a_ = a;
  return (b * a_);
};

/// \relates  Linop<std::complex<T>>
///
/// \brief These are used to define operator compositions, L1(L2). It is
/// defined as operator multiplication, L1 * L2.
/// Based on,
/// \f{align}
/// \left[ \sum_{i=0}^{m}\,b_i(y)\,D^i \right] \,
/// \left[ \sum_{j=0}^{n}\,a_i(y)\,D^i \right]  \;=\;
/// \sum_{i = 0}^{m}\,b_i(y)\,
/// \sum_{j = 0}^{n}\sum_{k = 0}^{i}\,C^i_k\,a_j^k(y)\,Dj+i-k
/// \f}
/// Implemented by recursion, not by the formula above.
template <class T> Linop<T> operator*(Linop<T> L1, Linop<T> L2) {
  Linop<T> out;
  out.set(L1.n + L2.n);
  if (L1.NCC == 0 && L1.NCC == 0) {
    for (int i = 0; i < L1.n + L2.n + 1; i++) {
      out.coef[i] = 0.0;
    }
    for (int i = 0; i < L1.n + 1; i++) {
      out += L1.coef[i] * diff(L2, L1.n - i);
    }
  } else {
    out.ncc(L1.n + L2.n);
    for (int i = 0; i < L1.n + L2.n + 1; i++) {
      out.coefFun[i] = 0.0;
    }
    if (L1.NCC == 0) {
      for (int i = 0; i < L1.n + 1; i++) {
        out += L1.coef[i] * diff(L2, L1.n - i);
      }
    } else {
      for (int i = 0; i < L1.n + 1; i++) {
        out += L1.coefFun[i] * diff(L2, L1.n - i);
      }
    }
  }
  return out;
};

template <class T>
Linop<std::complex<T> > operator*(Linop<std::complex<T> > L1, Linop<T> L2) {
  Linop<std::complex<T> > out;
  out.set(L1.n + L2.n);
  if (L1.NCC == 0 && L1.NCC == 0) {
    for (int i = 0; i < L1.n + L2.n + 1; i++) {
      out.coef[i] = 0.0;
    }
    for (int i = 0; i < L1.n + 1; i++) {
      out += L1.coef[i] * diff(L2, L1.n - i);
    }
  } else {
    out.ncc(L1.n + L2.n);
    for (int i = 0; i < L1.n + L2.n + 1; i++) {
      out.coefFun[i] = 0.0;
    }
    if (L1.NCC == 0) {
      for (int i = 0; i < L1.n + 1; i++) {
        out += L1.coef[i] * diff(L2, L1.n - i);
      }
    } else {
      for (int i = 0; i < L1.n + 1; i++) {
        out += L1.coefFun[i] * diff(L2, L1.n - i);
      }
    }
  }
  return out;
};

template <class T> Linop<T> operator*(Linop<T> L1, Linop<std::complex<T> > L2) {
  Linop<std::complex<T> > out;
  out.set(L1.n + L2.n);
  if (L1.NCC == 0 && L1.NCC == 0) {
    for (int i = 0; i < L1.n + L2.n + 1; i++) {
      out.coef[i] = 0.0;
    }
    for (int i = 0; i < L1.n + 1; i++) {
      out += L1.coef[i] * diff(L2, L1.n - i);
    }
  } else {
    out.ncc(L1.n + L2.n);
    for (int i = 0; i < L1.n + L2.n + 1; i++) {
      out.coefFun[i] = 0.0;
    }
    if (L1.NCC == 0) {
      for (int i = 0; i < L1.n + 1; i++) {
        out += L1.coef[i] * diff(L2, L1.n - i);
      }
    } else {
      for (int i = 0; i < L1.n + 1; i++) {
        out += L1.coefFun[i] * diff(L2, L1.n - i);
      }
    }
  }
  return out;
};

/// \brief Addition for two chebfuns. If both functions are not in same space,
/// then default evaluation will be in physical space.
template <class T> Chebfun<T> operator+(Chebfun<T> a, Chebfun<T> b) {
  Chebfun<T> out;
  if (a.dct_flag == b.dct_flag) {
    out.v = a.v + b.v;
    out.dct_flag = a.dct_flag;
    return out;
  } else {
    if (a.dct_flag == SIS_PHYS_SPACE) {
      b.c2p();
      out.v = a.v + b.v;
      out.dct_flag = SIS_PHYS_SPACE;
      b.p2c();
      return out;
    } else {
      a.c2p();
      out.v = a.v + b.v;
      out.dct_flag = SIS_PHYS_SPACE;
      a.p2c();
      return out;
    }
  }
};

template <class T>
Chebfun<std::complex<T> > operator+(Chebfun<T> a, Chebfun<std::complex<T> > b) {
  Chebfun<std::complex<T> > out, a_;
  a_ = a;
  return (a_ + b);
};

template <class T>
Chebfun<std::complex<T> > operator+(Chebfun<std::complex<T> > a, Chebfun<T> b) {
  Chebfun<std::complex<T> > out, b_;
  b_ = b;
  return (a + b_);
};

template <class T> Linop<T> pow(Linop<T> in, int a) {
  if (a < 0) {
    std::cout << "Cannot have a negative power for linop ..." << '\n';
    std::cout << "in " << __LINE__ << '\n';
    exit(1);
  } else if (a == 1) {
    return in;
  } else {
    return pow(in, a - 1) * in;
  }
};

template <class T> LinopMat<T> operator+(LinopMat<T> left, LinopMat<T> right) {
  LinopMat<T> temp;
  if (left.r != right.r || left.c != right.c) {
    std::cout << "Incompatible shapes of LinopMats for addition." << '\n';
    exit(1);
  }
  temp.resize(left.r, left.c);
  for (int i = 0; i < left.r; i++) {
    for (int j = 0; j < left.c; j++) {
      temp(i, j) = left(i, j) + right(i, j);
    }
  }
  return temp;
};

template <class T>
LinopMat<std::complex<T> > operator+(LinopMat<std::complex<T> > left,
                                     LinopMat<T> right) {
  LinopMat<std::complex<T> > temp;
  if (left.r != right.r || left.c != right.c) {
    std::cout << "Incompatible shapes of LinopMats for addition." << '\n';
    exit(1);
  }
  temp.resize(left.r, left.c);
  for (int i = 0; i < left.r; i++) {
    for (int j = 0; j < left.c; j++) {
      temp(i, j) = left(i, j) + right(i, j);
    }
  }
  return temp;
};

template <class T>
LinopMat<std::complex<T> > operator+(LinopMat<T> left,
                                     LinopMat<std::complex<T> > right) {
  LinopMat<std::complex<T> > temp;
  if (left.r != right.r || left.c != right.c) {
    std::cout << "Incompatible shapes of LinopMats for addition." << '\n';
    exit(1);
  }
  temp.resize(left.r, left.c);
  for (int i = 0; i < left.r; i++) {
    for (int j = 0; j < left.c; j++) {
      temp(i, j) = left(i, j) + right(i, j);
    }
  }
  return temp;
};

/// \brief Multiplies all linops with the constant.
template <class T> LinopMat<T> operator*(T a, LinopMat<T> right) {
  for (int i = 0; i < right.r; i++) {
    for (int j = 0; j < right.c; j++) {
      right(i, j) = a * right(i, j);
    }
  }
  return right;
};

/// \brief Multiplies all linops with the constant.
template <class T>
LinopMat<std::complex<T> > operator*(T a, LinopMat<std::complex<T> > right) {
  for (int i = 0; i < right.r; i++) {
    for (int j = 0; j < right.c; j++) {
      right(i, j) = a * right(i, j);
    }
  }
  return right;
};

/// \brief Multiplies all linops with the constant.
template <class T>
LinopMat<std::complex<T> > operator*(std::complex<T> a, LinopMat<T> right) {
  LinopMat<std::complex<T> > out(right.r, right.c);
  for (int i = 0; i < right.r; i++) {
    for (int j = 0; j < right.c; j++) {
      out(i, j) = a * right(i, j);
    }
  }
  return out;
};

/// \brief Multiplication of LinopMats.
template <class T> LinopMat<T> operator*(LinopMat<T> left, LinopMat<T> right) {
  if (left.c != right.r) {
    std::cout << "Incompatible shapes of LinopMats for multiplication." << '\n';
    exit(1);
  }

  LinopMat<T> out(left.r, right.c);
  out.setConstant(0.0);
  for (int i = 0; i < left.r; i++) {
    for (int j = 0; j < right.c; j++) {
      for (int k = 0; k < left.c; k++) {
        out(i, j) = out(i, j) + (left(i, k) * right(k, j));
      }
    }
  }
  return out;
};
/// \brief Multiplication of LinopMats.
template <class T>
LinopMat<std::complex<T> > operator*(LinopMat<std::complex<T> > left,
                                     LinopMat<T> right) {
  if (left.c != right.r) {
    std::cout << "Incompatible shapes of LinopMats for multiplication." << '\n';
    exit(1);
  }

  LinopMat<std::complex<T> > out(left.r, right.c);
  out.setConstant(0.0);
  for (int i = 0; i < left.r; i++) {
    for (int j = 0; j < right.c; j++) {
      for (int k = 0; k < left.c; k++) {
        out(i, j) = out(i, j) + (left(i, k) * right(k, j));
      }
    }
  }
  return out;
};
/// \brief Multiplication of LinopMats.
template <class T>
LinopMat<std::complex<T> > operator*(LinopMat<T> left,
                                     LinopMat<std::complex<T> > right) {
  if (left.c != right.r) {
    std::cout << "Incompatible shapes of LinopMats for multiplication." << '\n';
    exit(1);
  }

  LinopMat<std::complex<T> > out(left.r, right.c);
  out.setConstant(0.0);
  for (int i = 0; i < left.r; i++) {
    for (int j = 0; j < right.c; j++) {
      for (int k = 0; k < left.c; k++) {
        out(i, j) = out(i, j) + (left(i, k) * right(k, j));
      }
    }
  }
  return out;
};

/// \brief Subtraction of LinopMats.
template <class T> LinopMat<T> operator-(LinopMat<T> left, LinopMat<T> right) {
  LinopMat<T> temp;
  if (left.r != right.r || left.c != right.c) {
    std::cout << "Incompatible shapes of LinopMats for subtraction." << '\n';
    exit(1);
  }
  temp.resize(left.r, left.c);
  for (int i = 0; i < left.r; i++) {
    for (int j = 0; j < left.c; j++) {
      temp(i, j) = left(i, j) - right(i, j);
    }
  }
  return temp;
};

template <class T>
LinopMat<std::complex<T> > operator-(LinopMat<std::complex<T> > left,
                                     LinopMat<T> right) {
  LinopMat<std::complex<T> > temp;
  if (left.r != right.r || left.c != right.c) {
    std::cout << "Incompatible shapes of LinopMats for addition." << '\n';
    exit(1);
  }
  temp.resize(left.r, left.c);
  for (int i = 0; i < left.r; i++) {
    for (int j = 0; j < left.c; j++) {
      temp(i, j) = left(i, j) - right(i, j);
    }
  }
  return temp;
};

template <class T>
LinopMat<std::complex<T> > operator-(LinopMat<T> left,
                                     LinopMat<std::complex<T> > right) {
  LinopMat<std::complex<T> > temp;
  if (left.r != right.r || left.c != right.c) {
    std::cout << "Incompatible shapes of LinopMats for addition." << '\n';
    exit(1);
  }
  temp.resize(left.r, left.c);
  for (int i = 0; i < left.r; i++) {
    for (int j = 0; j < left.c; j++) {
      temp(i, j) = left(i, j) - right(i, j);
    }
  }
  return temp;
};

template <class T> LinopMat<T> operator/(LinopMat<T> left_, T right) {
  LinopMat<T> outMat;
  outMat.resize(left_.r, left_.c);
  for (int i = 0; i < left_.r; i++) {
    for (int j = 0; j < left_.c; j++) {
      outMat(i, j) = left_(i, j) / right;
    }
  }
  return outMat;
};

template <class T>
LinopMat<std::complex<T> > operator/(LinopMat<std::complex<T> > left_,
                                     T right) {
  LinopMat<std::complex<T> > outMat;
  outMat.resize(left_.r, left_.c);
  for (int i = 0; i < left_.r; i++) {
    for (int j = 0; j < left_.c; j++) {
      outMat(i, j) = left_(i, j) / right;
    }
  }
  return outMat;
};

/// \brief Linear equation solver.
template <class T>
ChebfunMat<std::complex<T> >
linSolve(const LinopMat<std::complex<T> > &Lmat_,
         const BcMat<std::complex<T> > &bcmat_,
         const ChebfunMat<std::complex<T> > &forc_) {
  LinopMat<std::complex<T> > Lmat = Lmat_;
  ChebfunMat<std::complex<T> > forc = forc_;
  BcMat<std::complex<T> > bcmat = bcmat_;

  Discretize<std::complex<T> > Dis;
  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> A0 =
      Dis.MatAppend(Lmat, bcmat);
  Eigen::ColPivHouseholderQR<
      Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> >
      qr;
  qr.compute(A0);
  if (!qr.isInvertible()) {
    std::cout << "Solution not possible. This system is singular." << '\n';
    exit(1);
  }
  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> force =
      forc.ChebfunMat2EigenMat();
  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> force2(
      A0.rows(), force.cols());
  std::cout << "force: " << size(force2) << '\n' << std::flush;
  for (int i = 0; i < force2.cols(); i++) {
    force2.block(0, i, force.rows(), 1) = force.col(i);
    force2.block(force.rows(), i, bcmat.vals.rows(), 1) = bcmat.vals;
  }
  std::cout << "A0: " << size(A0) << '\n' << std::flush;

  // Solve:
  force2 = qr.inverse() * force2;

  // Parse it from form (a0 ... an C0 C1) to (a0 ... an), append constants to
  // to the series
  MatGen<std::complex<T> > Mat;
  int row_counter = 0;
  for (int i = 0; i < forc.c; i++) {
    int n = Dis.highest_each_column[i];
    Mat.compute(n);
    force.block(i * (N + 1), 0, N + 1, forc.c) =
        Mat.mats2[n] *
        force2.block(row_counter, 0, N + 1 + n, forc.c);
    row_counter += N + 1 + n;
  }


  std::cout << "force.size: " << size(force) << '\n';
  ChebfunMat<std::complex<T> > out(forc.r, forc.c);
  // Assignment automatically assings Eigen matrix to a ChebfunMat.
  out = force;
  // This is raw assingn in the previous line. Specify that everything is in
  // Cheb space. Then convert to physical space for convinience:
  for (int i = 0; i < out.r; i++) {
    for (int j = 0; j < out.c; j++) {
      out(i, j).dct_flag = SIS_CHEB_SPACE;
      out(i, j).c2p();
    }
  }
  return out;
}

} // namespace sis
/*/// \brief Defining addition of ChebfunMats.
template <class T>
 ChebfunMat<T> operator+( ChebfunMat<T> a,  ChebfunMat<T> b){
   ChebfunMat<T> out;
  if (a.r != b.r || a.c != b.c){
    std::cout << "addition of vector of ChebfunMats not possible"
    << ". In " << __LINE__ << '\n';
    exit(1);
    out.resize(a.r,a.c);
    for (int i = 0; i < out.r; i++){
      for (int j = 0; j < out.c; j++){
        out(i,j) = a(i,j) + b(i,j);
      }
    }
  }
};

/// \brief Defining addition of vector of ChebfunMats. Used for eigenvectors,
/// Singular values, adjoints etc.
template <class T>
std::vector< ChebfunMat<T> > operator+(std::vector< ChebfunMat<T> > a,
                      std::vector< ChebfunMat<T> > b){
  if (a.size() != b.size()){
    std::cout << "addition of vector of ChebfunMats not possible"
    << ". In " << __LINE__ << '\n';
    exit(1);
  }
  std::vector<sis::ChebfunMat<T> > out;
  out.resize(a.size());
  for (int i = 0; i < a.size(); i++ ){
    if (a[i].r != b[i].r || a[i].c != b[i].c){
      std::cout << "addition of vector of ChebfunMats not possible"
      << ". In " << __LINE__ << '\n';
      exit(1);
    } else {
      out[i] = a[i] + b[i];
    }
  }
  return out;
};*/

/// \brief This is used to define a multiplying a variable to a complex type.
/// C++ only allows multiplying a complex type with another complex type, so
/// we define this as a work around.
// template <class T>
// std::complex<T> operator*(T a, std::complex<T> b) {
//  return std::complex<T>(a, 0.0) * b;
//}
