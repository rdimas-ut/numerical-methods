#include <iostream>
#include <cmath>
#include <fstream>
#include <string.h>
#include "iterativeLA.h"
#include "matrix.h"

//#define MONITOR // output history of l2-error
//define MONITOR2 // output history of solution and l2-error
//#define PLOT

#ifdef PLOT
  #include "matlabPlot.h"
#endif
using namespace std;

ila_state jacobi(const Matrix& A, const Vector& b, Vector& x,
	     int& maxIter, double tol) {
  // CHECK DATA

  int n = A.n(0);
  if(A.n(1) != n || b.n() != n || x.n() != n) return ILA_BAD_DATA;
  if(tol <= 0) return ILA_BAD_DATA;
  if(maxIter <= 0) maxIter = 1;

  for(int i=0; i<n; i++) {
    if(A(i,i) == 0) return ILA_BAD_DIAGONAL;
  }

#ifdef PLOT
    Vector res(maxIter);
    ofstream fout("Jacobi_plot.m");
    Vector r(n);
    matVecMult(A,x,r);
    r -= b;
    r *=(-1);
    double norm = scDot(r,r);
    res(0) = sqrt( norm );
#endif

  // APPLY JACOBI
  Vector xOld(x);
  for(int iter=0; iter<maxIter; iter++) {

    // Get new x
    for(int i=0; i<n; i++) {
      double sum = 0;
      for(int j=0; j<n; j++) {
	if(j==i) continue;
	sum += A(i,j)*xOld(j);
      }
      x(i) = ( -sum + b(i) ) / A(i,i);
    }

    // Check error tolerance
    xOld -= x;
    double l2error = l2norm(xOld) / (l2norm(x)+1);
#ifdef MONITOR
    cout << "Iter " << iter+1 << ", l2-error " << l2error << endl;
#endif
#ifdef MONITOR2
    cout << "Iter " << iter+1 << ", approx. solution: "
         << x << ", l2-error " << l2error << endl;
#endif
#ifdef PLOT
  Vector r(n);
  matVecMult(A,x,r);
  r -= b;
  r *=(-1);
  double norm = scDot(r,r);
  res(iter+1) = sqrt(norm);
#endif
    if( l2error <= tol) {
      maxIter = iter+1;
#ifdef PLOT
  res.resize(maxIter+1);
  matlabPlot(fout,res,"Convergence of Jacobi method",
                    "iter","L2-norm of residual");
  fout.close();
#endif
      return ILA_SUCCESS;
    }
    xOld = x;
  }

  return ILA_WONT_STOP;
}

ila_state gauss_seidel(const Matrix& A, const Vector& b, Vector& x,
	     int& maxIter, double tol) {
  // CHECK DATA
  int n = A.n(0);
  if(A.n(1) != n || b.n() != n || x.n() != n) return ILA_BAD_DATA;
  if(tol <= 0) return ILA_BAD_DATA;
  if(maxIter <= 0) maxIter = 1;

  // RECORD X0 RESIDUAL AND CREATE OFSTREAM
  #ifdef PLOT
      Vector res(maxIter);
      ofstream fout("Gauss_Seidel_plot.m");
      Vector r(n);
      matVecMult(A,x,r);
      r -= b;
      r *=(-1);
      double norm = scDot(r,r);
      res(0) = sqrt(norm);
  #endif

  // APPLY GAUSS_SEIDEL

  Vector xOld(x);
  for(int iter=0; iter<maxIter; iter++) {
    // Get new x
    for(int i=0; i<n; i++) {
      double sum = 0;
      for(int j=0; j<n; j++) {
	if(j==i) continue;
	sum += A(i,j)*x(j);
      }
      x(i) = ( -sum + b(i) ) / A(i,i);
    }
    // Check error tolerance
    xOld -= x;
    double l2error = l2norm(xOld) / (l2norm(x)+1);
#ifdef MONITOR
    cout << "Iter " << iter+1 << ", l2-error " << l2error << endl;
#endif
#ifdef MONITOR2
    cout << "Iter " << iter+1 << ", approx. solution: "
         << x << ", l2-error " << l2error << endl;
#endif
#ifdef PLOT
  Vector r(n);
  matVecMult(A,x,r);
  r -= b;
  r *=(-1);
  double norm = scDot(r,r);
  res(iter+1) = sqrt(norm);
#endif
    if( l2error <= tol) {
      maxIter = iter+1;
#ifdef PLOT
  res.resize(maxIter+1);
  matlabPlot(fout,res,"Convergence of Gauss Seidel method",
                    "iter","L2-norm of residual");
  fout.close();
#endif
      return ILA_SUCCESS;
    }
    xOld = x;
  }
  return ILA_WONT_STOP;
}

ila_state sor(const Matrix& A, const Vector& b, Vector& x,
	     int& maxIter, double tol, double omega) {

  // CHECK DATA
  int n = A.n(0);
  if(A.n(1) != n || b.n() != n || x.n() != n) return ILA_BAD_DATA;
  if(tol <= 0) return ILA_BAD_DATA;
  if(maxIter <= 0) maxIter = 1;
  if((omega < 1) or (omega > 2)){return ILA_BAD_DATA;}

  // RECORD X0 RESIDUAL AND CREATE OFSTREAM
#ifdef PLOT
      string omega_s = to_string(omega);
      Vector res(maxIter);
      ofstream fout("SOR_plot_"+omega_s+".m");
      Vector r(n);
      matVecMult(A,x,r);
      r -= b;
      r *=(-1);
      double norm = scDot(r,r);
      res(0) = sqrt(norm);
#endif

  // APPLY SOR

  Vector xOld(x);
  for(int iter=0; iter<maxIter; iter++) {

    // Get new x
    for(int i=0; i<n; i++) {
      double x_i = 0;
      double sum = 0;
      for(int j=0; j<n; j++) {
	if(j==i) continue;
	sum += A(i,j)*x(j);
      }
      x_i = ( -sum + b(i) ) / A(i,i);
      x(i) = ((1 - omega)*xOld(i)) + (omega*x_i);
    }

    // Check error tolerance
    xOld -= x;
    double l2error = l2norm(xOld) / (l2norm(x)+1);
#ifdef MONITOR
    cout << "Iter " << iter+1 << ", l2-error " << l2error << endl;
#endif
#ifdef MONITOR2
    cout << "Iter " << iter+1 << ", approx. solution: "
         << x << ", l2-error " << l2error << endl;
#endif
#ifdef PLOT
  Vector r(n);
  matVecMult(A,x,r);
  r -= b;
  r *=(-1);
  double norm = scDot(r,r);
  res(iter+1) = sqrt(norm);
#endif
    if( l2error <= tol) {
      maxIter = iter+1;
#ifdef PLOT
  res.resize(maxIter+1);
  matlabPlot(fout,res,"Convergence of SOR method",
                          "iter","L2-norm of residual");
  fout.close();
#endif
      return ILA_SUCCESS;
    }
    xOld = x;
  }
  return ILA_WONT_STOP;
}

ila_state conj_grad(const Matrix& A, const Vector& b, Vector& x,
	    	 int& maxIter, double tol) {

  // CHECK DATA

  int n = A.n(0);

  if(A.n(1) != n || b.n() != n || x.n() != n) return ILA_BAD_DATA;
  if(tol <= 0) return ILA_BAD_DATA;
  if(maxIter <= 0) maxIter = 1;

  for(int i=0; i<n-1; i++) {
    for(int j=i+1; j<n; j++) {
      if(A(i,j) != A(j,i) ) return ILA_NOT_SYMMETRIC;
    }
  }

  // INITIALIZE CONJUGATE GRADIENTS

  // Set initial residual r = b - Ax
  Vector r(n);

  matVecMult(A,x,r);
  r-=b; r*=(-1);

  double alpha = scDot(r,r);

#ifdef MONITOR
    cout << endl << "Iter 0, res l2-error " << sqrt(alpha) << endl;
#endif
#ifdef MONITOR2
    cout << endl << "Iter 0, guess: " << x
                 << ", res l2-error " << sqrt(alpha) << endl;
#endif
#ifdef PLOT
  Vector res(maxIter);
  ofstream fout("CG_plot.m");
  res(0) = sqrt(alpha);
#endif

  // Set initial search direction d = r
  Vector d(r); // Creates d and sets d = r

  double tolSq = tol*tol;

  // CONJUGATE GRADIENT LOOP

  for(int iter=0; iter<maxIter; iter++) {

    if(scDot(d,d) <= tolSq) {
      maxIter = iter+1;
#ifdef PLOT
      res.resize(maxIter+1);
      matlabPlot(fout,res,"Convergence of Conjugate Gradient method",
                          "iter","L2-norm of residual");
      fout.close();
#endif
      return ILA_SUCCESS;
    }

    // Set u = Ad
    Vector u(n);

    matVecMult(A,d,u);

    // Update x = x + td and r = r - tu
    double t = alpha / scDot(d,u);

    for(int i=0; i<n; i++) {
      x(i) += t*d(i);
      r(i) -= t*u(i);
    }

    // Get new search direction d = r + s*d;
    double beta = scDot(r,r);

#ifdef MONITOR
    cout << "Iter " << iter+1 << ", res l2-error " << sqrt(beta) << endl;
#endif
#ifdef MONITOR2
    cout << "Iter " << iter+1 << ", approx. solution: "
         << x << ", res l2-error " << sqrt(beta) << endl;
#endif
#ifdef PLOT
    res(iter+1) = sqrt(beta);
#endif

    if(beta <= tolSq) {
      maxIter = iter+1;
#ifdef PLOT
      res.resize(maxIter+1);
      matlabPlot(fout,res,"Convergence of Conjugate Gradient method",
                          "iter","L2-norm of residual");
      fout.close();
#endif
      return ILA_SUCCESS;
    }

    double s = beta / alpha;

    for(int i=0; i<n; i++) {
      d(i) = r(i) + s*d(i);
    }

    alpha = beta;
  }

  return ILA_WONT_STOP;
}

ila_state prec_conj_grad(ila_state (*prec)(const Matrix&, Vector&, double),
const Matrix& A, const Vector& b, Vector& x,
double omega, int& maxIter, double tol){

    // CHECK DATA

    int n = A.n(0);

    if(A.n(1) != n || b.n() != n || x.n() != n) return ILA_BAD_DATA;
    if(tol <= 0) return ILA_BAD_DATA;
    if(maxIter <= 0) maxIter = 1;

    for(int i=0; i<n-1; i++) {
      for(int j=i+1; j<n; j++) {
        if(A(i,j) != A(j,i) ) return ILA_NOT_SYMMETRIC;
      }
    }

    // INITIALIZE CONJUGATE GRADIENTS

    // Set initial residual r = b - Ax
    Vector r(n);

    matVecMult(A,x,r);
    r-=b; r*=(-1);

    double alpha = scDot(r,r);

  #ifdef MONITOR
      cout << endl << "Iter 0, res l2-error " << sqrt(alpha) << endl;
  #endif
  #ifdef MONITOR2
      cout << endl << "Iter 0, guess: " << x
                   << ", res l2-error " << sqrt(alpha) << endl;
  #endif
 #ifdef PLOT
    string title = "null";
    Vector res(maxIter);
    if (omega == 0.0){
      title = "CG_plot_Jacobi.m";
    } else if (omega == 1.0 ){
      title = "CG_plot_Gauss-Seidel.m";
    } else {
      string omega_s = to_string(omega);
      title = "CG_plot_SOR_"+omega_s+".m";
    }
    ofstream fout(title);
    res(0) = sqrt(alpha);
  #endif

    // Preconditions r and sets it to z
    Vector z(r);

    prec(A,z,omega);
    Vector d(z); // Creates d and sets d = z

    // sets alpha as r^{T}z
    alpha = scDot(r,z);

    double tolSq = tol*tol;

    // CONJUGATE GRADIENT LOOP

    for(int iter=0; iter<maxIter; iter++) {
      if(scDot(d,d) <= tolSq) { //Here I don't know why d, d so right now its r, r
        maxIter = iter+1;
  #ifdef PLOT
        res.resize(maxIter+1);
        matlabPlot(fout,res,"Convergence of Preconditioned Conjugate Gradient method",
                            "iter","L2-norm of residual");
        fout.close();
  #endif
        return ILA_SUCCESS;
      }

      // Set u = Ad
      Vector u(n);

      matVecMult(A,d,u);

      // Update x = x + td and r = r - tu
      double t = alpha / scDot(d,u);

      for(int i=0; i<n; i++) {
        x(i) += t*d(i);
        r(i) -= t*u(i);
      }

      // Gets new z
      Vector z(r);
      prec(A,z,omega);

      // Get new search direction d = r + s*d;
      double beta = scDot(r,z);
      double res_norm = scDot(r,r);

  #ifdef MONITOR
      cout << "Iter " << iter+1 << ", res l2-error " << sqrt(res_norm) << endl;
  #endif
  #ifdef MONITOR2
      cout << "Iter " << iter+1 << ", approx. solution: "
           << x << ", res l2-error " << sqrt(res_norm) << endl;
  #endif
  #ifdef PLOT
      res(iter+1) = sqrt(res_norm);
  #endif

      if(beta <= tolSq) {
        maxIter = iter+1;
  #ifdef PLOT
        res.resize(maxIter+1);
        matlabPlot(fout,res,"Convergence of Conjugate Gradient method",
                            "iter","L2-norm of residual");
        fout.close();
  #endif
        return ILA_SUCCESS;
      }

      double s = beta / alpha;

      for(int i=0; i<n; i++) {
        d(i) = z(i) + s*d(i);
      }
      alpha = beta;
    }

  return ILA_WONT_STOP;
}

ila_state jac_prec(const Matrix& A, Vector& x, double omega){
  int n = x.n();
  for(int i=0;i<n;i++){
    x(i) = x(i)/A(i,i);
  }
  return ILA_SUCCESS;
}

ila_state ssor_prec(const Matrix& A, Vector& x, double omega){
  int n = x.n();
  Vector c(n);
  c(0) = x(0);
  for(int i=1;i < n;i++){
    double c_entry = x(i);
    for(int j=0;j<i;j++){
      c_entry -= (A(i,j)*(1/A(j,j))*c(j)*omega);
    }
    c(i) = c_entry;
  }

  for(int i=n-1;i > -1;i--){
    double sums = 0;
    for(int j=n-1;j>i;j--){
      sums += A(i,j)*x(j);
    }
    sums *= omega;
    x(i) = (c(i) - sums)/A(i,i);
  }
  return ILA_SUCCESS;
}
