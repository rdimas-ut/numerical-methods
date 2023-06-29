///////////////////////////////////////////////////////////////////////////////
// Conjugate Gradients
//
// Solves Ax = b for x using  Conjugate Gradients with no preconditioner.
//
// INPUTS:
//   A        the matrix
//   b        the RHS vector
//   x        initial guess at the solution
//   maxIter  maximum number of iterations to take
//   tol      error tolerance
// OUTPUT:
//   x        the solution (if successful)
// RETURN:
//   state    an error status
///////////////////////////////////////////////////////////////////////////////

#include "matrix.h"
#include <cmath>
using namespace std;

//#define MONITOR  // output history of l2-error of residual
//#define MONITOR2 // output history of solution and l2-error

#define PLOT     // output history of l2-error of residual in Matlab format
#ifdef PLOT
  #include "matlabPlot.h"
#endif

enum state {SUCCESS, WONT_STOP, NOT_SYMMETRIC, BAD_DATA};

state conj_grad(const Matrix& A, const Vector& b, Vector& x,
	    	 int& maxIter, double tol) {

  // CHECK DATA

  int n = A.n(0);

  if(A.n(1) != n || b.n() != n || x.n() != n) return BAD_DATA;
  if(tol <= 0) return BAD_DATA;
  if(maxIter <= 0) maxIter = 1;

  for(int i=0; i<n-1; i++) {
    for(int j=i+1; j<n; j++) {
      if(A(i,j) != A(j,i) ) return NOT_SYMMETRIC;
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
      maxIter = iter;
#ifdef PLOT
      res.resize(maxIter+1);
      matlabPlot(fout,res,"Convergence of Conjugate Gradient method",
                          "iter","L2-norm of residual");
      fout.close();
#endif
      return SUCCESS;
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
      return SUCCESS;
    }

    double s = beta / alpha;

    for(int i=0; i<n; i++) {
      d(i) = r(i) + s*d(i);
    }

    alpha = beta;
  }

  return WONT_STOP;
}

// TEST PROGRAM

#include <iostream>
using namespace std;

int main() {
  int n;
  cout << "Enter size: " << flush;
  cin >> n;

  Matrix A(n,n);
  Vector x(n);
  Vector b(n);

  x = 0.0;

  cout << "Enter A by rows: " << flush;
  cin >> A;

  cout << "Enter b: " << flush;
  cin >> b;

  int maxIter;
  double tolerance;
  cout << "Enter maxIter and tolerance: " << flush;
  cin >> maxIter >> tolerance;

  state s = conj_grad(A,b,x,maxIter,tolerance);

  switch(s) {
  case WONT_STOP:
    cout << "ERROR: Exceeded maximum number of iterations." << endl;
    return 1;
  case NOT_SYMMETRIC:
    cout << "ERROR: Matrix is not symmetric." << endl;
    return 1;
  default:
    cout << "ERROR: Unspecified." << endl;
    return 1;
  case SUCCESS:
    cout << "The solution is:" << endl;
    cout << x << endl;

    Vector y(n);
    matVecMult(A,x,y);
    y -= b;
    cout << "The number of iterations is: " << maxIter << endl;
    cout << "The max-norm of residual is: " << maxNorm(y) << endl;
    cout << "The residual is: " << endl;
    cout << y << endl;
    return 0;
  }
}
