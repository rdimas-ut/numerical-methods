#include <iostream>
#include "iterativeNL.h"
#include "matrix.h"
#include "gaussElim.h"

using namespace std;

#define MONITOR

inl_state fixedpt(void g(const Vector&,Vector&), Vector& x,
	      double tol, int& iter) {
  if(iter < 1) iter = 1;
  if(tol <= 0) return INL_BAD_DATA;

  int maxIter = iter;
  Vector y(x);
  double error;
  for(iter = 1; iter <= maxIter; iter++) {
    g(y,x);      // x = g(y) is the new guess
    y -= x;      // y is now the change to x
    error = maxNorm(y);
#ifdef MONITOR
    cout << "Iter " << iter << ": x= " << x << ", err = " << error << endl;
#endif
    if (error <= tol) return INL_SUCCESS;
    y = x;
  }
  return INL_WONT_STOP;
}

inl_state newton(void f(const Vector&,Vector&), void df(const Vector&,Matrix&),
        Vector& x, double tol, int& iter){
  if(iter < 1) iter = 1;
  if(tol <= 0) return INL_BAD_DATA;

  int maxIter = iter;
  Vector y(x);
  Vector s(x.n());
  Matrix a(x.n(), x.n());
  Permutation p(x.n());
  ge_state g;
  double error;
  for(iter = 0; iter < maxIter; iter++) {
    df(y,a); // a is equal now to DF(x_k)
    f(x,s); // s is qual now to F(x_k)
    // catches initial guesses that are equal to roots
    if (iter == 0){
      double initial_error = maxNorm(s);
      if (initial_error <= tol) {
        cout << scientific << "Iter " << iter << ": x= " << x << ", err = " << error << endl;
        return INL_SUCCESS;
      }
    }
    if(maxNorm(s) == 0.0) return INL_BAD_ITERATE;
    s *= -1; // s is qual now to -F(x_k)
    g = solve(a, p, s); // solves the system DF(x_k)s=-F(x_k)
    switch(g){
      case GE_SUCCESS:
        break;
      case GE_SINGULAR:
        return INL_BAD_ITERATE;
        break;
      case GE_BADDATA:
        return INL_BAD_DATA;
        break;
    }
    x += s; // x is now new guess for x
    y -= x; // y is now the change to x
    error = maxNorm(y);
  #ifdef MONITOR
    cout << scientific << "Iter " << iter+1 << ": x= " << x << ", err = " << error << endl;
  #endif
    if (error <= tol) return INL_SUCCESS;
    y = x;
  }
  return INL_WONT_STOP;
}
