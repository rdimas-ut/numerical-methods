#include <iostream>
#include "iterativeNL.h"
#include "matrix.h"
#include "gaussElim.h"

using namespace std;

#define MONITOR

inl_state fixedpt(void g(const Vector&, Vector&), Vector& x, double tol,
                  int& iter) {
  if (iter < 1) iter = 1;
  if (tol <= 0) return INL_BAD_DATA;

  int maxIter = iter;
  Vector y(x);
  double error;
  for (iter = 1; iter <= maxIter; iter++) {
    g(y, x);  // x = g(y) is the new guess
    y -= x;   // y is now the change to x
    error = maxNorm(y);
#ifdef MONITOR
    cout << "Iter " << iter << ": x= " << x << ", err = " << error << endl;
#endif
    if (error <= tol) return INL_SUCCESS;
    y = x;
  }
  return INL_WONT_STOP;
}

inl_state newton(void f(const Vector&, Vector&),
                 void df(const Vector&, Matrix&), Vector& x, double tol,
                 int& iter) {
  if (iter < 1) iter = 1;
  if (tol <= 0) return INL_BAD_DATA;

  int maxIter = iter;
  Vector y(x);
  Vector s(x.n());
  Matrix a(x.n(), x.n());
  Permutation p(x.n());
  ge_state g;
  double error;

  for (iter = 0; iter < maxIter; iter++) {
    df(y, a);  // a is equal now to DF(x_k)
    f(x, s);   // s is qual now to F(x_k)
    // catches initial guesses that are equal to roots
    if (iter == 0) {
      double initial_error = maxNorm(s);
      if (initial_error <= tol) {
        cout << scientific << "Iter " << iter << ": x= " << x
             << ", err = " << error << endl;
        return INL_SUCCESS;
      }
    }
    if (maxNorm(s) == 0.0) return INL_BAD_ITERATE;
    s *= -1;             // s is qual now to -F(x_k)
    g = solve(a, p, s);  // solves the system DF(x_k)s=-F(x_k)
    switch (g) {
      case GE_SUCCESS:
        break;
      case GE_SINGULAR:
        return INL_BAD_ITERATE;
        break;
      case GE_BADDATA:
        return INL_BAD_DATA;
        break;
    }
    x += s;  // x is now new guess for x
    y -= x;  // y is now the change to x
    error = maxNorm(y);
#ifdef MONITOR
    cout << scientific << "Iter " << iter + 1 << ": x= " << x
         << ", err = " << error << endl;
#endif
    if (error <= tol) return INL_SUCCESS;
    y = x;
  }
  return INL_WONT_STOP;
}

inl_state broyden2(void f(const Vector&, Vector&),
                   void df(const Vector&, Matrix&), Vector& x, double tol,
                   int& iter, int init) {
  if (iter < 1) iter = 1;
  if (tol <= 0) return INL_BAD_DATA;
  int n = x.n();
  int maxIter = iter;
  double error; double factor;
  Vector r(n); Vector s(n);
  Vector l_delta(n); Vector u_delta(n);
  Matrix b(n, n); Matrix a(n, n); Matrix c(n, n);
  Permutation p(n);

  // Sets up B_{0} matrix as I or DF(x_{0}) depending on selection.
  if (init == 0) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        b(i, j) = 0;
      }
    }
    for (int i = 0; i < n; i++) {
      b(i, i) = 1;
    }
  } else {
    df(x, b);
    Matrix facA(b);
    luFactorize(facA, p);
    invFactoredMatrix(facA, p, b);
  }

  for (iter = 0; iter < maxIter; iter++) {
    f(x, r);
    error = maxNorm(r);  // r is qual now to F(x_i)
    // catches initial guesses that are equal to roots
    if (iter == 0 && error <= tol) {
      cout << scientific << "Iter " << iter << ": x= " << x
           << ", err = " << error << endl;
      return INL_SUCCESS;
    }
    // New x_{i+1}
    matVecMult(b, r, s);
    l_delta = x;
    x -= s;

    // Calculates l_delta
    l_delta -= x;
    l_delta *= -1;
    error = maxNorm(l_delta);
#ifdef MONITOR
    cout << scientific << "Iter " << iter + 1 << ": x= " << x
         << ", err = " << error << endl;
#endif
    if (error <= tol) return INL_SUCCESS;

    // r = B_{i}u_delta
    f(x, s);
    s -= r;
    matVecMult(b, s, r);

    // factor = l_delta^{T}B_{i}u_delta
    factor = scDot(r, l_delta);
    factor = (1 / factor);

    // r = l_delta - B_{i}u_delta
    r -= l_delta;
    r *= -1;

    // s = l_delta^{T}B_{i}
    // Note that s is a column vector
    for (int i = 0; i < x.n(); i++) {
      a(0, i) = l_delta(i);
    }
    matMatMult(a, b, c);
    for (int i = 0; i < x.n(); i++) {
      s(i) = c(0, i);
    }

    // c = factor*rs
    // outer product of r and s with a scaling factor
    c = outer_product(r, s);
    c *= factor;

    // b = b + c
    b += c;
  }
  return INL_WONT_STOP;
}

inl_state broyden1(void f(const Vector&, Vector&),
                   void df(const Vector&, Matrix&), Vector& x, double tol,
                   int& iter, int init) {
  if (iter < 1) iter = 1;
  if (tol <= 0) return INL_BAD_DATA;
  int n = x.n();
  int maxIter = iter;
  double error;double factor;
  Vector r(n); Vector s(n);
  Vector l_delta(n); Vector u_delta(n);
  Matrix a(n, n); Matrix c(n, n);

  // Sets up A_{0} matrix as I or DF(x_{0}) depending on selection.
  if (init == 0) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        a(i, j) = 0;
      }
    }
    for (int i = 0; i < n; i++) {
      a(i, i) = 1;
    }
  } else {
    df(x, a);
  }

  for (iter = 0; iter < maxIter; iter++) {
    Permutation p(n);
    f(x, r);
    error = maxNorm(r);  // r is qual now to F(x_i)
    // catches initial guesses that are equal to roots
    if (iter == 0 && error <= tol) {
      cout << scientific << "Iter " << iter << ": x= " << x
           << ", err = " << error << endl;
      return INL_SUCCESS;
    }
    // New x_{i+1}
    l_delta = x;
    matVecMult(a, x, s);
    s -= r;
    solve(a, p, s);
    x=s;
    // Calculates l_delta
    l_delta -= x;
    l_delta *= -1;
    error = maxNorm(l_delta);
#ifdef MONITOR
    cout << scientific << "Iter " << iter + 1 << ": x= " << x
         << ", err = " << error << endl;
#endif
    if (error <= tol) return INL_SUCCESS;

    // s = A_{i}l_delta
    matVecMult(a, l_delta, s);

    // factor = l_delta^{T}l_delta
    factor = scDot(l_delta, l_delta);
    factor = (1 / factor);

    // u_delta
    f(x,u_delta);
    u_delta -= r;

    // s = u_delta - A_{i}l_delta
    s -= u_delta; s *= -1;
    // c = factor*(u_delta - A_{i}l_delta)l_delta^{T}
    c = outer_product(s, l_delta);
    c *= factor;
    // a = a + c
    a += c;
  }
  return INL_WONT_STOP;
}
