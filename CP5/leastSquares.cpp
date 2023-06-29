#include <cmath>
#include "matrix.h"
#include "gaussElim.h"
#include "leastSquares.h"

using namespace std;

ls_state normalEquationsSetup(model t, Matrix& A, Vector& b, Vector& c){
  if (A.n(1) != 2){
    cout << "Bad data give!" << endl;
    return BAD_DATA;
  }
  b.resize(A.n(0));
  for (int i=0;i < A.n(0); i++) b(i) = A(i,1);

  //Setup for a specific model
  Matrix S(1,1);
  if (t == LINEAR){
    S.resize(A.n(0), 2);
    c.resize(2);
    for (int i = 0; i < A.n(0); i++){
      for (int j = 0; j < 2; j++){
        if (j == 1){ S(i,j) = 1;
        } else {
          S(i,j) = A(i,j); }
      }
    }
  } else if (t == QUAD) {
    S.resize(A.n(0), 3);
    c.resize(3);
    for (int i = 0; i < A.n(0); i++){
      for (int j = 0; j < 2; j++){
        if (j == 1){ S(i,2) = 1;
        } else {
          S(i,1) = A(i,j);
          S(i,0) = pow(A(i,j),2);  }
      }
    }
  } else if (t == CUBIC) {
    S.resize(A.n(0), 4);
    c.resize(4);
    for (int i = 0; i < A.n(0); i++){
      for (int j = 0; j < 2; j++){
        if (j == 1){ S(i,3) = 1;
        } else {
          S(i,2) = A(i,j);
          S(i,1) = pow(A(i,j),2);
          S(i,0) = pow(A(i,j),3);}
      }
    }
  } else if (t == EXP) {
    S.resize(A.n(0), 2);
    c.resize(2);
    for (int i = 0; i < A.n(0); i++){
      b(i) = log(b(i));
      for (int j = 0; j < 2; j++){
        if (j == 0){ S(i,j) = 1;
        } else {
          S(i,j) = A(i,0); }
      }
    }
  } else if (t == POWER) {
    S.resize(A.n(0), 2);
    c.resize(2);
    for (int i = 0; i < A.n(0); i++){
      b(i) = log(b(i));
      for (int j = 0; j < 2; j++){
        if (j == 0){ S(i,j) = 1;
        } else {
          S(i,j) = log(A(i,0)); }
      }
    }
  }
  A.resize(S.n(0), S.n(1));
  A = S;
  return SUCCESS;
}

ls_state incsSystemSolver(Matrix& A, Vector& x, Vector& b){
  if (A.n(0) != b.n()){
    cout << "Bad data give!" << endl;
    return BAD_DATA;
  }

  //Solves the linear system BSc = Bb
  Matrix B(A.n(1), A.n(0));
  Matrix C(A.n(1), A.n(1));
  Vector y(A.n(1));

  tranpose(A, B);
  matMatMult(B, A, C);
  matVecMult(B, b, y);
  Permutation P(C.n(0));
  solve(C, P, y);
  x = y;
  return SUCCESS;
}

ls_state normalEquationSolver(model t, Matrix& A, Vector& c){
  ls_state r;
  Vector b(1);
  Matrix B(A.n(0), A.n(1)); B = A;
  double se; double rmse;
  r = normalEquationsSetup(t, B, b, c);
  if (r == BAD_DATA) {return BAD_DATA;}
  r = incsSystemSolver(B, c, b);
  if (t == EXP || t == POWER){ c(0) = exp(c(0)); }

  Vector y(A.n(0));
  for (int i=0;i<A.n(0);i++) {
    b(i) = A(i, 1);
    if (t==LINEAR) {
      y(i) = (A(i, 0)*c(0)) + c(1);
    } else if (t==QUAD) {
      y(i) = (A(i, 0)*pow(c(0), 2)) + (A(i, 0)*c(1)) + c(2);
    } else if (t==CUBIC) {
      y(i) = (A(i, 0)*pow(c(0), 3)) + (A(i, 0)*pow(c(1), 2)) + (A(i, 2)*c(1)) + c(3);
    } else if (t==EXP) {
      y(i) = c(0)*exp(c(1)*A(i,0));
    } else if (t==POWER) {
      y(i) = c(0)*exp(c(1));
    }
  }
  for (int i=0;i<A.n(0);i++) { se += pow(y(i),2); }
  rmse = pow(se/A.n(0), .5); cout << endl;
  for (int i=0; i < c.n() - 1; i++) {
    cout << "c" << i << " = " << c(i) << ", ";}
  cout << "c" << c.n()-1 << " = " << c(c.n()-1) << endl;
  cout << "2-norm = " << pow(se, .5)
       << ", SE = " << se << ", RMSE = " << rmse << endl;

  return r;
}

ls_state pureIncsSolver(Matrix& A, Vector& x, Vector& b){
   ls_state r; r = incsSystemSolver(A, x, b);
   double se = 0; double rmse = 0;
   Vector y(b.n());
   matVecMult(A, x, y);
   y -= b; y *= -1;
   for (int i=0;i<A.n(0);i++) { se += pow(y(i),2); }
   rmse = pow(se/A.n(0), .5); cout << endl;
   for (int i=0; i < x.n()-1; i++) {
     cout << "x" << i << " = " << x(i) << ", ";}
   cout << "x" << x.n()-1 << " = " << x(x.n()-1) << endl;
   cout << "2-norm = " << pow(se, .5)
        << ", SE = " << se << ", RMSE = " << rmse << endl;
   return r;
}
