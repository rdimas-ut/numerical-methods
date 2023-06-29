#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include "iterativeNL.h"
using namespace std;

void f(const Vector& x, Vector& y) {
  y(0) = pow(x(0),2) + pow(x(1),2) - 1;
  y(1) = pow(x(0)-1,2) + pow(x(1),2) - 1;
}

void df(const Vector& x, Matrix& a) {
  a(0,0) = 2*x(0);
  a(0,1) = 2*x(1);
  a(1,0) = 2*(x(0) - 1);
  a(1,1) = 2*x(1);
}
void f2(const Vector& x, Vector& y) {
  y(0) = pow(x(0),2) + (4*pow(x(1),2)) - 4;
  y(1) = (4*pow(x(0)-1,2)) + pow(x(1),2) - 4;
}

void df2(const Vector& x, Matrix& a) {
  a(0,0) = 2*x(0);
  a(0,1) = 8*x(1);
  a(1,0) = 8*x(0);
  a(1,1) = 2*x(1);
}

void f3(const Vector& x, Vector& y) {
  y(0) = pow(x(0),2) - (4*pow(x(1),2)) - 4;
  y(1) = pow(x(0)-1,2) + pow(x(1),2) - 4;
}

void df3(const Vector& x, Matrix& a) {
  a(0,0) = 2*x(0);
  a(0,1) = 0 - 8*x(1);
  a(1,0) = 2*(x(0) - 1);
  a(1,1) = 2*x(1);
}

int main() {
  double tol;
  int maxIter, iter, selection;
  int init;
  int example;
  inl_state s;

  Vector x(2);
  Vector x0(x.n());

  // Input

  cout << "Enter guess at fixed point (" << x.n() << " numbers): " << flush;
  cin >> x0;
  cout << "Enter tolerance and max iterations: " << flush;
  cin  >> tol >> maxIter;
  cout.precision(11);
  cout << "Select solver from fixecpt(0), newton(1), broyden1(2), or broyden2(3) : " << flush;
  cin >> selection;
  cout << "Select init: " << flush;
  cin >> init;
  cout << "Select example: " << flush;
  cin >> example;
  cout << "\n";

  // Solve
  iter = maxIter;
  x = x0;
  clock_t t = clock();
  switch(selection) {
    case 0:
      s = fixedpt(f, x, tol, iter);
      break;
    case 1:
      s = newton(f3, df3 , x, tol, iter);
      break;
    case 2:
      switch(example) {
        case 1:
          s = broyden1(f, df, x, tol, iter, init);
          break;
        case 2:
          s = broyden1(f2, df2, x, tol, iter, init);
          break;
        case 3:
          s = broyden1(f3, df3, x, tol, iter, init);
          break;
      }
      break;
    case 3:
      switch(example) {
        case 1:
          s = broyden2(f, df, x, tol, iter, init);
          break;
        case 2:
          s = broyden2(f2, df2, x, tol, iter, init);
          break;
        case 3:
          s = broyden2(f3, df3, x, tol, iter, init);
          break;
      }
    break;
  }

  // Report results

  switch(s) {
  case INL_SUCCESS: {
    t = clock() - t;
    cout << "\nIn " << iter << " iterations, the fixed point is:\n";
    cout << x << endl;
    cout << "Total Run time = " << ((float)t)/CLOCKS_PER_SEC << flush;
    return 0;
    break;
  }
  case INL_WONT_STOP:
    cout << "ERROR: Failed to converge in " << maxIter << " iterations!"
         << endl;
    break;
  case INL_BAD_DATA:
    cout << "ERROR: Data is bad!" << endl;
    break;
  case INL_BAD_ITERATE:
    cout << "ERROR: Bad iterate!" << endl;
    break;
  default:
    cout << "ERROR: Coding error!" << endl;
  }

  return 1;
}
