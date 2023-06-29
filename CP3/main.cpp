#include <iostream>
#include <iomanip>
#include <cmath>
#include "iterativeNL.h"
using namespace std;

void f(const Vector& x, Vector& y) {
  y(0) = 6*pow(x(0),3) + x(0)*x(1) - 3*pow(x(1), 3) - 4;
  y(1) = pow(x(0), 2) - 18*x(0)*pow(x(1), 2) + 16*pow(x(1),3) + 1;
}

void df(const Vector& x, Matrix& a) {
  a(0,0) = 18*pow(x(0),2) + x(1);
  a(0,1) = x(0) - 9*pow(x(1),2);
  a(1,0) = 2*x(0) - 18*pow(x(1),2);
  a(1,1) = -36*x(0)*x(1) + 48*pow(x(1), 2);
}

int main() {
  double tol;
  int maxIter, iter, selection;
  inl_state s;

  Vector x(2);
  Vector x0(x.n());

  // Input

  cout << "Enter guess at fixed point (" << x.n() << " numbers): " << flush;
  cin >> x0;
  cout << "Enter tolerance and max iterations: " << flush;
  cin  >> tol >> maxIter;
  cout.precision(11);
  cout << "Select solver from fixecpt(0) or newton(1): " << flush;
  cin >> selection;
  cout << "\n";

  // Solve

  switch(selection) {
    case 0:
      iter = maxIter;
      x = x0;
      s = fixedpt(f, x, tol, iter);
      break;
    case 1:
      iter = maxIter;
      x = x0;
      s = newton(f, df , x, tol, iter);
      break;
  }

  // Report results

  switch(s) {
  case INL_SUCCESS: {
    cout << "\nIn " << iter << " iterations, the fixed point is:\n";
    cout << x << endl;
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
