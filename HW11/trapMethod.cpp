///////////////////////////////////////////////////////////////////////////////
// Trapezoidal Method for Solving an ODE
//
// y' = f(t,y), a < t < b,
// y(a) = y0
//
// Inputs:
//   f         The name of the function "f" in the ODE
//   y         Space for the solution (with correct size on input)
//   y0        Initial condition
//   a         Initial time
//   b         Final time
// Outputs:
//   y         The solution (at equally spaced points)
///////////////////////////////////////////////////////////////////////////////
#include "trapMethod.h"
#include "matlabPlot.h"
using namespace std;

void trapMethod(double f(double,double), vector<double>& y, double y0, double a, double b) {
  int n = y.size()-1;
  double h = (b-a)/n;

  double t = a;
  y[0] = y0;

  for(int i=0; i<n; i++) {
    double f1 = f(t,y[i]);
    double yEuler = y[i] + h*f1;

    t += h;
    y[i+1] = yEuler;
  }
}

void trapMethod(double f(double,double), vector<double>& x, vector<double>& y, double y0, double a, double b) {
  int n = y.size()-1;
  double h = (b-a)/n;

  double t = a;
  y[0] = y0;
  x[0] = t;
  for(int i=0; i<n; i++) {
    double f1 = f(t,y[i]);
    double yEuler = y[i] + h*f1;

    t += h;

    y[i+1] = yEuler;
    x[i+1] = t;
  }
}
