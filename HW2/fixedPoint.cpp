///////////////////////////////////////////////////////////////////////////////
// FIXED POINT (PICARD) ITERATION METHOD
//
// Solves the problem
//   g(x) = x
// using fixed point iteration
//
// state fixedPoint(double& x, double tolerance, int maxIteration, int debug)
//
// Inputs:
//   x             The initial guess at the fixed point.
//   tolerance     The convergence tolerance (must be > 0).
//   maxIteration  The maximum number of iterations that can be taken.
//   debug         Boolean for printing out information on every iteration.
// Outputs:
//   x             The  fixed point.
// Return:
//   state         An error status code.
//     SUCCESS     Sucessful termination.
//     WONT_STOP   Error: Exceeded maximum number of iterations.
//
// We assume the function g is given:
//   double g(double);
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <cmath>
using namespace std;

enum state {SUCCESS=0, WONT_STOP, BAD_DATA, BAD_ITERATE};

double g(double x) {
  return pow(3-log(x),1.0/2);
}

state fixedPoint(double& x, double tolerance, int maxIteration, int debug) {
  if(debug) cout << "Iter " << 0 << ": x = " << x << endl;

  for(int iteration = 1; iteration <= maxIteration; iteration++) {
    double gx = g(x);
    double dx = x - gx;
    x = gx;

    if(debug) cout << "Iter " << iteration << ": x = " << x << endl;

    // Check error tolerance
    if(fabs(dx) <= tolerance*(fabs(x)+1)) return SUCCESS;
  }

  return WONT_STOP;
}

int main() {
  double x,tol;
  int maxIter;
  int debug;

  // Input

  cout << "Enter x: " << flush;
  cin >> x;
  cout << "Enter tolerance and maxIteration: " << flush;
  cin >> tol >> maxIter;
  cout << "Monitor iterations? (1/0): " << flush;
  cin >> debug;

  // Solve for fixed point

  state s = fixedPoint(x, tol, maxIter, debug);

  // Report results

  switch(s) {
  case SUCCESS: {
    int prec = (int) (log10(1.0/tol) + log10(fabs(x))) + 1;
    cout.precision(prec);
    cout << "The fixed point is " << x << endl;
    cout << "  g(" << x << ") = " << g(x) << endl;
    return 0;
  }
  case WONT_STOP:
    cout << "ERROR: Failed to converge in " << maxIter << " iterations!" << endl;
    return 1;
  default:
    cout << "ERROR: Coding error!" << endl;
    return 1;
  }
}
