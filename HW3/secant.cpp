///////////////////////////////////////////////////////////////////////////////
// SECANT METHOD
//
// state secant(double p0, double& x,
//              double tolerance, int maxIteration, int debug)
//
// Inputs:
//   p0             The initial guess at the solution.
//   x              The second guess at the solution.
//   tolerance      The convergence tolerance (must be > 0).
//   maxIteration   The maximum number of iterations that can be taken.
//   debug          Boolean to set debugging output.
// Outputs:
//   x              The solution.
// Return:
//   state          An error status code.
//     SUCCESS      Sucessful termination.
//     WONT_STOP    Error: Exceeded maximum number of iterations.
//     BAD_ITERATE  Error: The function had a vanishing derivative.
//
//  Remark: We assume we are given the function
//   f              The name of the function for which a root is sought.
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
using namespace std;
static int prec;

enum state {SUCCESS=0, WONT_STOP, BAD_DATA, BAD_ITERATE};

double f(double x) { return 5*x + 10; }

state secant(double p0, double& x, double tolerance, int maxIteration, int debug) {
  if(debug) {
    cout << "Guess 0: x = " << p0 << endl;
    cout << "Guess 1: x = " << x << endl;
  }

  double fp0 = f(p0);
  for(int iteration = 1; iteration <= maxIteration; iteration++) {
    double fx = f(x);
    if(fx == fp0) return BAD_ITERATE;

    double dx = -fx*(x-p0)/(fx - fp0);
    p0 = x; fp0 = fx;
    x += dx;

    if(debug) {
      cout << "Iter " << iteration << ": x = " << fixed << x
           << ", dx = " << fixed << dx << endl;
    }

    if(fabs(dx) <= tolerance) return SUCCESS;
  }

  return WONT_STOP;
}

int main() {
  double root,root0,tol;
  int maxIter,debug;

  // Input

  cout << "Enter 2 guesses at root: ";
  cin >> root0 >> root;
  cout << "Enter tolerance, max iteration, and debug flag: ";
  cin >> tol >> maxIter >> debug;

  // Solve

  prec = (int) (log10(1.0/tol));
  cout.precision(prec);
  state s = secant(root0, root, tol, maxIter, debug);

  // Report results

  switch(s) {
  case SUCCESS: {
    prec = (int) (log10(1.0/tol) + log10(fabs(root))) + 1;
    cout << "The root is " << fixed << root << endl;
    return 0;
  }
  case WONT_STOP:
    cout << "ERROR: Failed to converge in " << maxIter << " iterations!"
         << endl;
    break;
  case BAD_ITERATE:
    cout << "ERROR: Obtained a vanishing derivative!" << endl;
    break;
  default:
    cout << "ERROR: Coding error!" << endl;
  }

  return 1;
}
