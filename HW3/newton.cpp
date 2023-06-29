///////////////////////////////////////////////////////////////////////////////
// NEWTON'S METHOD
//
// state newton(double& x, double tolerance, int maxIteration)
//
// Inputs:
//   x              The initial guess at the solution.
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
//  Remark: We assume we are given the two functions
//   f              The name of the function for which a root is sought.
//   df             The name of the derivative of the function.  The derivative
//                  of f must be computed by hand and coded correctly as df,
//                  or newton will not work!
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
using namespace std;
static int prec;

enum state {SUCCESS=0, WONT_STOP, BAD_DATA, BAD_ITERATE};
// ((pow(x,2)+1)/(x+1))-(3*x), (2*pow(x,2) + (2*x) - 1)/(pow(x+1,2))-3
double f(double x) { return  (2*exp(x-1) - pow(x,2) - 1); }
double df(double x) { return 2*exp(x-1) - 2*x; }

state newton(double& x, double tolerance, int maxIteration, int debug) {
  if(debug) {
    cout << "Guess: x = " << fixed << x << endl;
  }
    double r;
    double t;
    r=1;
  for(int iteration = 1; iteration <= maxIteration; iteration++) {
    double dfx = df(x);
    if(dfx == 0.0) return BAD_ITERATE;

    double dx = -f(x)/dfx;
    
    
    
    x += dx;
      
      if(fabs(dx) < tolerance) {
          cout << "Iter " << iteration << ": x = " << fixed << x
          << ", dx = " << fixed << dx << ", e = " << abs(x-r)
          << endl;
         return SUCCESS;
          
      }
      
    if(debug and (iteration==1)) {
        cout << "Iter " << iteration << ": x = " << fixed << x
        << ", dx = " << fixed << dx << ", e = " << abs(x-r)
        << endl;
      }

    if(debug and (iteration>1)) {
        cout << "Iter " << iteration << ": x = " << fixed << x
        << ", dx = " << fixed << dx << ", e = " << abs(x-r)
        << ", e ratio = " << abs(x-r)/t << endl;
    }
    
    t = abs(double (x-r));
    
  }

  return WONT_STOP;
}

int main() {
  double root,tol;
  int maxIter,debug;

  // Input

  cout << "Enter guess at root: ";
  cin >> root;
  cout << "Enter tolerance,  max iteration, and debug flag: ";
  cin >> tol >> maxIter >> debug;

  // Solve

  prec = (int) (log10(1.0/tol));
  cout.precision(prec);
  state s = newton(root, tol, maxIter, debug);

  // Report results

  switch(s) {
  case SUCCESS: {
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
