/******************************************************************************
We approximate
  - (p(x) u')' + q(x) u  = f(x),  a <= x <= b
  u(a) = u(b) = 0
using the Finite Element method with piecewise linears.  We evaluate the
integrals using the composite trapezoidal rule.

INPUTS
  alpha, beta  doubles ya and yb
  x       A vector of length n+2.  The grid points.
  p,q,f   The names of the functions in the ODE.
OUTPUTS
  u       A vector of length n+2.  The solution.
RETURN
  fem1d_state   The return status.
                LINFE_SUCCESS
                LINFE_BAD_DATA  if there is a problem with the data.
                LINFE_SINGULAR  if there is a problem solving the linear system
******************************************************************************/
#include "fem1d.h"

#define TEST

fem1d_state linFEM1d(double& alpha, double& beta, const Vector& x, Vector& u, func p, func q, func f) {
  int n = x.n();
  if(n != u.n()) return LINFE_BAD_DATA;
  n -= 2;

  u(0) = alpha;
  u(n+1) = beta;

  // Set up a tridiagonal system

  Matrix A(3,n+2);  // Note: A(.,0), A(.,n+1), A(0,1), and A(2,n) are not used.
  const int diag=1; // diagonal -- a_{i,i}
  const int supd=2; // superdiagonal -- a_{i,i+1}
  const int subd=0; // subdiagonal -- a_{i,i-1}

  double hL = x(1) - x(0), hR;
  if(hL <= 0) return LINFE_BAD_DATA;

  double pL = p(x(0)), pC = p(x(1)), pR;
  double qC, fC;

  for(int i=1; i<=n; i++) {
    hR = x(i+1) - x(i);
    if(hR <= 0) return LINFE_BAD_DATA;

    pR = p(x(i+1));
    qC = q(x(i));
    fC = f(x(i));

    A(diag,i)   = (pL + pC)/(2*hL) + (pR + pC)/(2*hR) + (hL + hR)*qC/2;
    A(subd,i+1) = -(pR + pC)/(2*hR);
    A(supd,i)   = -(pR + pC)/(2*hR);

    u(i) = (hL + hR)*fC/2;

    // Nonhomogenous Modification
    if(i == 1){
      u(i) -= -(alpha*((pR + pC)/(2*hR)));
    } else if (i == n) {
      u(i) -= -(beta*((pR + pC)/(2*hR)));
    }

    // Reset for next iteration
    hL = hR;
    pL = pC;
    pC = pR;
  }

  // Solve tridiagonal system using gaussian elimination w/o pivoting

  if(A(diag,1) == 0) return LINFE_SINGULAR;
  A(diag,1) = 1/A(diag,1);  // Store inverse of diagonal.  Eliminates to 1.
  A(supd,1) *= A(diag,1);
  u(1) *= A(diag,1);

  for(int i=2; i<=n; i++) {
    // Eliminate subd entry
    A(diag,i) -= A(subd,i)*A(supd,i-1);
    u(i) -= A(subd,i)*u(i-1);

    // rescale so diagonal is 1
    if(A(diag,i) == 0) return LINFE_SINGULAR;
    A(diag,i) = 1/A(diag,i);  // Store inverse of diagonal.  Eliminates to 1.
    A(supd,i) *= A(diag,i);
    u(i) *= A(diag,i);
  }

  for(int i=n-1; i>= 1; i--) {
    u(i) -= u(i+1)*A(supd,i);
  }

  return LINFE_SUCCESS;
}

#ifdef TEST
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include "matlabPlot.h"
using namespace std;

int test=0;

double p(double x) {
  switch(test) {

  case 1:
    return (-(9));

  case 2:
    return (-(pow(M_E,2*x)));

  default:
    cout << "ERROR: Undefined test number in p" << endl;
    return 1e+30;
    break;
  }
}

double q(double x) {
  switch(test) {

  case 1:
    return (pow(M_PI, 2));

  case 2:
    return (-(3*pow(M_E, 2*x)));

  default:
    cout << "ERROR: Undefined test number in q" << endl;
    return 1e+30;
    break;
  }
}

double f(double x) {
  switch(test) {

  case 1: {
    return 0.0; }

  case 2:
    return 0.0;

  default:
    cout << "ERROR: Undefined test number in f" << endl;
    return 1e+30;
    break;
  }
}

double y_tr(double x) {
  switch(test) {

  case 1:
    return ( (3*(sin((M_PI*x)/(3)))) - (cos((M_PI*x)/(3))) );

  case 2:
    return ( pow(M_E, 3 - (3*x)) );

  default:
    cout << "ERROR: Undefined test number in y_tr" << endl;
    return 1e+30;
    break;
  }
}

void interval(double& a, double& b) {
  switch(test) {

  case 1:
    a = 0.0; b = 1.5;
    break;

  case 2:
    a = 0.0; b = 1.0;
    break;

  default:
    cout << "ERROR: Undefined test number in interval" << endl;
    break;
  }
}

void boundary(double& alpha, double& beta) {
  switch(test) {

  case 1:
    alpha = -1.0; beta = 3.0;
    break;

  case 2:
    alpha = pow(M_E, 3); beta = 1.0;
    break;

  default:
    cout << "ERROR: Undefined test number in interval" << endl;
    break;
  }
}

int main() {

  int maxTest = 2;
  cout << "Enter test number: " << flush;
  cin >> test;
  if ((test<1)||(test>maxTest)) {
    cout << "ERROR: Invalid test number " << test
         << ". Must be between 1 and " << maxTest << "!" << endl;
    return 1;
  }

  int n;
  cout << "Enter number of points: ";
  cin >> n;
  if (n<3) {
    cout << "ERROR: Invalid number of points. Must be at least 3!" << endl;
    return 1;
  }

  Vector x(n);
  Vector y(n);
  double a, b, alpha, beta;
  interval(a,b);
  boundary(alpha, beta);
  double h = (b-a)/(n-1);
  for(int i=0; i<n; i++) x(i) = a+i*h;

  fem1d_state s = linFEM1d(alpha, beta, x, y, p, q, f);

  switch(s) {
  case  LINFE_BAD_DATA:
    cout << "ERROR: Bad data." << endl;
    return 1;
  case LINFE_SINGULAR:
    cout << "ERROR: Singular matrix encountered." << endl;
    return 1;
  default:
    break;
  }

  Vector yTrue(n);
  Vector error(n);
  int i;
  double rmsError = 0;
  for(i=0; i<n; i++) {
    yTrue(i) = y_tr(x(i));
    error(i) = yTrue(i) - y(i);
    rmsError += error(i)*error(i);
  }

  yTrue(0) = alpha; yTrue(n-1) = beta;
  rmsError = sqrt(rmsError/n);

  cout << "The solution is:\n";
  cout << y << endl;

  cout << "The true solution is:\n";
  cout << yTrue << endl;

  // cout << "The errors are:\n";
  // cout << error << endl;

  // cout << "The root mean square error is: " << rmsError << endl;


  // ofstream fout("numsol.m");
  // matlabPlot(fout, x, y, "FEM with p.w. linears and composite trapezoidal rule",
  //             "X", "Approximate Solution");
  // fout.close();
  //
  // fout.open("sol.m");
  // matlabPlot(fout, x, yTrue, "FEM with p.w. linears and composite trapezoidal rule",
  //             "X", "True Solution");
  // fout.close();
  //
  // fout.open("err.m");
  // matlabPlot(fout, x, error, "FEM with p.w. linears and composite trapezoidal rule",
  //            "X", "Error");
  // fout.close();

  double maximum_error = 0.0;
  for (int i = 0; i < n; i++){
    if (abs(error(i)) > maximum_error){
       maximum_error = abs(error(i));}
  }
  cout << "The rms error is:";
  cout << rmsError << endl;

  return 0;
}

#endif
