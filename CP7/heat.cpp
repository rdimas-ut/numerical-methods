/******************************************************************************
We solve the heat equation for u(x,t):

  u_t - c u_xx = f(x,t),   a <= x <= b, 0 < t <= T
  u(a,t) = u(b,t) = 0
  u(x,y) = g(x)

using finite differences (forward euler, backward euler, and
crank-nicolson).  We use m equally spaced points in the x-direction,
and n equally spaced points in time.

INPUTS
  c       The coefficient in the PDE.  Note that c > 0 is required.
  f,g     The names of the functions in the PDE.
  a,b     The limits of x (a <= x <= b).
  T       The final time (0 < t <= T).
OUTPUTS
  u       A matrix (2-D vector) of size (m+1) X (n+1).  The solution.
RETURN
  heat_state  The return status.
          HEAT_SUCCESS
          HEAT_BAD_DATA
******************************************************************************/
#include "heat.h"
#include <math.h>

// FORWARD EULER

heat_state heatFE(Matrix& u, double c, func2 f, func g, 
                  double a, double b, double T) {
  int i,j;
  double x,t;

  if(c <= 0) return HEAT_BAD_DATA;

  int m = u.n(0)-1; // number of grid subintervals in x
  int n = u.n(1)-1; // number of time steps

  if(m <= 0 || n < 0) return HEAT_BAD_DATA;

  double h = (b-a)/m; // x point spacing
  double k = T/n;   // time step size

  // SET INITIAL AND BOUNDARY VALUES

  x = a;
  for(i=1; i<m; i++) {
    x += h;
    u(i,0) = g(x);
  }

  for(j=0; j<=n; j++) {
    u(0,j) = 0;
    u(m,j) = 0;
  }

  // SET UP LINEAR SYSTEM

  double sigma = c*k/h/h;  // sub and superdiagonals
  double tau = 1 - 2*sigma; // diagonals

  // LOOP OVER TIME INTERVALS

  t = 0;
  for(int j=1; j<=n; j++) {
    x = a;
    for(i=1; i<m; i++) {
      x += h;
      u(i,j) = tau*u(i,j-1) + sigma*( u(i+1,j-1) + u(i-1,j-1) ) + k*f(x,t);
    }
    t += k;
  }

  return HEAT_SUCCESS;
}

// BACKWARD EULER

heat_state heatBE(Matrix& u, double c, func2 f, func g, 
                  double a, double b, double T) {
  int i,j;
  double x,t;

  if(c <= 0) return HEAT_BAD_DATA;

  int m = u.n(0)-1; // number of grid subintervals in x
  int n = u.n(1)-1; // number of time steps

  if(m <= 0 || n < 0) return HEAT_BAD_DATA;

  double h = (b-a)/m; // x point spacing
  double k = T/n;   // time step size

  // SET INITIAL AND BOUNDARY VALUES

  x = a;
  for(i=1; i<m; i++) {
    x += h;
    u(i,0) = g(x);
  }

  for(j=0; j<=n; j++) {
    u(0,j) = 0;
    u(m,j) = 0;
  }

  // SET UP TRIDIAGONAL LINEAR SYSTEM MATRIX

  Matrix A(3,m);  // Note: A(.,0), A(0,1), and A(2,m-1) are not used.
  const int subd=0; // subdiagonal -- a_{i,i-1}
  const int diag=1; // diagonal -- a_{i,i}
  const int supd=2; // superdiagonal -- a_{i,i+1}

  double sigma = c*k/h/h;  // negative of sub and superdiagonals
  double tau = 1 + 2*sigma; // diagonals

  for(int i=1; i<m; i++) {
    A(subd,i) = -sigma;
    A(diag,i) =  tau;
    A(supd,i) = -sigma;
  }

  // FACTOR MATRIX

  A(diag,1) = 1/A(diag,1);  // Store inverse of diagonal.  Eliminates to 1.
  A(supd,1) *= A(diag,1);

  for(int i=2; i<m; i++) {
    // Eliminate subd entry
    A(diag,i) -= A(subd,i)*A(supd,i-1);

    // rescale so diagonal is 1
    A(diag,i) = 1/A(diag,i);  // Store inverse of diagonal.  Eliminates to 1.
    A(supd,i) *= A(diag,i);
  }

  // LOOP OVER TIME INTERVALS AND SOLVE

  t = 0;
  for(int j=1; j<=n; j++) {
    t += k;

    // Set up RHS
    x = a;
    for(i=1; i<m; i++) {
      x += h;
      u(i,j) = u(i,j-1) + k*f(x,t);
    }

    // Solve tridiagonal system using gaussian elimination w/o pivoting
    u(1,j) *= A(diag,1);
    for(i=2; i<m; i++) {
      u(i,j) -= A(subd,i)*u(i-1,j);
      u(i,j) *= A(diag,i);
    }

   for(int i=m-2; i>= 1; i--) {
      u(i,j) -= u(i+1,j)*A(supd,i);
    }
  }

  return HEAT_SUCCESS;
}

// CRANK_NICOLSON

heat_state heatCN(Matrix& u, double c, func2 f, func g, 
                  double a, double b, double T) {
  int i,j;
  double x,t;

  if(c <= 0) return HEAT_BAD_DATA;

  int m = u.n(0)-1; // number of grid subintervals in x
  int n = u.n(1)-1; // number of time steps

  if(m <= 0 || n < 0) return HEAT_BAD_DATA;

  double h = (b-a)/m; // x point spacing
  double k = T/n;   // time step size

  // SET INITIAL AND BOUNDARY VALUES

  x = a;
  for(i=1; i<m; i++) {
    x += h;
    u(i,0) = g(x);
  }

  for(j=0; j<=n; j++) {
    u(0,j) = 0;
    u(m,j) = 0;
  }

  // SET UP TRIDIAGONAL LINEAR SYSTEM MATRIX

  Matrix A(3,m);  // Note: A(.,0), A(0,1), and A(2,m-1) are not used.
  const int subd=0; // subdiagonal -- a_{i,i-1}
  const int diag=1; // diagonal -- a_{i,i}
  const int supd=2; // superdiagonal -- a_{i,i+1}

  double sigma = c*k/h/h;  // twice negative of sub and superdiagonals
  double tau = 1 + sigma; // diagonals

  for(int i=1; i<m; i++) {
    A(subd,i) = -sigma/2;
    A(diag,i) =  tau;
    A(supd,i) = -sigma/2;
  }

  // FACTOR MATRIX

  A(diag,1) = 1/A(diag,1);  // Store inverse of diagonal.  Eliminates to 1.
  A(supd,1) *= A(diag,1);

  for(int i=2; i<m; i++) {
    // Eliminate subd entry
    A(diag,i) -= A(subd,i)*A(supd,i-1);

    // rescale so diagonal is 1
    A(diag,i) = 1/A(diag,i);  // Store inverse of diagonal.  Eliminates to 1.
    A(supd,i) *= A(diag,i);
  }

  // LOOP OVER TIME INTERVALS AND SOLVE

  t = 0;
  for(int j=1; j<=n; j++) {
    t += k;

    // Set up RHS
    x = a;
    for(i=1; i<m; i++) {
      x += h;
      u(i,j) = (1 - sigma)*u(i,j-1)
	+ sigma*( u(i+1,j-1) + u(i-1,j-1) )/2 + k*f(x,t-k/2);
    }

    // Solve tridiagonal system using gaussian elimination w/o pivoting
    u(1,j) *= A(diag,1);
    for(i=2; i<m; i++) {
      u(i,j) -= A(subd,i)*u(i-1,j);
      u(i,j) *= A(diag,i);
    }

   for(int i=m-2; i>= 1; i--) {
      u(i,j) -= u(i+1,j)*A(supd,i);
    }
  }

  return HEAT_SUCCESS;
}
