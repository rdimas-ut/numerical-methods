#ifndef HEAT_INCLUDED
#define HEAT_INCLUDED
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
  u       A matrix (2-D vector) of size m X (n+1).  The solution.
RETURN
  heat_state   The return status.
            HEAT_SUCCESS
            HEAT_BAD_DATA
******************************************************************************/
#include "matrix.h"

enum heat_state {HEAT_SUCCESS, HEAT_BAD_DATA};

typedef double func(double);
typedef double func2(double,double);

// FORWARD EULER
heat_state heatFE(Matrix& u, double c, func2 f, func g, 
                  double a, double b, double T);

// BACKWARD EULER
heat_state heatBE(Matrix& u, double c, func2 f, func g, 
                  double a, double b, double T);

//CRANK_NICOLSON
heat_state heatCN(Matrix& u, double c, func2 f, func g, 
                  double a, double b, double T);

#endif
