#ifndef FEM1D_INCLUDED
#define FEM1D_INCLUDED
/******************************************************************************
We approximate
  - (p(x) u')' + q(x) u  = f(x),  a <= x <= b
  u(a) = u(b) = 0
using the Finite Element method with piecewise linears.  We evaluate the
integrals using the composite trapezoidal rule.

INPUTS
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

#include "matrix.h"
#include "matlabPlot.h"

enum fem1d_state {LINFE_SUCCESS, LINFE_BAD_DATA, LINFE_SINGULAR};

typedef double func(double);

fem1d_state linFEM1d(double& alpha, double& beta, const Vector& x, Vector& u, func p, func q, func f);

#endif
