#ifndef leastSquares_Included
#define leastSquares_Included

#include "matrix.h"
#include "gaussElim.h"

enum ls_state {SUCCESS, BAD_DATA};
enum model{LINEAR, QUAD, CUBIC, EXP, POWER};

ls_state incsSystemSolver(Matrix& A, Vector& x, Vector& b); // Solves inconsistent systems
ls_state normalEquationsSetup(model t, Matrix& A, Vector& b, Vector& c); //Setups the inconsistent system for normalEquations
ls_state normalEquationSolver(model t, Matrix& A, Vector& c); //Setups and solve inconsistent syste for normalEquations
ls_state pureIncsSolver(Matrix& A, Vector& x, Vector& b); // Solves inc systems and produceceses 2-norm, se, and rmse.

#endif
