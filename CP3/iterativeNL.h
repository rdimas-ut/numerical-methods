#ifndef IterativeNL_Included
#define IterativeNL_Included

#include "matrix.h"
#include "gaussElim.h"

enum inl_state {INL_SUCCESS=0, INL_WONT_STOP, INL_BAD_DATA, INL_BAD_ITERATE};

inl_state fixedpt(void g(const Vector&,Vector&),
        Vector& x, double tol, int& iter);

inl_state newton(void f(const Vector&,Vector&), void df(const Vector&,Matrix&),
        Vector& x, double tol, int& iter);

#endif
