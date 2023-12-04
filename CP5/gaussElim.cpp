#include "gaussElim.h"
#include <cmath>

// #define TEST
//#define LARGE_TEST
ge_state invFactoredMatrix(const Matrix& a, const Permutation& p, Matrix& invA){
  int n = a.n(0);
  Matrix perM(n,n);
  invA = a;
  invA *= -1;

  //Sets up invA
  for (int i=0;i<n;i++){
    for (int j=0;j<n;j++){
      if (i < j){ invA(i,j) = 0;}
      if (i == j){ invA(i,j) = 1;}
    }
  }
  // Finds L^{-1}
  for (int j = 1; j < n-1;j++){
    for (int i = 2; i < n;i++){
      for (int k = 0; k < j;k++){
        invA(i,k) += (invA(i-1,k)*invA(i,j));
        perM(i,k) = invA(i,k);
      }
    }
  }
  //Finds the product of U^{-1}L^{-1}
  for (int c = n-1; c > 0;c--){
    for (int r = c-1; r > -1;r--){
      double mult =  0 - (a(r,c)/a(c,c));
      for (int s = 0; s < n;s++){
        invA(r,s) += invA(c,s)*mult;
        perM(r,s) = invA(r,s);
      }
    }
  }
  for (int p = 0; p < n ;p++){
    for (int q = 0; q < n;q++){
      invA(p,q) = (invA(p,q)/a(p,p));
      perM(p,q) = invA(p,q);
      }
    }
  //Finds the product of U^{-1}L^{-1}P
  for (int x = 0; x < n;x++){
   for (int y = 0; y < n;y++){
    invA(y,x) = perM(y,p(x));
   }
  }
  return GE_SUCCESS;
}
static void swapRows(Matrix& a, int i, int j) {
  // E_i <--> E_j
  for(int k=0; k<a.n(1); k++) {
    double aa = a(i,k);
    a(i,k) = a(j,k);
    a(j,k) = aa;
  }
}

static void swap(Vector& v, int i, int j) {
  double vv = v(i);
  v(i) = v(j);
  v(j) = vv;
}

static void rowReplacement(Matrix& a, int i, int j) {
  // E_j --> E_j - (a(j,i)/a(i,i))E_i
  double f = a(j,i)/a(i,i);
  a(j,i) = f;
  for(int k=i+1; k<a.n(1); k++) {
    a(j,k) -= f*a(i,k);
  }
}

ge_state luFactorize(Matrix& a, Permutation& p) {
  int n = a.n(0);
  int i,j,k;

  if(a.n(1) != n || p.n() != n) return GE_BADDATA; // error (wrong sizes)

  // Set up permutation
  p.identity();

  // Determine scale factors for scaled partial pivoting
  Vector s(n);
  for(int i=0; i<n; i++) {
    s(i) = fabs(a(i,0));
    for(int j=1; j<n; j++) {
      if( s(i) < fabs(a(i,j)) ) s(i) = fabs(a(i,j));
    }
  }

  // Loop on Columns
  for(j=0; j<n; j++) {

    // Get nonzero pivot (use scaled partial pivoting)
    double pivot = fabs(a(j,j))/s(j);
    i = j;
    for(k=j+1; k<n; k++) {
      double q = fabs(a(k,j))/s(k);
      if(q > pivot) {
	pivot = q;
        i = k;
      }
    }
    if(pivot == 0) return GE_SINGULAR;
    if(i != j) {
      swapRows(a,i,j);
      p.swap(i,j);
      swap(s,i,j);
    }

    // Loop on rows
    for(i=j+1; i<n; i++) rowReplacement(a,j,i);
  }
  return GE_SUCCESS;
}

ge_state luSolve(const Matrix& a, const Permutation& p, Vector& x) {
  int n = a.n(0);
  int i,j,k;

  if(a.n(1) != n || p.n() != n || x.n() != n) return GE_BADDATA;

  // Apply permutation to x
  p.permute(x);

  // FORWARD SUBSTITUTION

  // Loop on columns and rows of a
  for(j=0; j<n; j++) {
    for(i=j+1; i<n; i++) {
      x(i) -= a(i,j)*x(j);
    }
  }

  // BACKWARD SUBSTITUTION

  for(i=n-1; i>=0; i--) {
    for(j=i+1; j<n; j++) {
      x(i) -= a(i,j)*x(j);
    }
    x(i) /= a(i,i);
  }

  return GE_SUCCESS;
}

ge_state solve(Matrix& a, Permutation& p, Vector& x) {
  ge_state s = luFactorize(a,p);
  if(s != GE_SUCCESS) return s;
  return luSolve(a,p,x);
}



#ifdef TEST
#include <iostream>
#include <cmath>
using namespace std;

int main(void) {
  Matrix a(3,3);
  Matrix aa(3,3);
  Vector x(3);
  Vector b(3);
  Permutation p(3);

  cout << "Enter a (3x3): ";
  cin >> a;

  cout << "Enter b: ";
  cin >> b;

  x = b;
  aa = a;

  ge_state s = solve(a,p,x);

  switch(s) {
  case GE_SINGULAR:
    cout << "ERROR: Singular matrix" << endl;
    return 1;
  case GE_BADDATA:
    cout << "ERROR: Bad data (wrong sizes)" << endl;
    return 1;
  case GE_SUCCESS:
    cout << "The solution is: ("
	 << x(0) << ", " << x(1) << ", " << x(2) << ")\n";

    cout << "The permutation is: ("
	 << p(0) << ", " << p(1) << ", " << p(2) << ")\n";
    cout << "The parity is: " << p.parity() << "\n";

    cout << "The factored matrix is: / "
	 << a(0,0) << ", " << a(0,1) << ", " << a(0,2) << " \\\n";
    cout << "                        | "
	 << a(1,0) << ", " << a(1,1) << ", " << a(1,2) << " |\n";
    cout << "                        \\ "
	 << a(2,0) << ", " << a(2,1) << ", " << a(2,2) << " /\n";

    // Compute aa*x - b (= 0?)
    Vector r(3);
    matVecMult(aa,x,r);
    r -= b;

    cout << "The norm of the error is: " << l2norm(r) << endl;

    return 0;
  }
  return 2;
}
#endif

#ifdef LARGE_TEST
#include <iostream>
#include <cmath>
#include <ctime>
using namespace std;

int main(void) {
  int n;
  cout << "Enter n: ";
  cin >> n;

  cout << "SETUP:" << endl;

  Matrix a(n,n);
  Vector x(n);
  Permutation p(n);

  for(int j=0; j<n; j++) {
    for(int i=0; i<n; i++) {
      a(i,j) = 0;
    }
    a(j,j) = 1;
    x(j) = j;
  }

  cout << "SOLVE:" << endl;

  clock_t t = clock();
  ge_state s = solve(a,p,x);

  switch(s) {
  case GE_SINGULAR:
    cout << "ERROR: Singular matrix" << endl;
    return 1;
  case GE_BADDATA:
    cout << "ERROR: Bad data (wrong sizes)" << endl;
    return 1;
  case GE_SUCCESS:
    t = clock() - t;
    cout << "FINISHED.  Total time = " << ((float)t)/CLOCKS_PER_SEC << endl;
    return 0;
  }
  return 2;
}
#endif
