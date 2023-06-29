#include <iostream>
#include <math.h>
#include "heat.h"
using namespace std;
#define _USE_MATH_DEFINES

/*
double uTrue(double x, double t) { return exp(-M_PI*M_PI*t)*sin(M_PI*x) + x*(1-x); }
double f(double x, double t) { return 2; }
double g(double x) { return uTrue(x,0); };
*/
double uTrue(double x, double t) { return exp(-M_PI*t)*sin(M_PI*x); }
double f(double x, double t) { return 0; }
double g(double x) { return uTrue(x,0); };

int main() {
  int n,m;
  double a, b, T;
  double c = (1/M_PI);

  cout << "Enter m n a b T: ";
  cin >> m >> n >> a >> b >> T;

  if ((a >= b) || (m < 1)) {
    cout << "ERROR: Invalid interval [" << a << "," << b
         << "] or number of space steps " << m << endl;
    return 1;
  }
  if ((T <= 0) || (n < 1)) {
    cout << "ERROR: Invalid final time " << T
         << " or number of timesteps " << n << endl;
    return 1;
  }
  Matrix uFE(m+1,n+1);
  Matrix uBE(m+1,n+1);
  Matrix uCN(m+1,n+1);

  heat_state s;
  s = heatFE(uFE, c, f, g, a, b, T);
  switch(s) {
  case HEAT_BAD_DATA:
    cout << "ERROR: Bad data heatFE" << endl;
    return 1;
  default:
    break;
  }

  s = heatBE(uBE, c, f, g, a, b, T);
  switch(s) {
  case HEAT_BAD_DATA:
    cout << "ERROR: Bad data heatBE" << endl;
    return 1;
  default:
    break;
  }

  s = heatCN(uCN, c, f, g, a, b, T);
  switch(s) {
  case HEAT_BAD_DATA:
    cout << "ERROR: Bad data heatCN" << endl;
    return 1;
  default:
    break;
  }

  double h = (b-a)/m;
  double k = T/n;
  double errorFE = 0, errorBE = 0, errorCN = 0;
  double x = a;
  for(int i=1; i<m; i++) {
    x += h;
    double t = 0;
    for(int j=1; j<=n; j++) {
      t += k;
      double uu = uTrue(x,t), error;
      error = fabs(uFE(i,j) - uu);  if(errorFE < error) errorFE = error;
      error = fabs(uBE(i,j) - uu);  if(errorBE < error) errorBE = error;
      error = fabs(uCN(i,j) - uu);  if(errorCN < error) errorCN = error;
    }
  }

  char ans;

  cout << "Print uFE? [y/n/f/a]: ";
  cin >> ans;
  if(ans == 'y') {
    cout << uFE << endl;
  } else if (ans == 'a'){
    cout << uFE(3, n) << endl;
  } else if(ans == 'f') {
    int i=0;
    while(1) {
      cout << uFE(i,n);
      i++;
      if(i>m) break;
      if(!(i%5)) { cout << "\n"; } else { cout << "  "; }
    }
    cout << "\n";
  }

  cout << "Print uBE? [y/n/f/a]: ";
  cin >> ans;
  if(ans == 'y') {
    cout << uBE << endl;
  } else if (ans == 'a'){
    cout << uBE(3, n) << endl;
  } else if(ans == 'f') {
    int i=0;
    while(1) {
      cout << uBE(i,n);
      i++;
      if(i>m) break;
      if(!(i%5)) { cout << "\n"; } else { cout << "  "; }
    }
    cout << "\n";
  }

  cout << "Print uCN? [y/n/f/a]: ";
  cin >> ans;
  if(ans == 'y') {
    cout << uCN << endl;
  } else if (ans == 'a'){
    cout << uCN(12, n) << endl;
  } else if(ans == 'f') {
    int i=0;
    while(1) {
      cout << uCN(i,n);
      i++;
      if(i>m) break;
      if(!(i%5)) { cout << "\n"; } else { cout << "  "; }
    }
    cout << "\n";
  }

  cout << "Error FE: " << errorFE << endl;
  cout << "Error BE: " << errorBE << endl;
  cout << "Error CN: " << errorCN << endl;

  return 0;
}
// 10 25 1 0 1 .25
