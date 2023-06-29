#include <iostream>
#include <cmath>
using namespace std;

double arcLength (double df(double), const double a, const double b, const int m){
  int two_m = 2*m;
  double x = a;
  double even, odd, initial, last, rvalue;
  double h = (b - a)/two_m;
  even = 0; odd = 0;
  initial = sqrt(1 + pow(df(x),2));
  x += h;

  for (int i = 1; i < two_m; i++){
    if (i % 2 == 0){
      even += sqrt(1 + pow(df(x),2));
      x += h;
    } else {
      odd += sqrt(1 + pow(df(x),2));
      x += h;
    }
  }
  last = sqrt(1 + pow(df(x),2));

  rvalue = (h/3)*(initial + last + (4*odd + 2*even));

  return rvalue;
}

double dfa(double x) { return 3*x*x; }
double dfb(double x) { return 1/(cos(x)*cos(x)); }
double dfc(double x) { return 1/(1+(x*x)); }

int main(){
  double a, b, re;
  int m;

  cout << "Enter values of a, b for function a: ";
  cin >> a >> b;
  cout << "Enter number of planes m: ";
  cin >> m;
  re = arcLength(dfa, a, b, m);
  cout << "ArcLength is: " << re << endl;

  cout << "Enter values of a, b for function b: ";
  cin >> a >> b;
  cout << "Enter number of planes m: ";
  cin >> m;
  re = arcLength(dfb, a, b, m);
  cout << "ArcLength is: " << re << endl;

  cout << "Enter values of a, b for function c: ";
  cin >> a >> b;
  cout << "Enter number of planes m: ";
  cin >> m;
  re = arcLength(dfc, a, b, m);
  cout << "ArcLength is: " << re << endl;

  return 0;
}
