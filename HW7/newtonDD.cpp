#include <vector>
using namespace std;

#define TEST


double newtonEval(double t, const vector<double>& coefs, const vector<double>& x) {
  int n = coefs.size();
  double value = coefs[n-1];
  for(int i=n-2; i>=0; i--) {
    value = value * (t - x[i]) + coefs[i];
  }
  return value;
}

// Set up divided difference coefficients
  int newtonDD(vector<double>& coefs, const vector<double>& x, const vector<double>& y) {
  int n = x.size();
  if(y.size() != n) return 1; // error

  coefs.resize(n);

  // DD level 0
  for(int i=0; i<n; i++) coefs[i] = y[i];

  // DD higher levels (bottom to top, overwrite lower entries as they are finished)
  for(int level=1; level<n; level++) {
    for(int i=n-1; i>=level; i--) {
      double dx = x[i] - x[i-level];
      if(dx==0) return 2;
      
      coefs[i] = ( coefs[i]-coefs[i-1] ) / dx;
    }
  }
  return 0;
}

#ifdef TEST
#define _USE_MATH_DEFINES
#include <iostream>
#include <string>
#include <math.h>

int main() {
    double pi_q = M_PI/4;

    vector<double> x(3);
    x[0] = 0; x[1] = M_PI/4; x[2] = M_PI/2;
    vector<double> y(3);
    y[0] = 1; y[1] = 1/(sqrt(2)); y[2] = 0;
    
    vector<double> r(3);
    r[0] = pi_q + (pi_q*cos(M_PI/6)); r[1] = pi_q + (pi_q*cos(M_PI/2)); r[2] = pi_q + (pi_q*cos((5*M_PI)/6));
    vector<double> t(3);
    t[0] = cos(r[0]); t[1] = cos(r[1]); t[2] = cos(r[2]);
    
    vector<double> coefs(3);
    vector<double> coefs2(3);
    
    newtonDD(coefs,x,y);
    newtonDD(coefs2,r,t);
    
    int n;
    cout << "Enter number of t: ";
    cin >> n;
    double a[n];
    double s[n];
    double o[n];
    cout << "Enter t's:";
    
    for (int i=0; i<n;i++){
        cin >> a[i];
    }
    for (int j=0; j<n;j++){
        double value = a[j];
        value = fmod(value, 2*M_PI);
        if (value < 0){value += (2*M_PI);}
        
        if (value < (M_PI/2)) {
            //regular value
            s[j] = newtonEval(value,coefs,x);
            o[j] = newtonEval(value,coefs2,r);
        } else if (value < M_PI) {
            //pi - value, neg
            value = M_PI - value;
            s[j] = - (newtonEval(value,coefs,x));
            o[j] = - (newtonEval(value,coefs2,r));
        } else if (value < (M_PI + (M_PI/2))) {
            //value - pi, neg
            value = value - M_PI;
            s[j] = - (newtonEval(value,coefs,x));
            o[j] = - (newtonEval(value,coefs2,r));
        } else {
            //2pi - value, pos
            value = (2*M_PI - value);
            s[j] = (newtonEval(value,coefs,x));
            o[j] = newtonEval(value,coefs2,r);
        }
        
    }
    cout << "x  " << "cos(x)   " << "cos1(x)   "<< "cos2(x)     " << "error1   " << "error2" <<"\n";
    for (int j=0; j<n;j++){
        cout << a[j] << ", " << cos(a[j])
        << ", " << s[j]<< ", " << o[j] << ", " << abs(cos(a[j]) - s[j])<< ", " <<abs(cos(a[j])-o[j]) << "\n";
    }
  return 0;
}
#endif
