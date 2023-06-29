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

    vector<double> x(3);
    x[0] = 0; x[1] = M_PI/4; x[2] = M_PI/2;
    vector<double> y(3);
    y[0] = 1; y[1] = 1/(sqrt(2)); y[2] = 0;
    vector<double> coefs(3);
    
    newtonDD(coefs,x,y);
    
    int n;
    cout << "Enter number of t: ";
    cin >> n;
    double a[n];
    double s[n];
    
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
        } else if (value < M_PI) {
            //pi - value, neg
            value = M_PI - value;
            s[j] = - (newtonEval(value,coefs,x));
        } else if (value < (M_PI + (M_PI/2))) {
            //value - pi, neg
            value = value - M_PI;
            s[j] = - (newtonEval(value,coefs,x));
        } else {
            //2pi - value, pos
            value = (2*M_PI - value);
            s[j] = (newtonEval(value,coefs,x));
        }
        
    }
    cout << "x  " << "cos(x)   " << "cos1(x)   "<< " error" <<"\n"; 
    for (int j=0; j<n;j++){
        
        cout << a[j] << ", " << cos(a[j])
        << ", " << s[j]<< ", " << abs(cos(a[j]) - s[j]) << "\n";
    }
    
  
  return 0;
}
#endif
