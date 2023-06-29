#include <iostream>
#include <vector>
#include <cmath>
using namespace std;


//// MISREAD CLAMPED SPLINE
int triDiagSolve(vector<double> diag, vector<double> superDiag, vector<double> subDiag, vector<double>& x){
  int n = diag.size();
  vector<double> cubicCol(n-1);
  double subf;
  double superf;
  double zero = 0;
  double up;
  bool cubic = (superDiag[0] != zero);
  // Lot of trouble error handlind any of this :(
  for (int i = 0; i < n; i++ ) {
    if (diag[i] == zero) {return 1;}
    if (i < n-1){
      subf = subDiag[i]/diag[i];
      if (isnan(subf) == false){
      subDiag[i] -= (diag[i]*subf); diag[i+1] -= (superDiag[i]*subf); x[i+1] -= (x[i]*subf);
      }
    }
    if ((i > 1) and (cubic == false)){
      superf = superDiag[i-1]/diag[i];
      if (isnan(superf) == false){
      superDiag[i-1] -= (diag[i]*superf); x[i-1] -= (x[i]*superf);
      }
    }

    if (cubic and (i > 0)) {
      superf = superDiag[i-1]/diag[i];
      if (i == (n - 1)){ double up = 1;} else {double up = superDiag[i];}
      if (isnan(superf) == false){
      superDiag[i-1] -= (diag[i]*superf); x[i-1] -= (x[i]*superf);
      cubicCol[i-1] = -(superDiag[i] * superf);
      }


      for (int r = i-2; r > -1; r--){
        superf = cubicCol[r]/diag[i];
        if (i == (n - 1)){ double up = 1;} else {double up = superDiag[i];}
        if (isnan(superf) == false){
        cubicCol[r] = -(up*superf);
        x[r] -= (x[i]*superf);
        }
      }
    }
    }
  for (int r = 0; r < n; r++) {
    x[r] /= diag[r];
    cout << x[r];
  }
  cout<< endl;
  return 0;
}

int main(){
  double zero = 0;
  int n;
  cout << "Enter n: ";
  cin >> n;
  vector<double> x(n), y(n), b(n), diag(n), subDiag(n-1), superDiag(n-1), cubic(2), e(n-1), d(n-1);
  string mode = "None";

  // Might wanna add error handling here

  cout << "Enter the x values: ";
  for (int i = 0; i < n; i++) { cin >> x[i]; }
  cout << "Enter the y values: ";
  for (int i = 0; i < n; i++) { cin >> y[i]; }

  while ((mode != "cubic") and (mode != "clamped")){
    cout << "Enter the spline mode, cubic or clamped: ";
    cin >> mode;
  }
  if (mode == "clamped") {
    cout << "Please enter v and w: ";
    for (int i = 0; i < 2;i++){cin >> cubic[i];}
  }

  //Add loop to
  for (int m = 0 ; m < n-2 ; m++){
    subDiag[m] = x[m+1] - x[m];
    diag[m+1] = (2*(x[m+1] - x[m]) + 2*(x[m+2] - x[m+1]));
    superDiag[m+1] = x[m+2] - x[m+1];
    b[m+1] = 3*(((y[m+2] - y[m+1])/(x[m+2] - x[m+1])) - ((y[m+1] - y[m])/(x[m+1] - x[m])));
  }

  if (mode == "clamped") {
    diag[0] = 2*(x[1] - x[0]);
    diag[n-1] = 2*(x[n-1] - x[n-2]);
    superDiag[0] = x[1] - x[0];
    subDiag[1] = x[2] - x[1];
    b[0] = 3*((y[1] - y[0])/(x[1] - x[0]) - cubic[0]);
    b[n-1] = 3*(cubic[1] - ((y[n-1] - y[n-2])/(x[n-1] - x[n-2])));
  } else {
    diag[0] = 1;
    diag[n-1] = 1;
    superDiag[0] = 0;
    subDiag[n-2] = 0;
    b[0] = 0;
    b[n-1] = 0;
  }
  int result;
  result = triDiagSolve(diag, superDiag, subDiag, b);

  for (int s = 0; s < n-1 ;s++) {
    d[s] = ((b[s+1] - b[s])/(3*(x[s+1] - x[s])));
    e[s] = (((y[s+1] - y[s])/(x[s+1] - x[s])) - (((x[s+1] - x[s])/3)*(2*b[s] + b[s+1])));
  }
  cout << d[0] << d[1] << endl;

  vector<string> ssStr(10);  ssStr[0] = "\u2080";ssStr[1] = "\u2081";ssStr[2] =  "\u2082";
                             ssStr[3] = "\u2083";ssStr[4] = "\u2084";ssStr[5] =  "\u2085";
                             ssStr[6] = "\u2086";ssStr[7] = "\u2087";ssStr[8] =  "\u2088";
                             ssStr[9] = "\u2089";

  vector<string> ssStrH(3);  ssStrH[0] = "\u00B9";ssStrH[1] = "\u00B2";ssStrH[2] =  "\u00B3";

  // Displays the equations w
  if (result == 0){
  for (int p = 0; p < n-1; p++) {
    string index = "";
    int zeros = to_string(p).size();
    string xstring;
    for (int i = 0; i < zeros; i++) {int t = (p / (int)pow(10,i)) % 10; index += ssStr[t+1];}
    cout << "S" << index << "(x) = " << y[p];
    for (int w = 0; w < 3; w++) {
        if (x[p] <= 0) {xstring = "(x + ";} else {xstring = "(x - "; }
        if (w == 0 ) {
          if (e[p] < 0){
            cout << " - " << (- e[p]);

          } else {
            cout << " + " <<  e[p];
          }
          cout << xstring << x[p] << ")";
        } else if (w==1){
          if (b[p] < 0){
            cout << " - " <<  (- b[p]);
          } else {
            cout << " + " << b[p];
          }
          cout << xstring << x[p] << ")" <<ssStrH[w];
        } else if (w==2) {
          if (d[p] < 0){
            cout << " - " <<  (- d[p]);

          } else {
            cout << " + " << d[p];
          }
          cout << xstring << x[p] << ")" <<ssStrH[w];
        }
      }
      cout << endl;
  }
} else {cout << "Error, zero pivot" << endl;}


  return 0;
}
