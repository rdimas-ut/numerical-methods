///////////////////////////////////////////////////////////////////////////////
//NAIVE GAUSSIAN
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
using namespace std;
double eps = 2e-52;

int elimination(double a[], double b[], int n) {
    for (int j=0; j < n-1; j++ ) {
        if(fabs(a[j*n + j]) < eps) return -1;
        for (int i = j+1; i < n;i++){
            double mult = a[i*n + j]/a[n*j + j];
            for (int k = 0; k < n;k++){
                a[i*n + k] = a[i*n + k] - (mult*a[j*n + k]);
            }
            b[i] = b[i] - (mult*b[j]);
        }
    }
    return 0;
}
int back_sub(double x[], double a[], double b[], int n){
    for (int i = n-1; i>-1;i--) {
        for (int j = i+1;j<n;j++){
            b[i] = b[i] - (a[i*n + j]*x[j]);
        }
        x[i] = b[i]/(a[i*n + i]);
    }
    return 0;
}
int main() {
    int n;
    cout << "Enter an n: ";
    cin >> n;

    double a[n*n];
    double b[n];
    double x[n];
    for (int i = 0;i<n;i++){x[i] = 1;}
    
    cout << "Enter a (3x3): ";
    cin >> a[0] >> a[1] >>a[2];
    cin >> a[3] >> a[4] >>a[5];
    cin >> a[6] >> a[7] >>a[8];
    cout << "Enter b: ";  cin >> b[0] >> b[1] >> b[2];
    
    int k = elimination(a, b, n);
    if (k == -1) {cout << "Singular Matrix" << "!";return 0;}
    
    back_sub(x, a, b, n);
    cout << "x=" <<x[0] << ", ";
    cout << "y=" <<x[1] << ", ";
    cout << "z=" <<x[2] << ", ";
    return 0;
}

