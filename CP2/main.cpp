#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#include "iterativeLA.h"

int main() {
  int n;
  cout << "Enter size of n by n Matrix A: " << flush;
  cin >> n;

  Matrix A(n,n);
  Vector x(n);
  Vector b(n);

  x = 0.;

  string yes_no = "empty";
  bool A_in_met = false;
  while (A_in_met == false) {
    cout << "Enter data by file? y/n : " << flush;
    cin >> yes_no;
    if (yes_no == "n") {
      cout << "Enter A by rows: " << flush;
      cin >> A;
      A_in_met = true;
    } else if (yes_no == "y") {
      string file_name = "null";
      cout << "Input file name: " << flush;
      cin >> file_name;
      ifstream matrix_files;
      matrix_files.open(file_name, ifstream::in);
      matrix_files >> A;
      matrix_files.close();
      A_in_met = true;

    }
  }
  yes_no = "empty";
  A_in_met = false;
  while (A_in_met == false) {
    cout << "Enter b by file? y/n : " << flush;
    cin >> yes_no;
    if (yes_no == "n") {
      cout << "Enter b: " << flush;
      cin >> b;
      A_in_met = true;
    } else if (yes_no == "y") {
      A_in_met = true;
      string file_name = "null";
      cout << "Input file name: " << flush;
      cin >> file_name;
      ifstream matrix_files;
      matrix_files.open(file_name, ifstream::in);
      matrix_files >> b;
      matrix_files.close();
      A_in_met = true;
    }
  }

  int maxIter;
  double tolerance;
  cout << "Enter maxIter and tolerance: " << flush;
  cin >> maxIter >> tolerance;

  ila_state s;

  string method = "Method";
  string precond = "Precond";
  bool valid_input = false;
  bool valid_precond = false;
  while (valid_input == false) {
    cout << "Please select Jacobi, Gauss-Seidel, SOR, or CCG: "
    << flush;
    cin >> method;
    if (method == "Jacobi"){
      valid_input = true;
      s = jacobi(A,b,x,maxIter,tolerance);
    } else if (method == "Gauss-Seidel") {
      valid_input = true;
      s = gauss_seidel(A,b,x,maxIter,tolerance);
    } else if (method == "SOR") {
      double omega = 0;
      cout << "Input the omega value:";
      cin >> omega;
      s = sor(A,b,x,maxIter,tolerance, omega);
      valid_input = true;
    } else if (method == "CCG") {
      valid_input = true;
      bool ccg_bool = false;
      string ccg_resp = "null";
      while (ccg_bool == false) {
        cout << "Do you wish to use a preconditioner? (y/n) : ";
        cin >> ccg_resp;
        if (ccg_resp == "y"){
          ccg_bool = true;
        } else if (ccg_resp == "n") {
          ccg_bool = true;
          valid_precond = true;
          s = conj_grad(A,b,x,maxIter,tolerance);
        }
      }
      while (valid_precond == false){
        cout << "Please select Jacobi, Gauss-Seidel, or SOR preconditioner: "
        << flush;
        cin >> precond;
      if (precond == "Jacobi"){
        valid_precond = true;
        s = prec_conj_grad(jac_prec,A,b,x,0.0,maxIter,tolerance);
      } else if (precond == "Gauss-Seidel") {
        valid_precond = true;
        s = prec_conj_grad(ssor_prec,A,b,x,1.0,maxIter,tolerance);
      } else if (precond == "SOR") {
        valid_precond = true;
        double omega = 0;
        cout << "Input the omega value:";
        cin >> omega;
        s = prec_conj_grad(ssor_prec,A,b,x,omega,maxIter,tolerance);
      }
    }
  }
  }

  int prec = (int) (log10(1.0/tolerance));
  cout.precision(prec);

  switch(s) {
  case ILA_WONT_STOP:
    cout << "ERROR: Exceeded maximum number of iterations." << endl;
    return 1;
  case ILA_BAD_DIAGONAL:
    cout << "ERROR: A diagonal entry of A was 0." << endl;
    return 1;
  default:
    cout << "ERROR: Unspecified." << endl;
    return 1;
  case ILA_SUCCESS:
    cout << "The solution is:" << endl;
    cout << x << endl;

    Vector y(n);
    matVecMult(A,x,y);
    y -= b;
    cout << "The number of iterations is: " << maxIter << endl;
    cout << "The max-norm of residual is: " << maxNorm(y) << endl;
    cout << "The residual is: " << endl;
    cout << y << endl;
    return 0;
  }
}
