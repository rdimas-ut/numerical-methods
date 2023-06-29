#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "leastSquares.h"
using namespace std;

int main() {
  // Input and options
  int n;int m;int options=-1;;
  ls_state s;
  model op [5] = {LINEAR, QUAD, CUBIC, EXP, POWER};
  model t=EXP;
  bool datafile = false;
  bool modeling = false;
  bool inconsistent = false;

  string ny1 = "ny"; string ny2 = "ny";
  //Selects mode: normal equation modeling or inconsistent systems
  while (!(modeling || inconsistent)) {
    if (!(ny1=="n"||ny1=="y")) {
      cout << "Modeling with Normal Equations? y/n : " << flush;
      cin >> ny1;
      modeling = (ny1=="y");
    }
    if (!(ny1=="n"||ny1=="y")) {
      continue;
    }

    if ((!(ny2=="n"||ny2=="y")) && (ny1=="n"||ny1=="y") && (!modeling)) {
      cout << "Solving inconsistent system? y/n : " << flush;
      cin >> ny2;
      inconsistent = (ny2=="y");
    }
    if (!(ny2=="n"||ny2=="y")) {
      continue;
    }

    if (inconsistent == false)
    {modeling = true; inconsistent=true;}
  }

  while (modeling && (!(options < 5 && options > -1))){
    cout << "Select model from, linear(0), quadratic(1), cubic(2), exp(3), and power(4): " << flush;
    cin >> options;
  }

  if (modeling && inconsistent)
  {cout << "Okay, goodbye!" << endl;return 0;}

  ny1 = "ny"; ny2 = "ny";

  //Input data for system or normal equations
  cout << "Enter size of n by m Matrix A: " << flush;
  cin >> n >> m;

  Matrix A(m,n);
  Vector c(n);
  Vector b(m);

  while (datafile == false) {
    cout << "Enter data by file? y/n : " << flush;
    cin >> ny1;
    if (modeling && ny1=="n") {
      cout << "Enter model data by rows: " << flush;
      cin >> A;
      datafile = true;
    } else if (inconsistent && ny1=="n") {
      cout << "Enter A by rows: " << flush;
      cin >> A;
      datafile = true;
    } else if (ny1=="y") {
      string file_name = "null";
      cout << "Input file name: " << flush;
      cin >> file_name;
      ifstream matrix_files;
      matrix_files.open(file_name, ifstream::in);
      matrix_files >> A;
      matrix_files.close();
      datafile = true;
    }
  }

  datafile = false;
  while ((!datafile) && inconsistent) {
    cout << "Enter b by file? y/n : " << flush;
    cin >> ny2;
    if (ny2 == "n") {
      cout << "Enter b: " << flush;
      cin >> b;
      datafile = true;
    } else if (ny2 == "y") {
      string file_name = "null";
      cout << "Input file name: " << flush;
      cin >> file_name;
      ifstream matrix_files;
      matrix_files.open(file_name, ifstream::in);
      matrix_files >> b;
      matrix_files.close();
      datafile = true;
    }
  }

  //Models data by normal equations or solves inconsistent system.
  if (modeling) {
    s = normalEquationSolver(op[options], A, c);
  } else if (inconsistent) {
    s = pureIncsSolver(A, c, b);
  }

  if (s == SUCCESS){
    return 0;
  } else {
    return 1;
  }
}
