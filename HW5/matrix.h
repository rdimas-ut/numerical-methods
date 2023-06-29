#ifndef Matrix_Included
#define Matrix_Included

#include <iostream>

class Vector {
 private:
  int my_size;
  double* array;

 public:
  Vector(int n) { my_size = n;  array = new double[my_size]; };
  ~Vector() { delete[] array; }

  int size(int=0) const { return my_size; }

  double operator() (int i) const { return array[i]; }
  double& operator() (int i) { return array[i]; }

  Vector& operator=(const Vector&);
  Vector& operator+=(const Vector&);
  Vector& operator-=(const Vector&);
  Vector& operator*=(double);        // multiplication by a scalar
  Vector& operator/=(double);        // division by a scalar
};

std::istream& operator>>(std::istream&, Vector&);
std::ostream& operator<<(std::ostream&, const Vector&);

class Matrix {
 private:
  int my_size[2];
  double* array;

 public:
  Matrix(int n0, int n1) { my_size[0] = n0; my_size[1] = n1;
    array = new double[my_size[0]*my_size[1]]; };
  ~Matrix() { delete[] array; }

  int size(int i) const { return my_size[i]; }

  double operator() (int i, int j) const { return array[i + my_size[0]*j]; }
  double& operator() (int i, int j) { return array[i + my_size[0]*j]; }

  Matrix& operator=(const Matrix&);
  Matrix& operator+=(const Matrix&);
  Matrix& operator-=(const Matrix&);
  Matrix& operator*=(double);         // multiplication by a scalar
  Matrix& operator/=(double);         // division by a scalar
};

std::istream& operator>>(std::istream&, Matrix&);
std::ostream& operator<<(std::ostream&, const Matrix&);

class Permutation {
 private:
  int my_size;
  int* array;
  int my_parity;

 public:
  Permutation(int n) { my_size = n;  array = new int[my_size]; identity(); };
  ~Permutation() { delete[] array; }

  int size(int=0) const { return my_size; }

  int operator() (int i) const { return array[i]; }

  void identity();
  void swap(int i, int j);
  double parity() const { return my_parity; }

  void permute(Vector& b) const;
};

#endif
