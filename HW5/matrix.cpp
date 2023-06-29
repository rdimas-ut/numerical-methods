#include "matrix.h"

// Vector /////////////////////////////////////////////////////////////////////
Vector& Vector::operator=(const Vector& v) {
  if(my_size != v.my_size)
    std::cerr << "Sizes differ in Vector::operator=(const Vector&)";
  for(int i=0; i<my_size; i++) array[i] = v.array[i];
  return *this;
}

Vector& Vector::operator+=(const Vector& v) {
  if(my_size != v.my_size)
    std::cerr << "Sizes differ in Vector::operator+=(const Vector&)";
  for(int i=0; i<my_size; i++) array[i] += v.array[i];
  return *this;
}

Vector& Vector::operator-=(const Vector& v) {
  if(my_size != v.my_size)
    std::cerr << "Sizes differ in Vector::operator-=(const Vector&)";
  for(int i=0; i<my_size; i++) array[i] -= v.array[i];
  return *this;
}

Vector& Vector::operator*=(double x) {
  for(int i=0; i<my_size; i++) array[i] *= x;
  return *this;
}

Vector& Vector::operator/=(double x) {
  for(int i=0; i<my_size; i++) array[i] /= x;
  return *this;
}

std::istream& operator>>(std::istream & vin, Vector& v) {
  for(int i=0; i<v.size(); i++) vin >> v(i);
  return vin;
}

std::ostream& operator<<(std::ostream & vout, const Vector& v) {
  vout << v(0);
  for(int i=1; i<v.size(); i++) vout << "  " << v(i);
  vout << "\n";
  return vout;
}

// Matrix /////////////////////////////////////////////////////////////////////

Matrix& Matrix::operator=(const Matrix& m) {
  if(my_size[0] != m.my_size[0] || my_size[1] != m.my_size[1])
    std::cerr << "Sizes differ in Matrix::operator=(const Matrix&)";
  for(int i=0; i<my_size[0]*my_size[1]; i++) array[i] = m.array[i];
  return *this;
}

Matrix& Matrix::operator+=(const Matrix& m) {
  if(my_size[0] != m.my_size[0] || my_size[1] != m.my_size[1])
    std::cerr << "Sizes differ in Matrix::operator+=(const Matrix&)";
  for(int i=0; i<my_size[0]*my_size[1]; i++) array[i] += m.array[i];
  return *this;
}

Matrix& Matrix::operator-=(const Matrix& m) {
  if(my_size[0] != m.my_size[0] || my_size[1] != m.my_size[1])
    std::cerr << "Sizes differ in Matrix::operator-=(const Matrix&)";
  for(int i=0; i<my_size[0]*my_size[1]; i++) array[i] -= m.array[i];
  return *this;
}

Matrix& Matrix::operator*=(double x) {
  for(int i=0; i<my_size[0]*my_size[1]; i++) array[i] *= x;
  return *this;
}

Matrix& Matrix::operator/=(double x) {
  for(int i=0; i<my_size[0]*my_size[1]; i++) array[i] /= x;
  return *this;
}

std::istream& operator>>(std::istream & matin, Matrix& a) {
  for(int i=0; i<a.size(0); i++)
  for(int j=0; j<a.size(1); j++) {
    matin >> a(i,j);
  }
  return matin;
}

std::ostream& operator<<(std::ostream & matout, const Matrix& a) {
  for(int i=0; i<a.size(0); i++) {
    matout << a(i,0);
    for(int j=1; j<a.size(1); j++) matout << "  " << a(i,j);
    matout << "\n";
  }
  return matout;
}

// Permutation ////////////////////////////////////////////////////////////////

void Permutation::identity() {
  for(int i=0; i<my_size; i++) array[i] = i;
  my_parity = 1;
}

void Permutation::swap(int i, int j) {
  int k = array[i];
  array[i] = array[j];
  array[j] = k;
  my_parity *= -1;
}

void Permutation::permute(Vector& b) const {
  if(my_size > b.size()) return;

  Vector c(b.size());

  for(int i=0; i<my_size; i++) c(i) = b(array[i]);
  for(int i=0; i<my_size; i++) b(i) = c(i);
}
