#pragma once
#include <cmath>

#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

class VecInt
{
private:
  size_t nn;
  int *v;

public:
  VecInt();
  VecInt(const size_t n);
  VecInt(const size_t n, const int &a);
  VecInt(const VecInt &a);
  size_t size();
  void print();
  VecInt &operator=(const VecInt &a);
  int &operator[](const size_t i);
  const int &operator[](const size_t i) const;
  ~VecInt();
};

class VecDoub
{
private:
  size_t nn;
  double *v;

public:
  VecDoub();
  VecDoub(const size_t n);
  VecDoub(const size_t n, const double &a);
  VecDoub(const VecDoub &a);
  size_t size();
  void resize(const size_t newn);
  void zero();
  void print();
  VecDoub &operator=(const VecDoub &a);
  double &operator[](const size_t i);
  const double &operator[](const size_t i) const;
  ~VecDoub();
};

static const bool UNIT = true;

class MatDoub
{
private:
  size_t nrows, ncols;
  double **v;

public:
  MatDoub();
  MatDoub(const size_t n, const size_t m);
  MatDoub(const size_t n, const size_t m, const bool is_unit);
  MatDoub(const size_t n, const size_t m, const double &a);
  MatDoub(const MatDoub &a);
  void print();
  size_t rows();
  size_t cols();
  void resize(const size_t newn, const size_t newm);
  void zero();
  MatDoub &operator=(const MatDoub &a);
  double *operator[](const size_t i);
  const double *operator[](const size_t i) const;
  ~MatDoub();
};

class Mat3DDoub
{
private:
  size_t nn, mm, kk;
  double ***v;

public:
  Mat3DDoub();
  Mat3DDoub(const size_t n, const size_t m, const size_t k);
  double **operator[](const size_t i);
  const double *const *operator[](const size_t i) const;
  size_t dim1();
  size_t dim2();
  size_t dim3();
  ~Mat3DDoub();
};