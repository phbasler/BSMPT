#include <BSMPT/utility/relaxation/data_structures.h>

VecInt::VecInt() : nn(0), v(nullptr)
{
}

VecInt::VecInt(const size_t n) : nn(n), v(n > 0 ? new int[n] : nullptr)
{
}

VecInt::VecInt(const size_t n, const int &a)
    : nn(n)
    , v(n > 0 ? new int[n] : nullptr)
{
  for (size_t i = 0; i < nn; i++)
    v[i] = a;
}

VecInt::VecInt(const VecInt &a) : nn(a.nn), v(nn > 0 ? new int[nn] : nullptr)
{
  for (size_t i = 0; i < nn; i++)
    v[i] = a[i];
}

void VecInt::print()
{
  for (size_t i = 0; i < nn; i++)
    std::cout << v[i] << "\t";
  std::cout << "\n";
}

size_t VecInt::size()
{
  return nn;
}

VecInt &VecInt::operator=(const VecInt &a)
{
  if (this != &a)
  {
    if (nn != a.nn)
    {
      if (v != nullptr) delete[] (v);
      nn = a.nn;
      v  = nn > 0 ? new int[nn] : nullptr;
    }
    for (size_t i = 0; i < nn; i++)
      v[i] = a[i];
  }
  return *this;
}

int &VecInt::operator[](const size_t i)
{
  return v[i];
}

const int &VecInt::operator[](const size_t i) const
{
  return v[i];
}

VecInt::~VecInt()
{
  if (v != nullptr) delete[] (v);
}
VecDoub::VecDoub() : nn(0), v(nullptr)
{
}

VecDoub::VecDoub(const size_t n) : nn(n), v(n > 0 ? new double[n] : nullptr)
{
}

VecDoub::VecDoub(const size_t n, const double &a)
    : nn(n)
    , v(n > 0 ? new double[n] : nullptr)
{
  for (size_t i = 0; i < nn; i++)
    v[i] = a;
}

VecDoub::VecDoub(const VecDoub &a)
    : nn(a.nn)
    , v(nn > 0 ? new double[nn] : nullptr)
{
  for (size_t i = 0; i < nn; i++)
    v[i] = a[i];
}

void VecDoub::print()
{
  for (size_t i = 0; i < nn; i++)
    std::cout << v[i] << "\t";
  std::cout << "\n";
}

size_t VecDoub::size()
{
  return nn;
}

void VecDoub::resize(const size_t newn)
{
  if (newn != nn)
  {
    if (v != nullptr)
    {
      delete[] (v);
    }
    nn = newn;
    v  = nn > 0 ? new double[nn] : nullptr;
  }
}

void VecDoub::zero()
{
  for (size_t i = 0; i < nn; i++)
    v[i] = 0.;
}

VecDoub &VecDoub::operator=(const VecDoub &a)
{
  if (this != &a)
  {
    if (nn != a.nn)
    {
      if (v != nullptr) delete[] (v);
      nn = a.nn;
      v  = nn > 0 ? new double[nn] : nullptr;
    }
    for (size_t i = 0; i < nn; i++)
      v[i] = a[i];
  }
  return *this;
}

double &VecDoub::operator[](const size_t i)
{
  return v[i];
}

const double &VecDoub::operator[](const size_t i) const
{
  return v[i];
}

VecDoub::~VecDoub()
{
  if (v != nullptr) delete[] (v);
}

MatDoub::MatDoub() : nrows(0), ncols(0), v(nullptr)
{
}

MatDoub::MatDoub(const size_t n, const size_t m)
    : nrows(n)
    , ncols(m)
    , v(n > 0 ? new double *[n] : nullptr)
{
  size_t i, nel = m * n;
  if (v) v[0] = nel > 0 ? new double[nel] : nullptr;
  for (i = 1; i < n; i++)
    v[i] = v[i - 1] + m;
}

MatDoub::MatDoub(const size_t n, const size_t m, const bool is_unit)
    : nrows(n)
    , ncols(m)
    , v(n > 0 ? new double *[n] : nullptr)
{
  (void)is_unit;
  size_t i, j, nel = m * n;
  if (v) v[0] = nel > 0 ? new double[nel] : nullptr;
  for (i = 1; i < n; i++)
    v[i] = v[i - 1] + m;
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      v[i][j] = i == j ? 1. : 0.;
}
MatDoub::MatDoub(const size_t n, const size_t m, const double &a)
    : nrows(n)
    , ncols(m)
    , v(n > 0 ? new double *[n] : nullptr)
{
  size_t i, j, nel = m * n;
  if (v) v[0] = nel > 0 ? new double[nel] : nullptr;
  for (i = 1; i < n; i++)
    v[i] = v[i - 1] + m;
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      v[i][j] = a;
}

MatDoub::MatDoub(const MatDoub &a)
    : nrows(a.nrows)
    , ncols(a.ncols)
    , v(nrows > 0 ? new double *[nrows] : nullptr)
{
  size_t i, j, nel = nrows * ncols;
  if (v) v[0] = nel > 0 ? new double[nel] : nullptr;
  for (i = 1; i < nrows; i++)
    v[i] = v[i - 1] + ncols;
  for (i = 0; i < nrows; i++)
    for (j = 0; j < ncols; j++)
      v[i][j] = a[i][j];
}

void MatDoub::print()
{
  size_t i, j;
  for (i = 0; i < ncols; i++)
  {
    std::cout << "------";
  }
  std::cout << "\n";
  for (i = 0; i < nrows; i++)
  {
    for (j = 0; j < ncols; j++)
    {
      std::cout << v[i][j] << std::setw(2) << " ";
    }
    std::cout << "\n";
  }
}

size_t MatDoub::rows()
{
  return nrows;
}

size_t MatDoub::cols()
{
  return ncols;
}

void MatDoub::resize(const size_t newn, const size_t newm)
{
  size_t i, nel;
  if (newn != nrows || newm != ncols)
  {
    if (v != nullptr)
    {
      delete[] (v[0]);
      delete[] (v);
    }
    nrows = newn;
    ncols = newm;
    v     = nrows > 0 ? new double *[nrows] : nullptr;
    nel   = nrows * ncols;
    if (v) v[0] = nel > 0 ? new double[nel] : nullptr;
    for (i = 1; i < nrows; i++)
      v[i] = v[i - 1] + ncols;
  }
}

void MatDoub::zero()
{
  for (size_t i = 0; i < nrows; i++)
    for (size_t j = 0; j < ncols; j++)
      v[i][j] = 0.;
}

MatDoub &MatDoub::operator=(const MatDoub &a)
{
  size_t i, j, nel;
  if (this != &a)
  {
    if (nrows != a.nrows || ncols != a.ncols)
    {
      if (v != nullptr)
      {
        delete[] (v[0]);
        delete[] (v);
      };
      nrows = a.nrows;
      ncols = a.ncols;
      v     = nrows > 0 ? new double *[nrows] : nullptr;
      nel   = nrows * ncols;
      if (v) v[0] = nel > 0 ? new double[nel] : nullptr;
      for (i = 1; i < nrows; i++)
        v[i] = v[i - 1] + ncols;
    }
    for (i = 0; i < nrows; i++)
      for (j = 0; j < ncols; j++)
        v[i][j] = a[i][j];
  }
  return *this;
}

/* MatDoub &MatDoub::operator=(MatDoub a)
{
  size_t i, j, nel;
  if (this != &a)
  {
    if (nrows != a.nrows || ncols != a.ncols)
    {
      if (v != nullptr)
      {
        delete[] (v[0]);
        delete[] (v);
      };
      nrows = a.nrows;
      ncols = a.ncols;
      v     = nrows > 0 ? new double *[nrows] : nullptr;
      nel   = nrows * ncols;
      if (v) v[0] = nel > 0 ? new double[nel] : nullptr;
      for (i = 1; i < nrows; i++)
        v[i] = v[i - 1] + ncols;
    }
    for (i = 0; i < nrows; i++)
      for (j = 0; j < ncols; j++)
        v[i][j] = a[i][j];
  }
  return *this;
} */
double *MatDoub::operator[](const size_t i)
{
  return v[i];
}

const double *MatDoub::operator[](const size_t i) const
{
  return v[i];
}

MatDoub::~MatDoub()
{
  if (v != nullptr)
  {
    delete[] v[0]; // free the contiguous data
    delete[] v;    // free row pointers
  }
}

Mat3DDoub::Mat3DDoub() : nn(0), mm(0), kk(0), v(nullptr)
{
}

Mat3DDoub::Mat3DDoub(const size_t n, const size_t m, const size_t k)
    : nn(n)
    , mm(m)
    , kk(k)
    , v(new double **[n])
{
  size_t i, j;
  v[0]    = new double *[n * m];
  v[0][0] = new double[n * m * k];
  for (j = 1; j < m; j++)
    v[0][j] = v[0][j - 1] + k;
  for (i = 1; i < n; i++)
  {
    v[i]    = v[i - 1] + m;
    v[i][0] = v[i - 1][0] + m * k;
    for (j = 1; j < m; j++)
      v[i][j] = v[i][j - 1] + k;
  }
}

double **Mat3DDoub::operator[](const size_t i)
{
  return v[i];
}

const double *const *Mat3DDoub::operator[](const size_t i) const
{
  return v[i];
}

size_t Mat3DDoub::dim1()
{
  return nn;
}
size_t Mat3DDoub::dim2()
{
  return mm;
}
size_t Mat3DDoub::dim3()
{
  return kk;
}

Mat3DDoub::~Mat3DDoub()
{
  if (v != nullptr)
  {
    delete[] (v[0][0]);
    delete[] (v[0]);
    delete[] (v);
  }
}