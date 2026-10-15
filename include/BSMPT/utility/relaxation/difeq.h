// https://github.com/NJdevPro/Numerical-Recipes

#pragma once

#include <BSMPT/utility/relaxation/data_structures.h>

struct Difeq
{
  virtual void smatrix(const int k,
                       const int k1,
                       const int k2,
                       const int jsf,
                       VecInt &indexv,
                       MatDoub &s,
                       MatDoub &y) = 0;
};
