// https://github.com/NJdevPro/Numerical-Recipes

#pragma once

#include <BSMPT/utility/relaxation/difeq.h>
#include <fstream>

class RelaxOde
{
private:
  MatDoub &y;
  const int ne, nb1, nb2, nenb1, ne2, m;
  VecInt kmax;
  VecDoub ermax;
  Mat3DDoub c;
  MatDoub s;
  bool VERBOSE = false;

public:
  double err;

  RelaxOde(const int itmax,
           const double conv,
           const double slowc,
           const VecDoub &scalv,
           VecInt &indexv,
           const int nbb,
           MatDoub &yy,
           Difeq &difeq)
      : y(yy)
      , ne(y.rows())
      , nb1(nbb)
      , nb2(ne - nb1)
      , nenb1(ne + nb1)
      , ne2(2 * ne)
      , m(y.cols())
      , kmax(ne)
      , ermax(ne)
      , c(ne, ne - nb1 + 1, m + 1)
      , s(ne, 2 * ne + 1, 0.)
  {
    int jv, k, nvars = ne * m;
    if (VERBOSE) std::cout << std::setw(8) << "Iter.";
    if (VERBOSE)
      std::cout << std::setw(10) << "Error" << std::setw(10) << "FAC"
                << std::endl;
    for (int it = 0; it < itmax; it++)
    {
      k = 0;
      difeq.smatrix(k, 0, m, ne2, indexv, s, y);
      pinvs(nb2, ne, ne, 0, 0);
      for (k = 1; k < m; k++)
      {
        int kp = k;
        difeq.smatrix(k, 0, m, ne2, indexv, s, y);
        red(ne, 0, nb1, nb1, ne, kp);
        pinvs(0, ne, nb1, 0, k);
      }
      k = m;
      difeq.smatrix(k, 0, m, ne2, indexv, s, y);
      red(nb2, ne, nenb1, nenb1, ne2, m);
      pinvs(0, nb2, nenb1, nb2, m);
      bksub(nb2, 0, m);
      err = 0.0;
      for (int j = 0; j < ne; j++)
      {
        jv          = indexv[j];
        double errj = 0.0, vmax = 0.0;
        int km = 0;
        for (k = 0; k < m; k++)
        {
          double vz = std::abs(c[jv][0][k]);
          if (vz > vmax)
          {
            vmax = vz;
            km   = k + 1;
          }
          errj += vz;
        }
        err += errj / scalv[j];
        ermax[j] = c[jv][0][km - 1] / scalv[j];
        kmax[j]  = km;
      }
      err /= nvars;
      double fac = (err > slowc ? slowc / err : 1.0);

      for (int j = 0; j < ne; j++)
      {
        jv = indexv[j];
        for (k = 0; k < m; k++)
          y[j][k] -= fac * c[jv][0][k];
      }
      if (VERBOSE) std::cout << std::setw(6) << it;
      if (VERBOSE) std::cout << std::setw(13) << err;
      if (VERBOSE) std::cout << std::setw(12) << fac << std::endl;
      (void)conv;
    }
  }

  /**
   * @brief pinvs diagonalizes the square subsection of s and stores unreduced
   * coefficients
   *
   * @param ie1
   * @param ie2
   * @param je1
   * @param jsf
   * @param jc1
   * @param k
   */
  void pinvs(const int ie1, int ie2, const int je1, const int jc1, const int k)
  {
    int jpiv = 0, jp = 0, jcoff = 0, irow = 0, ipiv = 0, icoff = 0;
    int je2 = 0;
    double pivinv, piv, big;
    const int iesize = ie2 - ie1;
    VecInt indxr(iesize);
    VecDoub pscl(iesize);
    je2 = je1 + iesize;
    for (int i = ie1; i < ie2; i++)
    {
      big = 0.0;
      for (int j = je1; j < je2; j++)
        if (std::abs(s[i][j]) > big) big = std::abs(s[i][j]);
      if (big == 0.0)
      {
        if (VERBOSE)
          std::cout << "Singular matrix in pinvs with k = " << k
                    << " with rows= (" << ie1 << ",\t" << ie2 << ")\n";

        if (VERBOSE)
          std::cout << "Singular matrix in pinvs with k = " << k
                    << " with columns= (" << je1 << ",\t" << je2 << ")\n";
        if (VERBOSE) std::cout << "\n\n";
        for (int ii = ie1; ii < ie2; ii++)
        {
          for (int jj = je1; jj < je2; jj++)
            if (VERBOSE) std::cout << s[ii][jj] << "\t";
          if (VERBOSE) std::cout << "\n";
        }
        if (VERBOSE) std::cout << "\n\n";
        throw("Singular matrix - row all 0, in pinvs");
      }
      pscl[i - ie1]  = 1.0 / big;
      indxr[i - ie1] = 0;
    }
    for (int id = 0; id < iesize; id++)
    {
      piv = 0.0;
      for (int i = ie1; i < ie2; i++)
      {
        if (indxr[i - ie1] == 0)
        {
          big = 0.0;
          for (int j = je1; j < je2; j++)
          {
            if (std::abs(s[i][j]) > big)
            {
              jp  = j;
              big = std::abs(s[i][j]);
            }
          }
          if (big * pscl[i - ie1] > piv)
          {
            ipiv = i;
            jpiv = jp;
            piv  = big * pscl[i - ie1];
          }
        }
      }
      if (piv == 0.0) throw("No pivot found in pinvs");
      if (s[ipiv][jpiv] == 0.0) throw("Singular matrix in routine pinvs");
      indxr[ipiv - ie1] = jpiv + 1;
      pivinv            = 1.0 / s[ipiv][jpiv];
      for (int j = je1; j <= ne2; j++)
        s[ipiv][j] *= pivinv;
      s[ipiv][jpiv] = 1.0;
      for (int i = ie1; i < ie2; i++)
      {
        if (indxr[i - ie1] != jpiv + 1)
        {
          if (s[i][jpiv] != 0.0)
          {
            double dum = s[i][jpiv];
            for (int j = je1; j <= ne2; j++)
              s[i][j] -= dum * s[ipiv][j];
            s[i][jpiv] = 0.0;
          }
        }
      }
    }
    jcoff = jc1 - je2;
    icoff = ie1 - je1;
    for (int i = ie1; i < ie2; i++)
    {
      irow = indxr[i - ie1] + icoff;
      for (int j = je2; j <= ne2; j++)
        c[irow - 1][j + jcoff][k] = s[i][j];
    }
  }

  void bksub(const int jf, const int k1, const int k2)
  {
    int nbf = ne - nb1, im = 1;
    for (int k = k2 - 1; k >= k1; k--)
    {
      if (k == k1) im = nbf + 1;
      int kp = k + 1;
      for (int j = 0; j < nbf; j++)
      {
        double xx = c[j][jf][kp];
        for (int i = im - 1; i < ne; i++)
          c[i][jf][k] -= c[i][j][k] * xx;
      }
    }
    for (int k = k1; k < k2; k++)
    {
      int kp = k + 1;
      for (int i = 0; i < nb1; i++)
        c[i][0][k] = c[i + nbf][jf][k];
      for (int i = 0; i < nbf; i++)
        c[i + nb1][0][k] = c[i][jf][kp];
    }
  }

  /**
   * @brief  red eliminates leading columns of the s matrix using results from
   * prior blocks
   *
   * @param iz1
   * @param iz2
   * @param jz1
   * @param jz2
   * @param jm1
   * @param jm2
   * @param jmf
   * @param ic1
   * @param jc1
   * @param jcf
   * @param kc
   */
  void red(const int iz2,
           const int jz1,
           const int jz2,
           const int jm1,
           const int jm2,
           const int kc)
  {
    int l, j, i;
    double vx;
    int loff = -jm1, ic = nb2;
    for (j = jz1; j < jz2; j++)
    {
      for (l = jm1; l < jm2; l++)
      {
        vx = c[ic][l + loff][kc - 1];
        for (i = 0; i < iz2; i++)
          s[i][l] -= s[i][j] * vx;
      }
      vx = c[ic][nb2][kc - 1];
      for (i = 0; i < iz2; i++)
        s[i][ne2] -= s[i][j] * vx;
      ic += 1;
    }
  }
};