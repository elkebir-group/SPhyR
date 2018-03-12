/*
 * columngenflipclustered.cpp
 *
 *  Created on: 10-mar-2018
 *      Author: M. El-Kebir
 */

#include "columngenflipclustered.h"

ColumnGenFlipClustered::ColumnGenFlipClustered(const Matrix& B,
                                               int k,
                                               bool lazy,
                                               double alpha,
                                               double beta,
                                               int l,
                                               const StlIntVector& z)
  : ColumnGenFlip(B, l, k, lazy, alpha, beta)
  , _l(l)
  , _z(z)
{
}

void ColumnGenFlipClustered::initActiveVariables()
{
  const int m = _B.getNrTaxa();
  const int n = _B.getNrCharacters();
  
  _nrActiveVariables = 0;
  _activeVariables = StlBool3Matrix(m, StlBoolMatrix(_n, StlBoolVector(_k + 2, false)));
  
  const double log_alpha = log(_alpha);
  const double log_1_minus_alpha = log(1 - _alpha);
  const double log_beta = log(_beta);
  const double log_1_minus_beta = log(1 - _beta);
  
  for (int p = 0; p < m; p++)
  {
    for (int f = 0; f < _l; ++f)
    {
      int count0 = 0;
      int count1 = 0;
      for (int c = 0; c < n; c++)
      {
        if (_z[c] == f)
        {
          if (_B.getEntry(p, c) == 0)
          {
            ++count0;
          }
          else if (_B.getEntry(p, c) == 1)
          {
            ++count1;
          }
        }
      }
      
      if (log_1_minus_beta * count0 + log_alpha * count1 > log_beta * count0 + log_1_minus_alpha * count1)
      {
        _activeVariables[p][f][0] = true;
        ++_nrActiveVariables;
      }
      else
      {
        _activeVariables[p][f][1] = true;
        ++_nrActiveVariables;
      }
    }
  }
  updateVariableBounds();
}

void ColumnGenFlipClustered::initObjective()
{
  const int m = _B.getNrTaxa();
  const int n = _B.getNrCharacters();
  
  const double log_alpha = log(_alpha);
  const double log_1_minus_alpha = log(1 - _alpha);
  const double log_beta = log(_beta);
  const double log_1_minus_beta = log(1 - _beta);
  
  for (int c = 0; c < n; ++c)
  {
    const int f = _z[c];
    for (int p = 0; p < m; p++)
    {
      int d_pc = _B.getEntry(p, c);
      assert(d_pc == 0 || d_pc == 1 || d_pc == -1);
      
      if (d_pc == 0)
      {
        for (int j = 0; j <= _k + 1; ++j)
        {
          if (j == 0)
          {
            _obj += log_1_minus_alpha * _A[p][f][j];
          }
          else if (j == 1)
          {
            _obj += log_beta * _A[p][f][j];
          }
          else
          {
            _obj += log_1_minus_alpha * _A[p][f][j];
          }
        }
      }
      else if (d_pc == 1)
      {
        for (int j = 0; j <= _k + 1; ++j)
        {
          if (j == 0)
          {
            _obj += log_alpha * _A[p][f][j];
          }
          else if (j == 1)
          {
            _obj += log_1_minus_beta * _A[p][f][j];
          }
          else
          {
            _obj += log_alpha * _A[p][f][j];
          }
        }
      }
    }
  }
  
  _obj *= 1000;
  _model.add(IloMaximize(_env, _obj));
}
