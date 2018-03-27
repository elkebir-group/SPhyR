/*
 * columngenflip.cpp
 *
 *  Created on: 10-mar-2018
 *      Author: M. El-Kebir
 */

#include "columngenflip.h"
#include <ilconcert/ilothread.h>
#include <lemon/time_measure.h>

ColumnGenFlip::ColumnGenFlip(const Matrix& B,
                             int k,
                             bool lazy,
                             double alpha,
                             double beta)
  : ColumnGen(B, k, lazy)
  , _alpha(alpha)
  , _beta(beta)
{
}

ColumnGenFlip::ColumnGenFlip(const Matrix& B,
                             int m,
                             int n,
                             int k,
                             bool lazy,
                             double alpha,
                             double beta)
  : ColumnGen(B, m, n, k, lazy)
  , _alpha(alpha)
  , _beta(beta)
{
}

void ColumnGenFlip::initFixedEntriesConstraints()
{
}

void ColumnGenFlip::initObjective()
{
  const double log_alpha = log(_alpha);
  const double log_1_minus_alpha = log(1 - _alpha);
  const double log_beta = log(_beta);
  const double log_1_minus_beta = log(1 - _beta);
  
  IloExpr x(_env);
  IloExpr y(_env);
  for (int p = 0; p < _m; p++)
  {
    for (int c = 0; c < _n; c++)
    {
      int b_pc = _B.getEntry(p, c);
      if (b_pc == 0)
      {
        for (int j = 0; j <= _k + 1; ++j)
        {
          if (j == 0)
          {
            _obj += log_1_minus_alpha * _A[p][c][j];
          }
          else if (j == 1)
          {
            _obj += log_beta * _A[p][c][j];
          }
          else
          {
            _obj += log_1_minus_alpha * _A[p][c][j];
          }
        }
      }
      else if (b_pc == 1)
      {
        for (int j = 0; j <= _k + 1; ++j)
        {
          if (j == 0)
          {
            _obj += log_alpha * _A[p][c][j];
          }
          else if (j == 1)
          {
            _obj += log_1_minus_beta * _A[p][c][j];
          }
          else
          {
            _obj += log_alpha * _A[p][c][j];
          }
        }
      }
    }
  }
  
  IloExpr lossSum(_env);
  double unit = std::max(log_alpha, std::max(log_beta, std::max(log_1_minus_alpha, log_1_minus_beta)));
  for (int p = 0; p < _m; ++p)
  {
    for (int c = 0; c < _n; ++c)
    {
      for (int i = 2; i <= _k + 1; ++i)
      {
        lossSum += _A[p][c][i] * pow(1./(_m * _n), _k + 2 - i);// (unit / (m * _n));
      }
    }
  }
  
  _obj += lossSum * unit;
  _obj = 1000 * _obj;
  _model.add(IloMaximize(_env, _obj));
}

void ColumnGenFlip::activate(int p, int c, int i)
{
  if (i == 0)
  {
    if (!_activeVariables[p][c][1])
    {
      _activeVariables[p][c][1] = true;
      _A[p][c][1].setUB(1);
      ++_nrActiveVariables;
    }
    for (int j = 2; j <= _k + 1; ++j)
    {
      if (_k >= 1 && !_activeVariables[p][c][j])
      {
        _activeVariables[p][c][j] = true;
        _A[p][c][j].setUB(1);
        ++_nrActiveVariables;
      }
    }
  }
  else if (i == 1)// && g_tol.nonZero(_alpha))
  {
    if (!_activeVariables[p][c][0])
    {
      _activeVariables[p][c][0] = true;
      _A[p][c][0].setUB(1);
      ++_nrActiveVariables;
    }
    for (int j = 2; j <= _k + 1; ++j)
    {
      if (!_activeVariables[p][c][j])
      {
        _activeVariables[p][c][j] = true;
        _A[p][c][j].setUB(1);
        ++_nrActiveVariables;
      }
    }
  }
}

