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
                             double alpha,
                             double beta)
  : ColumnGen(B, k)
  , _alpha(alpha)
  , _beta(beta)
{
}

void ColumnGenFlip::initFixedEntriesConstraints()
{
}

void ColumnGenFlip::initObjective()
{
  const int m = _B.getNrTaxa();
  const int n = _B.getNrCharacters();
  
  const double log_alpha = log(_alpha);
  const double log_1_minus_alpha = log(1 - _alpha);
  const double log_beta = log(_beta);
  const double log_1_minus_beta = log(1 - _beta);
  
  IloExpr obj(_env);
  IloExpr x(_env);
  IloExpr y(_env);
  for (int p = 0; p < m; p++)
  {
    for (int c = 0; c < n; c++)
    {
      int b_pc = _B.getEntry(p, c);
      if (b_pc == 0)
      {
        for (int j = 0; j <= _k + 1; ++j)
        {
          if (j == 0)
          {
            obj += log_1_minus_alpha * _A[p][c][j];
          }
          else if (j == 1)
          {
            obj += log_beta * _A[p][c][j];
          }
          else
          {
            obj += log_1_minus_alpha * _A[p][c][j];
          }
        }
      }
      else if (b_pc == 1)
      {
        for (int j = 0; j <= _k + 1; ++j)
        {
          if (j == 0)
          {
            obj += log_alpha * _A[p][c][j];
          }
          else if (j == 1)
          {
            obj += log_1_minus_beta * _A[p][c][j];
          }
          else
          {
            obj += log_alpha * _A[p][c][j];
          }
        }
      }
    }
  }
  
  _model.add(IloMaximize(_env, 1000 * obj));
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
    if (_k >= 1 && !_activeVariables[p][c][2])
    {
      _activeVariables[p][c][2] = true;
      _A[p][c][2].setUB(1);
      ++_nrActiveVariables;
    }
  }
  else if (i == 1)
  {
    if (!_activeVariables[p][c][0])
    {
      _activeVariables[p][c][0] = true;
      _A[p][c][0].setUB(1);
      ++_nrActiveVariables;
    }
    if (_k >= 1 && !_activeVariables[p][c][2])
    {
      _activeVariables[p][c][2] = true;
      _A[p][c][2].setUB(1);
      ++_nrActiveVariables;
    }
  }
  else if (i < _k + 1)
  {
    if (!_activeVariables[p][c][i + 1])
    {
      _activeVariables[p][c][i + 1] = true;
      _A[p][c][i + 1].setUB(1);
      ++_nrActiveVariables;
    }
  }
}

