/*
 * columngenflipclustered.cpp
 *
 *  Created on: 10-mar-2018
 *      Author: M. El-Kebir
 */

#include "columngenflipclustered.h"

ColumnGenFlipClustered::ColumnGenFlipClustered(const Matrix& B,
                                               const StlIntMatrix& multiplicities,
                                               double baseL,
                                               int k,
                                               bool lazy,
                                               double alpha,
                                               double beta,
                                               int t,
                                               const StlIntVector& zC,
                                               int s,
                                               const StlIntVector& zT)
  : ColumnGenFlip(B, s, t, k, lazy, alpha, beta)
  , _multiplicities(multiplicities)
  , _baseL(baseL)
  , _zT(zT)
  , _zC(zC)
  , _s(s)
  , _t(t)
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
  
  for (int h = 0; h < _s; h++)
  {
    for (int f = 0; f < _t; f++)
    {
      int count0 = 0; // negative
      int count1 = 0; // positive
      for (int p = 0; p < m; p++)
      {
        if (_zT[p] != h) continue;
          
        for (int c = 0; c < n; c++)
        {
          if (_zC[c] != f) continue;
          
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
      
//      if (g_tol.nonZero(_alpha))
      {
        if (log_1_minus_beta * count0 + log_alpha * count1 > log_beta * count0 + log_1_minus_alpha * count1)
        {
          _activeVariables[h][f][0] = true;
          ++_nrActiveVariables;

          for (int i = 2; i <= _k + 1; ++i)
          {
            _activeVariables[h][f][i] = true;
            ++_nrActiveVariables;
          }
        }
        else
        {
          _activeVariables[h][f][1] = true;
          ++_nrActiveVariables;
        }
      }
//      else
//      {
//        if (count1 > 0)
//        {
//          _activeVariables[h][f][1] = true;
//          ++_nrActiveVariables;
//        }
//        else
//        {
//          _activeVariables[h][f][0] = true;
//          ++_nrActiveVariables;
//
//          for (int i = 2; i <= _k + 1; ++i)
//          {
//            _activeVariables[h][f][i] = true;
//            ++_nrActiveVariables;
//          }
//        }
//      }
    }
  }
  updateVariableBounds();
}

void ColumnGenFlipClustered::initHotStart(const Matrix& E)
{
  IloNumVarArray startVar(_env);
  IloNumArray startVal(_env);
  
  for (int h = 0; h < _s; ++h)
  {
    for (int f = 0; f < _t; ++f)
    {
      for (int i = 0; i <= _k + 1; ++i)
      {
        startVar.add(_A[h][f][i]);
        startVal.add(i == E.getEntry(h, f) ? 1 : 0);
        
        if (E.getEntry(h, f) == i && !_activeVariables[h][f][i])
        {
          _activeVariables[h][f][i] = true;
          _A[h][f][i].setUB(1);
          ++_nrActiveVariables;
        }
      }
    }
  }
  
  _cplex.addMIPStart(startVar, startVal);
  startVar.end();
  startVal.end();
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
    const int f = _zC[c];
    for (int p = 0; p < m; p++)
    {
      const int h = _zT[p];
      int d_pc = _B.getEntry(p, c);
      assert(d_pc == 0 || d_pc == 1 || d_pc == -1);
      
      const int mult = _multiplicities[p][c];
      
      if (d_pc == 0)
      {
        for (int j = 0; j <= _k + 1; ++j)
        {
          if (j == 0)
          {
            _obj += mult * log_1_minus_beta * _A[h][f][j];
          }
          else if (j == 1)
          {
            _obj += mult * log_beta * _A[h][f][j];
          }
          else
          {
            _obj += mult * log_1_minus_beta * _A[h][f][j];
          }
        }
      }
      else if (d_pc == 1)// && g_tol.nonZero(_alpha))
      {
        for (int j = 0; j <= _k + 1; ++j)
        {
          if (j == 0)
          {
            _obj += mult * log_alpha * _A[h][f][j];
          }
          else if (j == 1)
          {
            _obj += mult * log_1_minus_alpha * _A[h][f][j];
          }
          else
          {
            _obj += mult * log_alpha * _A[h][f][j];
          }
        }
      }
    }
  }
  
  // TODO: minimize losses?
  IloExpr lossSum(_env);
  double unit = 0;
  unit = std::max(log_alpha, std::max(log_beta, std::max(log_1_minus_alpha, log_1_minus_beta)));
  
  for (int h = 0; h < _s; h++)
  {
    for (int f = 0; f < _t; f++)
    {
      for (int i = 2; i <= _k + 1; ++i)
      {
        lossSum += _A[h][f][i] * pow(1./(_s * _t), _k + 2 - i);
      }
    }
  }
  
  _obj += _baseL;
  _obj += lossSum * unit;
  _obj *= 1000;
  _model.add(IloMaximize(_env, _obj));
}
