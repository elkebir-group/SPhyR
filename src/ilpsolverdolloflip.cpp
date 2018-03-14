/*
 * ilpsolverdolloflip.cpp
 *
 *  Created on: 25-feb-2018
 *      Author: M. El-Kebir
 */

#include "ilpsolverdolloflip.h"

IlpSolverDolloFlip::IlpSolverDolloFlip(const Matrix& D,
                                       int k,
                                       double alpha,
                                       double beta)
  : IlpSolverDollo(D, k)
  , _alpha(alpha)
  , _beta(beta)
{
}

IlpSolverDolloFlip::IlpSolverDolloFlip(const Matrix& D,
                                       int k,
                                       int m,
                                       int n,
                                       double alpha,
                                       double beta)
  : IlpSolverDollo(D, k, m, n)
  , _alpha(alpha)
  , _beta(beta)
{
}

void IlpSolverDolloFlip::init()
{
  initVariables();
  initConstraints();
  initObjective();
}

void IlpSolverDolloFlip::initObjective()
{
  const int m = _D.getNrTaxa();
  const int n = _D.getNrCharacters();
  
  const double log_alpha = log(_alpha);
  const double log_1_minus_alpha = log(1 - _alpha);
  const double log_beta = log(_beta);
  const double log_1_minus_beta = log(1 - _beta);
//  const double logg_beta = log_beta - log_1_minus_beta;
//  const double logg_alpha = log_alpha - log_1_minus_alpha;
  
  IloExpr obj(_env);
  IloExpr x(_env);   // 0 (input) to 1 (output) FN
  IloExpr y(_env);   // 1 (input) to 0,2,3... (output) FP
  int X = 0;         // zeros in D
  int Y = 0;         // ones in D
  
  for (int p = 0; p < m; p++)
  {
    for (int c = 0; c < n; c++)
    {
      int d_pc = _D.getEntry(p, c);
      if (d_pc == 0)
      {
        ++X;
        for (int j = 0; j <= _k + 1; ++j)
        {
          if (j == 0)
          {
            obj += log_1_minus_alpha * _E[p][c][j];
          }
          else if (j == 1)
          {
            obj += log_beta * _E[p][c][j];
            x += _E[p][c][j];
          }
          else
          {
            obj += log_1_minus_alpha * _E[p][c][j];
          }
        }
      }
      else if (d_pc == 1)
      {
        ++Y;
        for (int j = 0; j <= _k + 1; ++j)
        {
          if (j == 0)
          {
            obj += log_alpha * _E[p][c][j];
            y += _E[p][c][j];
          }
          else if (j == 1)
          {
            obj += log_1_minus_beta * _E[p][c][j];
          }
          else
          {
            obj += log_alpha * _E[p][c][j];
            y += _E[p][c][j];
          }
        }
      }
    }
  }
//  obj.clear();
//  obj += logg_beta * x + logg_alpha * y + X * log_1_minus_beta + Y * log_1_minus_alpha;
  
  double unit = std::max(log_alpha, std::max(log_beta, std::max(log_1_minus_alpha, log_1_minus_beta)));
  for (int p = 0; p < m; ++p)
  {
    for (int c = 0; c < n; ++c)
    {
      for (int i = 2; i <= _k + 1; ++i)
      {
        obj += _E[p][c][i] * (unit / (m * n));
      }
    }
  }
  
  _model.add(IloMaximize(_env, 1000 * obj));
//  _model.add(IloMinimize(_env, obj));
}
