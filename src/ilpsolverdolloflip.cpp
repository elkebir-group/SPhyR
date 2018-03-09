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
                                       int n,
                                       double alpha,
                                       double beta)
  : IlpSolverDollo(D, k, n)
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
  
  IloExpr obj(_env);
  IloExpr x(_env);
  IloExpr y(_env);
  for (int p = 0; p < m; p++)
  {
    for (int c = 0; c < n; c++)
    {
      int d_pc = _D.getEntry(p, c);
      if (d_pc == 0)// || d_pc == -1)
      {
        for (int j = 0; j <= _k + 1; ++j)
        {
          if (j == 0)
          {
            obj += log_1_minus_alpha * _E[p][c][j];
          }
          else if (j == 1)
          {
            obj += log_beta * _E[p][c][j];
            y += _E[p][c][j];
          }
          else
          {
            obj += log_1_minus_alpha * _E[p][c][j];
          }
        }
      }
      else if (d_pc == 1)
      {
        for (int j = 0; j <= _k + 1; ++j)
        {
          if (j == 0)
          {
            obj += log_alpha * _E[p][c][j];
            x += _E[p][c][j];
          }
          else if (j == 1)
          {
            obj += log_1_minus_beta * _E[p][c][j];
          }
          else
          {
            obj += log_alpha * _E[p][c][j];
            x += _E[p][c][j];
          }
        }
      }
    }
  }
//  obj.clear();
//  double lambda = getLambda(_alpha, _beta);
//  obj += lambda * x;
//  obj += (1 - lambda) * y;
  
  double unit = std::max(log_alpha, std::max(log_beta, std::max(log_1_minus_alpha, log_1_minus_beta)));
//  std::cout << "Lambda = " << lambda << std::endl;
//  double unit = std::min(lambda, 1 - lambda);
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
  
  _model.add(IloMaximize(_env, obj));
//  _model.add(IloMinimize(_env, obj));
}
