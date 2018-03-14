/*
 * ilpsolverdolloflipcluster.cpp
 *
 *  Created on: 23-feb-2018
 *      Author: M. El-Kebir
 */

#include "ilpsolverdolloflipclustered.h"

IlpSolverDolloFlipClustered::IlpSolverDolloFlipClustered(const Matrix& D,
                                                         int k,
                                                         double alpha,
                                                         double beta,
                                                         int lC,
                                                         const StlIntVector& zC,
                                                         int lT,
                                                         const StlIntVector& zT)
  : IlpSolverDolloFlip(D, k, lT, lC, alpha, beta)
  , _zC(zC)
  , _zT(zT)
  , _lC(lC)
  , _lT(lT)
  , _bestCaseLikelihoods()
  , _worstCaseLikelihoods()
  , _L()
{
  const int m = _D.getNrTaxa();
  const int n = _D.getNrCharacters();
  
  const double log_alpha = log(_alpha);
  const double log_1_minus_alpha = log(1 - _alpha);
  const double log_beta = log(_beta);
  const double log_1_minus_beta = log(1 - _beta);
  
  _worstCaseLikelihoods = StlDoubleVector(n, 0);
  _bestCaseLikelihoods = StlDoubleVector(n, 0);
  for (int c = 0; c < n; ++c)
  {
    for (int p = 0; p < m; p++)
    {
      int d_pc = _D.getEntry(p, c);
      if (d_pc == -1)
      {
        continue;
      }
      
      if (d_pc == 0)
      {
        _worstCaseLikelihoods[c] += std::min(log_beta, log_1_minus_alpha);
        _bestCaseLikelihoods[c] += std::max(log_beta, log_1_minus_alpha);
      }
      else
      {
        _worstCaseLikelihoods[c] += std::min(log_alpha, log_1_minus_beta);
        _bestCaseLikelihoods[c] += std::max(log_alpha, log_1_minus_beta);
      }
    }
  }
}

void IlpSolverDolloFlipClustered::initVariables()
{
  IlpSolverDolloFlip::initVariables();
  
  const int n = _D.getNrCharacters();

  char buf[1024];
  
  _L = IloNumVarArray(_env, n);
  for (int c = 0; c < n; ++c)
  {
    snprintf(buf, 1024, "L_%d", c);
    _L[c] = IloNumVar(_env, _worstCaseLikelihoods[c], _bestCaseLikelihoods[c], buf);
  }
  
  _loss = IloBoolVarArray(_env, _lC);
  for (int f = 0; f < _lC; ++f)
  {
    snprintf(buf, 1024, "loss_%d", f);
    _loss[f] = IloBoolVar(_env, buf);
  }
}

void IlpSolverDolloFlipClustered::initMaxLoss(int maxL)
{
  
}

void IlpSolverDolloFlipClustered::initConstraints()
{
  IlpSolverDolloFlip::initConstraints();
  
  const int m = _D.getNrTaxa();
  const int n = _D.getNrCharacters();
  
  const double log_alpha = log(_alpha);
  const double log_1_minus_alpha = log(1 - _alpha);
  const double log_beta = log(_beta);
  const double log_1_minus_beta = log(1 - _beta);
  
  IloExpr sum(_env);

  // We cannot erase a mutation
  // TODO: is this necessary?
  for (int p = 0; p < m; p++)
  {
    for (int c = 0; c < n; ++c)
    {
      int d_pc = _D.getEntry(p, c);
      if (d_pc == 1)
      {
        for (int f = 0; f < _lC; ++f)
        {
          for (int i = 1; i <= _k + 1; ++i)
          {
            sum += _E[p][f][i];
          }
        }
        _model.add(sum >= 1);
        sum.clear();
      }
    }
  }
  
  for (int f = 0; f < _lC; ++f)
  {
    for (int i = 2; i <= _k + 1; ++i)
    {
      for (int p = 0; p < m; p++)
      {
        _model.add(_loss[f] >= _E[p][f][i]);
        sum += _E[p][f][i];
      }
      _model.add(_loss[f] <= sum);
      sum.clear();
    }
  }
  
  IloExpr obj(_env);
  for (int c = 0; c < n; ++c)
  {
    const int f = _zC[c];
    for (int p = 0; p < m; p++)
    {
      int d_pc = _D.getEntry(p, c);
      assert(d_pc == 0 || d_pc == 1 || d_pc == -1);
      
      if (d_pc == 0)
      {
        for (int j = 0; j <= _k + 1; ++j)
        {
          if (j == 0)
          {
            obj += log_1_minus_alpha * _E[p][f][j];
          }
          else if (j == 1)
          {
            obj += log_beta * _E[p][f][j];
          }
          else
          {
            obj += log_1_minus_alpha * _E[p][f][j];
          }
        }
      }
      else if (d_pc == 1)
      {
        for (int j = 0; j <= _k + 1; ++j)
        {
          if (j == 0)
          {
            obj += log_alpha * _E[p][f][j];
          }
          else if (j == 1)
          {
            obj += log_1_minus_beta * _E[p][f][j];
          }
          else
          {
            obj += log_alpha * _E[p][f][j];
          }
        }
      }
    }
    _model.add(_L[c] == obj);
    obj.clear();
  }
  
  // symmetry breaking, not sure how this works!!!
  IloExpr prevSum(_env);
  for (int f = 0; f < _lC; ++f)
  {
    for (int p = 0; p < m; p++)
    {
      for (int j = 0; j <= _k + 1; ++j)
      {
        if (j == 0)
        {
          sum += _E[p][f][j] * 1.0 / (m * m);
        }
        else if (j == 1)
        {
          sum += _E[p][f][j];
        }
        else
        {
          sum += _E[p][f][j] * 1.0 / m;
        }
      }
    }
    
    if (f > 0)
    {
      _model.add(prevSum >= sum);
    }
    
    std::swap(prevSum, sum);
    sum.clear();
  }
  
  //  // TODO: remove me
  //  _model.add(_E[0][0][1] == 1);
  //  _model.add(_E[1][0][1] == 1);
  //  _model.add(_E[2][0][2] == 1);
  //  _model.add(_E[3][0][1] == 1);
  //
  //  _model.add(_E[0][1][1] == 1);
  //  _model.add(_E[1][1][0] == 1);
  //  _model.add(_E[2][1][1] == 1);
  //  _model.add(_E[3][1][0] == 1);
  //
  //  _model.add(_E[0][2][1] == 1);
  //  _model.add(_E[1][2][0] == 1);
  //  _model.add(_E[2][2][0] == 1);
  //  _model.add(_E[3][2][0] == 1);
  //
  //  _model.add(_E[0][3][0] == 1);
  //  _model.add(_E[1][3][0] == 1);
  //  _model.add(_E[2][3][1] == 1);
  //  _model.add(_E[3][3][0] == 1);
  //
  //  _model.add(_phi[0][0] == 1);
  //  _model.add(_phi[1][1] == 1);
  //  _model.add(_phi[2][2] == 1);
  //  _model.add(_phi[3][3] == 1);
}

void IlpSolverDolloFlipClustered::initHotStart(const Matrix& E,
                                               const StlIntVector& z)
{
  const int m = _D.getNrTaxa();
  const int n = _D.getNrCharacters();
  
  assert(E.getNrTaxa() == m);
  assert(E.getNrCharacters() == _lC);
  assert(z.size() == n);
  
  IloNumVarArray startVar(_env);
  IloNumArray startVal(_env);
  
  // set _E
  for (int p = 0; p < m; ++p)
  {
    for (int f = 0; f < _lC; ++f)
    {
      for (int i = 0; i <= _k + 1; ++i)
      {
        startVar.add(_E[p][f][i]);
        startVal.add(i == E.getEntry(p, f) ? 1 : 0);
      }
    }
  }
  
  _cplex.addMIPStart(startVar, startVal);
  startVar.end();
  startVal.end();
}

void IlpSolverDolloFlipClustered::initObjective()
{
  const int m = _D.getNrTaxa();
  const int n = _D.getNrCharacters();
  
  const double log_alpha = log(_alpha);
  const double log_1_minus_alpha = log(1 - _alpha);
  const double log_beta = log(_beta);
  const double log_1_minus_beta = log(1 - _beta);
  
  IloExpr obj(_env);
  for (int c = 0; c < n; ++c)
  {
    obj += _L[c];
  }
  
  StlIntVector nrCharactersPerCluster(_lC, 0);
  for (int c = 0; c < n; ++c)
  {
    ++nrCharactersPerCluster[_zC[c]];
  }
  
  // minimize number of losses
  IloExpr sum(_env);
  for (int f = 0; f < _lC; ++f)
  {
      sum += nrCharactersPerCluster[f] * _loss[f] * 1. / _lC;
  }
  
  double unit = std::max(log_alpha, std::max(log_beta, std::max(log_1_minus_alpha, log_1_minus_beta)));
  _model.add(IloMaximize(_env, 1000 * obj + sum * unit));
}
