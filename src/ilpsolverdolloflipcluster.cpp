/*
 * ilpsolverdolloflipcluster.cpp
 *
 *  Created on: 23-feb-2018
 *      Author: M. El-Kebir
 */

#include "ilpsolverdolloflipcluster.h"

IlpSolverDolloFlipCluster::IlpSolverDolloFlipCluster(const Matrix& D,
                                                     int k,
                                                     double alpha,
                                                     double beta,
                                                     int l)
  : IlpSolverDolloFlip(D, k, l, alpha, beta)
  , _l(l)
  , _bestCaseLikelihoods()
  , _worstCaseLikelihoods()
  , _L()
  , _solZ(D.getNrCharacters(), 0)
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

void IlpSolverDolloFlipCluster::initVariables()
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
  
  _w = IloBoolVarMatrix(_env, n);
  for (int c = 0; c < n; ++c)
  {
    _w[c] = IloBoolVarArray(_env, _l);
    for (int f = 0; f < _l; ++f)
    {
      snprintf(buf, 1024, "w_%d_%d", c, f);
      _w[c][f] = IloBoolVar(_env, buf);
    }
  }
  
  // TODO: record loss of mutation, and minimize this
//  _loss = IloBoolVarMatrix(_env, );
}

void IlpSolverDolloFlipCluster::initConstraints()
{
  IlpSolverDolloFlip::initConstraints();
  
  const int m = _D.getNrTaxa();
  const int n = _D.getNrCharacters();
  
  const double log_alpha = log(_alpha);
  const double log_1_minus_alpha = log(1 - _alpha);
  const double log_beta = log(_beta);
  const double log_1_minus_beta = log(1 - _beta);
  
  IloExpr sum(_env);
  for (int c = 0; c < n; ++c)
  {
    for (int f = 0; f < _l; ++f)
    {
      sum += _w[c][f];
    }
    _model.add(sum == 1);
    sum.clear();
  }
  
  // We cannot erase a mutation
  for (int p = 0; p < m; p++)
  {
    for (int c = 0; c < n; ++c)
    {
      int d_pc = _D.getEntry(p, c);
      if (d_pc == 1)
      {
        for (int f = 0; f < _l; ++f)
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
  
  IloExpr obj(_env);
  for (int c = 0; c < n; ++c)
  {
    for (int f = 0; f < _l; ++f)
    {
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
      _model.add(_L[c] <= obj + (_bestCaseLikelihoods[c] - _worstCaseLikelihoods[c]) * (1 - _w[c][f]));
      obj.clear();
    }
  }
  
  // symmetry breaking
  IloExpr prevSum(_env);
  for (int f = 0; f < _l; ++f)
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

void IlpSolverDolloFlipCluster::initHotStart(const Matrix& E,
                                             const StlIntVector& z)
{
  const int m = _D.getNrTaxa();
  const int n = _D.getNrCharacters();
  
  assert(E.getNrTaxa() == m);
  assert(E.getNrCharacters() == _l);
  assert(z.size() == n);
  
  IloNumVarArray startVar(_env);
  IloNumArray startVal(_env);
  
  // set _E
  for (int p = 0; p < m; ++p)
  {
    for (int f = 0; f < _l; ++f)
    {
      for (int i = 0; i <= _k + 1; ++i)
      {
        startVar.add(_E[p][f][i]);
        startVal.add(i == E.getEntry(p, f) ? 1 : 0);
      }
    }
  }
  
  // set _w
  for (int c = 0; c < n; ++c)
  {
    for (int f = 0; f < _l; ++f)
    {
      startVar.add(_w[c][f]);
      startVal.add(f == z[c] ? 1 : 0);
    }
  }
  
  // set _L
  const double log_alpha = log(_alpha);
  const double log_1_minus_alpha = log(1 - _alpha);
  const double log_beta = log(_beta);
  const double log_1_minus_beta = log(1 - _beta);
  
  for (int c = 0; c < n; ++c)
  {
    double L = 0;
    for (int p = 0; p < m; p++)
    {
      int d_pc = _D.getEntry(p, c);
      int e_pf = E.getEntry(p, z[c]);
      assert(d_pc == 0 || d_pc == 1 || d_pc == -1);
      
      if (d_pc == 0)
      {
        if (e_pf == 0)
        {
          L += log_1_minus_alpha;
        }
        else if (e_pf == 1)
        {
          L += log_beta;
        }
        else
        {
          L += log_1_minus_alpha;
        }
      }
      else if (d_pc == 1)
      {
        if (e_pf == 0)
        {
          L += log_alpha;
        }
        else if (e_pf == 1)
        {
          L += log_1_minus_beta;
        }
        else
        {
          L += log_alpha;
        }
      }
    }
    startVar.add(_L[c]);
    startVal.add(L);
  }
  
  _cplex.addMIPStart(startVar, startVal);
  startVar.end();
  startVal.end();
}

void IlpSolverDolloFlipCluster::initClusters(const StlIntVector& z)
{
  assert(z.size() == _D.getNrCharacters());
  
  const int n = _D.getNrCharacters();
  for (int c = 0; c < n; ++c)
  {
    assert(0 <= z[c] && z[c] < _l);
    _model.add(_w[c][z[c]] == 1);
  }
}

void IlpSolverDolloFlipCluster::processSolution()
{
  IlpSolverDolloFlip::processSolution();
  
  const int n = _D.getNrCharacters();
  for (int c = 0; c < n; ++c)
  {
    for (int f = 0; f < _l; ++f)
    {
      if (g_tol.nonZero(_cplex.getValue(_w[c][f])))
      {
        _solZ[c] = f;
      }
    }
  }
}

void IlpSolverDolloFlipCluster::printSolution(std::ostream& out) const
{
  IlpSolverDolloFlip::printSolution(out);
  
  const int n = _D.getNrCharacters();
  
  out << std::endl;
  for (int c = 0; c < n; ++c)
  {
    out << "L[" << c << "] = " << _cplex.getValue(_L[c]);
    for (int f = 0; f < _l; ++f)
    {
      if (g_tol.nonZero(_cplex.getValue(_w[c][f])))
      {
        out << ", cluster " << f;
      }
    }
    out << std::endl;
  }
}

void IlpSolverDolloFlipCluster::initObjective()
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
  
  // minimize number of missing entries!
  IloExpr sum(_env);
  for (int p = 0; p < m; p++)
  {
    for (int f = 0; f < _l; ++f)
    {
      for (int j = _k + 1; j >= 1; --j)
      {
        sum += _E[p][f][j] * pow(1./(m *_l), _k + 2 - j);
      }
    }
  }
  
  double unit = std::max(log_alpha, std::max(log_beta, std::max(log_1_minus_alpha, log_1_minus_beta)));
  _model.add(IloMaximize(_env, 1000 * obj + sum * unit));
}
