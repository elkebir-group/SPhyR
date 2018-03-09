/*
 * ilpsolverdollo.cpp
 *
 *  Created on: 25-feb-2018
 *      Author: M. El-Kebir
 */

#include "ilpsolverdollo.h"
#include "dollocallback.h"

IlpSolverDollo::IlpSolverDollo(const Matrix& D,
                               int k)
  : _D(D)
  , _k(k)
  , _n(D.getNrCharacters())
  , _env()
  , _model(_env)
  , _cplex(_model)
  , _E(_env)
  , _solE(_D.getNrTaxa(), _n)
{
}

IlpSolverDollo::IlpSolverDollo(const Matrix& D,
                               int k,
                               int n)
  : _D(D)
  , _k(k)
  , _n(n)
  , _env()
  , _model(_env)
  , _cplex(_model)
  , _E(_env)
  , _solE(_D.getNrTaxa(), _n)
{
}

void IlpSolverDollo::init()
{
  initVariables();
  initConstraints();
  initFixedEntries();
  initObjective();
}

void IlpSolverDollo::initVariables()
{
  const int m = _D.getNrTaxa();
  
  char buf[1024];
  
  _E = IloBoolVar3Matrix(_env, m);
  for (int p = 0; p < m; p++)
  {
    _E[p] = IloBoolVarMatrix(_env, _n);
    for (int c = 0; c < _n; ++c)
    {
      _E[p][c] = IloBoolVarArray(_env, _k + 2);
      for (int i = 0; i < _k + 2; ++i)
      {
        snprintf(buf, 1024, "e_%d_%d_%d", p, c, i);
        _E[p][c][i] = IloBoolVar(_env, buf);
      }
    }
  }
}

void IlpSolverDollo::initFixedEntries()
{
  const int m = _D.getNrTaxa();

  IloExpr sum(_env);
  for (int p = 0; p < m; p++)
  {
    for (int c = 0; c < _n; c++)
    {
      int d_pc = _D.getEntry(p, c);
      if (d_pc == 1)
      {
        _model.add(_E[p][c][1] == 1);
      }
      else if (d_pc == 0)
      {
        _model.add(_E[p][c][1] == 0);
      }
    }
  }
}

void IlpSolverDollo::initConstraints()
{
  const int m = _D.getNrTaxa();
  
  IloExpr sum(_env);
  
  // Each entry has a unique state
  for (int p = 0; p < m; p++)
  {
    for (int c = 0; c < _n; c++)
    {
      for (int i = 0; i <= _k + 1; ++i)
      {
        sum += _E[p][c][i];
      }
      _model.add(sum == 1);
      sum.clear();
    }
  }
  
  // disable superfluous states and symmetry breaking
  IloExpr prevSum(_env);
  for (int c = 0; c < _n; c++)
  {
    for (int i = 2; i <= _k + 1; ++i)
    {
      for (int p = 0; p < m; p++)
      {
        sum += _E[p][c][i];
      }
      if (i > 2)
      {
        _model.add(prevSum >= sum);
      }
      
      std::swap(prevSum, sum);
      sum.clear();
    }
  }
  
  // Generalized three gamete condition
//  for (int p = 0; p < m; ++p)
//  {
//    for (int q = 0; q < m; ++q)
//    {
//      if (q == p) continue;
//      for (int r = 0; r < m; ++r)
//      {
//        if (r == p || r == q) continue;
//        for (int c = 0; c < _n; ++c)
//        {
//          for (int d = c + 1; d < _n; ++d)
//          {
//            _model.add((1 - _E[p][c][0]) + _E[q][c][0] + (1 - _E[r][c][0])
//                       + _E[p][d][0] + (1 - _E[q][d][0]) + (1 - _E[r][d][0]) <= 5);
//
//            for (int i = 2; i <= _k + 1; ++i)
//            {
//              _model.add(_E[p][c][i] + (1 - _E[q][c][i]) + _E[r][c][i]
//                         + _E[p][d][0] + (1 - _E[q][d][0]) + (1 - _E[r][d][0]) <= 5);
//            }
//
//            for (int j = 2; j <= _k + 1; ++j)
//            {
//              _model.add((1 - _E[p][c][0]) + _E[q][c][0] + (1 - _E[r][c][0])
//                         + (1 - _E[p][d][j]) + _E[q][d][j] + _E[r][d][j] <= 5);
//            }
//
//            for (int i = 2; i <= _k + 1; ++i)
//            {
//              for (int j = 2; j <= _k + 1; ++j)
//              {
//                _model.add(_E[p][c][i] + (1 - _E[q][c][i]) + _E[r][c][i]
//                           + (1 - _E[p][d][j]) + _E[q][d][j] + _E[r][d][j] <= 5);
//              }
//            }
//          }
//        }
//      }
//    }
//  }
}

void IlpSolverDollo::initObjective()
{
  const int m = _D.getNrTaxa();
  
  IloExpr obj(_env);
  for (int p = 0; p < m; ++p)
  {
    for (int c = 0; c < _n; ++c)
    {
      for (int i = 2; i <= _k + 1; ++i)
      {
        obj += _E[p][c][i];
      }
    }
  }
  
  _model.add(IloMinimize(_env, obj));
}

void IlpSolverDollo::solve(int timeLimit,
                           int memoryLimit,
                           int nrThreads)
{
  const int m = _D.getNrTaxa();
  
  _env.setOut(_env.getNullStream());
  _env.setError(_env.getNullStream());
  _env.setWarning(_env.getNullStream());
  
  _cplex.use(IloCplex::Callback(new (_env) DolloCallback<IloCplex::UserCutCallbackI>(_env, _E, m, _n, _k)));
  _cplex.use(IloCplex::Callback(new (_env) DolloCallback<IloCplex::LazyConstraintCallbackI>(_env, _E, m, _n, _k)));
  
  _cplex.setParam(IloCplex::ParallelMode, -1);
  if (nrThreads > 0)
  {
    _cplex.setParam(IloCplex::Threads, nrThreads);
  }
  if (timeLimit > 0)
  {
    _cplex.setParam(IloCplex::TiLim, timeLimit);
  }
  if (memoryLimit > 0)
  {
    _cplex.setParam(IloCplex::WorkMem, memoryLimit);
  }
  
//  _cplex.exportModel("/tmp/test.lp");
  if (_cplex.solve())
  {
    std::cout << "CPLEX: [" << _cplex.getObjValue() << " , " << _cplex.getObjValue() << "]" << std::endl;
    processSolution();
    std::cout << _solE << std::endl;
  }
}

void IlpSolverDollo::processSolution()
{
  const int m = _D.getNrTaxa();
  for (int p = 0; p < m; p++)
  {
    for (int c = 0; c < _n; ++c)
    {
      for (int i = 0; i <= _k + 1; ++i)
      {
        bool e_pci = g_tol.nonZero(_cplex.getValue(_E[p][c][i]));
        if (e_pci)
        {
          _solE.setEntry(p, c, i);
        }
      }
    }
  }
}

void IlpSolverDollo::printSolution(std::ostream& out) const
{
  const int m = _D.getNrTaxa();
  
  for (int p = 0; p < m; p++)
  {
    bool first = true;
    for (int c = 0; c < _n; ++c)
    {
      if (first)
        first = false;
      else
        out << " ";
      
      for (int i = 0; i <= _k + 1; ++i)
      {
        bool e_pci = g_tol.nonZero(_cplex.getValue(_E[p][c][i]));
        if (e_pci)
        {
          out << i;
          break;
        }
      }
    }
    out << std::endl;
  }
}
