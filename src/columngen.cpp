/*
 * columngen.cpp
 *
 *  Created on: 10-mar-2018
 *      Author: M. El-Kebir
 */

#include "columngen.h"
#include <ilconcert/ilothread.h>
#include <lemon/time_measure.h>

ColumnGen::ColumnGen(const Matrix& B,
                     int k)
  : _B(B)
  , _k(k)
  , _env()
  , _model(_env)
  , _cplex(_model)
  , _A(_env)
  , _vars(_env)
  , _activeVariables()
  , _nrActiveVariables(0)
  , _solA(_B.getNrTaxa(), _B.getNrCharacters())
{
}

void ColumnGen::init()
{
  initVariables();
  initConstraints();
  initFixedEntriesConstraints();
  initActiveVariables();
  initObjective();
}

void ColumnGen::initActiveVariables()
{
  const int m = _B.getNrTaxa();
  const int n = _B.getNrCharacters();
  
  _nrActiveVariables = 0;
  _activeVariables = StlBool3Matrix(m, StlBoolMatrix(n, StlBoolVector(_k + 2, false)));
  
  for (int p = 0; p < m; p++)
  {
    for (int c = 0; c < n; c++)
    {
      if (_B.getEntry(p, c) == 0)
      {
        _activeVariables[p][c][0] = true;
        ++_nrActiveVariables;
      }
      else
      {
        _activeVariables[p][c][1] = true;
        ++_nrActiveVariables;
      }
    }
  }
  updateVariableBounds();
}

void ColumnGen::initObjective()
{
  const int m = _B.getNrTaxa();
  const int n = _B.getNrCharacters();
  
  double factor = 1. / (m * n);
  
  IloExpr obj(_env);
  for (int p = 0; p < m; ++p)
  {
    for (int c = 0; c < n; ++c)
    {
      for (int i = 2; i <= _k + 1; ++i)
      {
        obj += pow(factor, _k + 1 - i) * _A[p][c][i];
      }
    }
  }
  
  _model.add(IloMinimize(_env, 1000 * obj));
}

void ColumnGen::updateVariableBounds()
{
  const int m = _B.getNrTaxa();
  const int n = _B.getNrCharacters();
  
  for (int p = 0; p < m; p++)
  {
    for (int c = 0; c < n; c++)
    {
      for (int i = 0; i <= _k + 1; ++i)
      {
        if (_activeVariables[p][c][i])
        {
          _A[p][c][i].setUB(1);
        }
        else
        {
          _A[p][c][i].setUB(0);
        }
      }
    }
  }
}

void ColumnGen::initFixedEntriesConstraints()
{
  const int m = _B.getNrTaxa();
  const int n = _B.getNrCharacters();
  
  for (int p = 0; p < m; p++)
  {
    for (int c = 0; c < n; c++)
    {
      int b_pc = _B.getEntry(p, c);
      if (b_pc == 1)
      {
        _model.add(_A[p][c][1] == 1);
      }
      else if (b_pc == 0)
      {
        _model.add(_A[p][c][1] == 0);
      }
    }
  }
}

void ColumnGen::initConstraints()
{
  const int m = _B.getNrTaxa();
  const int n = _B.getNrCharacters();
  
  IloExpr sum(_env);
  
  // Each entry has a unique state
  for (int p = 0; p < m; p++)
  {
    for (int c = 0; c < n; c++)
    {
      for (int i = 0; i <= _k + 1; ++i)
      {
        sum += _A[p][c][i];
      }
      _model.add(sum == 1);
      sum.clear();
    }
  }
  
  sum.end();
}

void ColumnGen::initVariables()
{
  const int m = _B.getNrTaxa();
  const int n = _B.getNrCharacters();
  
  char buf[1024];
  
  _A = IloBoolVar3Matrix(_env, m);
  _vars = IloBoolVarArray(_env, m * n * (_k + 2));
  for (int p = 0; p < m; p++)
  {
    _A[p] = IloBoolVarMatrix(_env, n);
    for (int c = 0; c < n; ++c)
    {
      _A[p][c] = IloBoolVarArray(_env, _k + 2);
      for (int i = 0; i < _k + 2; ++i)
      {
        snprintf(buf, 1024, "a_%d_%d_%d", p, c, i);
        _A[p][c][i] = IloBoolVar(_env, buf);
        _vars[getIndex(p, c, i)] = _A[p][c][i];
      }
    }
  }
}

void ColumnGen::activate(int p, int c, int i)
{
  assert(_B.getEntry(p, c) == 0);
  if (i == 0 && _k >= 1)
  {
    if (!_activeVariables[p][c][2])
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

bool ColumnGen::separate()
{
  const int m = _B.getNrTaxa();
  const int n = _B.getNrCharacters();
  
  IloNumArray vals = IloNumArray(_env, _vars.getSize());
  _cplex.getValues(vals, _vars);
  
  std::set<IntPair> activationSet;
  typedef std::set<int> IntSet;
  
  bool res = false;
  for (int c = 0; c < n; ++c)
  {
    for (int d = c + 1; d < n; ++d)
    {
      // condition 1
      for (int j = 2; j <= _k + 1; ++j)
      {
        for (int j_prime = 1; j_prime <= _k + 1; ++j_prime)
        {
          if (j_prime == j) continue;
          
          IntSet pSet;
          for (int p = 0; p < m; p++)
          {
            if (g_tol.nonZero(vals[getIndex(p, c, 1)]) && g_tol.nonZero(vals[getIndex(p, d, j_prime)]))
            {
              pSet.insert(p);
            }
          }
          
          IntSet qSet;
          for (int q = 0; q < m; q++)
          {
            if (g_tol.nonZero(vals[getIndex(q, c, 0)]) && g_tol.nonZero(vals[getIndex(q, d, j)]))
            {
              qSet.insert(q);
            }
          }
          
          IntSet rSet;
          for (int r = 0; r < m; r++)
          {
            if (g_tol.nonZero(vals[getIndex(r, c, 1)]) && g_tol.nonZero(vals[getIndex(r, d, j)]))
            {
              rSet.insert(r);
            }
          }
          
          for (int p : pSet)
          {
            for (int q : qSet)
            {
              for (int r : rSet)
              {
                activate(q, c, 0);
                if (j_prime != 1)
                {
                  activate(p, d, j_prime);
                }
                activate(q, d, j);
                activate(r, d, j);
                
                _model.add(_vars[getIndex(p, c, 1)] + _vars[getIndex(p, d, j_prime)] +
                           _vars[getIndex(q, c, 0)] + _vars[getIndex(q, d, j)] +
                           _vars[getIndex(r, c, 1)] + _vars[getIndex(r, d, j)] <= 5);
                
                res = true;
              }
            }
          }
        }
      }
      
      // condition 2
      for (int i = 2; i <= _k + 1; ++i)
      {
        for (int i_prime = 1; i_prime <= _k + 1; ++i_prime)
        {
          if (i == i_prime) continue;
          
          IntSet pSet;
          for (int p = 0; p < m; p++)
          {
            if (g_tol.nonZero(vals[getIndex(p, c, i)]) && g_tol.nonZero(vals[getIndex(p, d, 0)]))
            {
              pSet.insert(p);
            }
          }
          
          IntSet qSet;
          for (int q = 0; q < m; q++)
          {
            if (g_tol.nonZero(vals[getIndex(q, c, i_prime)]) && g_tol.nonZero(vals[getIndex(q, d, 1)]))
            {
              qSet.insert(q);
            }
          }
          
          IntSet rSet;
          for (int r = 0; r < m; r++)
          {
            if (g_tol.nonZero(vals[getIndex(r, c, i)]) && g_tol.nonZero(vals[getIndex(r, d, 1)]))
            {
              rSet.insert(r);
            }
          }
          
          for (int p : pSet)
          {
            for (int q : qSet)
            {
              for (int r : rSet)
              {
                activate(p, c, i);
                activate(p, d, 0);
                if (i_prime != 1)
                {
                  activate(q, c, i_prime);
                }
                activate(r, c, i);
                
                _model.add(_vars[getIndex(p, c, i)] + _vars[getIndex(p, d, 0)] +
                           _vars[getIndex(q, c, i_prime)] + _vars[getIndex(q, d, 1)] +
                           _vars[getIndex(r, c, i)] + _vars[getIndex(r, d, 1)] <= 5);
                
                res = true;
              }
            }
          }
        }
      }
      
      // condition 3
      for (int i = 2; i <= _k + 1; ++i)
      {
        for (int i_prime = 1; i_prime <= _k + 1; ++i_prime)
        {
          if (i == i_prime) continue;
          for (int j = 2; j <= _k + 1; ++j)
          {
            for (int j_prime = 1; j_prime <= _k + 1; ++j_prime)
            {
              if (j == j_prime) continue;
              
              IntSet pSet;
              for (int p = 0; p < m; p++)
              {
                if (g_tol.nonZero(vals[getIndex(p, c, i)]) && g_tol.nonZero(vals[getIndex(p, d, j_prime)]))
                {
                  pSet.insert(p);
                }
              }
              
              IntSet qSet;
              for (int q = 0; q < m; q++)
              {
                if (g_tol.nonZero(vals[getIndex(q, c, i_prime)]) && g_tol.nonZero(vals[getIndex(q, d, j)]))
                {
                  qSet.insert(q);
                }
              }
              
              IntSet rSet;
              for (int r = 0; r < m; r++)
              {
                if (g_tol.nonZero(vals[getIndex(r, c, i)]) && g_tol.nonZero(vals[getIndex(r, d, j)]))
                {
                  rSet.insert(r);
                }
              }
              
              for (int p : pSet)
              {
                for (int q : qSet)
                {
                  for (int r : rSet)
                  {
                    activate(p, c, i);
                    if (j_prime != 1)
                    {
                      activate(p, d, j_prime);
                    }
                    if (i_prime != 1)
                    {
                      activate(q, c, i_prime);
                    }
                    activate(r, c, i);
                    activate(r, d, j);
                    
                    _model.add(_vars[getIndex(p, c, i)] + _vars[getIndex(p, d, j_prime)] +
                               _vars[getIndex(q, c, i_prime)] + _vars[getIndex(q, d, j)] +
                               _vars[getIndex(r, c, i)] + _vars[getIndex(r, d, j)] <= 5);
                    
                    res = true;
                  }
                }
              }
            }
          }
        }
      }
      
      // condition 4
      for (int i = 1; i <= _k + 1; ++i)
      {
        for (int i_prime = 1; i_prime <= _k + 1; ++i_prime)
        {
          for (int j = 1; j <= _k + 1; ++j)
          {
            for (int j_prime = 1; j_prime <= _k + 1; ++j_prime)
            {
              IntSet pSet;
              for (int p = 0; p < m; p++)
              {
                if (g_tol.nonZero(vals[getIndex(p, c, i)]) && g_tol.nonZero(vals[getIndex(p, d, 0)]))
                {
                  pSet.insert(p);
                }
              }
              
              IntSet qSet;
              for (int q = 0; q < m; q++)
              {
                if (g_tol.nonZero(vals[getIndex(q, c, 0)]) && g_tol.nonZero(vals[getIndex(q, d, j)]))
                {
                  qSet.insert(q);
                }
              }
              
              IntSet rSet;
              for (int r = 0; r < m; r++)
              {
                if (g_tol.nonZero(vals[getIndex(r, c, i_prime)]) && g_tol.nonZero(vals[getIndex(r, d, j_prime)]))
                {
                  rSet.insert(r);
                }
              }
              
              for (int p : pSet)
              {
                for (int q : qSet)
                {
                  for (int r : rSet)
                  {
                    if (i != 1)
                    {
                      activate(p, c, i);
                    }
                    activate(p, d, 0);
                    activate(q, c, 0);
                    if (j != 1)
                    {
                      activate(q, d, j);
                    }
                    if (i_prime != 1)
                    {
                      activate(r, c, i_prime);
                    }
                    if (j_prime != 1)
                    {
                      activate(r, d, j_prime);
                    }
                
                    _model.add(_vars[getIndex(p, c, i)] + _vars[getIndex(p, d, 0)] +
                               _vars[getIndex(q, c, 0)] + _vars[getIndex(q, d, j)] +
                               _vars[getIndex(r, c, i_prime)] + _vars[getIndex(r, d, j_prime)] <= 5);
                    
                    res = true;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  return res;
}

void ColumnGen::processSolution()
{
  const int m = _B.getNrTaxa();
  const int n = _B.getNrCharacters();
  for (int p = 0; p < m; p++)
  {
    for (int c = 0; c < n; ++c)
    {
      for (int i = 0; i <= _k + 1; ++i)
      {
        bool a_pci = g_tol.nonZero(_cplex.getValue(_A[p][c][i]));
        if (a_pci)
        {
          _solA.setEntry(p, c, i);
        }
      }
    }
  }
}

bool ColumnGen::solve(int timeLimit,
                      int memoryLimit,
                      int nrThreads,
                      bool verbose)
{
  if (!verbose)
  {
    _env.setOut(_env.getNullStream());
    _env.setError(_env.getNullStream());
    _env.setWarning(_env.getNullStream());
    _cplex.setOut(_env.getNullStream());
    _cplex.setError(_env.getNullStream());
    _cplex.setWarning(_env.getNullStream());
  }
  else
  {
    _env.setOut(std::cerr);
    _env.setError(std::cerr);
    _env.setWarning(std::cerr);
    _cplex.setOut(std::cerr);
    _cplex.setError(std::cerr);
    _cplex.setWarning(std::cerr);
  }
  
//  _cplex.setParam(IloCplex::MIPEmphasis, IloCplex::MIPEmphasisFeasibility);
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
  
  lemon::Timer timer;
  
  int iteration = 1;
  bool res = false;
  while (true)
  {
    std::cerr << "Iteration " << iteration << " -- elapsed time " << timer.realTime() << " s" << std::endl;
    std::cerr << "Iteration " << iteration << " -- number of constraints: " << _cplex.getNrows() << std::endl;
    std::cerr << "Iteration " << iteration << " -- number of active variables: " << _nrActiveVariables << std::endl;
    
    if (timeLimit != -1 && timer.realTime() > timeLimit)
    {
      std::cerr << "Time limit exceeded" << std::endl;
      return false;
    }
    
    _cplex.solve();
    if (_cplex.getStatus() != IloAlgorithm::Optimal)
    {
      res = false;
      break;
    }
    
    if (!separate())
    {
      res = true;
      break;
    }
    
    ++iteration;
  }
  
  if (res)
  {
    processSolution();
    std::cerr << "CPLEX: [" << _cplex.getObjValue() << " , " << _cplex.getBestObjValue() << "]" << std::endl;
  }
  std::cerr << "Elapsed time: " << timer.realTime() << std::endl;
  
  
  return res;
}
