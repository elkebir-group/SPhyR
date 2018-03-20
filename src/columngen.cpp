/*
 * columngen.cpp
 *
 *  Created on: 10-mar-2018
 *      Author: M. El-Kebir
 */

#include "columngen.h"
#include <ilconcert/ilothread.h>
#include "dollocallback.h"

ColumnGen::ColumnGen(const Matrix& B,
                     int k,
                     bool lazy)
  : _B(B)
  , _m(B.getNrTaxa())
  , _n(B.getNrCharacters())
  , _k(k)
  , _lazy(lazy)
  , _env()
  , _model(_env)
  , _cplex(_model)
  , _A(_env)
  , _vars(_env)
  , _obj(_env)
  , _activeVariables()
  , _nrActiveVariables(0)
  , _nrConstraints(0)
  , _solA(_B.getNrTaxa(), _B.getNrCharacters())
{
}

ColumnGen::ColumnGen(const Matrix& B,
                     int m,
                     int n,
                     int k,
                     bool lazy)
  : _B(B)
  , _m(m)
  , _n(n)
  , _k(k)
  , _lazy(lazy)
  , _env()
  , _model(_env)
  , _cplex(_model)
  , _A(_env)
  , _vars(_env)
  , _obj(_env)
  , _activeVariables()
  , _nrActiveVariables(0)
  , _nrConstraints(0)
  , _solA(m, n)
{
}

void ColumnGen::init()
{
  initVariables();
  initConstraints();
  initFixedColumns();
  initFixedEntriesConstraints();
  initActiveVariables();
  initObjective();
}

void ColumnGen::writeActiveVariables(std::ostream& out) const
{
  for (int p = 0; p < _m; ++p)
  {
    for (int c = 0; c < _n; ++c)
    {
      if (c > 0)
      {
        out << " ";
      }
      for (int i = 0; i <= _k + 1; ++i)
      {
        if (_activeVariables[p][c][i])
        {
          out << i;
        }
      }
    }
    out << std::endl;
  }
  out << std::endl;
}

void ColumnGen::initActiveVariables()
{
  const int m = _B.getNrTaxa();
  const int n = _B.getNrCharacters();
  
  _nrActiveVariables = 0;
  _activeVariables = StlBool3Matrix(m, StlBoolMatrix(_n, StlBoolVector(_k + 2, false)));
  
  for (int p = 0; p < m; p++)
  {
    for (int c = 0; c < n; c++)
    {
      if (_B.getEntry(p, c) == 0)
      {
        _activeVariables[p][c][0] = true;
        ++_nrActiveVariables;
      }
      else if (_B.getEntry(p, c) == 1)
      {
        _activeVariables[p][c][1] = true;
        ++_nrActiveVariables;
      }
      else
      {
        _activeVariables[p][c][0] = true;
        ++_nrActiveVariables;
        
        _activeVariables[p][c][1] = true;
        ++_nrActiveVariables;
      }
    }
  }
  updateVariableBounds();
}

void ColumnGen::initObjective()
{
  double factor = 1. / (_m * _n);
  
  for (int p = 0; p < _m; ++p)
  {
    for (int c = 0; c < _n; ++c)
    {
      for (int i = 2; i <= _k + 1; ++i)
      {
        _obj += pow(factor, _k + 1 - i) * _A[p][c][i];
      }
    }
  }
  
  _obj *= 1000;
  
  _model.add(IloMinimize(_env, _obj));
}

void ColumnGen::updateVariableBounds()
{
  for (int p = 0; p < _m; p++)
  {
    for (int c = 0; c < _n; c++)
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

void ColumnGen::initFixedColumns()
{
  for (int c = 0; c < _n; c++)
  {
    int nrZeros = 0;
    int nrOnes = 0;
    int nrMissing = 0;
    
    for (int p = 0; p < _m; p++)
    {
      int b_pc = _B.getEntry(p, c);
      if (b_pc == 0)
      {
        ++nrZeros;
      }
      else if (b_pc == 1)
      {
        ++nrOnes;
      }
      else
      {
        ++nrMissing;
      }
    }
    
    if (nrOnes == 1)
    {
      for (int p = 0; p < _m; p++)
      {
        if (_B.getEntry(p, c) == 1)
        {
          _model.add(_A[p][c][1] == 1);
        }
        else
        {
          _model.add(_A[p][c][0] == 1);
        }
      }
    }
    else if (nrOnes + nrMissing == _m)
    {
      for (int p = 0; p < _m; p++)
      {
        _model.add(_A[p][c][1] == 1);
      }
    }
    else if (nrZeros + nrMissing == _m)
    {
      for (int p = 0; p < _m; p++)
      {
        _model.add(_A[p][c][0] == 1);
      }
    }
  }
}

void ColumnGen::initFixedEntriesConstraints()
{
  for (int p = 0; p < _m; p++)
  {
    for (int c = 0; c < _n; c++)
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
  IloExpr sum(_env);
  
  // Each entry has a unique state
  for (int p = 0; p < _m; p++)
  {
    for (int c = 0; c < _n; c++)
    {
      for (int i = 0; i <= _k + 1; ++i)
      {
        sum += _A[p][c][i];
      }
      _model.add(sum == 1);
      sum.clear();
    }
  }
  
  IloExpr prevSum(_env);
  for (int c = 0; c < _n; c++)
  {
    for (int i = 2; i <= _k + 1; ++i)
    {
      for (int p = 0; p < _m; p++)
      {
        sum += _A[p][c][i];
      }
      if (i > 2)
      {
        _model.add(prevSum >= sum);
      }
      
      std::swap(prevSum, sum);
      sum.clear();
    }
  }
  
  sum.end();
  prevSum.end();
}

void ColumnGen::initVariables()
{
  char buf[1024];
  
  _A = IloBoolVar3Matrix(_env, _m);
  _vars = IloBoolVarArray(_env, _m * _n * (_k + 2));
  for (int p = 0; p < _m; p++)
  {
    _A[p] = IloBoolVarMatrix(_env, _n);
    for (int c = 0; c < _n; ++c)
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
  if (_B.getEntry(p, c)) return;
  
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

int ColumnGen::separate()
{
  IloNumArray vals = IloNumArray(_env, _vars.getSize());
  _cplex.getValues(vals, _vars);
  
  ViolatedConstraintList constraints;
  
  std::set<IntPair> activationSet;
  typedef std::set<int> IntSet;
  
  bool res = false;
  for (int c = 0; c < _n; ++c)
  {
    for (int d = c + 1; d < _n; ++d)
    {
      // condition 1
      for (int i = 1; i <= _k + 1; ++i)
      {
        for (int i_prime = 1; i_prime <= _k + 1; ++i_prime)
        {
          for (int j = 1; j <= _k + 1; ++j)
          {
            for (int j_prime = 1; j_prime <= _k + 1; ++j_prime)
            {
              IntSet pSet;
              for (int p = 0; p < _m; p++)
              {
                if (g_tol.nonZero(vals[getIndex(p, c, i)]) && g_tol.nonZero(vals[getIndex(p, d, 0)]))
                {
                  pSet.insert(p);
                }
              }
              
              IntSet qSet;
              for (int q = 0; q < _m; q++)
              {
                if (g_tol.nonZero(vals[getIndex(q, c, 0)]) && g_tol.nonZero(vals[getIndex(q, d, j)]))
                {
                  qSet.insert(q);
                }
              }
              
              IntSet rSet;
              for (int r = 0; r < _m; r++)
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
                    activate(p, c, i);
                    activate(p, d, 0);
                    activate(q, c, 0);
                    activate(q, d, j);
                    activate(r, c, i_prime);
                    activate(r, d, j_prime);
                    
                    ViolatedConstraint constraint;
                    constraint[0] = Triple(p, c, i);
                    constraint[1] = Triple(p, d, 0);
                    constraint[2] = Triple(q, c, 0);
                    constraint[3] = Triple(q, d, j);
                    constraint[4] = Triple(r, c, i_prime);
                    constraint[5] = Triple(r, d, j_prime);
                    constraints.push_back(constraint);
                  }
                }
              }
            }
          }
        }
      }
      
      // condition 2
      for (int i = 1; i <= _k + 1; ++i)
      {
        for (int i_prime = 1; i_prime <= _k + 1; ++i_prime)
        {
          for (int j = 2; j <= _k + 1; ++j)
          {
            for (int j_prime = 1; j_prime <= _k + 1; ++j_prime)
            {
              if (j_prime == j) continue;
          
              IntSet pSet;
              for (int p = 0; p < _m; p++)
              {
                if (g_tol.nonZero(vals[getIndex(p, c, i)]) && g_tol.nonZero(vals[getIndex(p, d, j_prime)]))
                {
                  pSet.insert(p);
                }
              }
              
              IntSet qSet;
              for (int q = 0; q < _m; q++)
              {
                if (g_tol.nonZero(vals[getIndex(q, c, 0)]) && g_tol.nonZero(vals[getIndex(q, d, j)]))
                {
                  qSet.insert(q);
                }
              }
              
              IntSet rSet;
              for (int r = 0; r < _m; r++)
              {
                if (g_tol.nonZero(vals[getIndex(r, c, i_prime)]) && g_tol.nonZero(vals[getIndex(r, d, j)]))
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
                    activate(p, d, j_prime);
                    activate(q, c, 0);
                    activate(q, d, j);
                    activate(r, c, i_prime);
                    activate(r, d, j);
                    
                    ViolatedConstraint constraint;
                    constraint[0] = Triple(p, c, i);
                    constraint[1] = Triple(p, d, j_prime);
                    constraint[2] = Triple(q, c, 0);
                    constraint[3] = Triple(q, d, j);
                    constraint[4] = Triple(r, c, i_prime);
                    constraint[5] = Triple(r, d, j);
                    constraints.push_back(constraint);
                  }
                }
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
          if (i_prime == i) continue;
          for (int j = 1; j <= _k + 1; ++j)
          {
            for (int j_prime = 1; j_prime <= _k + 1; ++j_prime)
            {
              IntSet pSet;
              for (int p = 0; p < _m; p++)
              {
                if (g_tol.nonZero(vals[getIndex(p, c, i)]) && g_tol.nonZero(vals[getIndex(p, d, 0)]))
                {
                  pSet.insert(p);
                }
              }
              
              IntSet qSet;
              for (int q = 0; q < _m; q++)
              {
                if (g_tol.nonZero(vals[getIndex(q, c, i_prime)]) && g_tol.nonZero(vals[getIndex(q, d, j)]))
                {
                  qSet.insert(q);
                }
              }
              
              IntSet rSet;
              for (int r = 0; r < _m; r++)
              {
                if (g_tol.nonZero(vals[getIndex(r, c, i)]) && g_tol.nonZero(vals[getIndex(r, d, j_prime)]))
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
                    activate(q, c, i_prime);
                    activate(q, d, j);
                    activate(r, c, i);
                    activate(r, d, j_prime);
                    
                    ViolatedConstraint constraint;
                    constraint[0] = Triple(p, c, i);
                    constraint[1] = Triple(p, d, 0);
                    constraint[2] = Triple(q, c, i_prime);
                    constraint[3] = Triple(q, d, j);
                    constraint[4] = Triple(r, c, i);
                    constraint[5] = Triple(r, d, j_prime);
                    constraints.push_back(constraint);
                  }
                }
              }
            }
          }
        }
      }
      
      // condition 4
      for (int i = 2; i <= _k + 1; ++i)
      {
        for (int i_prime = 1; i_prime <= _k + 1; ++i_prime)
        {
          if (i_prime == i) continue;
          for (int j = 2; j <= _k + 1; ++j)
          {
            for (int j_prime = 1; j_prime <= _k + 1; ++j_prime)
            {
              if (j_prime == j) continue;
              
              IntSet pSet;
              for (int p = 0; p < _m; p++)
              {
                if (g_tol.nonZero(vals[getIndex(p, c, i)]) && g_tol.nonZero(vals[getIndex(p, d, j_prime)]))
                {
                  pSet.insert(p);
                }
              }
              
              IntSet qSet;
              for (int q = 0; q < _m; q++)
              {
                if (g_tol.nonZero(vals[getIndex(q, c, i_prime)]) && g_tol.nonZero(vals[getIndex(q, d, j)]))
                {
                  qSet.insert(q);
                }
              }
              
              IntSet rSet;
              for (int r = 0; r < _m; r++)
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
                    activate(p, d, j_prime);
                    activate(q, c, i_prime);
                    activate(q, d, j);
                    activate(r, c, i);
                    activate(r, d, j);
                    
                    ViolatedConstraint constraint;
                    constraint[0] = Triple(p, c, i);
                    constraint[1] = Triple(p, d, j_prime);
                    constraint[2] = Triple(q, c, i_prime);
                    constraint[3] = Triple(q, d, j);
                    constraint[4] = Triple(r, c, i);
                    constraint[5] = Triple(r, d, j);
                    constraints.push_back(constraint);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  IloConstraintArray modelConstraints(_env);
  IloConstraintArray lazyConstraints(_env);
  IloExpr sum(_env);
  for (const ViolatedConstraint& violatedConstraint : constraints)
  {
    for (int idx = 0; idx < 6; ++idx)
    {
      sum += _vars[getIndex(violatedConstraint[idx])];
    }
    if (!_lazy)
    {
      modelConstraints.add(sum <= 5);
    }
    else
    {
      lazyConstraints.add(sum <= 5);
    }
    sum.clear();
  }
  
  if (lazyConstraints.getSize() > 0)
  {
//    std::cerr << "Added " << lazyConstraints.getSize() << " lazy constraints" << std::endl;
    _cplex.addLazyConstraints(lazyConstraints);
    _cplex.addUserCuts(lazyConstraints);
    
    return lazyConstraints.getSize();
  }
  if (modelConstraints.getSize() > 0)
  {
//    std::cerr << "Added " << modelConstraints.getSize() << " model constraints" << std::endl;
    _model.add(modelConstraints);
    
    return modelConstraints.getSize();
  }
  
  return 0;
}

void ColumnGen::processSolution()
{
  for (int p = 0; p < _m; p++)
  {
    for (int c = 0; c < _n; ++c)
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
  
//  const int m = _B.getNrTaxa();
//  IloFastMutex mutex;
//  _cplex.use(IloCplex::Callback(new (_env) DolloCallback<IloCplex::UserCutCallbackI>(_env, _A, m, _n, _k, &mutex)));
//  _cplex.use(IloCplex::Callback(new (_env) DolloCallback<IloCplex::LazyConstraintCallbackI>(_env, _A, m, _n, _k, &mutex)));
  
  _nrConstraints = _cplex.getNrows();
  
  int iteration = 1;
  bool res = false;
  while (true)
  {
    std::cerr << "Step " << iteration << " -- elapsed time " << g_timer.realTime() << " s" << std::endl;
    std::cerr << "Step " << iteration << " -- number of constraints: " << _nrConstraints << std::endl;
    std::cerr << "Step " << iteration << " -- number of active variables: " << _nrActiveVariables << std::endl;
    
    if (timeLimit != -1 && g_timer.realTime() > timeLimit)
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
    
    int separatedConstraints = separate();
    _nrConstraints += separatedConstraints;
    std::cerr << "Step " << iteration << " -- introduced " << separatedConstraints << " constraints" << std::endl;
    if (separatedConstraints == 0)
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
  std::cerr << "Elapsed time: " << g_timer.realTime() << std::endl;
  
  
  return res;
}
