/*
 * dolloheuristic.cpp
 *
 *  Created on: 9-mar-2018
 *      Author: M. El-Kebir
 */

#include "dolloheuristic.h"
#include "matrix.h"

DolloHeuristic::DolloHeuristic(IloEnv env,
                               const IloBoolVar3Matrix& E,
                               const int m,
                               const int n,
                               const int k,
                               IloFastMutex* pMutex)
  : IloCplex::HeuristicCallbackI(env)
  , _m(m)
  , _n(n)
  , _k(k)
  , _vars()
  , _maxIterations(100)
  , _currentIterations(0)
  , _nodeId()
  , _pMutex(pMutex)
{
  _vars = IloBoolVarArray(env, _m * _n * (_k + 2));
  
  for (int p = 0; p < _m; ++p)
  {
    for (int c = 0; c < _n; ++c)
    {
      for (int i = 0; i <= _k + 1; ++i)
      {
        _vars[getIndex(p, c, i)] = E[p][c][i];
      }
    }
  }
}

void DolloHeuristic::main()
{
  IloNumArray vals = IloNumArray(getEnv(), _vars.getSize());
  
  _pMutex->lock();
  getValues(vals, _vars);
  _pMutex->unlock();
  
  Matrix A(_m, _n);
  
  std::map<IntPair, std::set<int> > options;
  
  for (int p = 0; p < _m; ++p)
  {
    for (int c = 0; c < _n; ++c)
    {
//      std::cout << "( ";
      for (int i = 0; i <= _k + 1; ++i)
      {
        double val = vals[getIndex(p, c, i)];
//        std::cout << val << " ";
        if (g_tol.less(val, 1.) && g_tol.nonZero(val))
        {
          options[IntPair(p, c)].insert(i);
        }
        else if (g_tol.nonZero(val))
        {
          A.setEntry(p, c, i);
        }
      }
//      std::cout << ") ";
    }
//    std::cout << std::endl;
  }
  
  std::cout << A;
  
  Matrix::ViolationList list;
  A.identifyViolations(_k, list);
  
  IntPairSet violationEntries;
  for (const Matrix::Violation& violation : list)
  {
    std::cout << "Condition " << violation._condition
              << " ; "<< "p = " << violation._p << " ; q = " << violation._q
              << " ; r = " << violation._r << " ; c = " << violation._c
              << " ; d = " << violation._d << std::endl;
    std::cout << A.getEntry(violation._p, violation._c)
              << " " << A.getEntry(violation._p, violation._d) << std::endl;
    std::cout << A.getEntry(violation._q, violation._c)
              << " " << A.getEntry(violation._q, violation._d) << std::endl;
    std::cout << A.getEntry(violation._r, violation._c)
              << " " << A.getEntry(violation._r, violation._d) << std::endl;
    std::cout << std::endl;
    violationEntries.insert(IntPair(violation._p, violation._c));
    violationEntries.insert(IntPair(violation._q, violation._c));
    violationEntries.insert(IntPair(violation._r, violation._c));
    violationEntries.insert(IntPair(violation._p, violation._d));
    violationEntries.insert(IntPair(violation._q, violation._d));
    violationEntries.insert(IntPair(violation._r, violation._d));
  }
  
  std::cout << "Total number of violations: " << list.size() << std::endl;
  std::cout << "Violated entries: " << violationEntries.size()
            << " / " << A.getNrTaxa() * A.getNrCharacters()
            << " = " << (double) violationEntries.size() / (A.getNrTaxa() * A.getNrCharacters()) << std::endl;
  
  
  return;
}

