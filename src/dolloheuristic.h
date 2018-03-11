/*
 * dolloheuristic.h
 *
 *  Created on: 9-mar-2018
 *      Author: M. El-Kebir
 */

#ifndef DOLLOHEURISTIC_H
#define DOLLOHEURISTIC_H

#include <ilcplex/ilocplex.h>
#include <ilconcert/ilothread.h>
#include "utils.h"

class DolloHeuristic : public IloCplex::HeuristicCallbackI
{
public:
  typedef IloArray<IloArray<IloBoolVarArray> > IloBoolVar3Matrix;
  
  DolloHeuristic(IloEnv env,
                 const IloBoolVar3Matrix& E,
                 const int m,
                 const int n,
                 const int k,
                 IloFastMutex* pMutex);
  
  IloCplex::CallbackI *duplicateCallback() const
  {
    return (new (getEnv()) DolloHeuristic(*this));
  }
  
  void main();
  
private:
  const int _m;
  const int _n;
  const int _k;
  IloBoolVarArray _vars;
  int _maxIterations;
  int _currentIterations;
  IloCplex::MIPCallbackI::NodeId _nodeId;
  IloFastMutex* _pMutex;
  
  int getIndex(int p, int c, int i) const
  {
    return (_n * (_k + 2)) * p + (_k + 2) * c + i;
  }
};

#endif // DOLLOHEURISTIC_H

