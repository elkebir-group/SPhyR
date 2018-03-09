/*
 * callback.h
 *
 *  Created on: 23-feb-2018
 *      Author: M. El-Kebir
 */

#ifndef CALLBACK_H
#define CALLBACK_H

#include <ilcplex/ilocplex.h>
#include "utils.h"

class UserCutCallback : public IloCplex::UserCutCallbackI
{
private:
  const int _m;
  const int _n;
  IloBoolVarArray _vars;
  
public:
  typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;
  
  UserCutCallback(IloEnv env,
                  const IloBoolVarMatrix& E,
                  const int m,
                  const int n)
    : UserCutCallbackI(env)
    , _m(m)
    , _n(n)
    , _vars()
  {
    _vars = IloBoolVarArray(env, _m * _n);
    for (int i = 0; i < _m; ++i)
    {
      for (int j = 0; j < _n; ++j)
      {
        _vars[i * _n + j] = E[i][j];
      }
    }
  }
  
  IloCplex::CallbackI *duplicateCallback() const
  {
    return (new (getEnv()) UserCutCallback(*this));
  }
  
  void main();
};

class LazyCutCallback : public IloCplex::LazyConstraintCallbackI
{
private:
  const int _m;
  const int _n;
  IloBoolVarArray _vars;
  
public:
  typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;
  
  LazyCutCallback(IloEnv env,
                  const IloBoolVarMatrix& E,
                  const int m,
                  const int n)
  : LazyConstraintCallbackI(env)
  , _m(m)
  , _n(n)
  , _vars()
  {
    _vars = IloBoolVarArray(env, _m * _n);
    for (int i = 0; i < _m; ++i)
    {
      for (int j = 0; j < _n; ++j)
      {
        _vars[i * _n + j] = E[i][j];
      }
    }
  }
  
  IloCplex::CallbackI *duplicateCallback() const
  {
    return (new (getEnv()) LazyCutCallback(*this));
  }
  
  void main();
};

#endif // CALLBACK_H
