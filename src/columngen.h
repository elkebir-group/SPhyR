/*
 * columngen.h
 *
 *  Created on: 10-mar-2018
 *      Author: M. El-Kebir
 */

#ifndef COLUMNGEN_H
#define COLUMNGEN_H

#include <ilcplex/ilocplex.h>
#include "matrix.h"

class ColumnGen
{
public:
  /// Constructor
  ///
  /// @param B Input matrix
  /// @param k Maximum number of losses
  ColumnGen(const Matrix& B,
            int k,
            bool lazy);
  
  /// Initialize solver
  virtual void init();
  
  /// Set upper bound on objective
  virtual void setObjUpperBound(double UB)
  {
    _model.add(_obj <= UB);
  }
  
  /// Return solution matrix
  const Matrix& getSolA() const
  {
    return _solA;
  }
  
  bool solve(int timeLimit,
             int memoryLimit,
             int nrThreads,
             bool verbose);
  
protected:
  ColumnGen(const Matrix& B,
            int n,
            int k,
            bool lazy);
  
  virtual void initActiveVariables();
  
  virtual void initVariables();
  
  virtual void initConstraints();
  
  virtual void initFixedColumns();
  
  virtual void initFixedEntriesConstraints();
  
  virtual void initObjective();
  
  void updateVariableBounds();
  
  virtual void activate(int p, int c, int i);
  
  void processSolution();
  
  int separate();
  
  typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;
  typedef IloArray<IloBoolVarMatrix> IloBoolVar3Matrix;
  
  typedef std::vector<StlBoolVector> StlBoolMatrix;
  typedef std::vector<StlBoolMatrix> StlBool3Matrix;
  
  int getIndex(int p, int c, int i) const
  {
    return (_n * (_k + 2)) * p + (_k + 2) * c + i;
  }
  
  struct Triple
  {
  public:
    Triple(int p, int c, int i)
      : _p(p)
      , _c(c)
      , _i(i)
    {
    }
    
    Triple()
      : _p(-1)
      , _c(-1)
      , _i(-1)
    {
    }
    
    int _p;
    int _c;
    int _i;
  };
  
  int getIndex(const Triple& triple) const
  {
    return getIndex(triple._p, triple._c, triple._i);
  }
  
  typedef std::array<Triple, 6> ViolatedConstraint;
  
  typedef std::list<ViolatedConstraint> ViolatedConstraintList;
  
protected:
  /// Input matrix
  const Matrix& _B;
  ///
  const int _n;
  /// Maximum number of losses
  const int _k;
  /// Lazy constraints
  const bool _lazy;
  /// Cplex environment
  IloEnv _env;
  /// Cplex model
  IloModel _model;
  /// Cplex solver
  IloCplex _cplex;
  /// _A[p][c][i] is the multi-state matrix
  IloBoolVar3Matrix _A;
  /// Flatten variable matrix _A
  IloBoolVarArray _vars;
  /// Objective function
  IloExpr _obj;
  /// Indicates which variables are active
  StlBool3Matrix _activeVariables;
  /// Number of active variables
  int _nrActiveVariables;
  /// Number of constraints
  int _nrConstraints;
  /// Solution matrix
  Matrix _solA;
};

#endif // COLUMNGEN_H
