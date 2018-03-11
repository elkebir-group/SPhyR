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
            int k);
  
  /// Initialize solver
  virtual void init();
  
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
  virtual void initActiveVariables();
  
  virtual void initVariables();
  
  virtual void initConstraints();
  
  virtual void initFixedEntriesConstraints();
  
  virtual void initObjective();
  
  void updateVariableBounds();
  
  virtual void activate(int p, int c, int i);
  
  void processSolution();
  
  bool separate();
  
  typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;
  typedef IloArray<IloBoolVarMatrix> IloBoolVar3Matrix;
  
  typedef std::vector<StlBoolVector> StlBoolMatrix;
  typedef std::vector<StlBoolMatrix> StlBool3Matrix;
  
  int getIndex(int p, int c, int i) const
  {
    const int n = _B.getNrCharacters();
    return (n * (_k + 2)) * p + (_k + 2) * c + i;
  }
  
protected:
  /// Input matrix
  const Matrix& _B;
  /// Maximum number of losses
  const int _k;
  /// Cplex environment
  IloEnv _env;
  /// Cplex model
  IloModel _model;
  /// Cplex solver
  IloCplex _cplex;
  /// _A[p][c][i] is the multi-state matrix
  IloBoolVar3Matrix _A;
  
  IloBoolVarArray _vars;
  ///
  StlBool3Matrix _activeVariables;
  ///
  int _nrActiveVariables;
  /// Solution matrix
  Matrix _solA;
};

#endif // COLUMNGEN_H
