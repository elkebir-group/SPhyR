/*
 * ilpsolverdollo.h
 *
 *  Created on: 25-feb-2018
 *      Author: M. El-Kebir
 */

#ifndef ILPSOLVERDOLLO_H
#define ILPSOLVERDOLLO_H

#include <ilcplex/ilocplex.h>
#include "matrix.h"

/// Dollo phylogeny without errors
class IlpSolverDollo
{
public:
  /// Constructor
  ///
  /// @param D Input matrix
  /// @param k Maximum number of losses
  IlpSolverDollo(const Matrix& D,
                 int k);
  
  /// Destructor
  virtual ~IlpSolverDollo()
  {
    _cplex.end();
    _model.end();
    _env.end();
  }
  
  /// Solve
  void solve(int timeLimit,
             int memoryLimit,
             int nrThreads);
  
  /// Initialize solver
  virtual void init();
  
  /// Write solution
  ///
  /// @param out Output stream
  virtual void printSolution(std::ostream& out) const;
  
  /// Return solution matrix
  const Matrix& getSolE() const
  {
    return _solE;
  }
  
protected:
  /// Constructor
  ///
  /// @param D Input matrix
  /// @param k Maximum number of losses
  /// @param n Number of characters
  IlpSolverDollo(const Matrix& D,
                 int k,
                 int n);
  
  /// Initialize variables
  virtual void initVariables();
  
  /// Initialize constraints
  virtual void initConstraints();
  
  /// Initialize contraints corresponding to fixed entries
  virtual void initFixedEntries();
  
  /// Initialize objective function (number of losses)
  virtual void initObjective();
  
  /// Extract cplex solution
  virtual void processSolution();
  
  typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;
  typedef IloArray<IloNumVarArray> IloNumVarMatrix;
  typedef std::map<IntPair, IloNumVar> IloNumVarMatrixMap;
  typedef IloArray<IloBoolVarMatrix> IloBoolVar3Matrix;

protected:
  /// Input matrix
  const Matrix& _D;
  /// Maximum number of losses
  const int _k;
  /// Number of characters
  const int _n;
  /// Cplex environment
  IloEnv _env;
  /// Cplex model
  IloModel _model;
  /// Cplex solver
  IloCplex _cplex;
  /// _E[p][c][i] is the multi-state matrix
  IloBoolVar3Matrix _E;
  /// Solution matrix
  Matrix _solE;
};

#endif // ILPSOLVERDOLLO_H
