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

/// This class provides a column generation approach for the k-DP problem
class ColumnGen
{
public:
  /// Constructor
  ///
  /// @param B Input matrix
  /// @param k Maximum number of losses per character
  /// @param lazy Introduce constraints into the lazy constraint pool
  ColumnGen(const Matrix& B,
            int k,
            bool lazy);
  
  /// Initialize solver
  virtual void init();
  
  /// Return solution matrix
  const Matrix& getSolA() const
  {
    return _solA;
  }
  
  /// Solve
  ///
  /// @param timeLimit Time limit in seconds
  /// @param memoryLimit Memory limit in megabytes
  /// @param nrThreads Number of threads the solver can use
  /// @param verbose Set to true to enable ILP solver output
  bool solve(int timeLimit,
             int memoryLimit,
             int nrThreads,
             bool verbose);
  
protected:
  /// Hidden constructor where output matrix dimensions may differ from input matrix
  ///
  /// @param B Input matrix
  /// @param m Number of taxa in output matrix
  /// @param n Number of character in output matrix
  /// @param k Maximum number of losses per character
  /// @param lazy Add constraints to the lazy constraint pool instead of the main constraint pool
  ColumnGen(const Matrix& B,
            int m,
            int n,
            int k,
            bool lazy);
  
  /// Initialize active variables (with non-empty domain)
  virtual void initActiveVariables();
  
  /// Initialize model variables
  virtual void initVariables();
  
  /// Initialize model constraints
  virtual void initConstraints();
  
  /// Initialize fixed columns
  virtual void initFixedColumns();
  
  /// Initialize fixed entries
  virtual void initFixedEntriesConstraints();
  
  /// Initialize objective function
  virtual void initObjective();
  
  /// Update variable bounds
  void updateVariableBounds();
  
  /// Activate variable
  ///
  /// @param p Taxon
  /// @param c Character
  /// @param i State
  virtual void activate(int p, int c, int i);
  
  /// Extract solution from ILP solver
  void processSolution();
  
  /// Identify violated constraints
  int separate();
  
  typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;
  typedef IloArray<IloBoolVarMatrix> IloBoolVar3Matrix;
  
  typedef std::vector<StlBoolVector> StlBoolMatrix;
  typedef std::vector<StlBoolMatrix> StlBool3Matrix;
  
  /// Construct 1D index from (p,c,i) triple
  ///
  /// @param p Taxon
  /// @param c Character
  /// @param i State
  int getIndex(int p, int c, int i) const
  {
    return (_n * (_k + 2)) * p + (_k + 2) * c + i;
  }
  
  /// Triple (p,c,i)
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
  
  /// Construct 1D index from (p,c,i) triple
  ///
  /// @param triple (p,c,i)
  int getIndex(const Triple& triple) const
  {
    return getIndex(triple._p, triple._c, triple._i);
  }
  
  /// Write activate variables to specified output stream
  ///
  /// @param out Output stream
  void writeActiveVariables(std::ostream& out) const;
  
  /// Forbidden submatrix
  typedef std::array<Triple, 6> ViolatedConstraint;
  
  /// List of forbidden submatrices
  typedef std::list<ViolatedConstraint> ViolatedConstraintList;
  
protected:
  /// Input matrix
  const Matrix& _B;
  /// Number of taxa
  const int _m;
  /// Number of characters
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
