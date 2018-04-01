/*
 * columngenflip.h
 *
 *  Created on: 10-mar-2018
 *      Author: M. El-Kebir
 */

#ifndef COLUMNGENFLIP_H
#define COLUMNGENFLIP_H

#include <ilcplex/ilocplex.h>
#include "matrix.h"
#include "columngen.h"

/// This class provides a column generation approach to the k-DPF problem
class ColumnGenFlip : public ColumnGen
{
public:
  /// Constructor
  ///
  /// @param B Input matrix
  /// @param k Maximum number of losses per character
  /// @param lazy Introduce constraints into the lazy constraint pool
  /// @param alpha False positive rate
  /// @param beta False negative rate
  ColumnGenFlip(const Matrix& B,
                int k,
                bool lazy,
                double alpha,
                double beta);
  
protected:
  /// Hidden constructor where output matrix dimensions may differ from input matrix
  ///
  /// @param B Input matrix
  /// @param m Number of taxa in output matrix
  /// @param n Number of character in output matrix
  /// @param k Maximum number of losses per character
  /// @param lazy Add constraints to the lazy constraint pool instead of the main constraint pool
  /// @param alpha False positive rate
  /// @param beta False negative rate
  ColumnGenFlip(const Matrix& B,
                int m,
                int n,
                int k,
                bool lazy,
                double alpha,
                double beta);
  
  /// Initialize fixed entries
  virtual void initFixedEntriesConstraints();
  
  /// Initialize objective function
  virtual void initObjective();
  
  /// Activate variable
  ///
  /// @param p Taxon
  /// @param c Character
  /// @param i State
  virtual void activate(int p, int c, int i);
  
protected:
  /// False positive rate
  const double _alpha;
  /// False negative rate
  const double _beta;
};

#endif // COLUMNGENFLIP_H
