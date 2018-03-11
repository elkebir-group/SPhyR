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

class ColumnGenFlip : public ColumnGen
{
public:
  /// Constructor
  ///
  /// @param B Input matrix
  /// @param k Maximum number of losses
  ColumnGenFlip(const Matrix& B,
                int k,
                double alpha,
                double beta);
protected:
  virtual void initFixedEntriesConstraints();
  
  virtual void initObjective();
  
  virtual void activate(int p, int c, int i);
  
protected:
  double _alpha;
  double _beta;
};

#endif // COLUMNGENFLIP_H
