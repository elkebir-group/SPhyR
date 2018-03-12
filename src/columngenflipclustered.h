/*
 * columngenflipclustered.h
 *
 *  Created on: 10-mar-2018
 *      Author: M. El-Kebir
 */

#ifndef COLUMNGENFLIPCLUSTERED_H
#define COLUMNGENFLIPCLUSTERED_H

#include "columngenflip.h"

class ColumnGenFlipClustered : public ColumnGenFlip
{
public:
  /// Constructor
  ///
  /// @param B Input matrix
  /// @param k Maximum number of losses
  ColumnGenFlipClustered(const Matrix& B,
                         int k,
                         bool lazy,
                         double alpha,
                         double beta,
                         int l,
                         const StlIntVector& z);
  
protected: 
  virtual void initObjective();
  
  virtual void initActiveVariables();
  
protected:
  const int _l;
  /// Cluster assignment
  const StlIntVector& _z;
};

#endif // COLUMNGENFLIPCLUSTERED_H
