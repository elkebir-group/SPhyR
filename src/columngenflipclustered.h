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
                         int lC,
                         const StlIntVector& zC,
                         int lT,
                         const StlIntVector& zT);
  
  void initHotStart(const Matrix& E);
  
protected: 
  virtual void initObjective();
  
  virtual void initActiveVariables();
  
  virtual void initFixedColumns()
  {
  }
  
protected:
  /// SNV (character) cluster assignment
  const StlIntVector& _zC;
  /// Cell (taxa) cluster assignment
  const StlIntVector& _zT;
  /// Number of SNV (character) clusters
  const double _lC;
  /// Number of cell (taxa) clusters
  const double _lT;
};

#endif // COLUMNGENFLIPCLUSTERED_H
