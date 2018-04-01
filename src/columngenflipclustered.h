/*
 * columngenflipclustered.h
 *
 *  Created on: 10-mar-2018
 *      Author: M. El-Kebir
 */

#ifndef COLUMNGENFLIPCLUSTERED_H
#define COLUMNGENFLIPCLUSTERED_H

#include "columngenflip.h"

/// This class provides a column generation approach for the k-DPFC subproblem where one is given the character and taxon clusters
class ColumnGenFlipClustered : public ColumnGenFlip
{
public:
  /// Constructor
  ///
  /// @param B Input matrix
  /// @param multiplicities The number of orginal entries that each entry b_pi corresponds to
  /// @param baseL Base log likelihood of removed characters and taxa of original input matrix
  /// @param k Maximum number of losses per character
  /// @param lazy Introduce constraints into the lazy constraint pool
  /// @param alpha False positive rate
  /// @param beta False negative rate
  /// @param t Number of character clusters
  /// @param zC Character cluster assignment
  /// @param s Number of taxon clusters
  /// @param zT Taxon cluster assignment
  ColumnGenFlipClustered(const Matrix& B,
                         const StlIntMatrix& multiplicities,
                         double baseL,
                         int k,
                         bool lazy,
                         double alpha,
                         double beta,
                         int t,
                         const StlIntVector& zC,
                         int s,
                         const StlIntVector& zT);
  
  void initHotStart(const Matrix& E);
  
protected:
  /// Initialize objective function
  virtual void initObjective();
  
  /// Initialize active variables
  virtual void initActiveVariables();
  
  /// Initialize fixed columns (there are none!)
  virtual void initFixedColumns()
  {
  }
  
protected:
  /// The number of orginal entries that each entry b_pi corresponds to
  const StlIntMatrix& _multiplicities;
  /// Base log likelihood of removed characters and taxa of original input matrix
  const double _baseL;
  /// Taxon cluster assignment
  const StlIntVector& _zT;
  /// Character cluster assignment
  const StlIntVector& _zC;
  /// Number of cell (taxa) clusters
  const double _s;
  /// Number of SNV (character) clusters
  const double _t;
};

#endif // COLUMNGENFLIPCLUSTERED_H
